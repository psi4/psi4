"""Pytest reporting for tests which depend on optional Psi4 addons."""

import json
from itertools import combinations_with_replacement
from pathlib import Path


_STATUS_ORDER = {"passed": 0, "skipped": 1, "failed": 2}


def _addons(item):
    """Return the addon names recorded by ``using``/``uusing`` markers."""
    # Each using()/uusing() call populates this cache while test modules are
    # imported. By collection completion it is the authoritative set of addon
    # marker names, avoiding any change to the established marker scheme.
    from addons import _using_cache

    return sorted(set(item.keywords).intersection(_using_cache))


def _cell_text(counts, *, diagonal=False):
    total = sum(counts.values())
    if not total:
        return "**〔—〕**" if diagonal else ""
    # A pass proves the addon is present and operational even if other tests in
    # the cell skip because a second addon is absent. Failures remain highest
    # priority; yellow is reserved for cells where every test skipped.
    icon = "🔴" if counts["failed"] else "🟢" if counts["passed"] else "🟡"
    details = " ".join(
        f"{counts[status]}{status[0].upper()}"
        for status in ("passed", "failed", "skipped")
        if counts[status]
    )
    text = f"{icon} {details}"
    return f"**〔{text}〕**" if diagonal else text


def _column_codes(addons):
    """Return stable one-character labels for a compact matrix."""
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    if len(addons) > len(alphabet):
        raise ValueError(f"Addon matrix supports at most {len(alphabet)} addons")
    return dict(zip(addons, alphabet))


def summarize(tests, known_addons=()):
    """Build per-addon and addon-pair result counts from test records."""
    addons = sorted(set(known_addons) | {addon for test in tests.values() for addon in test["addons"]})
    empty = lambda: {"registered": 0, "passed": 0, "failed": 0, "skipped": 0}
    totals = {addon: empty() for addon in addons}
    matrix = {a: {b: empty() for b in addons} for a in addons}

    for test in tests.values():
        status = test["status"]
        for addon in test["addons"]:
            totals[addon]["registered"] += 1
            totals[addon][status] += 1
        for a, b in combinations_with_replacement(test["addons"], 2):
            for row, col in {(a, b), (b, a)}:
                matrix[row][col]["registered"] += 1
                matrix[row][col][status] += 1
    return addons, totals, matrix


def render_markdown(tests, known_addons=()):
    addons, totals, matrix = summarize(tests, known_addons)
    lines = [
        "# Psi4 addon test report",
        "",
        "🟢 at least one passed · 🟡 all skipped · 🔴 at least one failed · **〔bold brackets〕** diagonal",
        "",
        "| Addon | Registered | Passed | Failed | Skipped |",
        "|---|---:|---:|---:|---:|",
    ]
    for addon in addons:
        row = totals[addon]
        lines.append(
            f"| {addon} | {row['registered']} | {row['passed']} | {row['failed']} | {row['skipped']} |"
        )

    lines.extend(["", "## Addon overlap matrix", ""])
    if addons:
        codes = _column_codes(addons)
        lines.append("| | " + " | ".join(codes[addon] for addon in addons) + " |")
        lines.append("|---|" + "---:|" * len(addons))
        for row_index, a in enumerate(addons):
            cells = [
                _cell_text(matrix[a][b], diagonal=col_index == row_index) if col_index <= row_index else ""
                for col_index, b in enumerate(addons)
            ]
            lines.append(f"| **{codes[a]} · {a}** | " + " | ".join(cells) + " |")
    else:
        lines.append("No addon tests were selected.")
    lines.append("")
    return "\n".join(lines)


class AddonReporter:
    def __init__(self, config, output):
        self.config = config
        self.output = Path(output)
        self.tests = {}
        self.known_addons = set()

    def pytest_runtest_logreport(self, report):
        properties = dict(report.user_properties)
        addons = properties.get("psi4_addons")
        self.known_addons.update(properties.get("psi4_known_addons", "").split(","))
        self.known_addons.discard("")
        if not addons:
            return
        current = self.tests.get(report.nodeid)
        status = "skipped" if report.skipped else "failed" if report.failed else "passed"
        # A passing setup must not mask a later skipped/failed call or teardown.
        if current is None or _STATUS_ORDER[status] > _STATUS_ORDER[current["status"]]:
            self.tests[report.nodeid] = {"addons": addons.split(","), "status": status}

    def pytest_sessionfinish(self):
        # Each invocation merges into the same report, which lets the ecosystem
        # workflow run geomeTRIC separately without fragmenting the totals.
        previous = {}
        known_addons = set(self.known_addons)
        if self.output.exists():
            try:
                old_payload = json.loads(self.output.read_text())
                previous = old_payload.get("tests", {})
                known_addons.update(old_payload.get("known_addons", ()))
            except (OSError, ValueError):
                pass
        # Non-xdist runs share the collection process with this reporter.
        try:
            from addons import _using_cache

            known_addons.update(_using_cache)
        except ImportError:
            pass
        previous.update(self.tests)
        known_addons = sorted(known_addons)
        addons, totals, matrix = summarize(previous, known_addons)
        payload = {"known_addons": known_addons, "addons": totals, "matrix": matrix, "tests": previous}
        self.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.output.with_suffix(self.output.suffix + ".tmp")
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        temporary.replace(self.output)
        self.output.with_suffix(".md").write_text(render_markdown(previous, known_addons))


def pytest_addoption(parser):
    group = parser.getgroup("psi4 addon reporting")
    group.addoption(
        "--addon-report",
        metavar="PATH",
        help="write addon test counts to PATH (JSON) and a sibling Markdown file",
    )


def pytest_configure(config):
    output = config.getoption("--addon-report")
    # xdist workers attach metadata to reports, but only the controller writes.
    if output and not hasattr(config, "workerinput"):
        config.pluginmanager.register(AddonReporter(config, output), "psi4-addon-reporter")


def pytest_collection_modifyitems(items):
    from addons import _using_cache

    known_addons = ",".join(sorted(_using_cache))
    for item in items:
        addons = _addons(item)
        if addons:
            item.user_properties.append(("psi4_addons", ",".join(addons)))
            # xdist's controller does not collect tests itself, so send the
            # complete addon universe along with reports from each worker.
            item.user_properties.append(("psi4_known_addons", known_addons))
