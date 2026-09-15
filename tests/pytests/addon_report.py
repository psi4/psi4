"""Pytest reporting for tests which depend on optional Psi4 addons."""

import argparse
import json
from itertools import combinations_with_replacement
from pathlib import Path


_STATUS_ORDER = {"passed": 0, "skipped": 1, "failed": 2}


def _report_addons(addons):
    """Normalize addon labels for reporting; Psi4 itself is implicit."""
    return sorted(set(addons) - {"psi4"})


def _addons(item):
    """Return the addon names recorded by ``using``/``uusing`` markers."""
    # Each using()/uusing() call populates this cache while test modules are
    # imported. By collection completion it is the authoritative set of addon
    # marker names, avoiding any change to the established marker scheme.
    from addons import _using_cache

    return _report_addons(set(item.keywords).intersection(_using_cache))


def _summary_status_icon(counts):
    """Status for the summary table."""
    if counts["failed"]:
        return "🔴"
    if counts["registered"] and counts["passed"] == counts["registered"]:
        return "🔵"
    if counts["passed"]:
        return "🟢"
    if counts["skipped"]:
        return "🟡"
    return "⚪"


def _count_text(counts):
    details = "+".join(
        f"{counts[status]}{status[0].upper()}"
        for status in ("passed", "failed", "skipped")
        if counts[status]
    )
    return f"{_summary_status_icon(counts)} {details or '0'}"


def summarize(tests, known_addons=()):
    """Build per-addon and addon-pair result counts from test records."""
    addons = _report_addons(set(known_addons) | {addon for test in tests.values() for addon in test["addons"]})
    empty = lambda: {"registered": 0, "passed": 0, "failed": 0, "skipped": 0}
    totals = {addon: empty() for addon in ["core", *addons, "psi4"]}
    matrix = {a: {b: empty() for b in addons} for a in addons}

    for test in tests.values():
        status = test["status"]
        test_addons = _report_addons(test["addons"])
        totals["psi4"]["registered"] += 1
        totals["psi4"][status] += 1
        if not test_addons:
            totals["core"]["registered"] += 1
            totals["core"][status] += 1
        for addon in test_addons:
            totals[addon]["registered"] += 1
            totals[addon][status] += 1
        for a, b in combinations_with_replacement(test_addons, 2):
            for row, col in {(a, b), (b, a)}:
                matrix[row][col]["registered"] += 1
                matrix[row][col][status] += 1
    return addons, totals, matrix


def render_markdown(tests, known_addons=()):
    addons, totals, matrix = summarize(tests, known_addons)
    lines = [
        "# Psi4 addon test report",
        "",
        "🔵 all registered passed · 🟢 at least one passed · 🟡 all skipped · 🔴 at least one failed · ⚪ no executed tests",
        "",
        "| Addon | Count | Also tested with |",
        "|---|:---|---|",
    ]
    for addon in ["core", *addons, "psi4"]:
        row = totals[addon]
        partners = ""
        if addon in matrix:
            partners = ", ".join(other for other in addons if other != addon and matrix[addon][other]["registered"])
        lines.append(f"| *{addon}* | {_count_text(row)} | {partners} |")

    lines.append("")
    return "\n".join(lines)


def render_combined_markdown(reports):
    """Render addon counts from multiple CI lanes as one table."""
    reports = sorted(reports, key=lambda report: report.get("label", ""))
    addons = sorted(
        {addon for report in reports for addon in report.get("addons", {}) if addon not in {"core", "psi4"}}
    )
    rows = ["core", *addons, "psi4"]
    lines = [
        "# Combined Psi4 addon test report",
        "",
        "🔵 all registered passed · 🟢 at least one passed · 🟡 all skipped · 🔴 at least one failed · ⚪ no executed tests · — report unavailable",
        "",
    ]
    if not reports:
        lines.extend(["No addon reports were available.", ""])
        return "\n".join(lines)

    partner_source = next(
        (report for report in reports if "linux" in (report.get("label") or "").lower() and report.get("tests")),
        next((report for report in reports if report.get("tests")), {}),
    )
    partners = {addon: set() for addon in addons}
    for test in partner_source.get("tests", {}).values():
        test_addons = _report_addons(test.get("addons", ()))
        for addon in test_addons:
            partners.setdefault(addon, set()).update(other for other in test_addons if other != addon)

    labels = [report.get("label") or f"Lane {index}" for index, report in enumerate(reports, 1)]
    lines.append("| Addon | Also tested with | " + " | ".join(labels) + " |")
    lines.append("|---|---|" + ":---|" * len(reports))
    for addon in rows:
        cells = []
        for report in reports:
            counts = report.get("addons", {}).get(addon)
            cells.append(_count_text(counts) if counts is not None else "—")
        also_tested_with = ", ".join(sorted(partners.get(addon, ())))
        lines.append(f"| *{addon}* | {also_tested_with} | " + " | ".join(cells) + " |")
    lines.append("")
    return "\n".join(lines)


def combine_reports(input_dir, output):
    reports = []
    for path in sorted(Path(input_dir).rglob("addon-report.json")):
        try:
            reports.append(json.loads(path.read_text(encoding="utf-8")))
        except (OSError, ValueError):
            continue
    Path(output).write_text(render_combined_markdown(reports), encoding="utf-8")


class AddonReporter:
    def __init__(self, config, output, label=None):
        self.config = config
        self.output = Path(output)
        self.label = label
        self.tests = {}
        self.known_addons = set()

    def pytest_runtest_logreport(self, report):
        properties = dict(report.user_properties)
        if "psi4_addons" not in properties:
            return
        addons = properties["psi4_addons"]
        self.known_addons.update(properties.get("psi4_known_addons", "").split(","))
        self.known_addons.discard("")
        current = self.tests.get(report.nodeid)
        status = "skipped" if report.skipped else "failed" if report.failed else "passed"
        # A passing setup must not mask a later skipped/failed call or teardown.
        if current is None or _STATUS_ORDER[status] > _STATUS_ORDER[current["status"]]:
            self.tests[report.nodeid] = {"addons": addons.split(",") if addons else [], "status": status}

    def pytest_sessionfinish(self):
        # Each invocation merges into the same report, which lets the ecosystem
        # workflow run geomeTRIC separately without fragmenting the totals.
        previous = {}
        previous_label = None
        known_addons = set(self.known_addons)
        if self.output.exists():
            try:
                old_payload = json.loads(self.output.read_text(encoding="utf-8"))
                previous = old_payload.get("tests", {})
                previous_label = old_payload.get("label")
                known_addons.update(old_payload.get("known_addons", ()))
            except (OSError, ValueError):
                pass
        # Non-xdist runs share the collection process with this reporter.
        try:
            from addons import _using_cache

            known_addons.update(_report_addons(_using_cache))
        except ImportError:
            pass
        previous.update(self.tests)
        known_addons = sorted(known_addons)
        _, totals, _ = summarize(previous, known_addons)
        payload = {
            "label": self.label or previous_label,
            "known_addons": known_addons,
            "addons": totals,
            "tests": previous,
        }
        self.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.output.with_suffix(self.output.suffix + ".tmp")
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        temporary.replace(self.output)
        self.output.with_suffix(".md").write_text(render_markdown(previous, known_addons), encoding="utf-8")


def pytest_addoption(parser):
    group = parser.getgroup("psi4 addon reporting")
    group.addoption(
        "--addon-report",
        metavar="PATH",
        help="write addon test counts to PATH (JSON) and a sibling Markdown file",
    )
    group.addoption(
        "--addon-report-label",
        metavar="LABEL",
        help="label this report when combining results from multiple CI lanes",
    )


def pytest_configure(config):
    output = config.getoption("--addon-report")
    # xdist workers attach metadata to reports, but only the controller writes.
    if output and not hasattr(config, "workerinput"):
        config.pluginmanager.register(
            AddonReporter(config, output, config.getoption("--addon-report-label")), "psi4-addon-reporter"
        )


def pytest_collection_modifyitems(items):
    from addons import _using_cache

    known_addons = ",".join(_report_addons(_using_cache))
    for item in items:
        addons = _addons(item)
        item.user_properties.append(("psi4_addons", ",".join(addons)))
        # xdist's controller does not collect tests itself, so send the
        # complete addon universe along with reports from each worker.
        item.user_properties.append(("psi4_known_addons", known_addons))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine Psi4 addon reports from CI lanes")
    parser.add_argument("--combine", metavar="DIRECTORY", required=True)
    parser.add_argument("--output", metavar="PATH", required=True)
    args = parser.parse_args()
    combine_reports(args.combine, args.output)
