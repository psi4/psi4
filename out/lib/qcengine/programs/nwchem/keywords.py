import collections
from typing import Any, Dict, Tuple


def format_keyword(keyword: str, val: Any, lop_off: bool = True, preserve_case: bool = False) -> Tuple[str, str]:
    """Function to reformat value `val` for `keyword` from python into nwchem-speak."""

    if preserve_case:
        key = keyword
    else:
        key = keyword.lower()

    # Transform string booleans into " "
    if val is True:
        return key, ""
    elif val is False:
        return "", ""

    # complete hack
    if keyword.upper() == "MEMORY":
        return keyword.lower(), f"{val} double"

    elif isinstance(val, list):
        text = " ".join([str(v) for v in val])
    elif isinstance(val, dict):
        text = []
        for k, v in val.items():
            merge = [k]
            merge.extend(str(v) if isinstance(v, (int, float)) else list(map(str, v)))
            text.append(" ".join(merge))
        text = " ".join(text)
    else:
        text = str(val)

    if lop_off:
        return key[7:], text
    else:
        return key, text


def format_keywords(keywords: Dict[str, Any]) -> str:
    """From NWCHEM-directed, non-default `keywords` dictionary, write a NWCHEM deck."""

    def rec_dd():
        return collections.defaultdict(rec_dd)

    grouped_options = rec_dd()

    for group_key, val in keywords.items():
        nesting = group_key.split("__")
        if len(nesting) == 1:
            key = nesting[0]
            grouped_options["aaaglobal"][key] = val
        elif len(nesting) == 2:
            g1, key = nesting
            grouped_options[g1][key] = val
        elif len(nesting) == 3:
            g1, g2, key = nesting
            grouped_options[g1][g2][key] = val
        else:
            raise ValueError("Nesting N!")

    grouped_lines = {}
    for group, opts in sorted(grouped_options.items()):
        lines = []
        group_level_lines = []
        for key, val in grouped_options[group].items():
            if isinstance(val, dict):
                g2_level_lines = []
                g2_level_lines.append(key.lower())
                for k2, v2 in val.items():
                    line2 = " ".join(format_keyword(k2, v2, lop_off=False))
                    g2_level_lines.append(line2)
                g2_level_lines = " ".join(g2_level_lines)
                lines.append(g2_level_lines)
            else:
                preserve_case = False
                try:
                    preserve_case = val.startswith("library ")
                except AttributeError:
                    pass
                line = " ".join(format_keyword(key, val, lop_off=False, preserve_case=preserve_case))
                if group.lower() == "basis" and any(
                    [word in line for word in ["spherical", "cartesian", "print", "noprint", "rel"]]
                ):
                    group_level_lines.append(line)
                else:
                    lines.append(line)
        if group == "aaaglobal":
            grouped_lines[group] = "\n".join(lines) + "\n"
        elif group.lower() == "set":
            grouped_lines[group] = "\n".join(f"set {l}" for l in lines) + "\n"
        else:
            grouped_lines[group] = (
                f"{group.lower()} " + " ".join(group_level_lines) + "\n  " + "\n  ".join(lines) + "\nend\n"
            )

    return "\n".join(grouped_lines.values()) + "\n"
