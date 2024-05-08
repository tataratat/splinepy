"""
A script for extra doc generation.
Feel free to extend!

Contents:
1. Create markdown table of show options.
"""

import os

import splinepy as spp

if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))

    # create md dir
    md_path = os.path.abspath(os.path.join(here, "..", "md"))
    os.makedirs(md_path, exist_ok=True)

    # 1. append show options to visualization
    with open(
        os.path.abspath(os.path.join(here, "../md/spline_plotting.md")), "a"
    ) as f:
        f.write("# List of show_options\n")
        derived = [
            spp.helpme.visualize.SplineShowOption,
            spp.helpme.visualize.MultipatchShowOption,
        ]
        for cls in derived:
            f.write(f"## {cls.__qualname__}\n\n")
            for option in cls._valid_options.values():
                t_str = str(option.allowed_types)
                t_str = (
                    t_str.replace("<class", "")
                    .replace("'", "")
                    .replace(">", "")
                )
                f.write(
                    f"<details><summary><strong>{option.key}"
                    "</strong></summary><p>\n"
                )
                f.write(f"\n{option.description}  \n")
                f.write(f"- _allowed types_: {t_str}  \n")
                f.write(f"- _default_: {option.default}  \n")
                f.write("</p></details> \n\n")
