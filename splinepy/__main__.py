import argparse

import splinepy
from splinepy._version import __version__


def show(args):
    fname = args.input_file
    print(fname)
    from splinepy import io

    loaded = io.load(fname)

    show_vedo = splinepy.show(loaded, interactive=not args.non_interactive)
    if args.output_file:
        show_vedo.screenshot(args.output_file)


def entry():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="sub-command help", dest="command")
    parser.add_argument(
        "-v",
        "--version",
        action="store_true",
        help="Show version of splinepy.",
    )

    parser_plot = subparsers.add_parser("show", help="Show the given spline.")
    parser_plot.add_argument(
        "-i", "--input-file", type=str, help="Input file name."
    )
    parser_plot.add_argument(
        "-o",
        "--output-file",
        type=str,
        help="Export graphic to file.",
        required=False,
    )
    parser_plot.add_argument(
        "-e",
        "--non-interactive",
        action="store_false",
        help="Show graphic in interactive window.",
    )

    args = parser.parse_args()

    if args.version:
        print(f"splinepy version {__version__}")
        return

    if args.command == "show":
        show(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    entry()
