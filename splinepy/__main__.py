import argparse

from splinepy._version import __version__


def convert(input_file_name, input_type, output_file_name, output_type):
    pass


def convert_cli(args):
    pass


def entry():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="sub-command help", dest="command")
    parser.add_argument(
        "-v",
        "--version",
        action="store_true",
        help="Show version of splinepy.",
    )

    parser_plot = subparsers.add_parser("plot", help="Plot the given spline.")
    parser_plot.add_argument(
        "-i", "--input-file", type=str, help="Input file name."
    )
    parser_plot.add_argument(
        "-o",
        "--output-file",
        type=str,
        help="Input file name.",
        required=False,
    )

    parser_convert = subparsers.add_parser(
        "convert", help="Convert the given spline."
    )
    parser_convert.add_argument(
        "-i", "--input-file", type=str, help="Input file name."
    )
    parser_convert.add_argument(
        "-o", "--output-file", type=str, help="Input file name."
    )
    parser_convert.add_argument(
        "-it",
        "--input-type",
        type=str,
        help="Input file type.",
        required=False,
    )
    parser_convert.add_argument(
        "-ot",
        "--output-type",
        type=str,
        help="Output file type.",
        required=False,
    )

    args = parser.parse_args()

    if args.version:
        print(f"splinepy version {__version__}")
        return

    if args.command == "plot":
        print("Plotting:", args)
    elif args.command == "convert":
        print("Converting:", args)


if __name__ == "__main__":
    entry()
