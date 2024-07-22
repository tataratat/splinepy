import argparse
import contextlib
from importlib import import_module

import splinepy
from splinepy._version import __version__


def convert(args):
    """Commandline interface functionality keyword convert.

    Convert splines from one format to another. See the io module for more information
    on what formats are supported for import and export.

    Args:
        args (namespace): Args namespace from argparse containing the options for
        convert.
    """
    fname = args.input_file
    from splinepy import io

    loaded = io.load(fname)
    if not isinstance(loaded, list):
        loaded = [loaded]
    if args.output_file:
        io.export(args.output_file, loaded)


def show(args):
    """Commandline interface functionality keyword show.

    Can load and show the file that is given in the args. Some options on what is shown
    can also be set.

    Args:
        args (namespace): Args namespace from argparse containing the options for show.
    """
    fname = args.input_file
    from splinepy import io

    loaded = io.load(fname)
    interactive = args.interactive or not args.output_file
    close = not bool(args.output_file)
    if interactive and not close:
        print(
            "Interactive mode is on and output file is set. "
            "Output file will be ignored. Instead, use S (Shift+S) in the "
            "interactive window to save a screenshot "
            "(saved as screenshot.png.)"
        )
    if not isinstance(loaded, list):
        loaded = [loaded]
    for loaded_object in loaded:
        loaded_object.show_options["knots"] = args.no_knot_vectors
        loaded_object.show_options["control_points"] = args.no_control_points
        loaded_object.show_options["control_point_ids"] = args.no_cp_ids
        loaded_object.show_options["control_mesh"] = args.no_control_mesh
        if args.resolution:
            loaded_object.show_options["resolutions"] = int(args.resolution)
    show_vedo = splinepy.show(loaded, interactive=interactive, close=close)
    if args.output_file:
        show_vedo.screenshot(args.output_file)


def entry():
    """Entry point for splinepy commandline interface.

    Currently only supports version and show commands.

    Try out `splinpy -h` for more information.
    """
    print()
    print(
        "                    %%\\ %%\\                                         "
    )
    print(
        "                    %% |\\__|                                        "
    )
    print(
        " %%%%%%%\\  %%%%%%\\  %% |%%\\ %%%%%%%\\   %%%%%%\\   %%%%%%\\  %%\\ "
        "  %%\\ "
    )
    print(
        "%%  _____|%%  __%%\\ %% |%% |%%  __%%\\ %%  __%%\\ %%  __%%\\ %% |  %% |"
    )
    print(
        "\\%%%%%%\\  %% /  %% |%% |%% |%% |  %% |%%%%%%%% |%% /  %% |%% |  %% |"
    )
    print(
        " \\____%%\\ %% |  %% |%% |%% |%% |  %% |%%   ____|%% |  %% |%% |  %% |"
    )
    print(
        "%%%%%%%  |%%%%%%%  |%% |%% |%% |  %% |\\%%%%%%%\\ %%%%%%%  |\\%%%%%%% |"
    )
    print(
        "\\_______/ %%  ____/ \\__|\\__|\\__|  \\__| \\_______|%%  ____/  \\____%% |"
    )
    print(
        "          %% |                                  %% |      %%\\   %% |"
    )
    print(
        "          %% |                                  %% |      \\%%%%%%  |"
    )
    print(
        "          \\__|                                  \\__|       \\______/"
    )
    print()
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
        "-i", "--input-file", type=str, help="Input file name.", required=True
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
        "--interactive",
        action="store_true",
        help="Show graphic in interactive window even if save to file is on.",
    )
    parser_plot.add_argument(
        "-c",
        "--no-control_points",
        action="store_false",
        help="Do not show control points.",
    )
    parser_plot.add_argument(
        "-k",
        "--no-knot-vectors",
        action="store_false",
        help="Do not show knot vectors.",
    )
    parser_plot.add_argument(
        "-n",
        "--no-cp-ids",
        action="store_false",
        help="Do not show control point ids.",
    )
    parser_plot.add_argument(
        "-m",
        "--no-control-mesh",
        action="store_false",
        help="Do not show control mesh.",
    )
    parser_plot.add_argument(
        "-r",
        "--resolution",
        action="store",
        help="Resolution if the spline. Number of interpolated points.",
    )
    parser_convert = subparsers.add_parser(
        "convert", help="Convert the given spline."
    )
    parser_convert.add_argument(
        "-i", "--input-file", type=str, help="Input file name.", required=True
    )
    parser_convert.add_argument(
        "-o",
        "--output-file",
        type=str,
        help="Export graphic to file.",
        required=True,
    )
    args = parser.parse_args()

    if args.version:
        dependent_versions = [
            "numpy",
            # "funi", # has currently no __version__ attribute
            "gustaf",
            "vedo",
            "scipy",
            "napf",
            "meshio",
            "k3d",
        ]
        print(f"splinepy version: {__version__}")
        for lib in dependent_versions:
            version = "Not Installed"
            with contextlib.suppress(ImportError):
                module = import_module(lib)
                version = module.__version__
            print(f"{lib} version: {version}")
        return

    if args.command == "show":
        if args.input_file is None:
            parser_plot.print_help()
            return
        show(args)
    elif args.command == "convert":
        if args.input_file is None or args.output_file is None:
            parser_convert.print_help()
            return
        convert(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    entry()
