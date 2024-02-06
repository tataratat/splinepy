import argparse
import contextlib
from importlib import import_module

import splinepy
from splinepy._version import __version__


def show(args):
    fname = args.input_file
    from splinepy import io

    loaded = io.load(fname)
    interactive = args.interactive
    close = not bool(args.output_file)
    if interactive and not close:
        print(
            "Interactive mode is on and output file is set. "
            "Output file will be ignored. Instead, use S (Shift+S) in the "
            "interactive window to save a screenshot "
            "(saved as screenshot.png.)"
        )
    show_vedo = splinepy.show(loaded, interactive=interactive, close=close)
    if args.output_file:
        show_vedo.screenshot(args.output_file)


def entry():
    print("")
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
    print("")
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
        "--interactive",
        action="store_true",
        help="Show graphic in interactive window even if save to file is on.",
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
        show(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    entry()
