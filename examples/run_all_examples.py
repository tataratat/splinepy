"""Call all *.py files in the working directory.
"""

import glob
import importlib
import itertools
import os
import pathlib
import sys

import gustaf

if __name__ == "__main__":
    gustaf.settings.DEFAULT_OFFSCREEN = True
    files_not_completed = []
    examples_root_dir = pathlib.Path(__file__).parent.absolute()
    current_dir = pathlib.Path(os.getcwd()).absolute()
    # find all example files
    for file_ in itertools.chain(
        glob.glob(str(examples_root_dir / "*.py")),
        glob.glob(str(examples_root_dir / "iga/*.py")),
    ):
        # ignore this file
        if "run_all_examples.py" in file_:
            continue
        # i like pathlib so i convert the path into it
        file = pathlib.Path(file_)
        print(
            f"file: {file}",
            f"file.parent: {file.parent}",
            f"current_dir: {current_dir}",
            f"file.stem: {file.stem}",
        )
        # change directory to the file's directory so that imports are
        # relative to the file's directory
        os.chdir(file.parent)
        print(f"Calling {file}")
        # import and call the file
        try:
            # add the file's directory to the path so importlib can
            # import it
            sys.path.insert(0, f"{file.parent}{os.sep}")
            example = importlib.import_module(file.stem)
        except Exception as e:
            print(f"Failed to import {file}.")
            print(f"Error: {e}")
            files_not_completed.append(file)
            continue
        try:
            # run example
            example.example()
        except Exception as e:
            print(f"Failed to call {file}.")
            print(f"Error: {e}")
            files_not_completed.append(file)
            continue
        # change directory back to the original directory
        os.chdir(current_dir)

    # print out all files that have failed and give system exit code
    if len(files_not_completed) > 0:
        print(
            f"Failed to call {len(files_not_completed)} files: "
            f"{files_not_completed}."
        )
        sys.exit(1)
    print("All files completed successfully.")
    sys.exit(0)
