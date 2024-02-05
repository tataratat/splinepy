"""Call all *.py files in the working directory.
"""

import glob
import itertools
import pathlib
import subprocess
import sys

if __name__ == "__main__":
    files_not_completed = []
    examples_root_dir = pathlib.Path(__file__).parent.absolute()
    for file in itertools.chain(
        glob.glob(str(examples_root_dir / "*.py")),
        glob.glob(str(examples_root_dir / "iga/*.py")),
    ):
        if "run_all_examples.py" in file:
            continue
        print(f"Calling {file}")
        proc_return = subprocess.run([sys.executable, file], check=False)
        if proc_return.returncode != 0:
            files_not_completed.append(file)
    if len(files_not_completed) > 0:
        print(
            f"Failed to call {len(files_not_completed)} files: "
            f"{files_not_completed}."
        )
        sys.exit(1)

    print("All files completed successfully.")
    sys.exit(0)
