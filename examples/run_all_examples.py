"""Call all *.py files in the working directory.
"""
import glob
import subprocess
import sys

if __name__ == "__main__":
    files_not_completed = []

    for file in glob.glob("*.py"):
        if file == "run_all_examples.py":
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
