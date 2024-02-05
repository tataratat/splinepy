"""Handles the processing of markdown files.

This script is used to process the markdown files in the project that are also
reused in the documentation. The script will replace the relative links in the
markdown files with relative links that are correct when the documentation is
built.

The paths change since the documentation is built from the docs folder and not
from the root of the project.

Author: Clemens Fricke
"""

import os
import pathlib
import re

# Path to this file.
file_path = os.path.abspath(os.path.dirname(__file__))
original_cwd = os.getcwd()
repo_root = str(pathlib.Path(__file__).resolve()).split("docs")[0]
os.chdir(repo_root)


def get_markdown_links(line: str) -> str:
    """Get the markdown links from a string.

    Args:
        line (str): text.

    Returns:
        str: Markdown links.
    """
    possible = re.findall(r"\[(.*?)\]\((.*?)\)", line)
    return possible if possible else ""


def get_github_path_from(link):
    """Substitute the path to the github repository.

    This will expand the link to the github repository. This is used to create
    pages that are independent of the documentation. Like for the long
    description on PyPI. Normally we try to use relative links in the files
    to guarantee that the documentation is also available offline. But for
    the long description on PyPI we need to use absolute links to an online
    repository like the github raw files.

    Args:
        link (_type_): _description_

    Returns:
        _type_: _description_
    """
    return str(pathlib.Path(link).resolve()).replace(
        repo_root, "https://raw.githubusercontent.com/tataratat/splinepy/main/"
    )


def get_common_parent(path1: str, path2: str) -> str:
    """Get the common parent of two paths.

    Args:
        path1 (str): Path 1.
        path2 (str): Path 2.

    Returns:
        str: Common parent path.
    """
    path1 = pathlib.Path(path1).resolve().parent
    path2 = pathlib.Path(path2).resolve()
    steps_back = 0
    while path1.is_relative_to(path2) is False:
        path2 = path2.parent
        steps_back += 1
    return path2, steps_back


# Folder to save the processed markdown files to.
folder_to_save_to = os.path.join(repo_root, "docs/md/")

# List of markdown files that are used in the documentation.
markdown_files = [
    pathlib.Path("README.md").resolve(),
    pathlib.Path("CONTRIBUTING.md").resolve(),
    pathlib.Path("docs/markdown/spline_intro.md").resolve(),
    pathlib.Path("docs/markdown/spline_plotting.md").resolve(),
]

# List of explicit link substitutions.
link_substitutions = {
    "docs/markdown/spline_plotting.md": "spline_intro.html#creating-and-plotting-splines"  # noqa E501
}


def process_file(
    file: str, relative_links: bool = True, return_content: bool = False
):
    """Process a markdown file.

    This function will process a markdown file. It will replace all relative
    links with links that are correct when the documentation is built.
    If relative_links is True, the links will be relative to the documentation
    folder. If relative_links is False, the links will be absolute links to
    the github repository.

    Args:
        file (str): Path to the markdown file.
        relative_links (bool, optional):
            Generate relative links. Defaults to True.
        return_content (bool, optional):
            Return the content instead of saving it. Defaults to False.

    Returns:
        str: Content of the markdown file. Only if return_content is True.
    """
    # read in the content of the markdown file
    with open(file) as f:
        content = f.read()
    os.chdir(os.path.dirname(file)) if os.path.dirname(file) else None
    # get all links from the markdown file
    links = get_markdown_links(content)

    for item in links:
        if item[1].startswith(("http", "#")):  # skip http links and anchors
            if "badge" in item[1]:
                continue
            content = content.replace(
                f"[{item[0]}]({item[1]})",
                f"<a href='{item[1]}'>{item[0]}</a>",
            )
            continue
        elif item[1] in link_substitutions:
            if relative_links:
                content = content.replace(
                    f"[{item[0]}]({item[1]})",
                    f"<a href='{link_substitutions[item[1]]}'>{item[0]}</a>",
                )
                continue
            else:
                content = content.replace(
                    f"[{item[0]}]({item[1]})",
                    "See documentation for examples.",
                )
        elif not relative_links:  # generate links to github repo
            new_path = get_github_path_from(pathlib.Path(item[1]).resolve())
        else:  # generate relative links
            common_sub_path, steps_back = get_common_parent(
                item[1], folder_to_save_to
            )
            new_path = "../" * steps_back + str(
                pathlib.Path(item[1]).resolve().relative_to(common_sub_path)
            )
        content = content.replace(item[1], str(new_path))

    os.chdir(original_cwd)

    if return_content:
        return content

    with open(
        os.path.join(folder_to_save_to, os.path.basename(file)), "w"
    ) as f:
        f.write(content)


if __name__ == "__main__":
    os.chdir(repo_root)
    os.makedirs(folder_to_save_to, exist_ok=True)
    # Process all markdown files
    for file in markdown_files:
        process_file(file)
