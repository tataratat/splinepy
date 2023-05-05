"""
io utils.
"""

import os


def make_meaningful(line, comment="#"):
    """
    Trys to turn a str input a meaningful line. Returns False if it is not
    possible

    Parameters
    -----------
    line: str
    comment: str
      str that defined comment section. "#" is default

    Returns
    --------

    """
    line = line.strip()
    if len(line) == 0 or line.startswith(comment):
        return False
    else:
        return line


def next_line(f, comment="#"):
    """
    Strips and returns next meaningful line from the opened, readable file.
    It tries 513 times before it says it is the end. If it is end, it
    will return None.
    This is to prevent never ending while loop.

    Parameters
    -----------
    f: _io.TextIoWrapper
    comment: str
      A symbol that defines comment section. "#" is default.

    Returns
    --------
    next_meaningful_str: str
    """
    line = f.readline().strip()
    counter = 0
    while line.startswith(comment) or len(line) == 0:
        counter += 1
        if counter >= 513:
            return None
        line = f.readline().strip()

    return line


def form_lines(*args):
    """
    Formulate a string, taking each args as a line.

    Parameters
    -----------
    *args: str

    Returns
    --------
    line_separated_str: str
    """
    line_separated_str = ""
    for a in args:
        line_separated_str += a + "\n"

    return line_separated_str


def abs_fname(fname):
    """
    Checks if fname is absolute. If not, returns abspath. Tilde safe.

    Parameters
    ----------
    fname: str

    Returns
    --------
    abs_fname: str
      Maybe same to fname, maybe not.
    """
    if os.path.isabs(fname):
        pass

    elif "~" in fname:
        fname = os.path.expanduser(fname)

    else:
        fname = os.path.abspath(fname)

    return fname


def expand_tabs(fname, overwrite=True, tab_expand=2):
    """Replaces tabs in a text file with spaces

    Parameters
    ----------
    fname : string
      Filename
    overwrite : bool
      Inplace modification, otherwise new file with prefix copy is created
    tab_expand : int
      Determines how many spaces to be set for each tab

    Returns
    -------
    None
    """
    import os.path as op

    # Safe guard in case an absolut path is handed to the function
    if overwrite:
        out_name = fname
    else:
        dir, file = op.split(fname)
        out_name = dir + "copy_" + file

    with open(fname) as inputFile:
        file_contents = inputFile.read()
    file_contents.replace("\t", " " * tab_expand)
    with open(out_name, "w") as exportFile:
        exportFile.write(file_contents)


def dict_to_spline(spline_dictionary):
    """Create a list of splines from a list of dictionaries

    Parameters
    ----------
    spline_dictionary : list
      List of dictionaries

    Returns
    -------
    spline_lits : list
      List of splines in called format (NAME_TO_TYPE)
    """
    from splinepy.settings import NAME_TO_TYPE
    from splinepy.spline import Spline

    spline_list = [Spline(**spd) for spd in spline_dictionary]
    return [NAME_TO_TYPE[spl.name](spline=spl) for spl in spline_list]
