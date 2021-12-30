"""
io utils.
"""

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
    Formulate a string, taking each *args as a line.

    Parameters
    -----------
    *args: *str

    Returns
    --------
    line_separated_str: str
    """
    line_separated_str = ""
    for a in args:
        line_separated_str += a + "\n"

    return line_separated_str


