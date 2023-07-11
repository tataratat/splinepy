"""
irit io.
"""
import re

import numpy as np

from splinepy import splinepy_core
from splinepy.io import ioutils

# Keywords that indicate new spline (irit is case insensitive, all lower case)
SPLINE_KEY_WORD_REGEX = r"\s+(bspline|bezier)\s+"


def load(fname, expand_tabs=True, strip_comments=False):
    """
    Read spline in `.itd` form.

    Parameters
    -----------
    fname: str
      Path to the irit file to read in.
    expand_tabs: bool
      Replace all tabs in the irit file. Defaults to True.
    strip_comments: bool
      (BETA) tries to identify comments and strips them from file

    Returns
    --------
    splines: list
      Spline Type defined in NAME_TO_TYPE
    """
    # Expand taps to spaces
    if expand_tabs:
        ioutils.expand_tabs(fname)

    # Delete lines not containing brackets in the beginning
    def extract_relevant_text(lines):
        if strip_comments:
            return " ".join(
                [
                    line.lower().strip()
                    for line in lines
                    if re.match(r".*[\[\]].*", line) is not None
                ]
            )
        else:
            return " ".join([line.lower().strip() for line in lines])

    with open(fname) as f:
        # Strip commentaries
        lines = f.readlines()
        relevant_text = extract_relevant_text(lines)

        # At this point we ignore all non-spline related data and do a keyword
        # based search for relevant information
        spline_strings = re.split(SPLINE_KEY_WORD_REGEX, relevant_text)

        # Init list
        spline_list = []

        # Loop over text file data
        for type, data in zip(spline_strings[1::2], spline_strings[2::2]):
            # Dictionary of current spline
            spline = {}

            # Extract data from data
            data = re.split(r"\]\s*\]", data)[0] + "]"

            # Split dimensions and degrees from data
            dim_and_deg = re.split(r"(e\d|p\d)", data)

            # must contain ctps/degs rational/poly dimension control points
            if len(dim_and_deg) != 3:
                raise ValueError(
                    "Unsupported irit format, could not identify spline from "
                    "following string data set:" + type + data
                )

            # Dim is stored in second string written (E/P)+DIM
            dim = int(dim_and_deg[1][1:])
            is_rational = dim_and_deg[1][0] == "p"

            # Extract degrees (and control point dimensions)
            deg_info = [
                int(s) for s in re.split(r"\s+", dim_and_deg[0]) if len(s) > 0
            ]

            # IRIT uses orders instead of degrees
            if type == "bezier":
                orders = np.array(deg_info)
                n_ctps = np.prod(orders)
            else:
                n_ctps = np.prod(deg_info[: (len(deg_info) // 2)])
                orders = np.array(deg_info[(len(deg_info) // 2) :])
            spline["degrees"] = orders - 1

            # Paradim is specified from degrees
            para_dim = orders.size

            # Extract control point and kv information
            ctps_and_kvs = re.findall(r"\[(.*?)\]", dim_and_deg[2])
            ctps = " ".join(
                [
                    s
                    for s in ctps_and_kvs
                    if re.match(r"\s*\d.*", s) is not None
                ]
            )

            if is_rational:
                ctps = np.fromstring(str(ctps), dtype=float, sep=" ").reshape(
                    -1, dim + 1
                )
                spline["weights"] = ctps[:, 0]
                spline["control_points"] = ctps[:, 1:] / ctps[:, 0:1]
            else:
                spline["control_points"] = np.fromstring(
                    str(ctps), dtype=float, sep=" "
                ).reshape(-1, dim)

            if spline["control_points"].shape[0] != n_ctps:
                raise ValueError(
                    "Could not find sufficient number of control points in "
                    "spline:" + type + data
                )

            # Extract knot-vectors
            if type == "bspline":
                # Extract knot-vectors
                knot_v_strings = [
                    s.split("kv")[1]
                    for s in ctps_and_kvs
                    if re.search(r"\s*kv.*", s) is not None
                ]
                if len(knot_v_strings) != para_dim:
                    raise ValueError(
                        "Could not find enough knot vectors in bspline string:"
                        + type
                        + data
                    )
                spline["knot_vectors"] = [
                    np.fromstring(kv, dtype=float, sep=" ")
                    for kv in knot_v_strings
                ]

            # Append dictionary to list
            spline_list.append(spline)

        return ioutils.dict_to_spline(spline_list)


def export(fname, splines):
    """
    Save splines as `.itd`.

    Parameters
    -----------
    fname: str
    splines: list

    Returns
    --------
    None
    """
    return splinepy_core.export_irit(fname, splines)
