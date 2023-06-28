"""
irit io.
"""
import re

import numpy as np
from nltk import Tree

from splinepy import splinepy_core
from splinepy.io import ioutils

KEYWORD_PATTERN = (
    r"(CURVE|SURFACE|TRIVAR|MULTIVAR)\s(BSPLINE|BEZIER)([\s\dEP]+)"
)
VARTYPE_N_ENTRIES = {"CURVE": 1, "SURFACE": 2, "TRIVAR": 3}


def load(fname, expand_tabs=True):
    """
    Read spline in `.itd` form.

    Parameters
    -----------
    fname: str

      Path to the irit file to read in.
    expand_tabs: bool
      Replace all tabs in the irit file. Defaults to True.

    Returns
    --------
    splines: list
      Spline Type defined in NAME_TO_TYPE
    """
    if expand_tabs:
        ioutils.expand_tabs(fname)

    def extract_relevant_text(lines):
        # Delete lines not containing brackets
        return " ".join(
            [line.strip() for line in lines if "[" in line or "]" in line]
        )

    def get_parentheses_indices(text):
        parenthesis_tuples = [
            (m.start(0), m.group()) for m in re.finditer(r"(\[|\])", text)
        ]
        parenthesis_indices = np.array([tpl[0] for tpl in parenthesis_tuples])
        some_dict = {"[": 1, "]": -1}
        parenthesis_list = [some_dict[tpl[1]] for tpl in parenthesis_tuples]
        parenthesis_levels = np.cumsum(np.array(parenthesis_list))
        return parenthesis_indices, parenthesis_levels

    def find_first_bigger_value(indices, p_indices):
        index_repeat_matrix = np.repeat(
            indices.reshape(-1, 1), len(p_indices), axis=1
        )
        return np.argmax(index_repeat_matrix < p_indices, axis=1)

    def find_relevant_parentheses_indices(
        keyword_bigger_indices, keyword_levels, p_levels, p_indices
    ):
        relevant_parentheses_indices = []

        for keyword_index, keyword_level in zip(
            keyword_bigger_indices, keyword_levels
        ):
            spline_next_parenthesis_index = keyword_index + np.argmax(
                p_levels[keyword_index:] < keyword_level
            )
            spline_before_parenthesis_index = (
                1 + np.where(p_levels[:keyword_index] < keyword_level)[0][-1]
            )
            relevant_parentheses_indices.append(
                (
                    p_indices[spline_before_parenthesis_index],
                    p_indices[spline_next_parenthesis_index],
                )
            )

        return relevant_parentheses_indices

    def extract_splineinfo(vartype, splinetype, splineinfo):
        spline_dict = {}
        # Extract splineinfo
        relevant_entries = re.findall(r"[E|P]*\d+", splineinfo)
        entry_ints = []
        for entry in relevant_entries:
            if entry.startswith("E") or entry.startswith("P"):
                is_projected = entry.startswith("P")
            else:
                entry_ints.append(int(entry))

        entry_ints = np.array(entry_ints)
        n_entries_vartype = VARTYPE_N_ENTRIES[vartype]

        if len(entry_ints) == n_entries_vartype:
            spline_dict["degrees"] = np.array(entry_ints) - 1
            n_control_pts = entry_ints
        else:
            half = len(entry_ints) // 2
            n_control_pts = entry_ints[:half]
            spline_dict["degrees"] = np.array(entry_ints[half:]) - 1

        is_bspline = splinetype == "BSPLINE"

        return (
            spline_dict,
            is_bspline,
            n_control_pts,
            n_entries_vartype,
            is_projected,
        )

    with open(fname) as f:
        lines = f.readlines()
        relevant_text = extract_relevant_text(
            lines
        )  # One string, deleted lines which don't contain brackets
        # Get indices in text of brackets and their level in a tree
        p_indices, p_levels = get_parentheses_indices(relevant_text)
        # Find relevant keywords, entry might be: ('CURVE', 'BEZIER', '4 E3')
        keyword_matches = re.findall(KEYWORD_PATTERN, relevant_text)
        keyword_indices = np.array(
            [m.start(0) for m in re.finditer(KEYWORD_PATTERN, relevant_text)]
        )
        keyword_bigger_indices = find_first_bigger_value(
            keyword_indices, p_indices
        )
        keyword_levels = p_levels[keyword_bigger_indices - 1]

        relevant_parentheses_indices = find_relevant_parentheses_indices(
            keyword_bigger_indices, keyword_levels, p_levels, p_indices
        )

        spline_dict_list = []

        for (start, end), (vartype, splinetype, splineinfo) in zip(
            relevant_parentheses_indices, keyword_matches
        ):
            (
                spline_dict,
                is_bspline,
                n_control_pts,
                n_entries_vartype,
                is_projected,
            ) = extract_splineinfo(vartype, splinetype, splineinfo)
            tree = Tree.fromstring(
                relevant_text[start : end + 1],
                brackets="[]",
                node_pattern=r"[KV\s]*[-\d.][-\d .e]+",
            )
            control_points = []
            entry_texts = [
                subtree.label().strip()
                for subtree in tree.subtrees(
                    lambda x: type(x) != str and x.label() != ""
                )
            ]
            if (
                len(entry_texts)
                != np.prod(n_control_pts) + is_bspline * n_entries_vartype
            ):
                raise ValueError(
                    "Entries in text block does not match the required number."
                )
            knot_vectors = []
            control_points = []
            weights = []
            for textentry in entry_texts:
                if textentry.startswith("KV"):
                    relevant_entry = textentry[2:].strip()
                    knot_vector = np.fromstring(relevant_entry, sep=" ")
                    knot_vectors.append(knot_vector)
                else:
                    entries_extracted = np.fromstring(textentry, sep=" ")
                    if is_projected:
                        weight = entries_extracted[0]
                        weights.append(weight)
                        control_points.append(entries_extracted[1:] / weight)
                    else:
                        control_points.append(entries_extracted)

            spline_dict["control_points"] = np.array(control_points)
            if spline_dict["control_points"].shape[0] != np.prod(
                n_control_pts
            ):
                raise ValueError(
                    "Number of control points does not match entries from file"
                )

            if len(weights) != 0:
                spline_dict["weights"] = np.array(weights)

            if len(knot_vectors) != 0:
                spline_dict["knot_vectors"] = knot_vectors

            spline_dict_list.append(spline_dict)

    return ioutils.dict_to_spline(spline_dict_list)


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
