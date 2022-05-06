"""
Export splines in custom json format
"""

import base64
import json
import logging
import numpy as np


def load(fname):
    """
    Loads splines from json file. Returns dict-splines

    Parameters
    -----------
    fname: str

    Returns
    --------
    spline_list: list
    """
    # Import data from file into dict format
    jsonbz = json.load(open(fname, "r"))

    spline_list = []
    base64encoding = jsonbz["Base64Encoding"]

    # Init return type
    return_dict = {
        "Bezier" : [],
        "NURBS" : [],
        "BSpline" : []
    }

    for jbz in jsonbz["SplineList"]:

        # Convert Base64 str into np array
        if base64encoding:
            jbz["control_points"] = np.frombuffer(
              base64.b64decode(
                jbz["control_points"].encode('ascii')
              ),
              dtype=np.float64
              ).reshape((-1, jbz["dim"]))

        # Transform data into Splines
        if jbz["SplineType"] == "Bezier":
            return_dict["Bezier"].append(
              {
                "control_points" : jbz["control_points"],
                "degrees" : jbz["degrees"]
              }
            )
        elif jbz["SplineType"] == "NURBS":
            return_dict["NURBS"].append(
              {
                "control_points" : jbz["control_points"],
                "degrees" : jbz["degrees"],
                "weights" : jbz["weights"],
                "knot_vectors" : jbz["knot_vectors"]
              }
            )
        elif jbz["SplineType"] == "BSpline":
            return_dict["BSpline"].append(
              {
                "control_points" : jbz["control_points"],
                "degrees" : jbz["degrees"],
                "knot_vectors" : jbz["knot_vectors"]
              }
            )
        else :
            logging.warning("Found unknown spline-type: " + str(jbz))

    logging.debug("Imported " + str(len(spline_list)) + " splines from file.")

    return return_dict

def export(fname, spline_list, list_name=None, base64encoding=False):
    """
    Exports a list of arbitrary splines in json-format

    Parameters
    ----------
    fname : str
      Export Filename
    spline_list: list
      List of arbitrary Spline-Types
    list_name: str
      Default is None and "SplineGroup" will be assigned. Used to define name.
    base64encoding: bool
      Default is False. If True, encodes float type spline properties before
      saving.

    Returns
    -------
    output_dict : dict
      Dictornay data written into file

    """
    # Determine Spline Group
    if list_name == None:
        list_name = "SplineGroup"

    # Ensure spline_list is actually a list
    if not isinstance(spline_list, list):
        spline_list = [spline_list]

    # Number of splines in group
    n_splines = len(spline_list)

    output_dict = {
      "Name" : list_name,
      "NumberOfSplines" : n_splines,
      "Base64Encoding" : base64encoding
    }

    # Append all splines to a dictionary
    spline_dictonary_list = []
    for i_spline in spline_list:
        if not type(i_spline).__qualname__ in ["Bezier", "BSpline", "NURBS"]:
            logging.warning(
                f"Skipping unknown spline-type `{type(i_spline).__qualname__}`"
            )
            continue

        # Create Dictionary
        i_spline_dict = i_spline.todict(tolist=True)
        i_spline_dict["SplineType"] = type(i_spline).__qualname__
        # Append para_dim and dim
        i_spline_dict["dim"] = i_spline.dim
        i_spline_dict["para_dim"] = i_spline.para_dim
        spline_dictonary_list.append(i_spline_dict)

        if base64encoding:
            i_spline_dict["control_points"] = base64.b64encode(
                np.array(i_spline_dict["control_points"])
            ).decode('utf-8')

    output_dict["SplineList"] = spline_dictonary_list

    with open(fname, 'w') as file_pointer:
        file_pointer.write(json.dumps(output_dict, indent=4))

    return output_dict
