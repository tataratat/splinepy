"""
Export splines in custom json format
"""

import base64
import json

import numpy as np

from splinepy import settings
from splinepy.utils.log import debug


def load(fname):
    """
    Loads splines from json file. Returns a list of splines

    Parameters
    -----------
    fname: str

    Returns
    --------
    spline_list: list
    """
    # Make Required Properties locally accessible
    from splinepy.spline import RequiredProperties

    # Import data from file into dict format
    with open(fname) as f:
        jsonbz = json.load(f)

    spline_list = []
    base64encoding = jsonbz["Base64Encoding"]

    for jbz in jsonbz["SplineList"]:
        # Convert Base64 str into np array
        if base64encoding:
            jbz["control_points"] = np.frombuffer(
                base64.b64decode(jbz["control_points"].encode("ascii")),
                dtype=np.float64,
            ).reshape((-1, jbz["dim"]))
            if jbz.get("weights") is not None:
                jbz["weights"] = np.frombuffer(
                    base64.b64decode(jbz["weights"].encode("ascii")),
                    dtype=np.float64,
                ).reshape((-1, 1))

        # Declare Type and get required properties
        req_props = RequiredProperties.of(jbz["SplineType"])
        data_dict = {}
        for prop in req_props:
            data_dict[prop] = jbz[prop]
        spline_list.append(
            settings.NAME_TO_TYPE[jbz["SplineType"]](**data_dict)
        )

    debug("Imported " + str(len(spline_list)) + " splines from file.")

    return spline_list


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
    if list_name is None:
        list_name = "SplineGroup"

    # Ensure spline_list is actually a list
    if not isinstance(spline_list, list):
        spline_list = [spline_list]

    # Number of splines in group
    n_splines = len(spline_list)

    output_dict = {
        "Name": list_name,
        "NumberOfSplines": n_splines,
        "Base64Encoding": base64encoding,
    }

    # Append all splines to a dictionary
    from splinepy.spline import Spline

    spline_dictonary_list = []
    for i_spline in spline_list:
        if not issubclass(type(i_spline), Spline):
            raise ValueError("Unsupported type in list")

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
            ).decode("utf-8")
            if "weights" in i_spline.required_properties:
                i_spline_dict["weights"] = base64.b64encode(
                    np.array(i_spline_dict["weights"])
                ).decode("utf-8")

    output_dict["SplineList"] = spline_dictonary_list

    with open(fname, "w") as file_pointer:
        file_pointer.write(json.dumps(output_dict, indent=4))

    return output_dict
