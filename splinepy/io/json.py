"""
Export splines in custom json format
"""

import json
import logging

# from splinepy import NURBS, BSpline, Bezier
# 
# 
# def load(fname):
#     # Import data from file into dict format
#     jsonbz = json.load(open(fname, "r"))
# 
#     spline_list = []
#     for jbz in jsonbz["SplineList"]:
#         if jbz["SplineType"] == "Bezier":
#             spline_list.append(
#                 splinepy.Bezier(
#                     degrees=jbz["Degrees"],
#                     control_points=jbz["ControlPoints"]
#                 )
#             )
#         elif jbz["SplineType"] == "NURBS":
#             spline_list.append(
#                 NURBS(
#                     degrees=jbz["Degrees"],
#                     control_points=jbz["ControlPoints"],
#                     knot_vector=jbz["knot_vector"],
#                     weights=jbz["weights"]
#                 )
#             )
#         elif jbz["SplineType"] == "BSpline":
#             spline_list.append(
#                 BSpline(
#                     degrees=jbz["Degrees"],
#                     control_points=jbz["ControlPoints"],
#                     knot_vector=jbz["knot_vector"]
#                 )
#             )
#         else :
#             logging.warning("Found unknown spline-type: " + str(jbz))
#     logging.debug("Imported " + str(len(spline_list)) + " splines from file.")
#     return spline_list

def export_splines(fname, spline_list, list_name=None):
  """
  Exports a list of arbitrary splines in json-format

  Parameters
  ----------
  fname : str
    Export Filename
  spline_list: list
    List of arbitrary Spline-Types

  Returns
  -------
  output_dict : dict
    Dictornay data written into file

  """
  if list_name == None:
    list_name = "SplineGroup"

  # Number of splines in group
  n_splines = len(spline_list)

  output_dict = {
    "Name" : list_name,
    "NumberOfSplines" : n_splines
  }

  # Append all splines to a dictionary
  spline_dictonary_list = []
  for i_spline in spline_list:
      # Create Dictionary
      i_spline_dict = i_spline.todict()
      # Append para_dim and dim
      i_spline_dict["dim"] = i_spline.dim
      i_spline_dict["para_dim"] = i_spline.para_dim
      i_spline_dict["SplineType"] = type(i_spline).__qualname__
      spline_dictonary_list.append(i_spline_dict)

  output_dict["SplineList"] = spline_dictonary_list

  with open(fname, 'w') as file_pointer:
      file_pointer.write(json.dumps(output_dict, indent=4))
  return output_dict


def write_json(spline, fname):
  """
  Pass Information to export_splines routin

  Parameters
  ----------
  spline : spline
    arbitrary spline-type
  fname : str
    filename

  Returns
  -------
  None
  """
  export_splines(fname, [spline])

