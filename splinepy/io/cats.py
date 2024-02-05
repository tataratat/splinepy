"""
xml io.
Keyword follows RWTH CATS's spline format.
Possibly will be turned into pure python io.
"""

import xml.etree.ElementTree as _ET
from sys import version as _python_version

import numpy as _np

from splinepy.io import ioutils as _ioutils
from splinepy.utils.log import debug as _debug
from splinepy.utils.log import warning as _warning

# List of spline keywords (bundled here in case they change - many unused)
CATS_XML_KEY_WORDS = {
    "spline_type": "SplineType",  # useless keyword (always one)
    "patch": "SplineEntry",  # key for new patch in List
    "n_patches": "NumberOfSplines",  # number of patches
    "para_dim": "splDim",  # Parametric dimension
    "dim": "spaceDim",  # Spatial dimension
    "field_dim": "numOfCntrlPntVars",  # Number of control point variables
    "n_ctps": "numCntrlPnts",  # Number of control points
    "n_element_vars": "numOfEleVars",  # Number of element variables
    "periodic": "closed",  # Is the spline closed? No[0] or yes[1] for each
    # parameter direction
    "coeff_names": "cntrlPntVarNames",  # Names of coefficients
    "control_points": "cntrlPntVars",  # Actual control points
    "conforming_ctps": "connectedCntrlPnts",  # Interfaces and conforming
    # control points
    "weights": "wght",  # Weights (every spline is a nurbs)
    "degrees": "deg",  # Degrees vector
    "knot_vectors": "kntVecs",  # Vector of knot vectors
    "knot_vector": "knotVec",  # individual knot vector
    "spline_list": "SplineList",  # Global list of patches
    "element_variables": "elmentVars",  # Variables defined on the element
}


def load(fname):
    """
    Read spline in `.xml` form.
    RWTH CATS spline format.

    Parameters
    -----------
    fname: str

    Returns
    --------
    splines: Multipatch object
      Spline Type defined in NAME_TO_TYPE
    """

    def _read_spline(xml_element):
        spline_dict = {}
        if not xml_element.tag.startswith(CATS_XML_KEY_WORDS["patch"]):
            _debug(
                f"Found unexpected keyword {xml_element.tag}, which will be "
                f"ignored"
            )
        # Read required information from Header in xml element

        dim = int(xml_element.attrib.get(CATS_XML_KEY_WORDS["dim"], -1))
        para_dim = int(
            xml_element.attrib.get(CATS_XML_KEY_WORDS["para_dim"], -1)
        )
        ctps_dim = int(
            xml_element.attrib.get(CATS_XML_KEY_WORDS["field_dim"], -1)
        )
        n_ctps = int(xml_element.attrib.get(CATS_XML_KEY_WORDS["n_ctps"], -1))
        is_periodic = (
            int(xml_element.attrib.get(CATS_XML_KEY_WORDS["periodic"], 0)) == 1
        )
        if -1 in [dim, para_dim, ctps_dim, n_ctps]:
            raise ValueError(
                f"Not enough information provided for xml-element "
                f"{xml_element}"
            )
        if ctps_dim > dim:
            _debug(
                "Spline seems to contain field information, but field "
                "retrieval is not implemented. All non-geometry info will be "
                "ignored"
            )
        if is_periodic:
            raise ValueError(
                "Periodic splines are not (yet) supported. Try saving your "
                "spline by repeating first and last control points."
            )
        for info in xml_element:
            # We ignore most keywords at the moment
            if CATS_XML_KEY_WORDS["coeff_names"] in info.tag:
                _debug(
                    f"Will use the following control point vars as control "
                    f"point coordinates : {info.text.split()[:dim]}"
                )
            # Control points
            elif CATS_XML_KEY_WORDS["control_points"] in info.tag:
                spline_dict["control_points"] = _np.fromstring(
                    info.text.replace("\n", " "), sep=" "
                ).reshape(-1, ctps_dim)[:, :dim]

            # degrees
            elif CATS_XML_KEY_WORDS["degrees"] in info.tag:
                spline_dict["degrees"] = _np.fromstring(
                    info.text.replace("\n", " "), sep=" ", dtype=_np.int64
                )
                if spline_dict["degrees"].size != para_dim:
                    raise ValueError(
                        f"Conflicting information provided. Got "
                        f"{spline_dict['degrees'].size} degrees, but "
                        f"parametric dimension is {para_dim}"
                    )

            # weights
            elif CATS_XML_KEY_WORDS["weights"] in info.tag:
                spline_dict["weights"] = _np.fromstring(
                    info.text.replace("\n", " "), sep=" "
                )

            # knot vectors
            elif CATS_XML_KEY_WORDS["knot_vectors"] in info.tag:
                spline_dict["knot_vectors"] = []
                for child_info in info:
                    if CATS_XML_KEY_WORDS["knot_vector"] not in child_info.tag:
                        _debug("Redundant item in knot_vectors block of xml")
                    spline_dict["knot_vectors"].append(
                        _np.fromstring(
                            child_info.text.replace("\n", " "), sep=" "
                        )
                    )

            # All other keywords will be ignored for the moment
            else:
                _debug(f"Ignoring keyword {info.tag}")

        return spline_dict

    # Parse  XML file
    _debug(f"Parsing xml-file '{fname}' ...")
    root = _ET.parse(fname).getroot()
    _debug("XML-file parsed start conversion")
    list_of_splines = []
    if root.tag.startswith(CATS_XML_KEY_WORDS["spline_list"]):
        if root.attrib.get(CATS_XML_KEY_WORDS["spline_type"], "1") != "1":
            raise ValueError(
                f"Unknown SplineType "
                f"{root.attrib[CATS_XML_KEY_WORDS['spline_type']]}"
            )
        n_patches = int(root.attrib.get(CATS_XML_KEY_WORDS["n_patches"], 0))
        for patch_element in root:
            list_of_splines.append(_read_spline(patch_element))
    else:
        _debug(f"Unused xml-keyword {root.tag}")
        return []

    _debug(f"Found a total of {len(list_of_splines)} " f"BSplines and NURBS")
    if len(list_of_splines) != n_patches:
        raise ValueError(
            f"Found {len(list_of_splines)} splines, but expected to find "
            f"{n_patches} patches"
        )

    # Return splines
    return _ioutils.dict_to_spline(list_of_splines)


def export(fname, spline_list, indent=True):
    """
    Save spline as `.xml`.

    Parameters
    -----------
    fname: str
      Export file name
    multipatch: Multipatch / list
      Splines to be exported
    indent: bool
      Option for pretty printing

    Returns
    --------
    None
    """
    from splinepy.multipatch import Multipatch as _Multipatch
    from splinepy.spline import Spline as _Spline

    # First transform spline-data into a multipatch-data if required
    if issubclass(type(spline_list), _Spline):
        # Transform to list
        spline_list = [spline_list]

    if isinstance(spline_list, _Multipatch):
        # @todo: current implementation only exports geometry, export fields
        if spline_list.fields:
            _warning("Fields are not supported yet")
        # Transform to multipatch
        spline_list = spline_list.patches

    if not isinstance(spline_list, list):
        raise ValueError(
            "export Function expects list for multipatch argument"
        )

    # Start a list of a list of patches
    spline_list_element = _ET.Element(
        CATS_XML_KEY_WORDS["spline_list"],
        SplineType=str(1),
        **{CATS_XML_KEY_WORDS["n_patches"]: str(len(spline_list))},
    )

    new_line_char = "\n" if indent else " "

    # All Splines (patches) are written into the spline list as entries
    for spline in spline_list:
        # Convert to non-bezier type (might make unnecessary copy)
        patch = spline.nurbs if spline.is_rational else spline.bspline

        # Write spline header
        patch_element = _ET.SubElement(
            spline_list_element,
            CATS_XML_KEY_WORDS["patch"],
            **{
                CATS_XML_KEY_WORDS["para_dim"]: str(patch.para_dim),
                CATS_XML_KEY_WORDS["dim"]: str(patch.dim),
                CATS_XML_KEY_WORDS["field_dim"]: str(patch.dim),
                CATS_XML_KEY_WORDS["n_ctps"]: str(patch.cps.shape[0]),
                CATS_XML_KEY_WORDS["n_element_vars"]: str(0),
                CATS_XML_KEY_WORDS["periodic"]: str(0),
            },
        )

        # Variable names (@todo fields)
        control_point_var_names = _ET.SubElement(
            patch_element,
            CATS_XML_KEY_WORDS["coeff_names"],
        )

        # These are usually hard coded variable names in xns or feafa
        # x y z t seems to only matter. see issue #257
        cp_varnames = [
            "x",
            "y",
            "z",
            "t",
            *[f"x{i}" for i in range(4, patch.dim)],
        ]
        control_point_var_names.text = " ".join(cp_varnames[: patch.dim])

        # Control point variables
        control_points_elements = _ET.SubElement(
            patch_element,
            CATS_XML_KEY_WORDS["control_points"],
        )
        control_points_elements.text = new_line_char.join(
            " ".join(str(xx) for xx in x) for x in patch.cps
        )

        # degrees
        degrees_elements = _ET.SubElement(
            patch_element,
            CATS_XML_KEY_WORDS["degrees"],
        )
        degrees_elements.text = " ".join(str(deg) for deg in patch.degrees)

        # knot-vectors
        knot_vectors_elements = _ET.SubElement(
            patch_element,
            CATS_XML_KEY_WORDS["knot_vectors"],
        )
        for kv in patch.knot_vectors:
            knot_vector_element = _ET.SubElement(
                knot_vectors_elements,
                CATS_XML_KEY_WORDS["knot_vector"],
            )
            knot_vector_element.text = " ".join(str(k) for k in kv)

        # weights if rational
        if patch.is_rational:
            weights_elements = _ET.SubElement(
                patch_element,
                CATS_XML_KEY_WORDS["weights"],
            )
            weights_elements.text = " ".join(
                str(w) for w in patch.weights.ravel()
            )

    # Beautify and export
    if int(_python_version.split(".")[1]) >= 9 and indent:
        # Pretty printing xml with indent only exists in version > 3.9
        _ET.indent(spline_list_element)

    elif int(_python_version.split(".")[1]) < 9 and indent:
        _debug(
            "Indented xml output is only supported from > python3.9.",
            "Output will not be indented.",
            f"Current python version: {_python_version}",
        )

    file_content = _ET.tostring(spline_list_element)
    with open(fname, "wb") as f:
        f.write(file_content)
