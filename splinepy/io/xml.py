"""
xml io.
Keyword follows RWTH CATS's spline format.
Possiblely will be turned into pure python io.
"""
import xml.etree.ElementTree as ET
from sys import version as python_version

from splinepy import splinepy_core
from splinepy.io import ioutils
from splinepy.multipatch import Multipatch
from splinepy.utils.log import debug, warning

# List of spline keywords (bundled here in case they change - many unused)
CATS_XML_KEY_WORDS = {
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
    return ioutils.dict_to_spline(splinepy_core.read_xml(fname))


def export(fname, multipatch, indent=True):
    """
    Save spline as `.xml`.

    Parameters
    -----------
    fname: str
      Export file name
    multipatch: Multipatch
      Splines to be exported
    indent: bool
      Option for pretty printing

    Returns
    --------
    None
    """
    from splinepy.spline import Spline

    # First transform spline-data into a multipatch-data if required
    if issubclass(type(multipatch), Spline):
        # Transform to list
        multipatch = [multipatch]

    if isinstance(multipatch, list):
        # Transform to multipatch
        multipatch = Multipatch(splines=multipatch)
        warning(
            "No interfaces (i.e. inter-face connectivity) information is"
            " given. Try converting list to Multipatch object. Interfaces"
            "will be computed on the fly as a result."
        )

    if not isinstance(multipatch, Multipatch):
        raise ValueError(
            "export Function expects list for multipatch argument"
        )

    # All entries in gismo xml file have a unique ID
    xml_data = ET.Element("xml")

    # Start a list of a list of patches
    spline_list_element = ET.SubElement(
        xml_data,
        CATS_XML_KEY_WORDS["spline_list"],
        SplineType=str(1),
        **{CATS_XML_KEY_WORDS["n_patches"]: str(len(multipatch.splines))},
    )

    # @todo: current implementation only exports geometry, export fields
    if multipatch.fields:
        warning("Fields are not supported yet")

    # All Splines (patches) are written into the spline list as entries
    for spline in multipatch.splines:
        # Convert to non-bezier type (might make unnecessary copy)
        if spline.is_rational:
            patch = spline.nurbs
        else:
            patch = spline.bspline

        # Write spline header

        patch_element = ET.SubElement(
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
        control_point_var_names = ET.SubElement(
            patch_element,
            CATS_XML_KEY_WORDS["coeff_names"],
        )
        # These don't really matter so let's just assign numbers
        control_point_var_names.text = " ".join(
            "x" + str(i) for i in range(patch.dim)
        )

        # Control point variables
        control_points_elements = ET.SubElement(
            patch_element,
            CATS_XML_KEY_WORDS["control_points"],
        )
        control_points_elements.text = "\n".join(" ".join(
            str(xx) for xx in x
        ) for x in patch.cps )

        # degrees
        degrees_elements = ET.SubElement(
            patch_element,
            CATS_XML_KEY_WORDS["degrees"],
        )
        degrees_elements.text = " ".join(
            str(deg) for deg in patch.degrees
        )

        # knot-vectors
        knot_vectors_elements = ET.SubElement(
            patch_element,
            CATS_XML_KEY_WORDS["knot_vectors"],
        )
        for kv in patch.knot_vectors:
            knot_vector_element = ET.SubElement(
                knot_vectors_elements,
                CATS_XML_KEY_WORDS["knot_vector"],
            )
            knot_vector_element.text = " ".join(str(k) for k in kv)
        
        # weights if rational
        if patch.is_rational:
            weights_elements = ET.SubElement(
                patch_element,
                CATS_XML_KEY_WORDS["weights"],
            )
            weights_elements.text = " ".join(
                str(w) for w in patch.weights.ravel()
            )
            

    # Beautify and export
    if int(python_version.split(".")[1]) >= 9 and indent:
        # Pretty printing xml with indent only exists in version > 3.9
        ET.indent(xml_data)

    elif int(python_version.split(".")[1]) < 9 and indent:
        debug(
            "Indented xml ouput is only supported from > python3.9.",
            "Output will not be indented.",
            f"Current python version: {python_version}",
        )

    file_content = ET.tostring(xml_data)
    with open(fname, "wb") as f:
        f.write(file_content)
