import xml.etree.ElementTree as ET
from sys import version as python_version

import numpy as np

from splinepy.utils.log import debug, warning


def export(fname, multipatch=None, indent=True):
    """Export as gismo readable xml geometry file
    Use gismo-specific xml-keywords to export (list of) splines. All Bezier
    patches are excluded as their respective non uniform counterpart

    Parameters
    ----------
    fname : string
      name of the output file
    spline_list : multipatch (preferred)
      (list of) Spline-Types in splinepy format
    indent: bool
      Appends white spaces using xml.etree.ElementTree.indent, if possible.

    Returns
    -------
    None
    """
    from splinepy import NURBS, BSpline, Multipatch
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

    if not isinstance(fname, str):
        raise ValueError("fname argument must be string")

    # All entries in gismo xml file have a unique ID
    xml_data = ET.Element("xml")
    index_offset = int(100)

    # First export Multipatch information
    multipatch_element = ET.SubElement(
        xml_data, "MultiPatch", parDim=str(multipatch.para_dim), id=str(0)
    )
    patch_range = ET.SubElement(multipatch_element, "patches", type="id_range")
    patch_range.text = (
        f"{index_offset} " f"{len(multipatch.splines) - 1 + index_offset}"
    )

    warning(
        "Gismo Export only works for structured meshes and no checks are "
        "performed whether this assumption is valid"
    )

    interface_data = ET.SubElement(multipatch_element, "interfaces")

    con_spline_id, con_face_id = np.where(multipatch.interfaces >= 0)
    if con_spline_id.size == 0:
        warning("No inter-face connections were found.")
    else:
        # First we identify all interfaces between patches
        is_lower_id = (
            multipatch.interfaces[con_spline_id, con_face_id] < con_spline_id
        )
        # Divide them into start points of the connection and end point of the
        # connection based on the index of the spline and the index of the
        # spline the boundary is connected with, i.e., the connection between
        # spline 3 and 7 will always point from spline 3 to 7 discarding the
        # inverse connection
        con_spline_id_start = con_spline_id[is_lower_id]
        con_face_id_start = con_face_id[is_lower_id]
        con_spline_id_end = con_spline_id[~is_lower_id]
        con_face_id_end = con_face_id[~is_lower_id]
        # Order them from start to end
        start_order = np.lexsort(
            (
                con_spline_id_start,
                multipatch.interfaces[con_spline_id_start, con_face_id_start],
            )
        ).reshape(-1, 1)
        end_order = np.lexsort(
            (
                multipatch.interfaces[con_spline_id_end, con_face_id_end],
                con_spline_id_end,
            )
        ).reshape(-1, 1)
        # HERE WE SHOULD IDENTIFY THE ORIENTATION
        # write to file
        # Reminder: Face enumeration starts at 1 in gismo (i.e. requires an
        # offset of 1)
        interface_array = np.hstack(
            (
                con_spline_id_start[start_order] + index_offset,
                con_face_id_start[start_order] + 1,
                con_spline_id_end[end_order] + index_offset,
                con_face_id_end[end_order] + 1,
                # This is the orientation:
                np.repeat(
                    np.arange(multipatch.para_dim, dtype=np.int32).reshape(
                        1, -1
                    ),
                    start_order.size,
                    axis=0,
                ),
                np.ones(
                    (start_order.size, multipatch.para_dim), dtype=np.int32
                ),
            )
        )
        interface_data.text = "\n".join(
            [" ".join([str(xx) for xx in x]) for x in interface_array]
        )

    # Start extracting all boundaries
    # Boundaries are stored as patch-id (global) face-id (local)
    boundary_spline, boundary_face = np.where(multipatch.interfaces < 0)
    if boundary_spline.size > 0:
        boundary_data = ET.SubElement(multipatch_element, "boundary")
        boundary_data.text = "\n".join(
            [
                str(sid + index_offset) + " " + str(bid + 1)
                for (sid, bid) in zip(boundary_spline, boundary_face)
            ]
        )
    ###
    # Boundary conditions
    ###
    # If there are multiple boundary IDs write them into boundary conditions
    # and fill with dummy values
    # Note here, that the boundary conditions are always stored in the local
    # enumeration system of them multipatch
    boundary_condition_list = multipatch.boundaries
    if len(boundary_condition_list) > 1:
        bcs_data = ET.SubElement(
            xml_data, "boundaryConditions", multipatch=str(0), id=str(1)
        )
        bcs_data.insert(
            0, ET.Comment(text="Please fill boundary conditions here")
        )
        for bc_data_i in boundary_condition_list:
            bc = ET.SubElement(bcs_data, "bc", type="Dirichlet", unknown="0")
            bc.text = "\n".join(
                [
                    str(sid) + " " + str(bid + 1)
                    for (sid, bid) in zip(bc_data_i[0], bc_data_i[1])
                ]
            )

    ###
    # Individual spline data
    ###
    for id, spline in enumerate(multipatch.splines):
        if not issubclass(type(spline), Spline):
            raise ValueError(
                "One of the splines handed to export is not a valid spline"
                " representation"
            )

        # Transform bezier types, as they are not supported in gismo
        if spline.name.startswith("Bezier"):
            type_name = "BSpline"
            spline = BSpline(
                **spline.todict(),
                knot_vectors=[
                    [0] * (a + 1) + [1] * (a + 1) for a in spline.degrees
                ],
            )
        elif spline.name.startswith("RationalBezier"):
            type_name = "Nurbs"
            spline = NURBS(
                **spline.todict(),
                knot_vectors=[
                    [0] * (a + 1) + [1] * (a + 1) for a in spline.degrees
                ],
            )
        elif spline.name.startswith("BSpline"):
            type_name = "BSpline"
        elif spline.name.startswith("NURBS"):
            type_name = "Nurbs"

        # Start element definition
        spline_element = ET.SubElement(
            xml_data,
            "Geometry",
            type="Tensor" + type_name + str(spline.para_dim),
            id=str(id + index_offset),
        )

        # Define Basis functions
        if "weights" in spline.required_properties:
            spline_basis_base = ET.SubElement(
                spline_element,
                "Basis",
                type="Tensor" + type_name + "Basis" + str(spline.para_dim),
            )
            spline_basis = ET.SubElement(
                spline_basis_base,
                "Basis",
                type="TensorBSplineBasis" + str(spline.para_dim),
            )
        else:
            spline_basis = ET.SubElement(
                spline_element,
                "Basis",
                type="Tensor" + type_name + "Basis" + str(spline.para_dim),
            )

        for i_para in range(spline.para_dim):
            basis_fun = ET.SubElement(
                spline_basis, "Basis", type="BSplineBasis", index=str(i_para)
            )
            knot_vector = ET.SubElement(
                basis_fun, "KnotVector", degree=str(spline.degrees[i_para])
            )
            knot_vector.text = " ".join(
                [str(k) for k in spline.knot_vectors[i_para]]
            )
        if "weights" in spline.required_properties:
            # Add weights
            weights = ET.SubElement(
                spline_basis_base,
                "weights",
            )
            weights.text = "\n".join(
                [str(w) for w in spline.weights.flatten()]
            )

        coords = ET.SubElement(
            spline_element,
            "coefs",
            geoDim=str(spline.dim),
        )
        coords.text = "\n".join(
            [" ".join([str(xx) for xx in x]) for x in spline.control_points]
        )

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


def load(fname):
    """Read gismo-keyword specified xml file

    Parameters
    ----------
    fname : str
      filename of the gismo xml

    Returns
    -------
    spline_dic_list : list<dict>
      list of dictionary types, to be put into the spline constructors
    """

    # Auxiliary function to unravel
    def retrieve_from_basis_(ETelement, SPdict):
        if "nurbs" in ETelement.attrib["type"].lower():
            # Recursive call for knot_vector
            for child in ETelement:
                if "weights" in child.tag:
                    SPdict["weights"] = np.fromstring(
                        child.text.replace("\n", " "), sep=" "
                    ).reshape(-1, dim)
                if "basis" in child.tag.lower():
                    retrieve_from_basis_(child, SPdict)
        elif "bspline" in ETelement.attrib["type"].lower():
            para_dim = int(ETelement.attrib["type"].lower().split("basis")[-1])
            knotvector = [None] * para_dim
            degrees = [None] * para_dim
            for child in ETelement:
                if "bsplinebasis" not in child.attrib["type"].lower():
                    raise ValueError("Something went wrong")
                if "knotvector" in child[0].tag.lower():
                    id = int(child.attrib["index"])
                    kv = np.fromstring(
                        child[0].text.replace("\n", " "), sep=" "
                    )
                    knotvector[id] = kv
                    degrees[id] = int(child[0].attrib["degree"])
            SPdict["knot_vectors"] = knotvector
            SPdict["degrees"] = degrees

    # Parse XML file
    debug(f"Parsing xml-file '{fname}' ...")
    root = ET.parse(fname).getroot()
    debug("XML-file parsed start conversion")

    # Init return value
    return_dict = {"BSpline": [], "NURBS": []}

    # Splines start with the keyword Geometry
    for child in root:
        if not child.tag.startswith("Geometry"):
            debug(
                f"Found unsupported keyword {child.tag}, which will be"
                " ignored"
            )
        spline_dict = {}
        debug(f"Found new spline in xml file with id {child.attrib.get('id')}")
        for info in child:
            if "basis" in info.tag.lower():
                retrieve_from_basis_(info, spline_dict)
            elif info.tag.startswith("coef"):
                dim = int(info.attrib.get("geoDim"))
                spline_dict["control_points"] = np.fromstring(
                    info.text.replace("\n", " "), sep=" "
                ).reshape(-1, dim)
        if spline_dict.get("weights") is None:
            return_dict["BSpline"].append(spline_dict)
        else:
            return_dict["NURBS"].append(spline_dict)

    debug(
        f"Found a total of {len(return_dict['BSpline'])} "
        f"BSplines and {len(return_dict['NURBS'])} NURBS"
    )

    return return_dict
