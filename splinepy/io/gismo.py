import copy
import xml.etree.ElementTree as ET
from sys import version as python_version

import numpy as np

from splinepy import settings
from splinepy.utils.log import debug, warning


def _spline_to_ET(root, multipatch, index_offset, fields_only=False):
    """
    Write spline data to xml element in gismo format

    Parameters
    ----------
    root : ElementTree.Subelement
      branch in element tree to which the spline info is to be added
    multipatch : Multipatch
      Multipatch containing the information requested
    index_offset : int
      index_offset for ids in xml export
      (All geometries are assigned individual ids that must be unique in the
      xml file)
    fields_only : bool (False)
      If set, exports only the fields associated to the multipatch data

    Returns
    -------
    None
    """
    from splinepy import NURBS, BSpline
    from splinepy.spline import Spline

    if fields_only and (len(multipatch.fields) == 0):
        return

    if fields_only:
        supports = np.array(
            [
                (i, j, 0)
                for j, field in enumerate(multipatch.fields)
                for i, v in enumerate(field)
                if v is not None
            ],
            dtype=np.int64,
        )

        # Very unintuitive solution to counting the support ids #indextrick
        indices = np.argsort(supports[:, 0], kind="stable")
        counter = np.arange(supports.shape[0])
        bincount = np.bincount(supports[:, 0])
        # Get number of occurences for previous element to correct index shift
        index_shift = np.cumsum(np.hstack([[0], bincount[:-1]]))
        counter -= np.repeat(index_shift, bincount)
        supports[indices, 2] = counter

        # Write Matrix
        design_v_support = ET.SubElement(
            root,
            "Matrix",
            rows=str(supports.shape[0]),
            cols=str(supports.shape[1]),
            id=str(10),
        )
        design_v_support.text = "\n".join(
            [" ".join([str(xx) for xx in x]) for x in supports]
        )

    for id, spline in enumerate(multipatch.splines):
        if fields_only:
            # Check supports
            support = supports[supports[:, 0] == id, 1]
            coefs = np.hstack(
                [multipatch.fields[j][id].control_points for j in support]
            )
            if "weights" in spline.required_properties:
                weights = np.hstack(
                    [multipatch.fields[j][id].weights for j in support]
                )
        else:
            coefs = spline.control_points
            if "weights" in spline.required_properties:
                weights = spline.weights

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
            root,
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
            weights_element = ET.SubElement(
                spline_basis_base,
                "weights",
            )
            weights_element.text = "\n".join(
                [" ".join([str(ww) for ww in w]) for w in weights]
            )
        coords = ET.SubElement(
            spline_element,
            "coefs",
            geoDim=str(coefs.shape[1]),
        )
        coords.text = "\n".join(
            [" ".join([str(xx) for xx in x]) for x in coefs]
        )


def export(
    fname,
    multipatch=None,
    indent=True,
    labeled_boundaries=True,
    options=None,
    export_fields=False,
):
    """Export as gismo readable xml geometry file
    Use gismo-specific xml-keywords to export (list of) splines. All Bezier
    patches are excluded as their respective non uniform counterpart

    Parameters
    ----------
    fname : string
      name of the output file
    spline_list : Multipatch (preferred)
      (list of) Spline-Types in splinepy format
    indent: bool
      Appends white spaces using xml.etree.ElementTree.indent, if possible.
    labeled_boundaries : bool
      Writes boundaries with labels into the MultiPatch part of the XML
    options : list
      List of dictionaries, that specify model related options, like boundary
      conditions, assembler options, etc.. The dictionaries must have the
      following keys, 'tag'->string, 'text'->string (optional),
      'attributes'->dictionary (optional), 'children'->list in the same format
      (optional)
    export_fields : bool
      Export fields to seperate files ending with field<id>.xml, e.g.,
      filename.xml.field1.xml
    collapse_fields : cool
      Only valid if export_fields is True, writes all fields in the same file
      by collapsing the control points, requires conformity (not checked)

    Returns
    -------
    None
    """
    from splinepy import Multipatch
    from splinepy.settings import NTHREADS, TOLERANCE
    from splinepy.spline import Spline
    from splinepy.splinepy_core import orientations

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

    interface_data = ET.SubElement(multipatch_element, "interfaces")

    con_spline_id, con_face_id = np.where(multipatch.interfaces >= 0)
    if con_spline_id.size == 0:
        warning("No inter-face connections were found.")
    else:
        # First we identify all interfaces between patches
        is_lower_id = (
            multipatch.interfaces[con_spline_id, con_face_id] > con_spline_id
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

        # Reorder arrays
        con_spline_id_start = con_spline_id_start[start_order]
        con_face_id_start = con_face_id_start[start_order]
        con_spline_id_end = con_spline_id_end[end_order]
        con_face_id_end = con_face_id_end[end_order]

        # Identify Orientation
        (
            axis_mapping,
            axis_orientation,
        ) = orientations(
            multipatch.splines,
            con_spline_id_start,
            con_face_id_start,
            con_spline_id_end,
            con_face_id_end,
            TOLERANCE,
            NTHREADS,
        )

        # write to file
        # Reminder: Face enumeration starts at 1 in gismo (i.e. requires an
        # offset of 1)
        interface_array = np.hstack(
            (
                con_spline_id_start + index_offset,
                con_face_id_start + 1,
                con_spline_id_end + index_offset,
                con_face_id_end + 1,
                # This is the orientation:
                axis_mapping,
                # Convert bool into -1 and 1 for gismo
                axis_orientation.astype(np.int32),
            )
        )
        interface_data.text = "\n".join(
            [" ".join([str(xx) for xx in x]) for x in interface_array]
        )

    if labeled_boundaries:
        # Export geometric boundaries and set a label
        for boundary_id, face_id_list in enumerate(
            multipatch.boundaries, start=1
        ):
            boundary_data = ET.SubElement(
                multipatch_element,
                "boundary",
                name="BID" + str(boundary_id),
            )
            boundary_data.text = "\n".join(
                [
                    str(patch_id + index_offset) + " " + str(local_face_id + 1)
                    for (patch_id, local_face_id) in zip(*face_id_list)
                ]
            )
    else:
        # This is valid in the old format (1. list all boundaries, 2. set BCs)
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
        # If there are multiple boundary IDs write them into boundary
        # conditions and fill with dummy values
        # Note here, that the boundary conditions are always stored in the
        # local enumeration system of them multipatch
        boundary_condition_list = multipatch.boundaries
        if len(boundary_condition_list) > 1:
            bcs_data = ET.SubElement(
                xml_data,
                "boundaryConditions",
                multipatch=str(len(multipatch.splines)),
                id=str(1),
            )
            bcs_data.insert(
                0, ET.Comment(text="Please fill boundary conditions here")
            )
            for bc_data_i in boundary_condition_list:
                bc = ET.SubElement(
                    bcs_data, "bc", type="Dirichlet", unknown="0"
                )
                bc.text = "\n".join(
                    [
                        str(sid) + " " + str(bid + 1)
                        for (sid, bid) in zip(bc_data_i[0], bc_data_i[1])
                    ]
                )

    ###
    # Individual spline data
    ###
    # Export fields first, as all necessary information is already available
    if export_fields:
        field_xml = copy.deepcopy(xml_data)
        _spline_to_ET(field_xml, multipatch, index_offset, fields_only=True)
        if int(python_version.split(".")[1]) >= 9 and indent:
            ET.indent(field_xml)
        file_content = ET.tostring(field_xml)
        with open(fname + ".fields.xml", "wb") as f:
            f.write(file_content)

    _spline_to_ET(xml_data, multipatch, index_offset)

    # Add addtional options to the xml file
    if options is not None:
        # Verify that the list stored in the correct format
        if not isinstance(options, list):
            options = [options]

        def _apply_options(ETelement, options_list):
            for gismo_dictionary in options_list:
                name = gismo_dictionary.get("tag", None)
                if name is None and not isinstance(name, str):
                    raise ValueError(
                        "Gismo option in unsupported format, tag must be set, "
                        "please check out export documentation"
                    )
                attributes = gismo_dictionary.get("attributes", dict())
                option_text = gismo_dictionary.get("text", None)
                optional_data = ET.SubElement(
                    ETelement,
                    name,
                    **attributes,
                )
                if option_text is not None:
                    optional_data.text = option_text
                if gismo_dictionary.get("children", None) is not None:
                    _apply_options(
                        optional_data,
                        options_list=gismo_dictionary["children"],
                    )

        _apply_options(xml_data, options)

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


def load(fname, load_options=True):
    """Read gismo-keyword specified xml file

    Parameters
    ----------
    fname : str
      filename of the gismo xml
    load_options : bool
      Retrieve additional options (else - multipatch/geometry only)

    Returns
    -------
    spline_dic_list : Multipatch
      Multipatch object with list of splines in NAME_TO_TYPE-type
    options : list
      List of additional options, like boundary conditions and functions
    """
    from splinepy.multipatch import Multipatch

    # Auxiliary function to unravel
    def retrieve_from_basis_(ETelement, SPdict):
        if "nurbs" in ETelement.attrib["type"].lower():
            # Recursive call for knot_vector
            for child in ETelement:
                if "weights" in child.tag:
                    SPdict["weights"] = np.fromstring(
                        child.text.replace("\n", " "), sep=" "
                    )
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

    # Auxiliary function to store additional information in dictionary
    def make_dictionary(ETelement):
        new_dictionary = {}
        new_dictionary["tag"] = ETelement.tag
        new_dictionary["attributes"] = ETelement.attrib
        new_dictionary["text"] = ETelement.text
        new_dictionary["children"] = []
        for element in ETelement.findall("./*"):
            new_dictionary["children"].append(make_dictionary(element))
        return new_dictionary

    # Parse XML file
    debug(f"Parsing xml-file '{fname}' ...")
    root = ET.parse(fname).getroot()
    debug("XML-file parsed start conversion")

    # Init return value
    list_of_splines = []
    list_of_options = []
    interface_array = None
    invalid_integer = -99919412
    # Splines start with the keyword Geometry
    for child in root:
        if child.tag.startswith("MultiPatch"):
            parametric_dimension = int(child.attrib["parDim"])
            # Patch ids
            patch_element = child.find("patches")
            if patch_element is None:
                debug("Unsupported format")
            if not patch_element.attrib.get("type") == "id_range":
                debug(f"Invalid patch type {patch_element.attrib.get('type')}")
            patch_range = np.fromstring(
                patch_element.text.replace("\n", " "), sep=" ", dtype=np.int64
            )
            offset = patch_range[0]
            n_splines = patch_range[1] - patch_range[0] + 1
            number_of_faces = parametric_dimension * 2

            # Extract interfaces and boundaries
            interfaces_element = child.find("interfaces")
            interface_array = (
                np.ones((n_splines, number_of_faces), dtype=np.int64)
                * invalid_integer
            )
            if interfaces_element is None or interfaces_element.text is None:
                # This can occur if only one patch is in a Multipatch
                debug("No connectivity found")
            else:
                interface_information = np.fromstring(
                    interfaces_element.text.replace("\n", " "),
                    sep=" ",
                    dtype=np.int64,
                ).reshape(-1, number_of_faces + 4)
                # Assign interfaces
                interface_array[
                    interface_information[:, 0] - offset,
                    interface_information[:, 1] - 1,
                ] = (
                    interface_information[:, 2] - offset
                )
                interface_array[
                    interface_information[:, 2] - offset,
                    interface_information[:, 3] - 1,
                ] = (
                    interface_information[:, 0] - offset
                )

            # Extract interfaces and boundaries
            boundary_elements = child.findall("boundary")
            if boundary_elements is None:
                debug("No boundary information found")
            else:
                for id, boundary_element in enumerate(
                    boundary_elements, start=1
                ):
                    if boundary_element.text is None:
                        continue
                    boundary_information = np.fromstring(
                        boundary_element.text.replace("\n", " "),
                        sep=" ",
                        dtype=np.int64,
                    ).reshape(-1, 2)
                    interface_array[
                        boundary_information[:, 0] - offset,
                        boundary_information[:, 1] - 1,
                    ] = -id

        elif child.tag.startswith("Geometry"):
            spline_dict = {}
            debug(f"Found new spline in xml with id {child.attrib.get('id')}")
            for info in child:
                if "basis" in info.tag.lower():
                    retrieve_from_basis_(info, spline_dict)
                elif info.tag.startswith("coef"):
                    dim = int(info.attrib.get("geoDim"))
                    spline_dict["control_points"] = np.fromstring(
                        info.text.replace("\n", " "), sep=" "
                    ).reshape(-1, dim)
            if spline_dict.get("weights") is None:
                list_of_splines.append(
                    settings.NAME_TO_TYPE["BSpline"](**spline_dict)
                )
            else:
                list_of_splines.append(
                    settings.NAME_TO_TYPE["NURBS"](**spline_dict)
                )
        else:
            if load_options:
                list_of_options.append(make_dictionary(child))
            else:
                debug(
                    f"Found unsupported keyword {child.tag}, which will be"
                    " ignored"
                )
                continue

    debug(f"Found a total of {len(list_of_splines)} " f"BSplines and NURBS")
    multipatch = Multipatch(list_of_splines)
    if interface_array is not None:
        if invalid_integer in interface_array:
            warning("Unmatched faces or insufficient information in xml file")
            invalid_splines = np.where(interface_array)[0][0] + offset
            warning(f"Check faces with ids {invalid_splines}")
            warning("Interface information ignored")

        multipatch.interfaces = interface_array

    if load_options:
        return multipatch, list_of_options
    else:
        return multipatch
