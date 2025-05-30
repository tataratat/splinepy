import base64 as _base64
import copy as _copy
import xml.etree.ElementTree as _ET
from sys import version as _python_version

import numpy as _np

from splinepy import settings as _settings
from splinepy.utils.data import enforce_contiguous as _enforce_contiguous
from splinepy.utils.log import debug as _debug
from splinepy.utils.log import warning as _warning


def _spline_to_ET(
    root,
    multipatch,
    index_offset,
    fields_only=False,
    as_base64=False,
    field_mask=None,
):
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
    as_base64 : bool (False)
      If set, coordinates, weights and knot-vectors are exported in B64
    field_mask : array-like
      Selection of active fields to be exported into the fields file

    Returns
    -------
    None
    """
    from splinepy.spline import Spline as _Spline

    if fields_only and (len(multipatch.fields) == 0):
        return

    if fields_only and field_mask is not None:
        if not isinstance(field_mask, (list, _np.ndarray)):
            raise ValueError("field_mask must be list of integers")
        field_mask = _np.unique(
            _enforce_contiguous(field_mask, dtype=_np.int64, asarray=True)
        )

        # Check if range is valid
        if (field_mask.min() < 0) or (
            field_mask.max() >= len(multipatch.fields)
        ):
            raise ValueError(
                "field_mask contains unsupported range, must be within "
                f"(0, {len(multipatch.fields)})"
            )

    def _array_to_text(array, is_matrix, as_b64):
        if as_b64:
            array = _enforce_contiguous(array, dtype=_np.float64, asarray=True)
            return _base64.b64encode(array).decode("ascii")
        elif is_matrix:
            return "\n".join(
                [" ".join([str(x) for x in row]) for row in array]
            )
        else:
            return " ".join(str(x) for x in array)

    if fields_only:
        if field_mask is None:
            field_mask = _np.arange(len(multipatch.fields))
        supports = _np.array(
            [
                (i, j, 0)
                for j, k in enumerate(field_mask)
                for i, v in enumerate(multipatch.fields[k].patches)
                if v is not None
            ],
            dtype=_np.int64,
        )

        # Sanity check (see error message)
        n_patches_in_support = _np.max(supports[:, 0]) + 1
        if _np.unique(supports[:, 0]).size != n_patches_in_support:
            raise ValueError(
                "The requested field(s) can not be exported, because the patch "
                "ids of the field support (i.e., where the fields are non-zero)"
                " is not contiguous. This can be the case, if a field is only "
                "non-zero on patches [0,1,4]. This case is currently not "
                "supported. If required, it is possible to (1) change the input"
                " routine and export a patch-mapping (former patch id to new "
                "patch id or (2) add zero-valued splines to fill the gaps. "
                "Please feel free to write an issue on github.com/taratat/"
                "splinepy for additional information."
            )

        # Very unintuitive solution to counting the support ids #indextrick
        indices = _np.argsort(supports[:, 0], kind="stable")
        counter = _np.arange(supports.shape[0])
        bincount = _np.bincount(supports[:, 0])
        # Get number of occurrences for previous element to correct index shift
        index_shift = _np.cumsum(_np.hstack([[0], bincount[:-1]]))
        counter -= _np.repeat(index_shift, bincount)
        supports[indices, 2] = counter

        # Write Matrix
        design_v_support = _ET.SubElement(
            root,
            "Matrix",
            cols=str(supports.shape[1]),
            id=str(10),
            rows=str(supports.shape[0]),
        )
        design_v_support.text = "\n".join(
            [" ".join([str(xx) for xx in x]) for x in supports]
        )
    else:
        n_patches_in_support = len(multipatch.patches)

    for id, spline in enumerate(multipatch.patches):
        if fields_only:
            # Check supports
            support = field_mask[supports[supports[:, 0] == id, 1]]
            if support.size == 0:
                continue
            coefs = _np.hstack(
                [
                    multipatch.fields[j].patches[id].control_points
                    for j in support
                ]
            )
            if "weights" in spline.required_properties:
                weights = _np.hstack(
                    [multipatch.fields[j][id].weights for j in support]
                )
        else:
            coefs = spline.control_points
            if "weights" in spline.required_properties:
                weights = spline.weights

        if not issubclass(type(spline), _Spline):
            raise ValueError(
                "One of the splines handed to export is not a valid spline"
                " representation"
            )

        # Transform bezier types, as they are not supported in gismo
        if spline.name.startswith("Bezier"):
            type_name = "BSpline"
            spline = spline.bspline  # noqa PLW2901
        elif spline.name.startswith("RationalBezier"):
            type_name = "Nurbs"
            spline = spline.nurbs  # noqa PLW2901
        elif spline.name.startswith("BSpline"):
            type_name = "BSpline"
        elif spline.name.startswith("NURBS"):
            type_name = "Nurbs"

        # Start element definition
        spline_element = _ET.SubElement(
            root,
            "Geometry",
            id=str(id + index_offset),
            type="Tensor" + type_name + str(spline.para_dim),
        )

        # Define Basis functions
        if "weights" in spline.required_properties:
            spline_basis_base = _ET.SubElement(
                spline_element,
                "Basis",
                type="Tensor" + type_name + "Basis" + str(spline.para_dim),
            )
            spline_basis = _ET.SubElement(
                spline_basis_base,
                "Basis",
                type="TensorBSplineBasis" + str(spline.para_dim),
            )
        else:
            spline_basis = _ET.SubElement(
                spline_element,
                "Basis",
                type="Tensor" + type_name + "Basis" + str(spline.para_dim),
            )

        # Set the format flag for output type
        format_flag = "B64Float64" if as_base64 else "ASCII"

        for i_para in range(spline.para_dim):
            basis_fun = _ET.SubElement(
                spline_basis,
                "Basis",
                index=str(i_para),
                type="BSplineBasis",
            )
            knot_vector = _ET.SubElement(
                basis_fun,
                "KnotVector",
                degree=str(spline.degrees[i_para]),
                format=format_flag,
            )
            knot_vector.text = _array_to_text(
                spline.knot_vectors[i_para], False, as_base64
            )

        if "weights" in spline.required_properties:
            # Add weights
            weights_element = _ET.SubElement(
                spline_basis_base,
                "weights",
                format=format_flag,
            )
            weights_element.text = _array_to_text(weights, True, as_base64)
        coords = _ET.SubElement(
            spline_element,
            "coefs",
            format=format_flag,
            geoDim=str(coefs.shape[1]),
        )
        coords.text = _array_to_text(coefs, True, as_base64)

    return n_patches_in_support


def export(
    fname,
    multipatch=None,
    indent=True,
    labeled_boundaries=True,
    options=None,
    export_fields=False,
    field_mask=None,
    as_base64=False,
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
      Export fields to separate file ending with fields.xml, e.g.,
      filename.xml.fields.xml. Only non-zero splines are exported to save
      memory
    field_mask : list
      Selection of active fields that are exported (to save memory and speed up
      export)
    as_base64 : bool (False)
      If set, coordinates, weights and knot-vectors are exported in B64
      encoding. Available in gismo, once gismo/gismo#634 merges.

    Returns
    -------
    None
    """
    from splinepy import Multipatch as _Multipatch
    from splinepy.spline import Spline as _Spline

    # First transform spline-data into a multipatch-data if required
    if issubclass(type(multipatch), _Spline):
        # Transform to list
        multipatch = [multipatch]

    if isinstance(multipatch, list):
        # Transform to multipatch
        multipatch = _Multipatch(splines=multipatch)
        _warning(
            "No interfaces (i.e. inter-face connectivity) information is"
            " given. Try converting list to Multipatch object. Interfaces"
            "will be computed on the fly as a result."
        )

    if not isinstance(multipatch, _Multipatch):
        raise ValueError(
            "export Function expects list for multipatch argument"
        )

    if not isinstance(fname, str):
        raise ValueError("fname argument must be string")

    # All entries in gismo xml file have a unique ID
    xml_data = _ET.Element("xml")
    index_offset = 100

    # licked it so its ours
    xml_data.insert(
        1,
        _ET.Comment(
            "generated by splinepy https://github.com/tataratat/splinepy"
        ),
    )

    # First export Multipatch information
    multipatch_element = _ET.SubElement(
        xml_data, "MultiPatch", id=str(0), parDim=str(multipatch.para_dim)
    )
    patch_range = _ET.SubElement(
        multipatch_element, "patches", type="id_range"
    )

    # Retrieve all interfaces (negative numbers refer to boundaries)
    interface_data = _ET.SubElement(multipatch_element, "interfaces")
    global_interface_id = _np.where(multipatch.interfaces.ravel() >= 0)[0]

    if global_interface_id.size == 0:
        _warning("No inter-face connections were found.")
    else:
        # Reminder: Face enumeration starts at 1 in gismo (i.e. requires an
        # offset of 1) Copy to not overwrite original interface orientations
        interface_array = multipatch.interface_orientations().copy()

        # Apply offsets
        interface_array[:, [0, 2]] += index_offset
        interface_array[:, [1, 3]] += 1  # Gismo starts counting at 1

        interface_data.text = "\n".join(
            [" ".join([str(xx) for xx in x]) for x in interface_array]
        )

    if labeled_boundaries:
        # Export geometric boundaries and set a label
        for boundary_id, face_id_list in enumerate(
            multipatch.boundaries, start=1
        ):
            boundary_data = _ET.SubElement(
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
        boundary_spline, boundary_face = _np.where(multipatch.interfaces < 0)
        if boundary_spline.size > 0:
            boundary_data = _ET.SubElement(multipatch_element, "boundary")
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
            bcs_data = _ET.SubElement(
                xml_data,
                "boundaryConditions",
                id=str(1),
                multipatch=str(len(multipatch.patches)),
            )
            bcs_data.insert(
                0, _ET.Comment(text="Please fill boundary conditions here")
            )
            for bc_data_i in boundary_condition_list:
                bc = _ET.SubElement(
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
        field_xml = _copy.deepcopy(xml_data)

        n_patches = _spline_to_ET(
            field_xml,
            multipatch,
            index_offset,
            fields_only=True,
            as_base64=as_base64,
            field_mask=field_mask,
        )
        field_xml.find("MultiPatch").find("patches").text = (
            f"{index_offset} " f"{n_patches - 1  + index_offset}"
        )
        if int(_python_version.split(".")[1]) >= 9 and indent:
            _ET.indent(field_xml)
        file_content = _ET.tostring(field_xml)
        with open(fname + ".fields.xml", "wb") as f:
            f.write(file_content)

    n_patches = _spline_to_ET(
        xml_data,
        multipatch,
        index_offset,
        as_base64=as_base64,
    )

    if n_patches != len(multipatch.patches):
        raise RuntimeError("Help - some patches were not recognised")

    patch_range.text = f"{index_offset} " f"{n_patches - 1 + index_offset}"

    # Add additional options to the xml file
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
                attributes = gismo_dictionary.get("attributes", {})
                option_text = gismo_dictionary.get("text", None)
                comment_text = gismo_dictionary.get("comment", None)
                if comment_text is not None:
                    xml_data.append(_ET.Comment(comment_text))
                optional_data = _ET.SubElement(
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

    if int(_python_version.split(".")[1]) >= 9 and indent:
        # Pretty printing xml with indent only exists in version > 3.9
        _ET.indent(xml_data)

    elif int(_python_version.split(".")[1]) < 9 and indent:
        _debug(
            "Indented xml output is only supported from > python3.9.",
            "Output will not be indented.",
            f"Current python version: {_python_version}",
        )

    file_content = _ET.tostring(xml_data)
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
    from splinepy.multipatch import Multipatch as _Multipatch

    def _matrix_from_node(xml_node):
        """
        Parameters
        ----------
        xml_node : XMLElementTree

        Returns
        -------
        array : np.array
          numpy array in unspecified type, one-dimensional
        """
        type_from_keyword = {
            "b64float64": _np.float64,
            "b64float32": _np.float32,
            "b64uint16": _np.uint16,
            "b64uint32": _np.uint32,
            "b64uint64": _np.uint64,
            "b64int16": _np.int16,
            "b64int32": _np.int32,
            "b64int64": _np.int64,
        }

        # Check format
        format_flag = xml_node.attrib.get("format", "ASCII")
        if format_flag.lower() == "ascii":
            return _np.fromstring(xml_node.text.replace("\n", " "), sep=" ")
        # Check for all excepted binary flags
        format = type_from_keyword.get(format_flag.lower(), None)

        if format is None:
            raise ValueError(
                "Unknown format in xml file: " + format_flag.lower()
            )

        return _np.frombuffer(
            _base64.b64decode(xml_node.text.strip().encode("ascii")),
            dtype=_np.float64,
        )

    # Auxiliary function to unravel
    def retrieve_from_basis_(ETelement, SPdict):
        if "nurbs" in ETelement.attrib["type"].lower():
            # Recursive call for knot_vector
            for child in ETelement:
                if "weights" in child.tag:
                    SPdict["weights"] = _matrix_from_node(child)
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
                    kv = _matrix_from_node(child[0])
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
    _debug(f"Parsing xml-file '{fname}' ...")
    root = _ET.parse(fname).getroot()
    _debug("XML-file parsed start conversion")

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
                _debug("Unsupported format")
            if patch_element.attrib.get("type") != "id_range":
                _debug(
                    f"Invalid patch type {patch_element.attrib.get('type')}"
                )
            patch_range = _matrix_from_node(patch_element).astype(_np.int64)
            offset = patch_range[0]
            n_splines = patch_range[1] - patch_range[0] + 1
            number_of_faces = parametric_dimension * 2

            # Extract interfaces and boundaries
            interfaces_element = child.find("interfaces")
            interface_array = (
                _np.ones((n_splines, number_of_faces), dtype=_np.int64)
                * invalid_integer
            )
            if interfaces_element is None or interfaces_element.text is None:
                # This can occur if only one patch is in a Multipatch
                _debug("No connectivity found")
            else:
                interface_information = (
                    _matrix_from_node(interfaces_element)
                    .astype(
                        dtype=_np.int64,
                    )
                    .reshape(-1, number_of_faces + 4)
                )
                # Assign interfaces
                interface_array[
                    interface_information[:, 0] - offset,
                    interface_information[:, 1] - 1,
                ] = (
                    (interface_information[:, 2] - offset) * number_of_faces
                    + interface_information[:, 3]
                    - 1
                )
                interface_array[
                    interface_information[:, 2] - offset,
                    interface_information[:, 3] - 1,
                ] = (
                    (interface_information[:, 0] - offset) * number_of_faces
                    + interface_information[:, 1]
                    - 1
                )

            # Extract interfaces and boundaries
            boundary_elements = child.findall("boundary")
            if boundary_elements is None:
                _debug("No boundary information found")
            else:
                for id, boundary_element in enumerate(
                    boundary_elements, start=1
                ):
                    if boundary_element.text is None:
                        continue
                    boundary_information = (
                        _matrix_from_node(boundary_element)
                        .astype(
                            dtype=_np.int64,
                        )
                        .reshape(-1, 2)
                    )
                    interface_array[
                        boundary_information[:, 0] - offset,
                        boundary_information[:, 1] - 1,
                    ] = -id

        elif child.tag.startswith("Geometry"):
            spline_dict = {}
            _debug(f"Found new spline in xml with id {child.attrib.get('id')}")
            for info in child:
                if "basis" in info.tag.lower():
                    retrieve_from_basis_(info, spline_dict)
                elif info.tag.startswith("coef"):
                    dim = int(info.attrib.get("geoDim"))
                    spline_dict["control_points"] = _matrix_from_node(
                        info
                    ).reshape(-1, dim)
            if spline_dict.get("weights") is None:
                list_of_splines.append(
                    _settings.NAME_TO_TYPE["BSpline"](**spline_dict)
                )
            else:
                list_of_splines.append(
                    _settings.NAME_TO_TYPE["NURBS"](**spline_dict)
                )
        elif load_options:
            list_of_options.append(make_dictionary(child))
        else:
            _debug(
                f"Found unsupported keyword {child.tag}, which will be"
                " ignored"
            )
            continue

    _debug(f"Found a total of {len(list_of_splines)} " f"BSplines and NURBS")
    multipatch = _Multipatch(list_of_splines)
    if interface_array is not None:
        if invalid_integer in interface_array:
            _warning("Unmatched faces or insufficient information in xml file")
            invalid_splines = _np.where(interface_array)[0][0] + offset
            _warning(f"Check faces with ids {invalid_splines}")
            _warning("Interface information ignored")

        multipatch.interfaces = interface_array

    if load_options:
        return multipatch, list_of_options
    else:
        return multipatch


def create_boundary_conditions_block(
    bc_id,
    dim,
    function_list,
    bc_list,
    unknown_id=0,
    multipatch_id=0,
    comment=None,
):
    """Create Python dictionary of boundary function's block to be used in
    gismo export function.

    Parameters
    ----------
    bc_id : int
        ID number of boundaryConditions-block
    dim : int
        Dimension of physical space
    function_list : list<str/tuple<str,str,...>>
        List of functions where entries are strings of the function of the BC.
        If tuple is given, it will be interpreted as component-wise functions
        (e.g. for vector valued functions)
    bc_list : list<tuple<str, str, int>>
        List of which boundary attains which type of boundary condition. First
        entry is the boundary name, second the type of condition (e.g.
        "Dirichlet"), third is on which variable/unknown it should be applied
    unknown_id : int (0)
        ID of the unknown/variable on which this boundaryConditions block should be
        applied
    multipatch_id : int (0)
        ID of the multipatch
    comment: str
        Comment above boundaryFunctions-block

    Returns
    -------
    bc_dict_for_xml : dict
      Dictionary of boundary conditions
    """
    bc_dict_for_xml = {}
    bc_dict_for_xml["tag"] = "boundaryConditions"
    bc_dict_for_xml["attributes"] = {
        "id": str(bc_id),
        "multipatch": str(multipatch_id),
    }

    # Define functions
    children_list = []
    for index, func_text in enumerate(function_list):
        function_dict = {
            "tag": "Function",
            "attributes": {
                "type": "FunctionExpr",
                "dim": str(dim),
                "index": str(index),
            },
        }
        if type(func_text) != tuple:
            function_dict["text"] = func_text
        # Else vector-valued function
        else:
            function_dict["attributes"]["c"] = str(dim)
            function_dict["children"] = [
                {
                    "tag": "c",
                    "attributes": {"index": str(component_index)},
                    "text": single_func_text,
                }
                for component_index, single_func_text in enumerate(func_text)
            ]

        children_list.append(function_dict)

    children_list += [
        {
            "tag": "bc",
            "attributes": {
                "name": bidname,
                "type": bctype,
                "unknown": str(unknown_id),
                # component omitted, might be sth for vector-valued unknowns
                "function": str(function_id),
            },
        }
        for bidname, bctype, function_id in bc_list
    ]

    bc_dict_for_xml["children"] = children_list

    if comment is not None:
        bc_dict_for_xml["comment"] = comment

    return bc_dict_for_xml


def create_function_block(dim, function_id, function_string, comment=None):
    """Create Python dictionary of custom function's block to be used in
    gismo export function.

    Parameters
    ----------
    dim: int
        Physical dimension of the problem
    id: int
        ID of function block
    function_string: str or tuple<str,str>
        Function as string. If tuple is given, function will be vector-valued.
    comment: str
        Comment above function block

    Returns
    -------
    function_dict: dict
        Dictionary to be used in export function
    """
    function_dict = {}
    function_dict["tag"] = "Function"
    function_dict["attributes"] = {
        "id": str(function_id),
        "type": "FunctionExpr",
        "dim": str(dim),
    }

    # Scalar-valued function
    if type(function_string) != tuple:
        function_dict["text"] = function_string
    # Vector-valued function
    else:
        function_dict["attributes"]["c"] = str(dim)
        function_dict["children"] = [
            {
                "tag": "c",
                "attributes": {"index": str(component_index)},
                "text": function_component_string,
            }
            for component_index, function_component_string in enumerate(
                function_string
            )
        ]

    if comment is not None:
        function_dict["comment"] = comment

    return function_dict


def create_assembly_options_block(
    options_id,
    dirichlet_strategy=11,
    dirichlet_values=101,
    interface_strategy=1,
    bda=2,
    bdb=1,
    bdo=0.333,
    qua=1,
    qub=1,
    qurule=1,
    comment=None,
):
    """Create Python dictionary of g+smo assembly options to be used in gismo
    export function

    Parameters
    ----------
    options_id: int
        ID of the OptionList-block
    dirichlet_strategy: int
        Method for enforcement of Dirichlet BCs [11..14]
    dirichlet_values: int
        Method for computation of Dirichlet DoF values [100..103]
    interface_strategy: int
        Method of treatment of patch interfaces [0..3]
    bda: int
        Estimated nonzeros per column of the matrix: bdA*deg + bdB
    bdb: int
        Estimated nonzeros per column of the matrix: bdA*deg + bdB
    bdo: int
        Overhead of sparse mem. allocation: (1+bdO)(bdA*deg + bdB) [0..1]
    qua: int
        Number of quadrature points: quA*deg + quB
    qub: int
        Number of quadrature points: quA*deg + quB
    qurule: int
        Quadrature rule [1:GaussLegendre,2:GaussLobatto]
    comment: str
        Comment above OptionList-block

    Returns
    -------
    assembly_dict: dict
        Dictionary to be used in export function
    """
    labels = [
        "DirichletStrategy",
        "DirichletValues",
        "InterfaceStrategy",
        "bdA",
        "bdB",
        "bdO",
        "quA",
        "quB",
        "quRule",
    ]
    descriptions = [
        "Method for enforcement of Dirichlet BCs [11..14]",
        "Method for computation of Dirichlet DoF values [100..103]",
        "Method of treatment of patch interfaces [0..3]",
        "Estimated nonzeros per column of the matrix: bdA*deg + bdB",
        "Estimated nonzeros per column of the matrix: bdA*deg + bdB",
        "Overhead of sparse mem. allocation: (1+bdO)(bdA*deg + bdB) [0..1]",
        "Number of quadrature points: quA*deg + quB",
        "Number of quadrature points: quA*deg + quB",
        "Quadrature rule [1:GaussLegendre,2:GaussLobatto]",
    ]
    values = [
        dirichlet_strategy,
        dirichlet_values,
        interface_strategy,
        bda,
        bdb,
        bdo,
        qua,
        qub,
        qurule,
    ]
    number_types = [
        "int",
        "int",
        "int",
        "real",
        "int",
        "real",
        "real",
        "int",
        "int",
    ]

    assembly_dict = {}
    assembly_dict["tag"] = "OptionList"
    assembly_dict["attributes"] = {"id": str(options_id)}
    assembly_dict["text"] = "\n    "

    children_list = []

    for label, description, value, number_type in zip(
        labels, descriptions, values, number_types
    ):
        children_list.append(
            {
                "tag": number_type,
                "attributes": {
                    "label": label,
                    "desc": description,
                    "value": str(value),
                },
            }
        )

    assembly_dict["children"] = children_list

    if comment is not None:
        assembly_dict["comment"] = comment

    return assembly_dict
