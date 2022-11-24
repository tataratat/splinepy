import xml.etree.ElementTree as ET
import numpy as np

from splinepy.utils.log import debug


def export(
    fname,
    spline_list,
):
    """Export as gismo readable xml geometry file
    Use gismo-specific xml-keywords to export (list of) splines. All Bezier
    patches are excluded as their respective non uniform counterpart

    Parameters
    ----------
    fname : string
      name of the output file
    spline_list : list<spline>, spline
      (list of) Spline-Types in splinepy format

    Returns
    -------
    None
    """
    from splinepy.spline import Spline
    from splinepy import NURBS, BSpline

    if issubclass(type(spline_list), Spline):
        # Transform to list
        spline_list = [spline_list]

    if not isinstance(spline_list, list):
        raise ValueError(
            "export Function expects list for spline_list argument")

    if not isinstance(fname, str):
        raise ValueError("fname argument must be string")

    xml_data = ET.Element('xml')
    for id, spline in enumerate(spline_list):
        if not issubclass(type(spline), Spline):
            raise ValueError(
                "One of the splines handed to export is not a valid spline"
                " representation"
            )

        # Transform bezier types, as they are not supported in gismo
        if spline.whatami.startswith("Bezier"):
            type_name = 'BSpline'
            spline = BSpline(
                    **spline.todict(),
                    knot_vectors=[
                            [0] * (a + 1) + [1] * (a + 1)
                            for a in spline.degrees
                    ]
            )
        elif spline.whatami.startswith("RationalBezier"):
            type_name = 'Nurbs'
            spline = NURBS(
                    **spline.todict(),
                    knot_vectors=[
                            [0] * (a + 1) + [1] * (a + 1)
                            for a in spline.degrees
                    ]
            )
        elif spline.whatami.startswith("BSpline"):
            type_name = 'BSpline'
        elif spline.whatami.startswith("NURBS"):
            type_name = 'Nurbs'

        # Start element definition
        spline_element = ET.SubElement(xml_data,
                                       'Geometry',
                                       type='Tensor' + type_name
                                       + str(spline.para_dim),
                                       id=str(id))

        # Define Basis functions
        if "weights" in spline.required_properties:
            spline_basis_base = ET.SubElement(
                spline_element,
                'Basis',
                type='Tensor' + type_name + 'Basis' + str(spline.para_dim),
            )
            spline_basis = ET.SubElement(
                spline_basis_base,
                'Basis',
                type='TensorBSplineBasis' + str(spline.para_dim),
            )
        else:
            spline_basis = ET.SubElement(
                spline_element,
                'Basis',
                type='Tensor' + type_name + 'Basis' + str(spline.para_dim),
            )


        for i_para in range(spline.para_dim):
            basis_fun = ET.SubElement(spline_basis,
                                      'Basis',
                                      type='BSplineBasis',
                                      index=str(i_para))
            knot_vector = ET.SubElement(basis_fun,
                                        'KnotVector',
                                        degree=str(spline.degrees[i_para]))
            knot_vector.text = ' '.join(
                [str(k) for k in spline.knot_vectors[i_para]])
        if "weights" in spline.required_properties:
            # Add weights
            weights = ET.SubElement(
                spline_basis_base,
                'weights',
            )
            weights.text = '\n'.join([str(w)
                                     for w in spline.weights.flatten()])

        coords = ET.SubElement(
            spline_element,
            'coefs',
            geoDim=str(spline.dim),
        )
        coords.text = '\n'.join(
            [' '.join([str(xx) for xx in x]) for x in spline.control_points])

    ET.indent(xml_data)
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
                        child.text.replace("\n", " "), sep=' ').reshape(-1,dim)
                if "basis" in child.tag.lower():
                    retrieve_from_basis_(child, SPdict)
        elif "bspline" in ETelement.attrib["type"].lower():
            para_dim = int(ETelement.attrib["type"].lower().split("basis")[-1])
            knotvector = [None] * para_dim
            degrees = [None] * para_dim
            for child in ETelement:
                if not "bsplinebasis" in child.attrib["type"].lower():
                    raise ValueError("Something went wrong")
                if "knotvector" in child[0].tag.lower():
                    id = int(child.attrib["index"])
                    kv = np.fromstring(
                        child[0].text.replace("\n", " "), sep=' ')
                    knotvector[id] = kv
                    degrees[id] = int(child[0].attrib["degree"])
            SPdict["knot_vectors"] = knotvector
            SPdict["degrees"] = degrees
    # Parse XML file
    debug(f"Parsing xml-file '{fname}' ...")
    root = ET.parse(fname).getroot()
    debug(f"XML-file parsed start conversion")

    # Init return value
    return_dict = {
        "BSpline" : [],
        "NURBS" : []
        }

    # Splines start with the keyword Geometry
    for child in root:
        if not child.tag.startswith("Geometry"):
            debug(
                f"Found unsupported keyword {child.tag}, which will be"
                " ignored"
            )
        spline_dict = {}
        debug(
            f"Found new spline in xml file with id {child.attrib.get('id')}"
        )
        for info in child:
            if "basis" in info.tag.lower():
                retrieve_from_basis_(info, spline_dict)
            elif info.tag.startswith("coef"):
                dim = int(info.attrib.get("geoDim"))
                spline_dict["control_points"]=np.fromstring(
                    info.text.replace("\n", " "), sep=' ').reshape(-1,dim)
        if spline_dict.get("weights") is None:
            return_dict["BSpline"].append(spline_dict)
        else:
            return_dict["NURBS"].append(spline_dict)

    debug(
        f"Found a total of {len(return_dict['BSpline'])} "
        f"BSplines and {len(return_dict['NURBS'])} NURBS"
    )

    return return_dict
