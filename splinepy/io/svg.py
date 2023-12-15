import xml.etree.ElementTree as ET
from sys import version as python_version

import numpy as np
from vedo.colors import get_color as _get_color

from splinepy.utils.log import debug


def _export_control_mesh(
    spline, svg_spline_element, box_min_x, box_max_y, box_margins
):
    """
    Export a spline's control mesh in svg format using polylines for the mesh
    lines, circles for the control points.

    Parameters
    ----------
    spline : Spline
      Providing the control points for the mesh
    svg_spline_element : xml.etree.ElementTree
      Parent element of the spline (most likely spline group), where the mesh
      is written into
    box_min_x : float/int
      minimum coordinate of the box, providing the necessary offsets for spline
      export
    box_max_y : float/int
      maximum y coordinate of the surrounding box (pictures use LHS coordinate
      system in top left corner)
    box_margins : float/int
      box margins are currently used to determine the text-offset, more options
      will be added soon
    """
    # First draw all control points
    if not spline.show_options.get("control_points", True):
        return

    # Relevant options:
    # - control_point_ids
    # - control_point_alpha
    # - control_point_c
    # - control_point_r

    # transform color into (relative) rgb using vedo helper function
    r, g, b = _get_color(spline.show_options.get("control_point_c", "red"))
    a = spline.show_options.get("control_point_alpha", 1.0)
    svg_mesh = ET.SubElement(
        svg_spline_element,
        "g",
        id="control_mesh",
    )
    svg_control_points = ET.SubElement(
        svg_mesh,
        "g",
        id="control_points",
        fill=f"rgba({100*r}%,{100*g}%,{100*b}%,{100*a}%)",
    )

    for ctp in spline.control_points:
        ET.SubElement(
            svg_control_points,
            "circle",
            cx=str(ctp[0] - box_min_x),
            cy=str(box_max_y - ctp[1]),
            r=str(spline.show_options.get("control_point_r", 0.02)),
        )
    if spline.show_options.get("control_point_ids", True):
        # Text Options
        svg_control_point_ids = ET.SubElement(
            svg_mesh,
            "g",
            id="control_point_ids",
        )

        # Set text options
        svg_control_point_ids.attrib["font-family"] = "sans-serif"
        svg_control_point_ids.attrib["font-size"] = ".1"
        svg_control_point_ids.attrib["text-anchor"] = "middle"
        svg_control_point_ids.attrib["font-family"] = "sans-serif"
        svg_control_point_ids.attrib["fill"] = "black"
        svg_control_point_ids.attrib["stroke"] = "none"

        for i, ctp in enumerate(spline.control_points):
            text_element = ET.SubElement(
                svg_control_point_ids,
                "text",
                x=str(ctp[0] - box_min_x),
                y=str(box_max_y - ctp[1]),
                dx=str(0.5 * box_margins),
                dy=str(0.5 * box_margins),
            )
            text_element.text = str(i)

    # Then connect the points using a control mesh
    if spline.show_options.get("control_mesh", True):
        # Retrieve options
        r, g, b = _get_color(spline.show_options.get("control_mesh_c", "red"))
        a = spline.show_options.get("control_mesh_alpha", 1.0)

        # Create a new group
        svg_mesh_polylines = ET.SubElement(
            svg_mesh,
            "g",
            id="mesh",
        )

        # Relevant options:
        # - control_mesh_c
        # - control_mesh_lw
        # - control_mesh_alpha (tbd)
        svg_mesh_polylines.attrib["fill"] = "none"
        svg_mesh_polylines.attrib[
            "stroke"
        ] = f"rgba({100*r}%,{100*g}%,{100*b}%,{100*a}%)"
        svg_mesh_polylines.attrib["stroke-width"] = spline.show_options.get(
            "control_mesh_lw", ".01"
        )

        # Draw the actual lines
        if spline.para_dim == 1:
            ET.SubElement(
                svg_mesh_polylines,
                "polyline",
                points=" ".join(
                    [
                        str(xx - box_min_x) + "," + str(box_max_y - xy)
                        for (xx, xy) in spline.cps
                    ]
                ),
            )
        if spline.para_dim == 2:
            mi = spline.multi_index
            cmr_xi, cmr_eta = spline.control_mesh_resolutions
            for i in range(cmr_eta):
                # for
                ET.SubElement(
                    svg_mesh_polylines,
                    "polyline",
                    points=" ".join(
                        [
                            str(xx - box_min_x) + "," + str(box_max_y - xy)
                            for (xx, xy) in spline.cps[mi[:, i]]
                        ]
                    ),
                )
            for i in range(cmr_xi):
                # for
                ET.SubElement(
                    svg_mesh_polylines,
                    "polyline",
                    points=" ".join(
                        [
                            str(xx - box_min_x) + "," + str(box_max_y - xy)
                            for (xx, xy) in spline.cps[mi[i, :]]
                        ]
                    ),
                )

    pass


def _export_spline(
    spline, svg_spline_element, box_min_x, box_max_y, tolerance
):
    """
    Export a spline's control mesh in svg format using polylines for the mesh
    lines, circles for the control points.

    In order to provide the necessary tools for knot lines and simultaneously
    fill the spline in 2D, surfaces will be plotted without outlines first,
    whereas
    Parameters
    ----------
    spline : Spline
      Providing the control points for the mesh
    svg_spline_element : xml.etree.ElementTree
      Parent element of the spline (most likely spline group), where the mesh
      is written into
    box_min_x : float/int
      minimum coordinate of the box, providing the necessary offsets for spline
      export
    box_max_y : float/int
      maximum y coordinate of the surrounding box (pictures use LHS coordinate
      system in top left corner)
    box_margins : float/int
      box margins are currently used to determine the text-offset, more options
      will be added soon
    """

    svg_paths = ET.SubElement(
        svg_spline_element,
        "g",
        id="spline_paths",
    )

    if spline.para_dim == 1:
        # Export is done in 2 stages
        # (1) Approximation (trying to preserve the continuity)
        if spline.degrees[0] <= 3 and not spline.is_rational:
            spline_copy = spline.copy()
            spline_copy.elevate_degrees([0] * (3 - spline.degrees[0]))
        else:
            # Create a cubic b-spline with knot multiplicities
            # min(3 - (degree - original_m), 1)
            # Approximate it
            # Evaluate both splines (maybe at greville aabscissae? -> For
            # rational splines more points should be tested)
            # If distance is too large -> refine the knots where the error is
            # too large and approximate again, repeat (not very efficient, but
            # should do the trick)
            raise NotImplementedError(
                f"I am working on it. Approximate with tolerance {tolerance}"
            )

        # Sanity check
        assert spline_copy.degrees[0] == 3
        assert not spline_copy.is_rational

        # Relevant options
        # - c
        # - alpha
        # - knots

        # Bezier Extract and join the paths
        r, g, b = _get_color(spline.show_options.get("c", "green"))
        a = spline.show_options.get("alpha", 1.0)
        splines = spline_copy.extract.beziers()
        svg_path = ET.SubElement(
            svg_paths,
            "path",
            # draw path
            d=" ".join(
                [
                    (
                        f"M {s.cps[0,0]-box_min_x},{box_max_y - s.cps[0,1]} C"
                        f" {s.cps[1,0]-box_min_x},{box_max_y - s.cps[1,1]}"
                        f" {s.cps[2,0]-box_min_x},{box_max_y - s.cps[2,1]}"
                        f" {s.cps[3,0]-box_min_x},{box_max_y - s.cps[3,1]}"
                    )
                    for s in splines
                ]
            ),
        )
        svg_path.attrib["fill"] = "none"
        svg_path.attrib["stroke"] = (
            f"rgba({100*r}%,{100*g}%" f",{100*b}%,{100*a}%)"
        )
        svg_path.attrib["stroke-width"] = str(
            spline.show_options.get("lw", 0.01)
        )

        # Plot knots
        if spline.show_options.get("knots", True):
            ukvs = spline.unique_knots[0]
            knot_projected = spline.evaluate(ukvs.reshape(-1, 1))
            for x, y in knot_projected:
                # Relevant options
                # - knot_alpha
                # - knot_c
                # - knot_lw
                r, g, b = _get_color(
                    spline.show_options.get("knot_alpha", "blue")
                )
                a = spline.show_options.get("knot_alpha", 1.0)
                dx = spline.show_options.get("knot_lw", 0.02)
                ET.SubElement(
                    svg_paths,
                    "rect",
                    fill=f"rgba({100*r}%,{100*g}%" f",{100*b}%,{100*a}%)",
                    x=str(x - box_min_x - 0.5 * dx),
                    y=str(box_max_y - y - 0.5 * dx),
                    height=str(dx),
                    width=str(dx),
                )

    if spline.para_dim == 2:
        raise NotImplementedError("I am working on this, sorry")

    pass


def export(fname, *splines, indent=True, box_margins=0.1):
    """
    Exports a number of splines into an svg plot

    Graphics are built to resemble vedo plots as close as possible by relying
    on spline show_options

    Parameters
    ----------
    fname : string
        name of the output file
    splines : Spline-Type
        Splines to be exported
    indent: bool
      Appends white spaces using xml.etree.ElementTree.indent, if possible.
      Makes it more human-readable
    box_margins : int, float
      Offset around the splines (relative units)

    Returns
    -------
    None
    """
    from splinepy.spline import Spline

    # Check arguments
    if not isinstance(fname, str):
        raise ValueError("fname argument must be string")
    if not isinstance(box_margins, (int, float)):
        raise ValueError("fname argument must be string")

    for spline in splines:
        if not isinstance(spline, Spline):
            raise ValueError("Unexpected spline type")
        if not 1 <= spline.para_dim <= 2:
            raise ValueError("Only lines and surfaces supported")
        if not spline.dim == 2:
            raise ValueError("Only 2D geometries supported")

    # Determine bounding box of all elements
    ctps_bounds = splines[0].control_point_bounds
    for spline in splines[1::]:
        ctps_bounds[0, :] = np.minimum(
            spline.control_point_bounds[0, :], ctps_bounds[0, :]
        )
        ctps_bounds[1, :] = np.maximum(
            spline.control_point_bounds[1, :], ctps_bounds[1, :]
        )

    # Initialize svg element
    svg_data = ET.Element("svg")
    box_size = (
        ctps_bounds[1, 0] - ctps_bounds[0, 0] + 2 * box_margins,
        ctps_bounds[1, 1] - ctps_bounds[0, 1] + 2 * box_margins,
    )
    svg_data.attrib["viewBox"] = f"0 0 {box_size[0]} {box_size[1]}"
    svg_data.attrib["height"] = str(400)
    svg_data.attrib["width"] = str(400 / box_size[1] * box_size[0])
    svg_data.attrib["style"] = "background-color:white"
    svg_data.attrib["xmlns"] = "http://www.w3.org/2000/svg"

    # Determine offsets
    box_min_x = ctps_bounds[0, 0] - box_margins
    box_max_y = ctps_bounds[1, 1] + box_margins

    # Write splines in svg file
    for i, spline in enumerate(splines):
        # Put every spline into a new dedicated group
        spline_group = ET.SubElement(svg_data, "g", id="spline" + str(i))
        _export_control_mesh(
            spline, spline_group, box_min_x, box_max_y, box_margins
        )
        _export_spline(spline, spline_group, box_min_x, box_max_y, 1e-5)

    # Dump into file
    if int(python_version.split(".")[1]) >= 9 and indent:
        # Pretty printing xml with indent only exists in version > 3.9
        ET.indent(svg_data)

    elif int(python_version.split(".")[1]) < 9 and indent:
        debug(
            "Indented xml output is only supported from > python3.9.",
            "Output will not be indented.",
            f"Current python version: {python_version}",
        )

    file_content = ET.tostring(svg_data)
    with open(fname, "wb") as f:
        f.write(file_content)
    pass
