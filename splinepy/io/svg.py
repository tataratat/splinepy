import xml.etree.ElementTree as _ET

import numpy as _np

try:
    from vedo import Plotter as _Plotter
    from vedo import color_map as _color_map
    from vedo.colors import get_color as _get_color
except ImportError as err:
    from gustaf.helpers.raise_if import ModuleImportRaiser

    _error_message_vedo_import = "Vedo library required for svg export"
    _Plotter = ModuleImportRaiser(_error_message_vedo_import, err)
    _color_map = ModuleImportRaiser(_error_message_vedo_import, err)
    _get_color = ModuleImportRaiser(_error_message_vedo_import, err)

from splinepy.utils.log import debug as _debug
from splinepy.utils.log import warning as _warning


def _rgb_2_hex(r, g, b):
    """
    Helper function to convert (r,g,b) to hex format

    Parameters
    ----------
    r : float [0,1]
      red value
    g : float [0,1]
      green value
    b : float [0,1]
      blue value

    Returns
    -------
    hex : string
    """
    return f"#{int(255 * r):02x}{int(255 * g):02x}{int(255 * b):02x}"


def _export_spline_field(spline, svg_element, box_min_x, box_max_y):
    """
    Export a spline's data-field as a pixel-bitmap. The pixel density is based
    on the sample resolution, (10 * resolution). The image is embedded in b64
    format

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

    Returns
    -------
    None
    """

    from splinepy.helpme.visualize import _process_scalar_field, _sample_spline

    def _write_png(bitmap, as_base64=True):
        """
        Write PNG data from a buffer.

        This function is taken from
        https://stackoverflow.com/questions/56564977/

        Nitesh Menon's answer with some comments added to it

        Replacing this by an external library like PIL or imageio is a @todo
        still.

        Parameters
        ----------
        bitmap : np.ndarray
          Picture with dimensions (width, height, 4) as RGBA image
        as_base64 : bool
          Write into readable base64 format

        Returns
        -------
        png : string
          Base64 encoded String with png data
        """
        import base64
        import struct  # Packing data
        import zlib  # Compression

        # Extract data
        buf = bytearray(bitmap)
        if (
            (bitmap.ndim != 3)
            or (bitmap.shape[2] != 4)
            or (bitmap.dtype is not _np.dtype("uint8"))
        ):
            raise ValueError("Was expecting bitmap in RGBA format")
        width = bitmap.shape[1]
        height = bitmap.shape[0]

        # Reverse the vertical line order and add null bytes at the start
        width_byte_4 = width * 4

        # Note: Contrary to the original post, the array is not concatenated in
        # reverse but in the original order (modified by jzwar)
        raw_data = b"".join(
            # Add a null byte and concatenate the raw data
            b"\x00" + buf[span : (span + width_byte_4)]
            # Iterate over the buffer
            for span in range(0, (height - 1) * width_byte_4 + 1, width_byte_4)
        )

        def png_pack(png_tag, data):
            """
            Pack PNG data into a chunk according to png standard.

            Parameters:
            png_tag : bytes
              PNG tag identifier.
            data : bytes
              Data to be packed.

            Returns
            -------
            bytes
              Packed PNG chunk.
            """
            chunk_head = png_tag + data  # Concatenate the tag and data
            return (
                # Pack the length of chunk
                struct.pack("!I", len(data))
                +
                # Concatenate the chunk header
                chunk_head
                +
                # Pack CRC-32 checksum
                struct.pack("!I", 0xFFFFFFFF & zlib.crc32(chunk_head))
            )

        pure_data = b"".join(
            [
                b"\x89PNG\r\n\x1a\n",  # PNG file signature
                png_pack(
                    b"IHDR", struct.pack("!2I5B", width, height, 8, 6, 0, 0, 0)
                ),  # PNG header chunk
                # PNG data chunk (compressed)
                png_pack(b"IDAT", zlib.compress(raw_data, 9)),
                png_pack(b"IEND", b""),  # PNG end chunk
            ]
        )

        # Return encoded
        if as_base64:
            return base64.b64encode(pure_data)
        else:
            return pure_data

    # Sample the spline
    resolution = spline.show_options.get("resolutions", 100)
    sampled_spline = _sample_spline(spline, resolution)
    data = spline.show_options.get("data", None)

    # Check
    if data is None:
        raise ValueError(
            "There is no data provided, although data plot is requested"
        )

    # Process the scalar field
    _process_scalar_field(spline, data, sampled_spline, res=resolution)
    # scalar bar distorts the dimensions of spline bitmap image
    sampled_spline.show_options["scalarbar"] = False
    # Set a resolution (This needs to be improved still)
    pixel_size = (resolution * 10, resolution * 10)

    # Create a Plotter
    plotter = _Plotter(
        shape=(1, 1),  # Only one field
        N=1,  # Number of things to plot
        size=pixel_size,  # Resolution, i.e. pixel density
        sharecam=True,
        offscreen=True,
        title="",
        bg=(255, 255, 255),
        axes=0,
    )
    plotter.show(sampled_spline.showable(), zoom="tightest")

    # Extract bounding box
    x_min, y_min, x_max, y_max = sampled_spline.bounds().ravel()

    # Extract the bitmap and transform to RGBA
    bitmap = plotter.screenshot(asarray=True)

    if bitmap.shape[2] == 3:
        alpha_layer = (
            _np.ones((bitmap.shape[0], bitmap.shape[1]), dtype=bitmap.dtype)
            * 255
        )
        alpha_layer[_np.all(bitmap == 255, axis=-1)] = 0
        bitmap = _np.concatenate(
            (bitmap, alpha_layer.reshape(*alpha_layer.shape, 1)), axis=-1
        )
    elif bitmap.shape[2] == 4:
        alpha_layer = bitmap[:, :, -1]
    else:
        raise RuntimeError("Failed to extract field bitmap")

    # Crop image
    pixel_layer = alpha_layer != 0
    has_pixels = _np.where(_np.any(pixel_layer, axis=1))[0]
    min_x_pixel = has_pixels.min()
    max_x_pixel = has_pixels.max()
    has_pixels = _np.where(_np.any(pixel_layer, axis=0))[0]
    min_y_pixel = has_pixels.min()
    max_y_pixel = has_pixels.max()
    bitmap = bitmap[min_x_pixel:max_x_pixel, min_y_pixel:max_y_pixel, :]

    # Write Bitmap into svg
    image_svg = _ET.SubElement(
        svg_element,
        "image",
        x=str(x_min - box_min_x),  # str(box_max_y - xy)
        y=str(box_max_y - y_max),
        width=str(x_max - x_min),
        height=str(y_max - y_min),
    )
    image_svg.attrib["xlink:href"] = "data:image/png;base64," + _write_png(
        bitmap=bitmap, as_base64=True
    ).decode("utf-8")


def _export_control_mesh(
    spline, svg_spline_element, box_min_x, box_max_y, **kwargs
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
    kwargs :
      Key word argument options for output specification
      `{"font-family","font-size","text-anchor","font-family","fill","stroke"}`

    Returns
    -------
    None
    """

    # Default behaviour in show is, no mesh and no ids, if no control points
    if not spline.show_options.get("control_points", True):
        return None

    svg_mesh = _ET.SubElement(
        svg_spline_element,
        "g",
        id="control_mesh",
    )

    # First draw mesh
    if spline.show_options.get("control_mesh", True):
        # Relevant options:
        # - control_mesh_c
        # - control_mesh_lw
        # - control_mesh_alpha (tbd)
        r, g, b = _get_color(spline.show_options.get("control_mesh_c", "red"))
        a = spline.show_options.get("control_mesh_alpha", 1.0)
        stroke_width = spline.show_options.get(
            "control_mesh_lw",
            _np.linalg.norm(_np.diff(spline.control_point_bounds)) * 0.001,
        )

        # Create a new group
        svg_mesh_polylines = _ET.SubElement(
            svg_mesh,
            "g",
            id="mesh",
            style=(
                f"fill:none;stroke:{_rgb_2_hex(r,g,b)};stroke-opacity:{a};"
                f"stroke-width:{stroke_width};stroke-linecap:round"
            ),
        )

        # Draw the actual lines
        if spline.para_dim == 1:
            _ET.SubElement(
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
                _ET.SubElement(
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
                _ET.SubElement(
                    svg_mesh_polylines,
                    "polyline",
                    points=" ".join(
                        [
                            str(xx - box_min_x) + "," + str(box_max_y - xy)
                            for (xx, xy) in spline.cps[mi[i, :]]
                        ]
                    ),
                )

    # Then control points
    if spline.show_options.get("control_points", True):

        # Relevant options:
        # - control_point_ids
        # - control_point_alpha
        # - control_point_c
        # - control_point_r

        # transform color into (relative) rgb using vedo helper function
        r, g, b = _get_color(spline.show_options.get("control_point_c", "red"))
        a = spline.show_options.get("control_point_alpha", 1.0)

        svg_control_points = _ET.SubElement(
            svg_mesh,
            "g",
            id="control_points",
            style=f"fill:{_rgb_2_hex(r, g, b)};fill-opacity:{a}",
        )

        for ctp in spline.control_points:
            _ET.SubElement(
                svg_control_points,
                "circle",
                cx=str(ctp[0] - box_min_x),
                cy=str(box_max_y - ctp[1]),
                r=str(spline.show_options.get("control_point_r", 0.02)),
            )

    # Lastly IDs
    if spline.show_options.get("control_point_ids", True):
        # Text Options
        svg_control_point_ids = _ET.SubElement(
            svg_mesh,
            "g",
            id="control_point_ids",
        )

        # Set text options
        svg_control_point_ids.attrib["font-family"] = kwargs.get(
            "font_family", "sans-serif"
        )
        svg_control_point_ids.attrib["font-size"] = str(
            kwargs.get("font_size", 0.1)
        )
        svg_control_point_ids.attrib["text-anchor"] = kwargs.get(
            "text_anchor", "middle"
        )
        svg_control_point_ids.attrib["fill"] = _rgb_2_hex(
            *_get_color(kwargs.get("text_color", "k"))
        )
        svg_control_point_ids.attrib["stroke"] = "none"

        # Text offset
        dx, dy = 0.0, 0.0

        for i, ctp in enumerate(spline.control_points):
            text_element = _ET.SubElement(
                svg_control_point_ids,
                "text",
                x=str(ctp[0] - box_min_x),
                y=str(box_max_y - ctp[1]),
                dx=str(dx),
                dy=str(dy),
            )
            text_element.text = str(i)


def _quiver_plot(
    spline, svg_spline_element, box_min_x, box_max_y, tolerance=0.0, **kwargs
):
    """
    Export a quiver plot using arrows as defined
    https://matplotlib.org/stable/_images/quiver_sizes.svg

    The arrows will be scaled using the `arrow_data_scale` flag

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
    tolerance : float
      tolerance as a cut-off length under which no arrow is plotted
    kwargs :
      Optional for arrow sizing,
      `{"arrow_width", "arrow_head_width", "arrow_head_length",
      "arrow_headaxis_length"}`

    Returns
    -------
    None
    """

    # Get arrow scales
    arrow_width = kwargs.get("arrow_width", 0.1)
    arrow_head_width = kwargs.get("arrow_head_width", 0.3)
    arrow_head_length = kwargs.get("arrow_head_length", 0.5)
    arrow_headaxis_length = kwargs.get("arrow_headaxis_length", 0.45)

    # Create a default arrow along x-axis with length 1 that will be scaled
    # down accordingly
    default_arrow_control_points_ = _np.array(
        [
            [0, -arrow_width],
            [1 - arrow_headaxis_length, -arrow_width],
            [1 - arrow_head_length, -arrow_head_width],
            [1, 0],
            [1 - arrow_head_length, arrow_head_width],
            [1 - arrow_headaxis_length, arrow_width],
            [0, arrow_width],
        ]
    ) * [1, 0.5]

    # Obtain Data to be plotted
    arrow_data_field = spline.show_options.get("arrow_data", None)
    if arrow_data_field is None:
        return None

    # Use Arrow Data to retrieve spline data (is already scaled)
    arrow_data = spline.extract.arrow_data(arrow_data_field)
    positions = arrow_data.vertices[arrow_data.edges[:, 0]]
    directions = arrow_data.vertices[arrow_data.edges[:, 1]] - positions

    # Discard all points that are below the tolerance
    arrow_length = _np.linalg.norm(directions, axis=1)
    over_tolerance = arrow_length > tolerance
    directions = directions[over_tolerance]
    positions = positions[over_tolerance]
    arrow_length = arrow_length[over_tolerance]

    # Create rotation matrices
    angles = _np.arctan2(directions[:, 1], directions[:, 0])
    rotation_matrices = _np.empty((angles.shape[0], 2, 2))
    rotation_matrices[:, 0, 0] = _np.cos(angles)
    rotation_matrices[:, 1, 0] = _np.sin(angles)
    rotation_matrices[:, 0, 1] = -rotation_matrices[:, 1, 0]
    rotation_matrices[:, 1, 1] = rotation_matrices[:, 0, 0]

    # Retrieve information on colors and values
    v_min = spline.show_options.get("vmin", _np.min(arrow_length))
    v_max = spline.show_options.get("vmax", _np.max(arrow_length))
    cmap_style = spline.show_options.get("cmap", "jet")
    colors = _color_map(arrow_length, name=cmap_style, vmin=v_min, vmax=v_max)

    # Map arrow control points
    arrow_control_points_ = _np.einsum(
        "nij,kj,n->nki",
        rotation_matrices,
        default_arrow_control_points_,
        arrow_length,
    )

    # Create a new group
    svg_quiver = _ET.SubElement(
        svg_spline_element,
        "g",
        id="quiver",
    )

    # Loop over arrows and write them into a group
    for i in range(arrow_length.shape[0]):
        _ET.SubElement(
            svg_quiver,
            "polyline",
            points=" ".join(
                [
                    str(xx - box_min_x) + "," + str(box_max_y - xy)
                    for (xx, xy) in (
                        arrow_control_points_[i, :, :] + positions[i, :]
                    )
                ]
            ),
            style=(f"fill:{_rgb_2_hex(*colors[i])};stroke:none"),
        )

    return None


def _export_spline(
    spline, svg_spline_element, box_min_x, box_max_y, tolerance=None, **kwargs
):
    """
    Export a spline in svg format cubic beziers as approximations (even for
    lower degree splines, for simplicity). Also exports the knot lines

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
    tolerance : float
      Splines that can not be represented exactly by cubic B-Splines (i.e.,
      rationals or high order splines) a tolerance is required for the
      approximation. Default uses an absolute deviation of 1% of the bounding
      box of the given spline
    kwargs:
      `{"linecap"}`

    Returns
    -------
    None
    """
    from splinepy.helpme.fit import curve as _fit_curve
    from splinepy.helpme.reparametrize import flip_axes as _flip_axes
    from splinepy.settings import TOLERANCE

    # Maximum number of refinements for approximation
    MAX_ITERATION = 4

    # Set tolerance for export to default if no user data
    if tolerance is None:
        tolerance = 0.01 * _np.linalg.norm(
            spline.control_point_bounds[0, :]
            - spline.control_point_bounds[1, :]
        )

    # Sanity check
    if spline.para_dim not in [1, 2]:
        raise ValueError("Unsupported spline dimension")

    # Approximation curve-wise
    def _approximate_curve(original_spline, tolerance):
        """
        Approximate a curve using third order B-Splines

        Continuity will be preserved

        original_spline : Spline
          Any type spline as a basis for approximation
        tolerance : double
          approximation tolerance (given a maximum number of iterations)

        Returns
        -------
        """
        if original_spline.degrees[0] <= 3 and not (
            original_spline.is_rational and original_spline.degrees[0] > 1
        ):
            spline_approximation = original_spline.copy()
            spline_approximation.elevate_degrees(
                [0] * (3 - original_spline.degrees[0])
            )
        else:
            # Use fit tool to approximate curve
            _debug(
                "SVG export only supports (up to) cubic polynomial splines --"
                " using approximation"
            )

            # Queries
            para_queries = original_spline.greville_abscissae(
                duplicate_tolerance=TOLERANCE
            )

            # For rational splines this might be insufficient
            if original_spline.is_rational and (
                original_spline.degrees[0] < 3
            ):
                para_queries = _np.sort(
                    _np.vstack(
                        (
                            para_queries,
                            _np.convolve(
                                para_queries.ravel(),
                                _np.ones(2) * 0.5,
                                "valid",
                            ).reshape(-1, 1),
                        )
                    ),
                    axis=0,
                )

            # Create knot-vector
            k_mult = original_spline.knot_multiplicities[0]
            k_mult[1:-1] = _np.maximum(
                1, 3 - original_spline.degrees[0] + k_mult[1:-1]
            )
            k_mult[0] = 4
            k_mult[-1] = 4
            new_knot_vector = _np.repeat(original_spline.unique_knots, k_mult)
            residual = 2 * tolerance
            spline_approximation, residual = _fit_curve(
                original_spline.evaluate(para_queries),
                degree=3,
                knot_vector=new_knot_vector,
                associated_queries=para_queries,
            )

            # Refine until approximation is found
            for _ in range(MAX_ITERATION):
                if residual < tolerance:
                    break
                # Loop until tolerance satisfied
                para_queries = _np.sort(
                    _np.vstack(
                        (
                            para_queries,
                            _np.convolve(
                                para_queries.ravel(),
                                _np.ones(2) * 0.5,
                                "valid",
                            ).reshape(-1, 1),
                        )
                    ),
                    axis=0,
                )

                # Insert a few knots
                spline_approximation.insert_knots(
                    0,
                    _np.convolve(
                        spline_approximation.unique_knots[0].ravel(),
                        _np.ones(2) * 0.5,
                        "valid",
                    ),
                )
                new_knot_vector = spline_approximation.knot_vectors[0]

                # Create better approximation
                spline_approximation, residual = _fit_curve(
                    original_spline.evaluate(para_queries),
                    degree=3,
                    knot_vector=new_knot_vector,
                    associated_queries=para_queries,
                )

            if residual > tolerance:
                _warning(
                    "Requested tolerance could not be reached within maximum"
                    " number of refinement steps"
                )

        # Sanity check
        if not (
            (spline_approximation.degrees[0] == 3)
            or (not spline_approximation.is_rational)
            or (
                original_spline.is_rational and original_spline.degrees[0] <= 1
            )
        ):
            raise RuntimeError(
                "Spline Approximation returned unexpected result that is "
                "non-compliant to svg standard. Expected cubic bezier spline."
            )

        return spline_approximation

    # Write the actual spline as the lowest layer
    svg_spline = _ET.SubElement(
        svg_spline_element,
        "g",
        id="spline_paths",
    )

    # Check if a field is to be plotted
    if spline.show_options.get("data", None) is not None:

        # spline.show()
        _export_spline_field(spline, svg_spline, box_min_x, box_max_y)

    else:
        # Export is done in 2 stages
        # (1) Approximation (trying to preserve the continuity)
        # (2) Export

        # First : Retrieve export options
        r, g, b = _get_color(spline.show_options.get("c", "green"))
        a = spline.show_options.get("alpha", 1.0)

        # check if the parametric dimension is 2, if the case, extract boundary
        if spline.para_dim == 1:
            # Approximation
            spline_copy = _approximate_curve(spline, tolerance)
            bezier_elements = spline_copy.extract.beziers()
            path_d = (
                f"M {bezier_elements[0].cps[0,0]-box_min_x},"
                f"{box_max_y - bezier_elements[0].cps[0,1]}"
            )
            path_d += " ".join(
                [
                    (
                        f" C {s.cps[1,0]-box_min_x},{box_max_y - s.cps[1,1]}"
                        f" {s.cps[2,0]-box_min_x},{box_max_y - s.cps[2,1]}"
                        f" {s.cps[3,0]-box_min_x},{box_max_y - s.cps[3,1]}"
                    )
                    for s in bezier_elements
                ]
            )
            _ET.SubElement(
                svg_spline,
                "path",
                # draw path
                d=path_d,
                style=(
                    f"fill:none;stroke:{_rgb_2_hex(r,g,b)};stroke-opacity:{a};"
                    f"stroke-width:{spline.show_options.get('lw', 0.01)};"
                    "stroke-linecap:round"
                ),
            )
        else:
            # Extract boundaries
            spline_boundaries = spline.extract.boundaries()

            # Flip boundaries on boundary 0 and 3
            _flip_axes(spline_boundaries[0], axes=[0], inplace=True)
            _flip_axes(spline_boundaries[3], axes=[0], inplace=True)

            bezier_elements = []
            for i in [2, 1, 3, 0]:
                spline_copy = _approximate_curve(
                    spline_boundaries[i], tolerance
                )
                bezier_elements += spline_copy.extract.beziers()

            path_d = (
                f"M {bezier_elements[0].cps[0,0]-box_min_x},"
                f"{box_max_y - bezier_elements[0].cps[0,1]}"
            )
            path_d += " ".join(
                [
                    (
                        f" C"
                        f" {s.cps[1,0]-box_min_x},{box_max_y - s.cps[1,1]}"
                        f" {s.cps[2,0]-box_min_x},{box_max_y - s.cps[2,1]}"
                        f" {s.cps[3,0]-box_min_x},{box_max_y - s.cps[3,1]}"
                    )
                    for s in bezier_elements
                ]
            )

            _ET.SubElement(
                svg_spline,
                "path",
                # draw path
                d=path_d,
                style=(
                    f"fill:{_rgb_2_hex(r,g,b)};fill-opacity:{a};stroke:none;"
                    f"stroke-linecap:{kwargs.get('linecap','round')}"
                ),
            )

        # Plot knots as rectangles
    if spline.show_options.get("knots", True):
        # Exports knots in separate group
        svg_knots = _ET.SubElement(
            svg_spline_element,
            "g",
            id="knots",
        )

        # Relevant options
        # - knot_alpha
        # - knot_c
        # - knot_lw
        r, g, b = _get_color(spline.show_options.get("knot_c", "blue"))
        a = spline.show_options.get("knot_alpha", 1.0)
        lw = spline.show_options.get("knot_lw", 0.02)
        if spline.para_dim == 1:
            ukvs = spline.unique_knots[0]
            knot_projected = spline.evaluate(ukvs.reshape(-1, 1))
            for x, y in knot_projected:

                _ET.SubElement(
                    svg_knots,
                    "rect",
                    x=str(x - box_min_x - 0.5 * lw),
                    y=str(box_max_y - y - 0.5 * lw),
                    height=str(lw),
                    width=str(lw),
                    style=(
                        f"fill:{_rgb_2_hex(r,g,b)};stroke:none;"
                        f"fill-opacity:{a};"
                    ),
                )

        else:

            # Extract knot lines as splines
            knot_lines = []
            for knot in spline.unique_knots[0]:
                knot_lines.append(spline.extract.spline(0, knot))
            for knot in spline.unique_knots[1]:
                knot_lines.append(spline.extract.spline(1, knot))

            # Export knots to paths
            for knot_line in knot_lines:
                # Approximation and export
                spline_copy = _approximate_curve(knot_line, tolerance)
                bezier_elements = spline_copy.extract.beziers()
                path_d = (
                    f"M {bezier_elements[0].cps[0,0]-box_min_x},"
                    f"{box_max_y - bezier_elements[0].cps[0,1]}"
                )
                path_d += " ".join(
                    [
                        (
                            f" C {s.cps[1,0]-box_min_x},"
                            f"{box_max_y - s.cps[1,1]}"
                            f" {s.cps[2,0]-box_min_x},{box_max_y - s.cps[2,1]}"
                            f" {s.cps[3,0]-box_min_x},{box_max_y - s.cps[3,1]}"
                        )
                        for s in bezier_elements
                    ]
                )
                _ET.SubElement(
                    svg_knots,
                    "path",
                    # draw path
                    d=path_d,
                    style=(
                        f"fill:none;stroke:{_rgb_2_hex(r,g,b)};"
                        f"stroke-opacity:{a};"
                        f"stroke-width:{lw};"
                        f"stroke-linecap:{kwargs.get('linecap','round')}"
                    ),
                )


def export(
    fname, *splines, indent=True, box_margins=0.1, tolerance=None, **kwargs
):
    """
    Exports a number of splines into an svg plot

    Graphics are built to resemble vedo plots as close as possible by relying
    on spline show_options

    Example
    -------
    .. code-block:: python

        spline = c.splinepy.Bezier(
            degrees=[2, 2],
            control_points=[
                [1, 0],
                [3, 1],
                [5, 0],  # First row
                [0, 2],
                [2, 3],
                [4, 2],  # Second row
                [1, 4],
                [3, 5],
                [5, 4],  # Third row
            ],
        )
        spline.show_options["c"] = (157, 157, 156)
        spline.show_options["knot_c"] = "red"
        spline.show_options["control_mesh"] = True
        spline.show_options["control_mesh_lw"] = 0.05
        spline.show_options["control_point_c"] = (0, 102, 157)
        spline.show_options["control_mesh_c"] = (0, 102, 157)
        spline.show_options["control_point_ids"] = False
        spline.show_options["control_point_r"] = 0.2
        splinepy.io.svg.export(
            "spline.svg",
            spline,
            box_margins=0.2,
            background_c=None,
        )

    results in

    .. raw:: html

        <object
            data="../../../tests/data/svg_mini_example.svg"
            type="image/svg+xml"
            align="center"
            >
        </object>


    Parameters
    ----------
    fname : string
        name of the output file
    splines : Spline-Type, list
        Splines to be exported
    indent: bool
      Appends white spaces using xml.etree.ElementTree.indent, if possible.
      Makes it more human-readable
    box_margins : int, float
      Offset around the splines (relative units)
    tolerance : float
      Splines that can not be represented exactly by cubic B-Splines (i.e.,
      rationals or high order splines) a tolerance is required for the
      approximation. Default uses an absolute deviation of 1% of the bounding
      box of the given spline
    kwargs:
      Specify more output options:

      - background_c (color, string, rgb, None) : background color
      - linecap ("string"): linecap option for end of lines and connections, see
        https://www.w3.org/TR/SVG2/
        for more information
      - font_family (string): Font for control point ids
      - font_size (float): Size of control point ids
      - text_anchor (string): position of the control point ids relative to
        point
      - text_color (color, string, rgb) : control point id color
      - arrow_width (float): Arrow options see
        https://matplotlib.org/stable/_images/quiver_sizes.svg
      - arrow_head_width (float) : Arrow options see
        https://matplotlib.org/stable/_images/quiver_sizes.svg
      - arrow_head_length (float) : Arrow options see
        https://matplotlib.org/stable/_images/quiver_sizes.svg
      - arrow_headaxis_length (float) : Arrow options see
        https://matplotlib.org/stable/_images/quiver_sizes.svg

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

    # Check if user passed list
    if (len(splines) == 1) and isinstance(splines[0], list):
        splines = splines[0]

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
        ctps_bounds[0, :] = _np.minimum(
            spline.control_point_bounds[0, :], ctps_bounds[0, :]
        )
        ctps_bounds[1, :] = _np.maximum(
            spline.control_point_bounds[1, :], ctps_bounds[1, :]
        )

    # Initialize svg element
    svg_data = _ET.Element("svg")

    # licked it so its ours
    svg_data.insert(
        1,
        _ET.Comment(
            "generated by splinepy https://github.com/tataratat/splinepy"
        ),
    )

    # Set the box margins and svg options globally
    box_size = (
        ctps_bounds[1, 0] - ctps_bounds[0, 0] + 2 * box_margins,
        ctps_bounds[1, 1] - ctps_bounds[0, 1] + 2 * box_margins,
    )
    svg_data.attrib["viewBox"] = f"0 0 {box_size[0]} {box_size[1]}"
    svg_data.attrib["height"] = str(400)
    svg_data.attrib["width"] = str(400 / box_size[1] * box_size[0])
    background_c = kwargs.get("background_c", "white")
    if background_c is not None:
        svg_data.attrib["style"] = f"background-color:{background_c}"
    svg_data.attrib["xmlns"] = "http://www.w3.org/2000/svg"
    svg_data.attrib["xmlns:xlink"] = "http://www.w3.org/1999/xlink"

    # Determine offsets
    box_min_x = ctps_bounds[0, 0] - box_margins
    box_max_y = ctps_bounds[1, 1] + box_margins

    # Write splines in svg file
    for i, spline in enumerate(splines):
        # Copy spline options from keyword arguments
        # (taken from visualize to ensure conformity)
        orig_show_options = None
        if kwargs:
            orig_show_options = spline.show_options
            spline._show_options = spline.__show_option__(spline)
            orig_show_options.copy_valid_options(spline.show_options)
            for key, value in kwargs.items():
                try:
                    spline.show_options[key] = value
                except BaseException:
                    continue

        # Put every spline into a new dedicated group
        spline_group = _ET.SubElement(svg_data, "g", id="spline" + str(i))
        _export_spline(
            spline,
            spline_group,
            box_min_x,
            box_max_y,
            tolerance=tolerance,
            **kwargs,
        )
        _quiver_plot(
            spline=spline,
            svg_spline_element=spline_group,
            box_min_x=box_min_x,
            box_max_y=box_max_y,
            **kwargs,
        )
        _export_control_mesh(
            spline, spline_group, box_min_x, box_max_y, **kwargs
        )

        # set original options back
        if orig_show_options is not None:
            spline._show_options = orig_show_options

    # Dump into file
    if indent and hasattr(_ET, "indent"):
        # Pretty printing xml with indent only exists in version > 3.9
        _ET.indent(svg_data)

    file_content = _ET.tostring(svg_data)
    with open(fname, "wb") as f:
        f.write(file_content)
    pass
