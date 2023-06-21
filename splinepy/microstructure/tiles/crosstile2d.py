import numpy as np

from splinepy.bezier import Bezier
from splinepy.microstructure.tiles.tilebase import TileBase


class CrossTile2D(TileBase):
    def __init__(self):
        """Simple crosstile with linear-quadratic branches and a trilinear
        center spline."""
        self._dim = 2
        self._evaluation_points = np.array(
            [
                [0.0, 0.5],
                [1.0, 0.5],
                [0.5, 0.0],
                [0.5, 1.0],
            ]
        )
        self._parameter_space_dimension = 1

    def closing_tile(
        self,
        parameters=None,
        closure=None,
        boundary_width=0.1,
        filling_height=0.5,
        **kwargs,
    ):
        """Create a closing tile to match with closed surface.

        Parameters
        ----------
        parameters : tuple(np.ndarray)
          radii of fitting cylinder at evaluation points
        closure : int
          parametric dimension that needs to be closed. Positiv values mean
          that minimum parametric dimension is requested. That means,
          i.e. -2 closes the tile at maximum z-coordinate.
          (must currently be either -2 or 2)
        boundary_width : float
          with of the boundary surronding branch
        filling_height : float
          portion of the height that is filled in parametric domain

        Returns
        -------
        list_of_splines : list
        """
        # Check parameters
        if closure is None:
            raise ValueError("No closing direction given")

        if parameters is None:
            self._logd("Tile request is not parametrized, setting default 0.2")
            parameters = tuple([np.ones(6) * 0.2])
        parameters = parameters[0]
        if not (np.all(parameters > 0) and np.all(parameters < 0.5)):
            raise ValueError("Thickness out of range (0, .5)")

        if not (0.0 < float(boundary_width) < 0.5):
            raise ValueError("Boundary Width is out of range")

        if not (0.0 < float(filling_height) < 1.0):
            raise ValueError("Filling must  be in (0,1)")

        # Constant auxiliary values
        inv_boundary_width = 1.0 - boundary_width
        inv_filling_height = 1.0 - filling_height
        ctps_mid_height_top = (1 + filling_height) * 0.5
        ctps_mid_height_bottom = 1.0 - ctps_mid_height_top
        v_one_half = 0.5
        v_one = 1.0
        v_zero = 0.0

        spline_list = []
        if closure == "x_min":
            # Minimum x position
            branch_thickness = parameters[1]

            block0_ctps = np.array(
                [
                    [v_zero, v_zero],
                    [filling_height, v_zero],
                    [v_zero, boundary_width],
                    [filling_height, boundary_width],
                ]
            )

            block1_ctps = np.array(
                [
                    [v_zero, boundary_width],
                    [filling_height, boundary_width],
                    [v_zero, inv_boundary_width],
                    [filling_height, inv_boundary_width],
                ]
            )

            block2_ctps = np.array(
                [
                    [v_zero, inv_boundary_width],
                    [filling_height, inv_boundary_width],
                    [v_zero, v_one],
                    [filling_height, v_one],
                ]
            )

            branch_ctps = np.array(
                [
                    [filling_height, boundary_width],
                    [ctps_mid_height_top, v_one_half - branch_thickness],
                    [v_one, v_one_half - branch_thickness],
                    [filling_height, inv_boundary_width],
                    [ctps_mid_height_top, v_one_half + branch_thickness],
                    [v_one, v_one_half + branch_thickness],
                ]
            )

            spline_list.append(
                Bezier(degrees=[1, 1], control_points=block0_ctps)
            )

            spline_list.append(
                Bezier(degrees=[1, 1], control_points=block1_ctps)
            )

            spline_list.append(
                Bezier(degrees=[1, 1], control_points=block2_ctps)
            )

            spline_list.append(
                Bezier(degrees=[2, 1], control_points=branch_ctps)
            )
            return spline_list
        elif closure == "x_max":
            # Maximum x position
            branch_thickness = parameters[0]

            block0_ctps = np.array(
                [
                    [inv_filling_height, v_zero],
                    [v_one, v_zero],
                    [inv_filling_height, boundary_width],
                    [v_one, boundary_width],
                ]
            )

            block1_ctps = np.array(
                [
                    [inv_filling_height, boundary_width],
                    [v_one, boundary_width],
                    [inv_filling_height, inv_boundary_width],
                    [v_one, inv_boundary_width],
                ]
            )

            block2_ctps = np.array(
                [
                    [inv_filling_height, inv_boundary_width],
                    [v_one, inv_boundary_width],
                    [inv_filling_height, v_one],
                    [v_one, v_one],
                ]
            )

            branch_ctps = np.array(
                [
                    [0, v_one_half - branch_thickness],
                    [ctps_mid_height_bottom, v_one_half - branch_thickness],
                    [inv_filling_height, boundary_width],
                    [v_zero, v_one_half + branch_thickness],
                    [ctps_mid_height_bottom, v_one_half + branch_thickness],
                    [inv_filling_height, inv_boundary_width],
                ]
            )

            spline_list.append(
                Bezier(degrees=[1, 1], control_points=block0_ctps)
            )

            spline_list.append(
                Bezier(degrees=[1, 1], control_points=block1_ctps)
            )

            spline_list.append(
                Bezier(degrees=[1, 1], control_points=block2_ctps)
            )

            spline_list.append(
                Bezier(degrees=[2, 1], control_points=branch_ctps)
            )
            return spline_list
        elif closure == "y_min":
            # Minimum y position
            branch_thickness = parameters[3]

            block0_ctps = np.array(
                [
                    [v_zero, v_zero],
                    [boundary_width, v_zero],
                    [v_zero, filling_height],
                    [boundary_width, filling_height],
                ]
            )

            block1_ctps = np.array(
                [
                    [boundary_width, v_zero],
                    [inv_boundary_width, v_zero],
                    [boundary_width, filling_height],
                    [inv_boundary_width, filling_height],
                ]
            )

            block2_ctps = np.array(
                [
                    [inv_boundary_width, v_zero],
                    [v_one, v_zero],
                    [inv_boundary_width, filling_height],
                    [v_one, filling_height],
                ]
            )

            branch_ctps = np.array(
                [
                    [boundary_width, filling_height],
                    [inv_boundary_width, filling_height],
                    [v_one_half - branch_thickness, ctps_mid_height_top],
                    [v_one_half + branch_thickness, ctps_mid_height_top],
                    [v_one_half - branch_thickness, v_one],
                    [v_one_half + branch_thickness, v_one],
                ]
            )

            spline_list.append(
                Bezier(degrees=[1, 1], control_points=block0_ctps)
            )

            spline_list.append(
                Bezier(degrees=[1, 1], control_points=block1_ctps)
            )

            spline_list.append(
                Bezier(degrees=[1, 1], control_points=block2_ctps)
            )

            spline_list.append(
                Bezier(degrees=[1, 2], control_points=branch_ctps)
            )
            return spline_list
        elif closure == "y_max":
            # Maximum y position
            branch_thickness = parameters[2]

            block0_ctps = np.array(
                [
                    [v_zero, inv_filling_height],
                    [boundary_width, inv_filling_height],
                    [v_zero, v_one],
                    [boundary_width, v_one],
                ]
            )

            block1_ctps = np.array(
                [
                    [boundary_width, inv_filling_height],
                    [inv_boundary_width, inv_filling_height],
                    [boundary_width, v_one],
                    [inv_boundary_width, v_one],
                ]
            )

            block2_ctps = np.array(
                [
                    [inv_boundary_width, inv_filling_height],
                    [v_one, inv_filling_height],
                    [inv_boundary_width, v_one],
                    [v_one, v_one],
                ]
            )

            branch_ctps = np.array(
                [
                    [v_one_half - branch_thickness, v_zero],
                    [v_one_half + branch_thickness, v_zero],
                    [v_one_half - branch_thickness, ctps_mid_height_bottom],
                    [v_one_half + branch_thickness, ctps_mid_height_bottom],
                    [boundary_width, inv_filling_height],
                    [inv_boundary_width, inv_filling_height],
                ]
            )

            spline_list.append(
                Bezier(degrees=[1, 1], control_points=block0_ctps)
            )

            spline_list.append(
                Bezier(degrees=[1, 1], control_points=block1_ctps)
            )

            spline_list.append(
                Bezier(degrees=[1, 1], control_points=block2_ctps)
            )

            spline_list.append(
                Bezier(degrees=[1, 2], control_points=branch_ctps)
            )
            return spline_list
        else:
            raise NotImplementedError(
                "Requested closing dimension is not supported"
            )

    def create_tile(self, parameters=None, center_expansion=1.0, **kwargs):
        """Create a microtile based on the parameters that describe the branch
        thicknesses.

        Thickness parameters are used to describe the inner radius of the
        outward facing branches

        Parameters
        ----------
        parameters : tuple(np.array)
            only first entry is used, defines the internal radii of the
            branches
        center_expansion : float
            thickness of center is expanded by a factor
        Returns
        -------
        microtile_list : list(splines)
        """

        if not isinstance(center_expansion, float):
            raise ValueError("Invalid Type")
        if not ((center_expansion > 0.5) and (center_expansion < 1.5)):
            raise ValueError("Center Expansion must be in (.5,1.5)")
        max_radius = min(0.5, (0.5 / center_expansion))
        # set to default if nothing is given
        if parameters is None:
            self._logd("Setting branch thickness to default 0.2")
            parameters = tuple([np.ones(4) * 0.2])
        [x_min_r, x_max_r, y_min_r, y_max_r] = parameters[0].tolist()
        for radius in [x_min_r, x_max_r, y_min_r, y_max_r]:
            if not isinstance(radius, float):
                raise ValueError("Invalid type")
            if not (radius > 0 and radius < max_radius):
                raise ValueError(
                    f"Radii must be in (0,{max_radius}) for "
                    f"center_expansion {center_expansion}"
                )

        # center radius
        center_r = (
            (x_min_r + x_max_r + y_min_r + y_max_r) / 4.0 * center_expansion
        )
        hd_center = 0.5 * (0.5 + center_r)
        one_half = 0.5

        # Init return value
        spline_list = []

        # Create the center-tile
        center_points = np.array(
            [
                [-center_r, -center_r],
                [center_r, -center_r],
                [-center_r, center_r],
                [center_r, center_r],
            ]
        ) + np.array([one_half, one_half])

        y_min_ctps = np.array(
            [
                [-y_min_r, -one_half],
                [y_min_r, -one_half],
                [-y_min_r, -hd_center],
                [y_min_r, -hd_center],
                [-center_r, -center_r],
                [center_r, -center_r],
            ]
        ) + np.array([one_half, one_half])

        y_max_ctps = np.array(
            [
                [-center_r, center_r],
                [center_r, center_r],
                [-y_max_r, hd_center],
                [y_max_r, hd_center],
                [-y_max_r, one_half],
                [y_max_r, one_half],
            ]
        ) + np.array([one_half, one_half])

        x_min_ctps = np.array(
            [
                [-one_half, -x_min_r],
                [-hd_center, -x_min_r],
                [-center_r, -center_r],
                [-one_half, x_min_r],
                [-hd_center, x_min_r],
                [-center_r, center_r],
            ]
        ) + np.array([one_half, one_half])

        x_max_ctps = np.array(
            [
                [center_r, -center_r],
                [hd_center, -x_max_r],
                [one_half, -x_max_r],
                [center_r, center_r],
                [hd_center, x_max_r],
                [one_half, x_max_r],
            ]
        ) + np.array([one_half, one_half])

        spline_list.append(
            Bezier(degrees=[1, 1], control_points=center_points)
        )

        spline_list.append(Bezier(degrees=[2, 1], control_points=x_min_ctps))

        spline_list.append(Bezier(degrees=[2, 1], control_points=x_max_ctps))

        spline_list.append(Bezier(degrees=[1, 2], control_points=y_min_ctps))

        spline_list.append(Bezier(degrees=[1, 2], control_points=y_max_ctps))

        return spline_list
