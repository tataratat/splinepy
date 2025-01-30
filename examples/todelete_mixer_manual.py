import numpy as np

import splinepy as sp

WDIR = ""
EPS = 1e-7
KNOT_INSERTIONS = 0
TILING = [2, 2, 4]
BOX_LENGTH = 2


def href(patch):
    for _ in range(KNOT_INSERTIONS):
        patch.uniform_refine([0, 1, 2])


class SulzerInverse:
    def __init__(self, macro_spline, tiling, thickness):
        self._macro_spline = macro_spline
        self._tiling = tiling
        self._bar_thickness = thickness

        self._microstructure = None
        self._material_microstructure = None
        self._mixer_length = (
            macro_spline.cps[:, 2].max() - macro_spline.cps[:2].min()
        )

        self.twisted_layers = None

    def set_macro_spline(self, spline):
        self._macro_spline = spline

    def set_tiling(self, tiling):
        self._tiling = tiling

    # def create_unit_tile(self):
    #     v_zero = 0.0
    #     v_one_half = 0.5
    #     v_one = 1.0

    #     # Auxiliary values
    #     triangle_height = 0.5 * (v_one - self._bar_thickness - self._bar_thickness)

    #     spline_list = []

    #     # Define control points of crossbars
    #     front_bottom = np.array(
    #         [
    #             [v_zero, v_zero, self._bar_thickness],
    #             [v_zero, v_zero, self._bar_thickness + triangle_height],
    #             [
    #                 v_zero,
    #                 triangle_height,
    #                 self._bar_thickness + triangle_height,
    #             ],
    #             [v_zero, v_zero, v_one - self._bar_thickness],
    #             [v_one_half, v_zero, self._bar_thickness],
    #             [v_one_half, v_zero, self._bar_thickness + triangle_height],
    #             [
    #                 v_one_half,
    #                 triangle_height,
    #                 self._bar_thickness + triangle_height,
    #             ],
    #             [v_one_half, v_zero, v_one - self._bar_thickness],
    #         ]
    #     )

    #     back_bottom = np.array(
    #         [
    #             [v_one_half, v_zero, self._bar_thickness],
    #             [v_one_half, v_zero, self._bar_thickness + triangle_height],
    #             [
    #                 v_one_half,
    #                 triangle_height,
    #                 self._bar_thickness + triangle_height,
    #             ],
    #             [v_one_half, v_zero, v_one - self._bar_thickness],
    #             [v_one, v_zero, self._bar_thickness],
    #             [v_one, v_zero, self._bar_thickness + triangle_height],
    #             [
    #                 v_one,
    #                 triangle_height,
    #                 self._bar_thickness + triangle_height,
    #             ],
    #             [v_one, v_zero, v_one - self._bar_thickness],
    #         ]
    #     )

    #     front_bottom_bar = np.array(
    #         [
    #             [v_zero, v_zero, v_one - self._bar_thickness],
    #             [v_zero, self._bar_thickness, v_one],
    #             [
    #                 v_zero,
    #                 triangle_height,
    #                 self._bar_thickness + triangle_height,
    #             ],
    #             [
    #                 v_zero,
    #                 self._bar_thickness + triangle_height,
    #                 v_one - triangle_height,
    #             ],
    #             [v_one_half, v_zero, v_one - self._bar_thickness],
    #             [v_one_half, self._bar_thickness, v_one],
    #             [
    #                 v_one_half,
    #                 triangle_height,
    #                 self._bar_thickness + triangle_height,
    #             ],
    #             [
    #                 v_one_half,
    #                 self._bar_thickness + triangle_height,
    #                 v_one - triangle_height,
    #             ],
    #         ]
    #     )

    #     front_right = np.array(
    #         [
    #             [v_zero, self._bar_thickness, v_one],
    #             [v_zero, self._bar_thickness + triangle_height, v_one],
    #             [
    #                 v_zero,
    #                 self._bar_thickness + triangle_height,
    #                 v_one - triangle_height,
    #             ],
    #             [v_zero, v_one - self._bar_thickness, v_one],
    #             [v_one_half, self._bar_thickness, v_one],
    #             [v_one_half, self._bar_thickness + triangle_height, v_one],
    #             [
    #                 v_one_half,
    #                 self._bar_thickness + triangle_height,
    #                 v_one - triangle_height,
    #             ],
    #             [v_one_half, v_one - self._bar_thickness, v_one],
    #         ]
    #     )

    #     back_right = np.array(
    #         [
    #             [v_one_half, self._bar_thickness, v_one],
    #             [v_one_half, self._bar_thickness + triangle_height, v_one],
    #             [
    #                 v_one_half,
    #                 self._bar_thickness + triangle_height,
    #                 v_one - triangle_height,
    #             ],
    #             [v_one_half, v_one - self._bar_thickness, v_one],
    #             [v_one, self._bar_thickness, v_one],
    #             [v_one, self._bar_thickness + triangle_height, v_one],
    #             [
    #                 v_one,
    #                 self._bar_thickness + triangle_height,
    #                 v_one - triangle_height,
    #             ],
    #             [v_one, v_one - self._bar_thickness, v_one],
    #         ]
    #     )

    #     front_top = np.array(
    #         [
    #             [
    #                 v_zero,
    #                 v_one - triangle_height,
    #                 self._bar_thickness + triangle_height,
    #             ],
    #             [v_zero, v_one, v_one - self._bar_thickness],
    #             [v_zero, v_one, self._bar_thickness],
    #             [v_zero, v_one, self._bar_thickness + triangle_height],
    #             [
    #                 v_one_half,
    #                 v_one - triangle_height,
    #                 self._bar_thickness + triangle_height,
    #             ],
    #             [v_one_half, v_one, v_one - self._bar_thickness],
    #             [v_one_half, v_one, self._bar_thickness],
    #             [v_one_half, v_one, self._bar_thickness + triangle_height],
    #         ]
    #     )

    #     back_top = np.array(
    #         [
    #             [
    #                 v_one_half,
    #                 v_one - triangle_height,
    #                 self._bar_thickness + triangle_height,
    #             ],
    #             [v_one_half, v_one, v_one - self._bar_thickness],
    #             [v_one_half, v_one, self._bar_thickness],
    #             [v_one_half, v_one, self._bar_thickness + triangle_height],
    #             [
    #                 v_one,
    #                 v_one - triangle_height,
    #                 self._bar_thickness + triangle_height,
    #             ],
    #             [v_one, v_one, v_one - self._bar_thickness],
    #             [v_one, v_one, self._bar_thickness],
    #             [v_one, v_one, self._bar_thickness + triangle_height],
    #         ]
    #     )

    #     front_top_bar = np.array(
    #         [
    #             [
    #                 v_zero,
    #                 self._bar_thickness + triangle_height,
    #                 triangle_height,
    #             ],
    #             [
    #                 v_zero,
    #                 v_one - triangle_height,
    #                 self._bar_thickness + triangle_height,
    #             ],
    #             [v_zero, v_one - self._bar_thickness, v_zero],
    #             [v_zero, v_one, self._bar_thickness],
    #             [
    #                 v_one_half,
    #                 self._bar_thickness + triangle_height,
    #                 triangle_height,
    #             ],
    #             [
    #                 v_one_half,
    #                 v_one - triangle_height,
    #                 self._bar_thickness + triangle_height,
    #             ],
    #             [v_one_half, v_one - self._bar_thickness, v_zero],
    #             [v_one_half, v_one, self._bar_thickness],
    #         ]
    #     )

    #     front_left = np.array(
    #         [
    #             [v_zero, self._bar_thickness, v_zero],
    #             [
    #                 v_zero,
    #                 self._bar_thickness + triangle_height,
    #                 triangle_height,
    #             ],
    #             [v_zero, self._bar_thickness + triangle_height, v_zero],
    #             [v_zero, v_one - self._bar_thickness, v_zero],
    #             [v_one_half, self._bar_thickness, v_zero],
    #             [
    #                 v_one_half,
    #                 self._bar_thickness + triangle_height,
    #                 triangle_height,
    #             ],
    #             [v_one_half, self._bar_thickness + triangle_height, v_zero],
    #             [v_one_half, v_one - self._bar_thickness, v_zero],
    #         ]
    #     )

    #     back_left = np.array(
    #         [
    #             [v_one_half, self._bar_thickness, v_zero],
    #             [
    #                 v_one_half,
    #                 self._bar_thickness + triangle_height,
    #                 triangle_height,
    #             ],
    #             [v_one_half, self._bar_thickness + triangle_height, v_zero],
    #             [v_one_half, v_one - self._bar_thickness, v_zero],
    #             [v_one, self._bar_thickness, v_zero],
    #             [
    #                 v_one,
    #                 self._bar_thickness + triangle_height,
    #                 triangle_height,
    #             ],
    #             [v_one, self._bar_thickness + triangle_height, v_zero],
    #             [v_one, v_one - self._bar_thickness, v_zero],
    #         ]
    #     )

    #     back_top_bar = np.array(
    #         [
    #             [
    #                 v_one_half,
    #                 v_one - triangle_height,
    #                 self._bar_thickness + triangle_height,
    #             ],
    #             [
    #                 v_one_half,
    #                 self._bar_thickness + triangle_height,
    #                 v_one - triangle_height,
    #             ],
    #             [v_one_half, v_one, v_one - self._bar_thickness],
    #             [v_one_half, v_one - self._bar_thickness, v_one],
    #             [
    #                 v_one,
    #                 v_one - triangle_height,
    #                 self._bar_thickness + triangle_height,
    #             ],
    #             [
    #                 v_one,
    #                 self._bar_thickness + triangle_height,
    #                 v_one - triangle_height,
    #             ],
    #             [v_one, v_one, v_one - self._bar_thickness],
    #             [v_one, v_one - self._bar_thickness, v_one],
    #         ]
    #     )

    #     back_bottom_bar = np.array(
    #         [
    #             [v_one_half, self._bar_thickness, v_zero],
    #             [v_one_half, v_zero, self._bar_thickness],
    #             [
    #                 v_one_half,
    #                 self._bar_thickness + triangle_height,
    #                 triangle_height,
    #             ],
    #             [
    #                 v_one_half,
    #                 triangle_height,
    #                 self._bar_thickness + triangle_height,
    #             ],
    #             [v_one, self._bar_thickness, v_zero],
    #             [v_one, v_zero, self._bar_thickness],
    #             [
    #                 v_one,
    #                 self._bar_thickness + triangle_height,
    #                 triangle_height,
    #             ],
    #             [
    #                 v_one,
    #                 triangle_height,
    #                 self._bar_thickness + triangle_height,
    #             ],
    #         ]
    #     )

    #     # Append all the patches into list
    #     for control_points in [
    #         front_bottom,
    #         front_bottom_bar,
    #         front_right,
    #         front_top,
    #         front_top_bar,
    #         front_left,
    #         back_bottom,
    #         back_right,
    #         back_top_bar,
    #         back_top,
    #         back_left,
    #         back_bottom_bar,
    #     ]:
    #         spline_list.append(
    #             sp.BSpline(
    #                 degrees=[1, 1, 1],
    #                 control_points=control_points,
    #                 knot_vectors=[[0.0, 0.0, 1.0, 1.0]]*3
    #             )
    #         )

    #     return spline_list

    def create_unit_tile(self):
        v_zero = 0.0
        v_one_half = 0.5
        v_one = 1.0

        # Auxiliary values
        triangle_height = 0.5 * (
            v_one - self._bar_thickness - self._bar_thickness
        )

        spline_list = []

        # Helper array for defining x-component of control points
        e = np.ones((4, 1))

        # Define control points of crossbars
        xs_front = np.vstack((v_zero * e, v_one_half * e))
        xs_back = np.vstack((v_one_half * e, v_one * e))

        # Define the y- and z-coordinates of the triangle on the bottom
        triangle_bottom_sw = np.array(
            [
                [v_zero, self._bar_thickness],
                [v_zero, v_one_half],
                [triangle_height * 0.5, v_one_half - triangle_height * 0.5],
                [triangle_height / 3, v_one_half],
            ]
        )

        triangle_bottom_se = np.array(
            [
                [v_zero, v_one_half],
                [v_zero, v_one - self._bar_thickness],
                [triangle_height / 3, v_one_half],
                [triangle_height * 0.5, v_one_half + triangle_height * 0.5],
            ]
        )
        triangle_bottom_n = np.array(
            [
                [triangle_height * 0.5, v_one_half - triangle_height * 0.5],
                [triangle_height / 3, v_one_half],
                [triangle_height, v_one_half],
                [triangle_height * 0.5, v_one_half + triangle_height * 0.5],
            ]
        )

        triangle_right_s = np.array(
            [
                [
                    v_one_half - triangle_height * 0.5,
                    v_one - triangle_height * 0.5,
                ],
                [self._bar_thickness, v_one],
                [v_one_half, v_one - triangle_height / 3],
                [v_one_half, v_one],
            ]
        )

        triangle_right_n = np.array(
            [
                [v_one_half, v_one - triangle_height / 3],
                [v_one_half, v_one],
                [
                    v_one_half + triangle_height * 0.5,
                    v_one - triangle_height * 0.5,
                ],
                [v_one - self._bar_thickness, v_one],
            ]
        )

        triangle_right_w = np.array(
            [
                [v_one_half, v_one - triangle_height],
                [
                    v_one_half - triangle_height * 0.5,
                    v_one - triangle_height * 0.5,
                ],
                [
                    v_one_half + triangle_height * 0.5,
                    v_one - triangle_height * 0.5,
                ],
                [v_one_half, v_one - triangle_height / 3],
            ]
        )

        triangle_top_nw, triangle_top_ne, triangle_top_s = (
            np.hstack(
                (v_one - cps[:, 0].reshape(-1, 1), cps[:, 1].reshape(-1, 1))
            )
            for cps in [
                triangle_bottom_sw,
                triangle_bottom_se,
                triangle_bottom_n,
            ]
        )

        triangle_left_s, triangle_left_n, triangle_left_e = (
            np.hstack(
                (cps[:, 0].reshape(-1, 1), v_one - cps[:, 1].reshape(-1, 1))
            )
            for cps in [triangle_right_s, triangle_right_n, triangle_right_w]
        )

        # Create patches which have front as well as back
        for yz_cps in [
            triangle_bottom_sw,
            triangle_bottom_se,
            triangle_bottom_n,
            triangle_top_nw,
            triangle_top_ne,
            triangle_top_s,
            triangle_right_s,
            triangle_right_n,
            triangle_right_w,
            triangle_left_s,
            triangle_left_n,
            triangle_left_e,
        ]:
            for xs in [xs_front, xs_back]:
                control_points = np.hstack((xs, np.tile(yz_cps, [2, 1])))
                spline_list.append(
                    sp.BSpline(
                        degrees=[1, 1, 1],
                        control_points=control_points,
                        knot_vectors=[[0.0, 0.0, 1.0, 1.0]] * 3,
                    )
                )

        # Create the bars
        front_bottom_bar_s = np.array(
            [
                [
                    triangle_height * 0.5,
                    v_one - self._bar_thickness - triangle_height * 0.5,
                ],
                [v_zero, v_one - self._bar_thickness],
                [
                    self._bar_thickness + triangle_height * 0.5,
                    v_one - triangle_height * 0.5,
                ],
                [self._bar_thickness, v_one],
            ]
        )

        front_bottom_bar_n = np.array(
            [
                [triangle_height, v_one_half],
                [
                    triangle_height * 0.5,
                    v_one - self._bar_thickness - triangle_height * 0.5,
                ],
                [v_one_half, v_one - triangle_height],
                [
                    self._bar_thickness + triangle_height * 0.5,
                    v_one - triangle_height * 0.5,
                ],
            ]
        )

        # Helper functions for control point translation
        def identity_function(cps):
            return cps

        def opposite_corner_function(cps):
            return v_one - cps

        # Append bar patches to spline list
        for yz_cps in [front_bottom_bar_s, front_bottom_bar_n]:
            # 1: create front bottom, 2: create front top
            for translation_function in [
                identity_function,
                opposite_corner_function,
            ]:
                yz_cps_translated = translation_function(yz_cps)
                control_points = np.hstack(
                    (xs_front, np.tile(yz_cps_translated, [2, 1]))
                )
                spline_list.append(
                    sp.BSpline(
                        degrees=[1, 1, 1],
                        control_points=control_points,
                        knot_vectors=[[0.0, 0.0, 1.0, 1.0]] * 3,
                    )
                )

        # Define bar control points
        back_bottom_bar_s = np.array(
            [
                [v_zero, self._bar_thickness],
                [triangle_height * 0.5, v_one_half - triangle_height * 0.5],
                [self._bar_thickness, v_zero],
                [v_one_half - triangle_height * 0.5, triangle_height * 0.5],
            ]
        )

        back_bottom_bar_n = np.array(
            [
                [triangle_height * 0.5, v_one_half - triangle_height * 0.5],
                [triangle_height, v_one_half],
                [v_one_half - triangle_height * 0.5, triangle_height * 0.5],
                [v_one_half, triangle_height],
            ]
        )

        # Add back bar patches to spline list
        for yz_cps in [back_bottom_bar_s, back_bottom_bar_n]:
            # 1: create back bottom, 2: create back top bar
            for translation_function in [
                identity_function,
                opposite_corner_function,
            ]:
                yz_cps_translated = translation_function(yz_cps)
                control_points = np.hstack(
                    (xs_back, np.tile(yz_cps_translated, [2, 1]))
                )
                spline_list.append(
                    sp.BSpline(
                        degrees=[1, 1, 1],
                        control_points=control_points,
                        knot_vectors=[[0.0, 0.0, 1.0, 1.0]] * 3,
                    )
                )

        return spline_list

    def compute_start_points(self):
        min_positions = [self._macro_spline.cps[:, i].min() for i in range(3)]
        max_positions = [self._macro_spline.cps[:, i].max() for i in range(3)]

        # Determine base point for each tile
        x_starts = np.linspace(
            min_positions[0], max_positions[0], self._tiling[0] + 1
        )[:-1]
        y_starts = np.linspace(
            min_positions[1], max_positions[1], self._tiling[1] + 1
        )[:-1]
        z_starts = np.linspace(
            min_positions[2], max_positions[2], self._tiling[2] + 1
        )[:-1]

        # Determine base points for each tile
        base_points = [
            array.ravel(order="F")
            for array in np.meshgrid(y_starts, x_starts, z_starts)
        ]
        self.base_points = np.vstack(base_points).T

        # Determine tiles dimensions
        self.dd = np.diff(
            np.vstack((min_positions, max_positions)), axis=0
        ).ravel() / np.array(self._tiling)

    def create_microstructure(self):
        self.compute_start_points()
        # Create base tile
        base_tile = self.create_unit_tile()
        for patch in base_tile:
            patch.cps *= self.dd

        all_tiles = []
        for base_point in self.base_points:
            tile_patches_base = base_tile.copy()
            for patch in tile_patches_base:
                new_patch = patch.copy()
                new_patch.cps += base_point
                href(new_patch)
                all_tiles.append(new_patch)

        self._microstructure = sp.Multipatch(splines=all_tiles)
        self._microstructure.determine_interfaces()

    def create_twisted_to_nontwisted_linkage(self, gap_length, cell_length):
        zero = 0.0
        half = 0.5 * cell_length
        one = cell_length

        spline_list = []

        xmin_ymin = np.array(
            [
                [gap_length, zero, zero],
                [half, zero, zero],
                [gap_length, half, zero],
                [half, half, zero],
                [zero, gap_length, gap_length],
                [half, gap_length, gap_length],
                [zero, half, gap_length],
                [half, half, gap_length],
            ]
        )

        xmax_ymin = np.array(
            [
                [half, zero, zero],
                [one - gap_length, zero, zero],
                [half, half, zero],
                [one - gap_length, half, zero],
                [half, gap_length, gap_length],
                [one, gap_length, gap_length],
                [half, half, gap_length],
                [one, half, gap_length],
            ]
        )

        xmin_ymax = np.array(
            [
                [gap_length, half, zero],
                [half, half, zero],
                [gap_length, one, zero],
                [half, one, zero],
                [zero, half, gap_length],
                [half, half, gap_length],
                [zero, one - gap_length, gap_length],
                [half, one - gap_length, gap_length],
            ]
        )

        xmax_ymax = np.array(
            [
                [half, half, zero],
                [one - gap_length, half, zero],
                [half, one, zero],
                [one - gap_length, one, zero],
                [half, half, gap_length],
                [one, half, gap_length],
                [half, one - gap_length, gap_length],
                [one, one - gap_length, gap_length],
            ]
        )

        for patch_cps in [xmin_ymin, xmax_ymin, xmin_ymax, xmax_ymax]:
            patch = sp.Bezier(
                degrees=[1, 1, 1], control_points=patch_cps
            ).bspline
            spline_list.append(patch)

        return spline_list

    def create_nontwisted_to_twisted_linkage(self, gap_length, cell_length):
        spline_list = self.create_twisted_to_nontwisted_linkage(
            gap_length, cell_length
        )

        for patch in spline_list:
            cps = np.array(patch.cps)
            patch.cps[:4, :2] = cps[4:, :2]
            patch.cps[4:, :2] = cps[:4, :2]

        return spline_list

    def create_twisted_microstructure(self, twisted_layers=None):
        assert isinstance(
            twisted_layers, list
        ), "twisted_layers must be a list"
        twisted_layers.sort()
        for layer in twisted_layers:
            assert layer < self._tiling[2], f"Layer {layer} exceeds the tiling"

        self.twisted_layers = twisted_layers

        self.compute_start_points()

        # Create single unit tile and size it for actual size in microstructure
        base_tile = self.create_unit_tile()
        for patch in base_tile:
            patch.cps *= self.dd

        # Twist tile
        twisted_tile = []
        for patch in base_tile:
            twisted_patch = patch.copy()
            cps = np.array(patch.cps)
            twisted_patch.cps[:, 0] = cps[:, 1]
            twisted_patch.cps[:, 1] = cps[:, 0]
            twisted_tile.append(twisted_patch)

        min_positions = [self._macro_spline.cps[:, i].min() for i in range(3)]
        max_positions = [self._macro_spline.cps[:, i].max() for i in range(3)]

        # Determine base point for each tile
        x_starts = np.linspace(
            min_positions[0], max_positions[0], self._tiling[0] + 1
        )[:-1]
        y_starts = np.linspace(
            min_positions[1], max_positions[1], self._tiling[1] + 1
        )[:-1]

        base_points = [
            array.ravel(order="F") for array in np.meshgrid(x_starts, y_starts)
        ]
        base_points.reverse()
        base_points = np.vstack(base_points).T
        base_points = np.hstack(
            (base_points, np.zeros((base_points.shape[0], 1)))
        )

        # Create one layer of normal and one of twisted tiles, not placing it yet
        layer_tiles = []
        twisted_layer_tiles = []
        for basis, tile_patches in zip(
            [base_tile, twisted_tile], [layer_tiles, twisted_layer_tiles]
        ):
            for base_point in base_points:
                tile_patches_base = basis.copy()
                for patch in tile_patches_base:
                    new_patch = patch.copy()
                    new_patch.cps += base_point
                    href(new_patch)
                    tile_patches.append(new_patch)

        all_tiles = []
        z_offset = 0.0
        dz = (max_positions[2] - min_positions[2]) / self._tiling[2]
        bar_thickness = self._bar_thickness / self._tiling[1]

        # Go through each layer in z-direction
        for i in range(self._tiling[2]):
            # Get layer
            if i in twisted_layers:
                tile_patches_base = twisted_layer_tiles.copy()
            else:
                tile_patches_base = layer_tiles.copy()
            # Place respective layer in respective z-position
            for patch in tile_patches_base:
                new_patch = patch.copy()
                new_patch.cps[:, 2] += z_offset
                all_tiles.append(new_patch)

            # Add to layer offset in z-direction
            z_offset += dz

            # Create transitions if necessary
            twist_to_nontwist_transition = (
                i in twisted_layers and i + 1 not in twisted_layers
            )
            nontwist_to_twist_transition = (
                i not in twisted_layers and i + 1 in twisted_layers
            )
            # Decide on which transition to take
            if twist_to_nontwist_transition:
                linkage_patches = self.create_twisted_to_nontwisted_linkage(
                    gap_length=bar_thickness, cell_length=dz
                )
            elif nontwist_to_twist_transition:
                linkage_patches = self.create_nontwisted_to_twisted_linkage(
                    gap_length=bar_thickness, cell_length=dz
                )
            # Add transition
            if twist_to_nontwist_transition or nontwist_to_twist_transition:
                for patch in linkage_patches:
                    patch.cps[:, 2] += z_offset
                for base_point in base_points:
                    patch_base = linkage_patches.copy()
                    for patch in patch_base:
                        new_patch = patch.copy()
                        new_patch.cps += base_point
                        href(new_patch)
                        all_tiles.append(new_patch)

                z_offset += bar_thickness
                self._mixer_length += bar_thickness

        self._microstructure = sp.Multipatch(splines=all_tiles)
        self._microstructure.determine_interfaces()

    def create_fore_and_afterrun(self, forelength=1.0):
        bar_thickness = self._bar_thickness / self._tiling[1]
        assert (
            self._bar_thickness < forelength / 2
        ), "Forelength must be bigger than 2*bar thickness"

        dx = 1 / (2 * self._tiling[0])
        dy = 1 / (2 * self._tiling[1])
        x_tiling_points = np.linspace(0, 1, self._tiling[0] * 2 + 1)
        y_tiling_points = np.linspace(0, 1, self._tiling[1] * 2 + 1)

        connection_box = sp.helpme.create.box(dx, dy, bar_thickness).bspline
        middle_connection_box = sp.helpme.create.box(
            dx, dy, forelength / 2 - bar_thickness
        ).bspline
        farthest_box = sp.helpme.create.box(dx, dy, forelength / 2).bspline

        forerun_patches = []
        afterrun_patches = []
        for dx_tile in x_tiling_points[:-1]:
            for iy, dy_tile in enumerate(y_tiling_points[:-1]):
                # Connection box
                new_connection = connection_box.copy()
                new_connection.cps += np.array(
                    [[dx_tile, dy_tile, -bar_thickness]]
                )
                if iy % 2 == 0:
                    new_connection.cps[(4, 5), 1] += bar_thickness
                else:
                    new_connection.cps[(6, 7), 1] -= bar_thickness
                href(new_connection)
                forerun_patches.append(new_connection)
                new_connection_end = connection_box.copy()
                new_connection_end.cps += np.array(
                    [[dx_tile, dy_tile, self._mixer_length]]
                )
                if iy % 2 == 0:
                    new_connection_end.cps[(0, 1), 1] += bar_thickness
                else:
                    new_connection_end.cps[(2, 3), 1] -= bar_thickness
                href(new_connection_end)
                afterrun_patches.append(new_connection_end)
                # Middle connection
                new_middle = middle_connection_box.copy()
                new_middle.cps += np.array(
                    [[dx_tile, dy_tile, -forelength / 2]]
                )
                href(new_middle)
                forerun_patches.append(new_middle)
                new_middle_end = middle_connection_box.copy()
                new_middle_end.cps += np.array(
                    [[dx_tile, dy_tile, self._mixer_length + bar_thickness]]
                )
                href(new_middle_end)
                afterrun_patches.append(new_middle_end)
                # Farthest box
                new_farthest = farthest_box.copy()
                new_farthest.cps += np.array([[dx_tile, dy_tile, -forelength]])
                href(new_farthest)
                forerun_patches.append(new_farthest)
                new_farthest_end = farthest_box.copy()
                new_farthest_end.cps += np.array(
                    [[dx_tile, dy_tile, self._mixer_length + forelength / 2]]
                )
                href(new_farthest_end)
                afterrun_patches.append(new_farthest_end)

        self._microstructure = sp.Multipatch(
            splines=self._microstructure.patches
            + forerun_patches
            + afterrun_patches
        )

        self._microstructure.determine_interfaces()

        def inlet_identifier(points):
            return points[:, 2] < -forelength + EPS

        def outlet_identifier(points):
            return points[:, 2] > self._mixer_length + forelength - EPS

        self._microstructure.boundary_from_function(
            inlet_identifier, boundary_id=2
        )
        self._microstructure.boundary_from_function(
            outlet_identifier, boundary_id=3
        )

    def show(self, show_material=False):
        object_to_show = (
            self._material_microstructure
            if show_material
            else self._microstructure
        )
        sp.show(object_to_show, control_points=False, use_saved=True)

    def export(self):
        import sys

        sys.path.insert(
            0, "/Users/markriegler/Documents/ttrrial/splinepy/splinepy/io/"
        )
        from gismo import AdditionalBlocks
        from gismo import export as gexport

        # Define extra blocks
        blocks = AdditionalBlocks()

        blocks.add_assembly_options(block_id=4, comment=" Assembly options ")

        blocks.add_boundary_conditions(
            block_id=2,
            dim=3,
            function_list=[
                ("0.0", "0.0", "0.0"),
                ("0.0", "0.0", "x*(1-x)*y*(1-y)"),
            ],
            bc_list=[("BID1", "Dirichlet", 0), ("BID2", "Dirichlet", 1)],
            unknown_id=1,
            comment=" Velocity boundary conditions ",
        )

        blocks.add_boundary_conditions(
            block_id=3,
            dim=3,
            function_list=["0"],
            bc_list=[("BID3", "Dirichlet", 0)],
            unknown_id=0,
            comment=" Pressure boundary conditions ",
        )

        filename = WDIR + f"mixer_3D_{self._tiling[0]}x{self._tiling[1]}b.xml"
        gexport(
            fname=filename,
            multipatch=self._microstructure,
            additional_blocks=blocks.to_list(),
        )
        print(f"Export to {filename}")

    def export_material_microstructure(self):
        filename = (
            WDIR + f"material_mixer_3D_{self._tiling[0]}x{self._tiling[1]}.xml"
        )
        sp.io.gismo.export(
            fname=filename, multipatch=self._material_microstructure
        )
        print(f"Export material MS to {filename}")

    def create_material_unit_tile(self):
        """Unit tile for the actual static mixer geometry. Not simulatable right now"""
        zero = 0.0
        half = 0.5
        one = 1.0
        thickness = self._bar_thickness

        spline_list = []

        # Define control points
        front_top = np.array(
            [
                [zero, zero, zero],
                [zero, one, one],
                [zero, thickness, zero],
                [zero, one, one - thickness],
                [half, zero, zero],
                [half, one, one],
                [half, thickness, zero],
                [half, one, one - thickness],
            ]
        )

        front_bottom = np.array(
            [
                [zero, zero, thickness],
                [zero, one - thickness, one],
                [zero, zero, zero],
                [zero, one, one],
                [half, zero, thickness],
                [half, one - thickness, one],
                [half, zero, zero],
                [half, one, one],
            ]
        )

        front_nw_corner = np.array(
            [
                [zero, one - thickness, zero],
                [zero, one - thickness / 2, thickness / 2],
                [zero, one, zero],
                [zero, one, thickness],
                [half, one - thickness, zero],
                [half, one - thickness / 2, thickness / 2],
                [half, one, zero],
                [half, one, thickness],
            ]
        )

        front_se_corner = np.array(
            [
                [zero, zero, one - thickness],
                [zero, zero, one],
                [zero, thickness / 2, one - thickness / 2],
                [zero, thickness, one],
                [half, zero, one - thickness],
                [half, zero, one],
                [half, thickness / 2, one - thickness / 2],
                [half, thickness, one],
            ]
        )

        back_bottom = np.array(
            [
                [half, zero, one - thickness],
                [half, zero, one],
                [half, one - thickness, zero],
                [half, one, zero],
                [one, zero, one - thickness],
                [one, zero, one],
                [one, one - thickness, zero],
                [one, one, zero],
            ]
        )

        back_top = np.array(
            [
                [half, zero, one],
                [half, thickness, one],
                [half, one, zero],
                [half, one, thickness],
                [one, zero, one],
                [one, thickness, one],
                [one, one, zero],
                [one, one, thickness],
            ]
        )

        back_sw_corner = np.array(
            [
                [half, zero, zero],
                [half, zero, thickness],
                [half, thickness, zero],
                [half, thickness / 2, thickness / 2],
                [one, zero, zero],
                [one, zero, thickness],
                [one, thickness, zero],
                [one, thickness / 2, thickness / 2],
            ]
        )

        back_ne_corner = np.array(
            [
                [half, one - thickness, one],
                [half, one, one],
                [half, one - thickness / 2, one - thickness / 2],
                [half, one, one - thickness],
                [one, one - thickness, one],
                [one, one, one],
                [one, one - thickness / 2, one - thickness / 2],
                [one, one, one - thickness],
            ]
        )

        for control_points in [
            front_top,
            front_bottom,
            front_nw_corner,
            front_se_corner,
            back_top,
            back_bottom,
            back_sw_corner,
            back_ne_corner,
        ]:
            spline_list.append(
                sp.BSpline(
                    degrees=[1, 1, 1],
                    control_points=control_points,
                    knot_vectors=[[0.0, 0.0, 1.0, 1.0]] * 3,
                )
            )

        return spline_list

    def create_twisted_to_nontwisted_material_linkage(
        self, gap_length, cell_length
    ):
        zero = 0.0
        one = cell_length

        spline_list = []

        xmin = np.array(
            [
                [zero, gap_length, zero],
                [zero, zero, gap_length],
                [gap_length / 2, gap_length / 2, gap_length / 2],
                [gap_length, zero, gap_length],
                [zero, one - gap_length, zero],
                [zero, one, gap_length],
                [gap_length / 2, one - gap_length / 2, gap_length / 2],
                [gap_length, one, gap_length],
            ]
        )

        xmax = np.array(
            [
                [one, gap_length, zero],
                [one, zero, gap_length],
                [one - gap_length / 2, gap_length / 2, gap_length / 2],
                [one - gap_length, zero, gap_length],
                [one, one - gap_length, zero],
                [one, one, gap_length],
                [one - gap_length / 2, one - gap_length / 2, gap_length / 2],
                [one - gap_length, one, gap_length],
            ]
        )

        ymin = np.array(
            [
                [zero, zero, zero],
                [zero, zero, gap_length / 2],
                [zero, gap_length, zero],
                [zero, zero, gap_length],
                [one, zero, zero],
                [one, zero, gap_length / 2],
                [one, gap_length, zero],
                [one, zero, gap_length],
            ]
        )

        ymax = np.array(
            [
                [zero, one, zero],
                [zero, one, gap_length / 2],
                [zero, one - gap_length, zero],
                [zero, one, gap_length],
                [one, one, zero],
                [one, one, gap_length / 2],
                [one, one - gap_length, zero],
                [one, one, gap_length],
            ]
        )

        for patch_cps in [xmin, xmax, ymin, ymax]:
            patch = sp.Bezier(
                degrees=[1, 1, 1], control_points=patch_cps
            ).bspline
            spline_list.append(patch)

        return spline_list

    def create_nontwisted_to_twisted_material_linkage(
        self, gap_length, cell_length
    ):
        spline_list = self.create_twisted_to_nontwisted_material_linkage(
            gap_length, cell_length
        )

        for patch in spline_list:
            cps = np.array(patch.cps)
            patch.cps[:4, :2] = cps[4:, :2]
            patch.cps[4:, :2] = cps[:4, :2]

        return spline_list

    def create_twisted_material_microstructure(self):
        assert (
            self.twisted_layers is not None
        ), "Simulation domain has not been created yet"

        # Create single unit tile and size it for actual size in microstructure
        base_tile = self.create_material_unit_tile()
        for patch in base_tile:
            patch.cps *= self.dd

        # Twist tile
        twisted_tile = []
        for patch in base_tile:
            twisted_patch = patch.copy()
            cps = np.array(patch.cps)
            twisted_patch.cps[:, 0] = cps[:, 1]
            twisted_patch.cps[:, 1] = cps[:, 0]
            twisted_tile.append(twisted_patch)

        min_positions = [self._macro_spline.cps[:, i].min() for i in range(3)]
        max_positions = [self._macro_spline.cps[:, i].max() for i in range(3)]

        # Determine base point for each tile
        x_starts = np.linspace(
            min_positions[0], max_positions[0], self._tiling[0] + 1
        )[:-1]
        y_starts = np.linspace(
            min_positions[1], max_positions[1], self._tiling[1] + 1
        )[:-1]

        base_points = [
            array.ravel(order="F") for array in np.meshgrid(x_starts, y_starts)
        ]
        base_points.reverse()
        base_points = np.vstack(base_points).T
        base_points = np.hstack(
            (base_points, np.zeros((base_points.shape[0], 1)))
        )

        # Create one layer of normal and one of twisted tiles, not placing it yet
        layer_tiles = []
        twisted_layer_tiles = []
        for basis, tile_patches in zip(
            [base_tile, twisted_tile], [layer_tiles, twisted_layer_tiles]
        ):
            for base_point in base_points:
                tile_patches_base = basis.copy()
                for patch in tile_patches_base:
                    new_patch = patch.copy()
                    new_patch.cps += base_point
                    href(new_patch)
                    tile_patches.append(new_patch)

        all_tiles = []
        z_offset = 0.0
        dz = (max_positions[2] - min_positions[2]) / self._tiling[2]
        bar_thickness = self._bar_thickness / self._tiling[1]

        # Go through each layer in z-direction
        for i in range(self._tiling[2]):
            # Get layer
            if i in self.twisted_layers:
                tile_patches_base = twisted_layer_tiles.copy()
            else:
                tile_patches_base = layer_tiles.copy()
            # Place respective layer in respective z-position
            for patch in tile_patches_base:
                new_patch = patch.copy()
                new_patch.cps[:, 2] += z_offset
                all_tiles.append(new_patch)

            # Add to layer offset in z-direction
            z_offset += dz

            # Create transitions if necessary
            twist_to_nontwist_transition = (
                i in self.twisted_layers and i + 1 not in self.twisted_layers
            )
            nontwist_to_twist_transition = (
                i not in self.twisted_layers and i + 1 in self.twisted_layers
            )
            # Decide on which transition to take
            if twist_to_nontwist_transition:
                linkage_patches = (
                    self.create_twisted_to_nontwisted_material_linkage(
                        gap_length=bar_thickness, cell_length=dz
                    )
                )
            elif nontwist_to_twist_transition:
                linkage_patches = (
                    self.create_nontwisted_to_twisted_material_linkage(
                        gap_length=bar_thickness, cell_length=dz
                    )
                )
            # Add transition
            if twist_to_nontwist_transition or nontwist_to_twist_transition:
                for patch in linkage_patches:
                    patch.cps[:, 2] += z_offset
                for base_point in base_points:
                    patch_base = linkage_patches.copy()
                    for patch in patch_base:
                        new_patch = patch.copy()
                        new_patch.cps += base_point
                        href(new_patch)
                        all_tiles.append(new_patch)

                z_offset += bar_thickness
                self._mixer_length += bar_thickness

        self._material_microstructure = sp.Multipatch(splines=all_tiles)
        self._material_microstructure.determine_interfaces()


if __name__ == "__main__":
    sulzer = SulzerInverse(
        macro_spline=sp.helpme.create.box(1, 1, BOX_LENGTH),
        tiling=TILING,
        thickness=0.1,
    )

    tile_splines = sulzer.create_unit_tile()

    sulzer.create_twisted_microstructure([2, 3])
    # sulzer.create_microstructure()
    sulzer.create_fore_and_afterrun()

    # sulzer.export()
    sulzer.show()

    # sulzer.create_twisted_material_microstructure()
    # sulzer.export_material_microstructure()
    # sulzer.show(show_material=True)
