try:
    from . import common as c
except BaseException:
    import common as c


class MultipatchTest(c.unittest.TestCase):
    def setUp(self) -> None:
        """ """
        # Spline that points to itself
        self._self_referencing_spline = c.splinepy.Bezier(
            degrees=[3, 1],
            control_points=[
                [0, 0],
                [2, 3],
                [-2, 3],
                [0, 0],
                [0, 1],
                [1, 2],
                [-1, 2],
                [0, 1],
            ],
        )

        # List of splinewith with this topology
        #
        #  o ---- o ---- o
        #  |   0  |   1  |
        #  o ---- o ---- o
        #  |   2  |
        #  o ---- o
        #
        self._rect_arc_1 = c.splinepy.Bezier(
            [2, 1], [[0, 0], [1, 0], [2, 0], [0, 1], [1, 2], [2, 1]]
        )
        self._rect_arc_2 = c.splinepy.Bezier(
            [1, 1], [[2, 0], [3, 0], [2, 1], [3, 1]]
        )
        self._rect_arc_3 = c.splinepy.Bezier(
            [2, 1], [[0, -1], [1, -2], [2, -1], [0, 0], [1, 0], [2, 0]]
        )
        self._list_of_splines = [
            self._rect_arc_1,
            self._rect_arc_2,
            self._rect_arc_3,
        ]

    def test_interfaces(self):
        """ """
        # init multipatch with multiple splines
        multipatch = c.splinepy.Multipatch()
        multipatch.patches = self._list_of_splines

        # init multipatch with single spline
        single_p_multipatch = c.splinepy.Multipatch()
        single_p_multipatch.patches = [self._self_referencing_spline]

        # Determine connectivities
        single_p_multipatch.determine_interfaces()
        multipatch.determine_interfaces()
        self.assertTrue(
            (
                single_p_multipatch.interfaces
                == c.np.array([[1, 0, -1, -1]], dtype=int)
            ).all()
        )
        self.assertTrue(
            (
                multipatch.interfaces
                == c.np.array(
                    # Global Face IDs and their connectivity
                    # [0  1   2   3]  [4   5   6   7]   [8   9  10  11]
                    [[-1, 4, 11, -1], [1, -1, -1, -1], [-1, -1, -1, 2]],
                    dtype=int,
                )
            ).all()
        )

    def test_boundaries(self):
        """ """
        # init multipatch with multiple splines
        multipatch = c.splinepy.Multipatch()
        multipatch.patches = self._list_of_splines
        multipatch.determine_interfaces()

        # Using a function
        def west_side(points):
            return points[:, 0] < 0.1

        multipatch.boundary_from_function(west_side)
        self.assertTrue(
            (
                multipatch.interfaces
                == c.np.array(
                    # Global Face IDs and their connectivity
                    # [0  1   2   3]  [4   5   6   7]   [8   9  10  11]
                    [[-2, 4, 11, -1], [1, -1, -1, -1], [-2, -1, -1, 2]],
                    dtype=int,
                )
            ).all()
        )

        multipatch.set_boundary([0, 1], [3, 3], 3)
        multipatch.set_boundary([1], [1])
        self.assertTrue(
            (
                multipatch.interfaces
                == c.np.array(
                    # Global Face IDs and their connectivity
                    # [0  1   2   3]  [4   5   6   7]   [8   9  10  11]
                    [[-2, 4, 11, -3], [1, -4, -1, -3], [-2, -1, -1, 2]],
                    dtype=int,
                )
            ).all()
        )

        # Delete all boundaries and determine new ones based on continuity
        multipatch.boundaries_from_continuity()
        self.assertTrue(
            (
                multipatch.interfaces
                == c.np.array(
                    # Global Face IDs and their connectivity
                    # [0  1   2   3]  [4   5   6   7]   [8   9  10  11]
                    [[-1, 4, 11, -2], [1, -3, -4, -5], [-1, -6, -7, 2]],
                    dtype=int,
                )
            ).all()
        )

    def test_interfaces_and_boundaries(self):
        # 2 --- 3 1 --- 0
        # |  1  | |  3  |
        # 0 --- 1 3 --- 2
        # 3 --- 2 0 --- 2
        # |  2  | |  4  |
        # 1 --- 0 1 --- 3
        #
        # with 1 and three flipped in the 3rd axis

        # Create some splines in a random parametric swaps
        b1 = c.splinepy.Bezier(
            degrees=[1, 1, 1],
            control_points=[
                [-2, 0, 3],
                [0, 0, 3],
                [-2, 1, 3],
                [0, 1, 3],
                [-2, 0, 0],
                [0, 0, 0],
                [-2, 1, 0],
                [0, 1, 0],
            ],
        )
        b2 = c.splinepy.Bezier(
            degrees=[1, 1, 1],
            control_points=[
                [0, -1, 0],
                [-2, -1, 0],
                [0, 0, 0],
                [-2, 0, 0],
                [0, -1, 3],
                [-2, -1, 3],
                [0, 0, 3],
                [-2, 0, 3],
            ],
        )
        b3 = c.splinepy.Bezier(
            degrees=[1, 1, 1],
            control_points=[
                [1, 1, 3],
                [0, 1, 3],
                [1, 0, 3],
                [0, 0, 3],
                [1, 1, 0],
                [0, 1, 0],
                [1, 0, 0],
                [0, 0, 0],
            ],
        )
        b4 = c.splinepy.Bezier(
            degrees=[1, 1, 1],
            control_points=[
                [0, 0, 0],
                [0, -1, 0],
                [1, 0, 0],
                [1, -1, 0],
                [0, 0, 3],
                [0, -1, 3],
                [1, 0, 3],
                [1, -1, 3],
            ],
        )

        # Multipatch
        multipatch = c.splinepy.Multipatch([b1, b2, b3, b4])
        multipatch.boundaries_from_continuity()
        self.assertTrue(
            (
                multipatch.interfaces
                == c.np.array(
                    [
                        [-1, 13, 9, -2, -3, -4],
                        [20, -1, -5, 2, -4, -3],
                        [-6, 1, -2, 18, -3, -4],
                        [15, -5, 6, -6, -4, -3],
                    ],
                    dtype=int,
                )
            ).all()
        )


if __name__ == "__main__":
    c.unittest.main()
