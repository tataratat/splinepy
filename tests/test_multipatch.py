try:
    from . import common as c
except BaseException:
    import common as c


class MultipatchTest(c.unittest.TestCase):

    def setUp(self) -> None:
        # Spline that points to itself
        self._self_referencing_spline = c.splinepy.Bezier(
                degrees=[3, 1],
                control_points=[
                        [0, 0], [2, 3], [-2, 3], [0, 0], [0, 1], [1, 2],
                        [-1, 2], [0, 1]
                ]
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
                self._rect_arc_1, self._rect_arc_2, self._rect_arc_3
        ]

        return

    def test_connectivity(self):
        """
        """
        # init multipatch with multiple splines
        multipatch = c.splinepy.Multipatch()
        multipatch.splines = self._list_of_splines

        # init multipatch with single spline
        single_p_multipatch = c.splinepy.Multipatch()
        single_p_multipatch.splines = [self._self_referencing_spline]

        # Determine connectivities
        single_p_multipatch.determine_connectivity()
        multipatch.determine_connectivity()
        self.assertTrue(
                (
                        single_p_multipatch.connectivity
                        == c.np.array([[0, 0, -1, -1]], dtype=int)
                ).all()
        )
        self.assertTrue(
                (
                        multipatch.connectivity == c.np.array(
                                [
                                        [-1, 1, 2, -1], [0, -1, -1, -1],
                                        [-1, -1, -1, 0]
                                ],
                                dtype=int
                        )
                ).all()
        )


if __name__ == "__main__":
    c.unittest.main()
