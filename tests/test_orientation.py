try:
    from . import common as c
except BaseException:
    import common as c


class orientationTest(c.unittest.TestCase):
    """Test Orientation between adjacent spline patches
         2--3
         |C |
         0--1
    0--1 2--3 1--3
    |D | |A | | B|
    2--3 0--1 0--2
         1--0
    |    |E |
    O-   3--2
    """

    def test_orientation(self):
        # Init Splines to be tested
        a_s = c.splinepy.Bezier(
            degrees=[1, 1], control_points=[[1, 1], [2, 1], [1, 2], [2, 2]]
        )
        b_s = c.splinepy.Bezier(
            degrees=[1, 1], control_points=[[2, 1], [2, 2], [3, 1], [3, 2]]
        )
        c_s = c.splinepy.Bezier(
            degrees=[1, 1], control_points=[[1, 2], [2, 2], [1, 3], [2, 3]]
        )
        d_s = c.splinepy.Bezier(
            degrees=[1, 1], control_points=[[0, 2], [1, 2], [0, 1], [1, 1]]
        )
        e_s = c.splinepy.Bezier(
            degrees=[1, 1], control_points=[[2, 1], [1, 1], [2, 0], [1, 0]]
        )

        # Right Arm
        (
            axis_mapping,
            axis_orientation,
        ) = c.splinepy.splinepy_core.get_orientation(a_s, 1, b_s, 2, 0.001)
        self.assertTrue(c.np.all([1, 0] == axis_mapping))
        self.assertTrue(c.np.all([True, True] == axis_orientation))
        # Upper Arm
        (
            axis_mapping,
            axis_orientation,
        ) = c.splinepy.splinepy_core.get_orientation(a_s, 3, c_s, 2, 0.001)
        self.assertTrue(c.np.all([0, 1] == axis_mapping))
        self.assertTrue(c.np.all([True, True] == axis_orientation))
        # Left Arm
        (
            axis_mapping,
            axis_orientation,
        ) = c.splinepy.splinepy_core.get_orientation(a_s, 0, d_s, 1, 0.001)
        self.assertTrue(c.np.all([0, 1] == axis_mapping))
        self.assertTrue(c.np.all([True, False] == axis_orientation))
        # Lower Arm
        (
            axis_mapping,
            axis_orientation,
        ) = c.splinepy.splinepy_core.get_orientation(a_s, 2, e_s, 2, 0.001)
        self.assertTrue(c.np.all([0, 1] == axis_mapping))
        self.assertTrue(c.np.all([False, False] == axis_orientation))


if __name__ == "__main__":
    c.unittest.main()
