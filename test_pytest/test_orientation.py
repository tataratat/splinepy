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

        # Provide connectivity data
        spline_list = [a_s, b_s, c_s, d_s, e_s]
        start_ids = [0, 0, 0, 0]
        start_face_ids = [1, 3, 0, 2]
        neighbor_ids = [1, 2, 3, 4]
        neighbor_face_ids = [2, 2, 1, 2]

        # Determine orientation single thread
        (
            axis_mapping,
            axis_orientation,
        ) = c.splinepy.splinepy_core.orientations(
            spline_list,
            start_ids,
            start_face_ids,
            neighbor_ids,
            neighbor_face_ids,
            0.00001,
            1,
        )

        # Check results
        expected_mappings = [[1, 0], [0, 1], [0, 1], [0, 1]]
        expected_orientations = [
            [True, True],
            [True, True],
            [True, False],
            [False, False],
        ]
        self.assertTrue(c.np.all(expected_mappings == axis_mapping))
        self.assertTrue(c.np.all(expected_orientations == axis_orientation))

        # Repeat with multithread execution
        (
            axis_mapping,
            axis_orientation,
        ) = c.splinepy.splinepy_core.orientations(
            spline_list,
            start_ids,
            start_face_ids,
            neighbor_ids,
            neighbor_face_ids,
            0.00001,
            3,
        )
        self.assertTrue(c.np.all(expected_mappings == axis_mapping))
        self.assertTrue(c.np.all(expected_orientations == axis_orientation))


if __name__ == "__main__":
    c.unittest.main()
