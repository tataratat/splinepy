try:
    from . import common as c
except BaseException:
    import common as c


class extractBoundariesTest(c.unittest.TestCase):
    def testExtraction2D(self):
        """Create a Spline, extract boundaries and check validity"""
        # uniform
        bez_el0 = c.splinepy.Bezier(
            degrees=[1, 1], control_points=[[0, 0], [1, 0], [0, 1], [1, 1]]
        )
        boundary_0 = c.splinepy.Bezier(
            degrees=[1], control_points=[[0, 0], [0, 1]]
        )
        boundary_3 = c.splinepy.Bezier(
            degrees=[1], control_points=[[0, 1], [1, 1]]
        )
        list_of_boundaries = bez_el0.extract_boundaries([0, 3])
        self.assertTrue(c.are_splines_equal(list_of_boundaries[0], boundary_0))
        self.assertTrue(c.are_splines_equal(list_of_boundaries[1], boundary_3))

        # non-uniform
        bez_el0 = c.splinepy.NURBS(
            knot_vectors=[[0, 0, 0.5, 1, 1], [0, 0, 1, 1]],
            weights=[1.0, 0.5, 1, 1, 0.5, 1],
            degrees=[1, 1],
            control_points=[
                [0, 0],
                [0.5, 0.5],
                [1, 0],
                [0, 1],
                [0.5, 1.5],
                [1, 1],
            ],
        )
        boundary_2 = c.splinepy.NURBS(
            knot_vectors=[[0, 0, 0.5, 1, 1]],
            weights=[1.0, 0.5, 1],
            degrees=[1],
            control_points=[[0, 0], [0.5, 0.5], [1, 0]],
        )
        list_of_boundaries = bez_el0.extract_boundaries([2])
        self.assertTrue(c.are_splines_equal(list_of_boundaries[0], boundary_2))

    def testExtraction3D(self):
        """Create a Spline, extract boundaries and check validity"""
        # uniform
        bez_el0 = c.splinepy.Bezier(
            degrees=[1, 1, 1],
            control_points=[
                [0, 0, 0],
                [1, 0, 0],
                [0, 1, 0],
                [1, 1, 0],
                [0, 0, 1],
                [1, 0, 1],
                [0, 1, 1],
                [1, 1, 1],
            ],
        )
        boundary_0 = c.splinepy.Bezier(
            degrees=[1, 1],
            control_points=[
                [0, 0, 0],
                [0, 1, 0],
                [0, 0, 1],
                [0, 1, 1],
            ],
        )
        list_of_boundaries = bez_el0.extract_boundaries([0])
        self.assertTrue(c.are_splines_equal(list_of_boundaries[0], boundary_0))


if __name__ == "__main__":
    c.unittest.main()
