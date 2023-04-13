try:
    from . import common as c
except BaseException:
    import common as c

class InplaceModificationTest(c.unittest.TestCase):

    def test_inplace_change_knot_vectors(self):
        """test inplace change of knot_vectors"""
        # let's test 3D splines
        dim = 3
        box_data = c.nd_box(dim)
        n = c.splinepy.NURBS(**box_data)
        box_data.pop("weights")
        b = c.splinepy.BSpline(**box_data)

        # elevate degrees and insert some knots
        knot_insert_dims = [1,2]
        for s in (n, b):
            s.elevate_degrees([0,0,1,2])
            for kid in knot_insert_dims:
                s.insert_knots(kid, c.np.linspace(.1, .9, 9))

        qres = 2
        raster_query = c.raster([[0] * dim, [1] * dim], [qres] * dim)

        for s in (n, b):
            # is this valid spline?
            assert c.np.allclose(raster_query, s.evaluate(raster_query))

            # modify knots and corresponding queries, so that
            # evaluated points are same as raster_query
            modified_query = raster_query.copy()
            factor = 2
            for kid in knot_insert_dims:
                s.knot_vectors[kid] *= factor
                modified_query[:, kid] *= factor

                # check modified flag
                assert s.knot_vectors[kid]._modified
                # full modified check
                assert c.splinepy.spline.is_modified(s)

            # evaluation check
            assert c.np.allclose(raster_query, s.evaluate(modified_query))