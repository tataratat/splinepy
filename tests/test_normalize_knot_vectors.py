import copy

try:
    from . import common as c
except BaseException:
    import common as c


class NormalizeKnotVectorsTest(c.unittest.TestCase):
    def test_bspline_normalize_knot_vectors(self):
        """ """
        # make "iga book" bspline
        bspline = c.splinepy.BSpline(**c.b2P2D)

        ref = copy.deepcopy(bspline.knot_vectors)

        # manipulate - this does not update cpp spline, but it doesn't matter
        for kv in bspline.knot_vectors:
            for i in range(len(kv)):
                kv[i] += 37  # apply some offset
                kv[i] *= 1.7956  # apply some scaling

        # normalize
        bspline.normalize_knot_vectors()

        for i, (ref_kv, kv) in enumerate(zip(ref, bspline.knot_vectors)):
            assert c.np.allclose(
                ref_kv, kv
            ), f"{i}. para dim failed to normalize"

    def test_nurbs_normalize_knot_vectors(self):
        """ """
        nurbs = c.splinepy.NURBS(**c.n2P2D)

        ref = copy.deepcopy(nurbs.knot_vectors)

        # manipulate - this does not update cpp spline, but it doesn't matter
        for kv in nurbs.knot_vectors:
            for i in range(len(kv)):
                kv[i] += 37  # apply some offset
                kv[i] *= 1.7956  # apply some scaling

        # normalize
        nurbs.normalize_knot_vectors()

        for i, (ref_kv, kv) in enumerate(zip(ref, nurbs.knot_vectors)):
            assert c.np.allclose(
                ref_kv, kv
            ), f"{i}. para dim failed to normalize"


if __name__ == "__main__":
    c.unittest.main()
