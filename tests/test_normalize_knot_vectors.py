import numpy as np


def test_bspline_normalize_knot_vectors(bspline_2p2d):
    """ """
    # make "iga book" bspline
    bspline = bspline_2p2d

    ref = [kv.numpy() for kv in bspline.knot_vectors]

    for kv in bspline.knot_vectors:
        kv[:] = np.multiply(np.add(kv, 37), 1.78934).tolist()

    # normalize
    bspline.normalize_knot_vectors()

    for i, (ref_kv, kv) in enumerate(zip(ref, bspline.knot_vectors)):
        assert np.allclose(ref_kv, kv), f"{i}. para dim failed to normalize"


def test_nurbs_normalize_knot_vectors(nurbs_2p2d):
    """ """
    nurbs = nurbs_2p2d

    ref = [kv.numpy() for kv in nurbs.knot_vectors]

    # manipulate - this does not update cpp spline, but it doesn't matter
    for kv in nurbs.knot_vectors:
        kv[:] = np.multiply(np.add(kv, 37), 1.78934).tolist()

    # normalize
    nurbs.normalize_knot_vectors()

    for i, (ref_kv, kv) in enumerate(zip(ref, nurbs.knot_vectors)):
        assert np.allclose(ref_kv, kv), f"{i}. para dim failed to normalize"
