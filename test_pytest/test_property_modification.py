import numpy as np
import pytest

import splinepy

# fixtures for test_new_control_point_pointers_creation
all_2p2d_splines = (
    "rational_bezier_2p2d",
    "bezier_2p2d",
    "bspline_2p2d",
    "nurbs_2p2d",
)


# often called function; no testfunction
def cps_are_synced(spl):
    assert np.allclose(
        spl.cps, spl.current_core_properties()["control_points"]
    )


@pytest.mark.parametrize("splinetype", all_2p2d_splines)
def test_inplace_change_degrees_para(splinetype, request):
    """inplace change of degrees should not be allowed if core spline is
    initialized"""
    spline = request.getfixturevalue(splinetype)

    spline_dict = spline.todict()
    s = spline.__class__()
    s.degrees = spline_dict["degrees"]
    s.degrees += 1

    # this shoundn't be fine
    s._new_core(**spline_dict)
    with pytest.raises(ValueError):
        s.degrees += 1


def test_inplace_change_knot_vectors(nd_box, raster):
    """test inplace change of knot_vectors"""
    # let's test 3D splines
    dim = 3
    box_data = nd_box(dim)
    n = splinepy.NURBS(**box_data)
    box_data.pop("weights")
    b = splinepy.BSpline(**box_data)

    # elevate degrees and insert some knots
    knot_insert_dims = [1, 2]
    for s in (n, b):
        s.elevate_degrees([0, 0, 1, 2])
        for kid in knot_insert_dims:
            s.insert_knots(kid, np.linspace(0.1, 0.9, 9))

    qres = 5
    raster_query = raster([[0] * dim, [1] * dim], [qres] * dim)

    for s in (n, b):
        # is this valid spline?
        assert np.allclose(raster_query, s.evaluate(raster_query))

        # modify knots and corresponding queries, so that
        # evaluated points are same as raster_query
        modified_query = raster_query.copy()
        factor = 2
        for kid in knot_insert_dims:
            s.knot_vectors[kid][:] = np.multiply(
                s.knot_vectors[kid], factor
            ).tolist()
            modified_query[:, kid] *= factor

        # evaluation check
        assert np.allclose(raster_query, s.evaluate(modified_query))


def test_inplace_change_control_points(nd_box):
    """test inplace changes of control points"""
    dim = 3
    res = [3] * dim
    box_data = nd_box(dim)

    n = splinepy.NURBS(**box_data)
    weights = box_data.pop("weights")
    b = splinepy.BSpline(**box_data)
    box_data.pop("knot_vectors")
    z = splinepy.Bezier(**box_data)
    r = splinepy.RationalBezier(**box_data, weights=weights)

    for s in (z, r, b, n):
        # init
        orig = s.copy()

        # check if we are good to start
        assert np.allclose(
            orig.sample(res), s.sample(res)
        ), f"{s.whatami} failed"

        # hold reference to control points
        s_cps = s.control_points
        # modify cps
        s.control_points /= 2
        assert np.allclose(
            orig.sample(res) / 2, s.sample(res)
        ), f"{s.whatami} failed"
        # check if they are still the same array
        assert s_cps is s.control_points, f"{s.whatami} failed"

        # modify ws
        if s.is_rational:
            s_ws = s.weights
            s.weights *= 0.5
            assert s_ws is s.weights, f"{s.whatami} failed"


def test_inplace_change_weights(nurbs_2p2d):
    """test inplace change of weights by comparing quarter circle"""
    n_q_circle = nurbs_2p2d.todict()
    res = [4] * 2
    # init rational
    n = splinepy.NURBS(**n_q_circle)
    kvs = n_q_circle.pop("knot_vectors")
    r = splinepy.RationalBezier(**n_q_circle)
    # init non-rational
    n_q_circle.pop("weights")
    b = splinepy.BSpline(**n_q_circle, knot_vectors=kvs)
    z = splinepy.Bezier(**n_q_circle)

    # make sure they aren't the same
    assert not np.allclose(b.sample(res), n.sample(res))
    assert not np.allclose(z.sample(res), r.sample(res))

    # modify weights of rational splines
    for rs in (n, r):
        # set them all to 1
        rs.weights[abs(rs.weights - 1.0) > 1e-10] = 1.0

    # now, they should be the same
    assert np.allclose(b.sample(res), n.sample(res))
    assert np.allclose(z.sample(res), r.sample(res))
    assert np.allclose(b.sample(res), r.sample(res))
    assert np.allclose(z.sample(res), n.sample(res))


def test_physical_space_array(nurbs_2p2d):
    """tests if physical space array is syncing correctly."""
    n = splinepy.NURBS(**nurbs_2p2d.todict())

    # all setitem cases
    # 1.
    n.cps[0] = 100.0
    cps_are_synced(n)

    n.cps[np.array(1)] = [0.99, 0.11]  # here, weight is not 1
    cps_are_synced(n)

    # 2.
    n.cps[np.array([1, 4, 3])] = 123
    cps_are_synced(n)

    b_mask = np.ones(len(n.cps), dtype=bool)
    b_mask[[0, 4]] = False
    n.cps[b_mask] = -1235.5
    cps_are_synced(n)

    # 3.
    n.cps[[2, 0, 5]] = 299354.0
    cps_are_synced(n)

    n.cps[b_mask.tolist()] = 975.5
    cps_are_synced(n)

    # 4.
    n.cps[3, :2] = [-10, -20]
    cps_are_synced(n)

    n.cps[5, [1, 0]] = [1.722, 24.29]
    cps_are_synced(n)

    n.cps[np.array(2), 1] = 7
    cps_are_synced(n)

    # 5.
    n.cps[[0, 5, 2], 0] = 92314.111234
    cps_are_synced(n)
    n.cps[b_mask, :2] = n.cps[[2, 5, 0, 1], -2:]
    cps_are_synced(n)
    n.cps[b_mask.tolist(), 0] = n.cps[b_mask[::-1], 1]
    cps_are_synced(n)

    # 6.
    n.cps[:1, :2] = 160.4
    cps_are_synced(n)
    n.cps[2:-1, :2] = 78883.2
    n.cps[-1:-4:-1, :1] = -1687
    cps_are_synced(n)

    # 7.
    n.cps[...] = 1237792.4634527
    n.cps[..., :1] = -62.5
    cps_are_synced(n)

    # 8. ufunc
    np.add(n.cps, n.cps, out=n.cps)
    cps_are_synced(n)
    np.multiply(n.cps, n.cps, out=n.cps)
    cps_are_synced(n)

    # 9. if you have ideas for use case not listed above, please add!

    # now, child arrays
    carr = n.cps[0]
    assert carr.base is n.cps

    carr += [1324235.1, -234.5]
    assert np.allclose(carr, n.cps[0])
    cps_are_synced(n)

    # copy should return a normal np.ndarray
    ndarr = n.cps.copy()
    assert isinstance(ndarr, np.ndarray)
    assert not isinstance(ndarr, type(n.cps))
    assert not isinstance(ndarr, splinepy.utils.data.PhysicalSpaceArray)


@pytest.mark.parametrize("splinetype", all_2p2d_splines)
def test_new_control_point_pointers_creation(splinetype, request):
    """Test cp pointers are newly created if they change.
    _sync_source_ptr compares len of cp pointers and PhysicalSpaceArray
    before syncing.
    """
    # get spline
    spline = request.getfixturevalue(splinetype)

    # initial check
    spline.cps._sync_source_ptr()
    cps_are_synced(spline)

    # check after degree elevation
    spline.elevate_degrees([0, 1])
    spline.cps._sync_source_ptr()
    cps_are_synced(spline)

    # and after knot vectors and reduce degree
    if spline.has_knot_vectors:
        spline.insert_knots(0, [0.2, 0.5, 0.7])
        spline.cps._sync_source_ptr()
        cps_are_synced(spline)

        spline.insert_knots(1, [0.2, 0.5, 0.7])
        spline.cps._sync_source_ptr()
        cps_are_synced(spline)

        spline.remove_knots(0, [0.7])
        spline.cps._sync_source_ptr()
        cps_are_synced(spline)

        # reduce degree only available for bsplines
        spline.reduce_degrees([0])
        spline.cps._sync_source_ptr()
        cps_are_synced(spline)
