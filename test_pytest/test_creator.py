import numpy as np
import pytest

import splinepy

# frequently used fixtures
all_splinetypes = (
    "rational_bezier_2p2d",
    "bezier_2p2d",
    "bspline_2p2d",
    "nurbs_2p2d",
)


def test_create_embedded(request):
    """
    Test embedding for (rational) splines
    """
    # make a couple of 2D splines
    bspline = request.getfixturevalue("bspline_2p2d")
    nurbs = request.getfixturevalue("nurbs_2p2d")
    embedded_bspline = bspline.create.embedded(3)
    assert np.allclose(
        bspline.cps[:, :2],
        embedded_bspline.cps[:, :2],
    )
    embedded_nurbs = nurbs.create.embedded(3)
    assert np.allclose(
        nurbs.cps[:, :2],
        embedded_nurbs.cps[:, :2],
    )


# Test Extrusion routines
@pytest.mark.parametrize("splinetype", all_splinetypes)
def test_create_extrude(splinetype, np_rng, request):
    """
    Test extrusion for different input types and arguments
    """
    # create 2D spline
    spline = request.getfixturevalue(splinetype)

    # Expect Failure - not a spline
    with pytest.raises(NotImplementedError):
        splinepy.helpme.create.extruded([4])

    if splinetype == "bspline_2p2d":
        # Expect Failure - no axis given
        with pytest.raises(ValueError):
            spline.create.extruded()
        # Expect Failure - axis wrong format
        with pytest.raises(ValueError):
            spline.create.extruded(extrusion_vector=[1])

    # Create a random axis
    axis = np_rng.random(3)
    x, y, z = np_rng.random(3)

    # Test results
    assert np.allclose(
        spline.create.extruded(extrusion_vector=axis).evaluate([[x, y, z]]),
        np.hstack((spline.evaluate([[x, y]]), np.zeros((1, 1)))) + z * axis,
    )


# Test Revolution Routine
@pytest.mark.parametrize("splinetype", all_splinetypes)
def test_create_revolution(splinetype, np_rng, request):
    """
    Test revolution routines for different input types and arguments
    """
    # Make 2D spline
    spline = request.getfixturevalue(splinetype)

    # Make some lines
    bezier_line = splinepy.Bezier(control_points=[[1, 0], [2, 1]], degrees=[1])
    nurbs_line = bezier_line.nurbs

    # Make a cuboid
    cuboid = splinepy.Bezier(
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
        degrees=[1, 1, 1],
    )

    # Expect Failure - not a spline
    with pytest.raises(NotImplementedError):
        splinepy.helpme.create.revolved([4])
    # Expect Failure - No rotation axis
    with pytest.raises(ValueError):
        cuboid.create.revolved()

    if splinetype == "bspline_2p2d":
        # Expect Failure - axis wrong format
        with pytest.raises(ValueError):
            spline.create.revolved(axis=[1])
        # Expect Failure - axis too small
        with pytest.raises(ValueError):
            spline.create.revolved(axis=[0, 0, 1e-18])

    # Revolve always around z-axis
    # init rotation matrix
    r_angle = np_rng.random()
    r_center = np.array([1, 0])
    cc, ss = np.cos(r_angle), np.sin(r_angle)
    R = np.array([[cc, -ss, 0], [ss, cc, 0], [0, 0, 1]])
    R2 = np.array([[cc, -ss], [ss, cc]])

    # Test 3D revolutions for bodies
    dim_bumped_cps = np.zeros((spline.control_points.shape[0], 1))
    ref_sol = np.matmul(
        np.hstack((spline.control_points, dim_bumped_cps)), R.T
    )

    revolved_cps = spline.create.revolved(
        axis=[0, 0, 1],
        center=[0, 0, 0],
        angle=r_angle,
        degree=False,
    ).control_points

    assert np.allclose(
        revolved_cps[-(len(ref_sol)) :, :], ref_sol
    ), f"{spline.whatami} failed revolution"

    # Test 2D Revolutions of lines
    for spline_g in (bezier_line, nurbs_line):
        assert np.allclose(
            spline_g.create.revolved(
                angle=r_angle, degree=False
            ).control_points[-2:, :],
            np.matmul(spline_g.control_points, R2.T),
        )

    # Test 2D Revolutions of lines around center
    for spline_g in (bezier_line, nurbs_line):
        assert np.allclose(
            spline_g.create.revolved(
                angle=r_angle, center=r_center, degree=False
            ).control_points[-2:, :],
            np.matmul(spline_g.control_points - r_center, R2.T) + r_center,
        ), f"{spline_g.whatami} failed revolution around center"
