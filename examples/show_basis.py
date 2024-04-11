import splinepy as spp

if __name__ == "__main__":
    # curve
    # define degrees
    ds1 = [2]

    # define knot vectors
    kvs1 = [[0.0, 0.0, 0.0, 1.0, 1.0, 1.0]]

    # define control points
    cps1 = [
        [0.0, 0.0],
        [-2.0, 4.0],
        [1.0, 3.0],
    ]

    # init bspline
    b = spp.BSpline(
        degrees=ds1,
        knot_vectors=kvs1,
        control_points=cps1,
    )

    # make richer
    b.elevate_degrees([0, 0])
    b.insert_knots(0, [0.1, 0.2, 0.4, 0.7, 0.8, 0.8, 0.8])

    # press `+` key to plot axes behind
    spp.show(
        ["BSpline", b],
        ["Basis", *b.extract.bases()],
        ["Basis (physical space)", *b.extract.bases(parametric_view=False)],
        resolutions=201,
    )

    # 2D Basis and 1D Basis
    # define degrees
    ds2 = [2, 2]

    # define knot vectors
    kvs2 = [
        [0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
        [0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
    ]

    # define control points
    cps2 = [
        [0, 0],
        [0, 1],
        [1, 1.5],
        [3, 1.5],
        [-1, 0],
        [-1, 2],
        [1, 4],
        [3, 4],
        [-2, 0],
        [-2, 2],
        [1, 5],
        [3, 5],
    ]

    # init bspline
    b = spp.BSpline(
        degrees=ds2,
        knot_vectors=kvs2,
        control_points=cps2,
    )

    b.elevate_degrees([0, 0, 1])
    b.insert_knots(0, [0.1, 0.6, 0.8, 0.8, 0.8])
    b.insert_knots(1, [0.4, 0.6, 0.8, 0.8, 0.8])

    spp.show(
        ["BSpline", b],
        ["Basis", *b.extract.bases()],
        ["Basis 11, 19, 43", *b.extract.bases([11, 19, 43])],
        ["Basis (physical space)", *b.extract.bases(parametric_view=False)],
        resolutions=101,
        knots=False,
    )
