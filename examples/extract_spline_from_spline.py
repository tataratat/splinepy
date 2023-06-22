import gustaf as gus

import splinepy

if __name__ == "__main__":
    # create 2p2d nurbs
    s_nurbs = splinepy.Bezier(
        degrees=[2, 1],
        control_points=[
            [0, 0],
            [0.5, -0.5],
            [1, 0],
            [0, 1],
            [0.5, 0.5],
            [1, 1],
        ],
    ).nurbs
    s_nurbs.insert_knots(0, [0.5])

    # extract line by specifying ranges
    line = s_nurbs.extract.spline({0: [0.4, 0.8], 1: 0.7})

    gus.show(
        ["Source spline", s_nurbs],
        ["extract.spline({0: [0.4, 0.8], 1: 0.7})", line, s_nurbs],
        ["Extracted spline", line],
    )
