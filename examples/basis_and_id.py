import splinepy
import numpy as np

if __name__ == "__main__":

    kv = [[0, 0, 0, 0.5, 1, 1, 1], [0, 0, 0, 1, 1, 1], [0, 0, 1, 1]]
    cp1 = np.array(
            [
                    [0, 0, 0],
                    [0, 1, 0],
                    [1, 1.5, 0],
                    [3, 1.5, 0],
                    [-1, 0, 0],
                    [-1, 2, 0],
                    [1, 4, 0],
                    [3, 4, 0],
                    [-2, 0, 0],
                    [-2, 2, 0],
                    [1, 5, 0],
                    [3, 5, 2],
            ]
    )

    cp1 = np.vstack((cp1, cp1 + [0, 0, 1]))

    # generate splines
    b = splinepy.BSpline(
            degrees=[2, 2, 1],
            knot_vectors=kv,
            control_points=cp1,
    )
    b.insert_knots(0, [.1, .2, .3, .4, .6, .7, .8, .9])
    b.insert_knots(1, [.1, .2, .3, .4, .5, .6, .7, .8, .9])
    b.insert_knots(2, [.1, .2, .3, .4, .5, .6, .7, .8, .9])

    n = b.nurbs
    # edit a weight so that the values aren't the same
    n.weights[0, 0] = .5  # this is a tall array

    # queries
    q = [
            [.01, .01, .01],
            #[.2, .3, .4],
            #[.3, .4, .5],
            #[.4, .5, .6],
            [.51, .61, .71],
    ]

    print("BSpline basis functions and support ids:")
    print(b.basis_functions(q))

    print("NURBS basis functions and support ids:")
    print(n.basis_functions(q))
