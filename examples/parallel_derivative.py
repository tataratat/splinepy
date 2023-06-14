"""
Parallel evaluation of spline, which requires to pickle spline.
"""

from functools import partial
from multiprocessing import Pool
from time import perf_counter as tic

import numpy as np

import splinepy

if __name__ == "__main__":
    # define a bspline
    ds = [2, 2]
    kvs = [
        [0, 0, 0, 0.5, 1, 1, 1],
        [0, 0, 0, 1, 1, 1],
    ]
    cps = [
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

    b = splinepy.BSpline(
        degrees=ds,
        knot_vectors=kvs,
        control_points=cps,
    )

    # elevate degree to make evaluation take a bit longer
    b.elevate_degrees(0)
    b.elevate_degrees(0)
    b.elevate_degrees(0)
    b.elevate_degrees(1)
    b.elevate_degrees(1)
    b.elevate_degrees(1)

    # let's try nurbs
    n = b.nurbs

    # workers
    workers = 4

    # form queries
    qs = [np.random.random((5000, 2)) for _ in range(workers)]
    der = [2, 3]

    # full multiprocessing time
    now = tic()
    pool = Pool(workers)
    res_mp = list(pool.imap(partial(n.derivative, orders=der), qs))
    pool.close()
    pool.join()
    mptime0 = tic() - now
    print(f"multiprocess (with setup) evaluation took {mptime0} seconds.")

    # only evaluation multiprocessing time
    pool = Pool(workers)
    now = tic()
    res_mp = list(pool.imap(partial(n.derivative, orders=der), qs))
    mptime1 = tic() - now
    print(f"multiprocess (without setup) evaluation took {mptime1} seconds.")
    pool.close()
    pool.join()

    # stack batched queries
    q = np.vstack(qs)

    # multithread eval
    now = tic()
    res_t = n.derivative(q, der, nthreads=4)
    ttime = tic() - now
    print(f"multithread evaluation took {ttime} seconds.")

    # serial eval
    now = tic()
    res_s = n.derivative(q, der)
    stime = tic() - now
    print(f"serial evaluation took {stime} seconds.")

    print(f"serial / multiprocess (with setup) = {stime / mptime0}")
    print(f"serial / multiprocess (without setup) = {stime / mptime1}")
    print(f"serial / threading = {stime / ttime}")

    # cool, but are they correct?
    assert np.allclose(
        np.vstack(res_mp), res_s
    ), "Something went wrong during evaluation"
    assert np.allclose(res_t, res_s), "Something went wrong during evaluation"
