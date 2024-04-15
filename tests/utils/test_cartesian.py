from itertools import product

import numpy as np

import splinepy


def test_cartesian_product(np_rng):
    """Test cartesian product test"""

    # define arrays for cartesian product
    # make different len, but not too long so that test
    # doesn't take too long.
    max_dim = 10
    pool = np.array([0.0, 1.1, 4.4, 8.8, 9.6])
    queries = [pool[: np_rng.integers(2, 5)] for _ in range(max_dim)]

    # in case they are all the same, let's change two entries
    queries[-1] = pool
    queries[-2] = pool[:3]

    def ref(arrays, reverse):
        """create reference answers using itertools.product"""
        if reverse:
            prod = np.array(list(product(*arrays[::-1])))
            return prod[:, ::-1]
        else:
            prod = np.array(list(product(*arrays)))
            return prod

    # test cartesian product
    for i in range(1, max_dim + 1):
        in_arr = queries[:i]
        assert np.allclose(
            ref(in_arr, False),
            splinepy.utils.data.cartesian_product(in_arr, False),
        ), f"{i}-dim failed with NO reverse"

    # test reversed dim - what's often used in lib
    for i in range(1, max_dim + 1):
        in_arr = queries[:i]
        assert np.allclose(
            ref(in_arr, True),
            splinepy.utils.data.cartesian_product(in_arr, True),
        ), f"{i}-dim failed with reverse"
