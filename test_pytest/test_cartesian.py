try:
    from . import common as c
except BaseException:
    import common as c

from itertools import product


class CartesianProductTest(c.unittest.TestCase):
    def test_cartesian_product(self):
        """
        Test cartesian product test
        """
        # define arrays for cartesian product
        # make different len, but not too long so that test
        # doesn't take too long.
        max_dim = 10
        pool = c.np.array([0.0, 1.1, 4.4, 8.8, 9.6])
        queries = [pool[: c.np.random.randint(2, 5)] for _ in range(max_dim)]

        # in case they are all the same, let's change two entries
        queries[-1] = pool
        queries[-2] = pool[:3]

        def ref(arrays, reverse):
            """create reference answers using itertools.product"""
            if reverse:
                prod = c.np.array(list(product(*arrays[::-1])))
                return prod[:, ::-1]
            else:
                prod = c.np.array(list(product(*arrays)))
                return prod

        # test cartesian product
        for i in range(1, max_dim + 1):
            in_arr = queries[:i]
            assert c.np.allclose(
                ref(in_arr, False),
                c.splinepy.utils.data.cartesian_product(in_arr, False),
            ), f"{i}-dim failed with NO reverse"

        # test reversed dim - what's often used in lib
        for i in range(1, max_dim + 1):
            in_arr = queries[:i]
            assert c.np.allclose(
                ref(in_arr, True),
                c.splinepy.utils.data.cartesian_product(in_arr, True),
            ), f"{i}-dim failed with reverse"


if __name__ == "__main__":
    c.unittest.main()
