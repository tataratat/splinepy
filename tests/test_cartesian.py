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
        # define some array
        a = c.np.array([0, 0.5, 1])

        def ref(arrays, reverse):
            """create reference answers using itertools.product"""

            prod = c.np.array(list(product(*arrays)))
            if reverse:
                return prod[:, ::-1]
            else:
                return prod

        # test cartesian product
        for d in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10):
            in_arr = [a] * d
            assert c.np.allclose(
                ref(in_arr, False),
                c.splinepy.utils.data.cartesian_product(in_arr, False),
            ), f"{d}-dim failed"

        # test reversed dim - what's often used in lib
        for d in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10):
            in_arr = [a] * d
            assert c.np.allclose(
                ref(in_arr, True),
                c.splinepy.utils.data.cartesian_product(in_arr, True),
            ), f"{d}-dim failed"


if __name__ == "__main__":
    c.unittest.main()
