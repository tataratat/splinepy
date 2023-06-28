try:
    from . import common as c
except BaseException:
    import common as c


class MultiIndexTest(c.unittest.TestCase):
    def test_multi_index(self):
        """
        Test MultiIndex using a round trip of np.ravel_multi_index
        """
        # 4d? can actually do anything
        for d in range(2, 10 + 1):
            resolutions = c.np.random.randint(2, 5, size=d)

            multi = c.splinepy.helpme.multi_index.MultiIndex(resolutions)

            # test individual id query
            n_q = 3  # number of queries
            query = [c.np.random.randint(0, r, size=n_q) for r in resolutions]

            # test index round trip
            # note reversed order - alternatively this should also be the same:
            # c.np.ravel_multi_index(query, resolutions, order="F")
            ref_raveled = c.np.ravel_multi_index(
                query[::-1], resolutions[::-1]
            )
            to_test = multi[tuple(query)]  # can't unpack inside []
            assert c.np.array_equal(ref_raveled, to_test)

            # test first and last slice that's orthogonal to the last dimension
            # first slice
            ref_raveled = c.np.arange(c.np.prod(resolutions[:-1]))
            to_test = multi[..., 0]
            assert c.np.array_equal(ref_raveled, to_test)

            # last slice
            upper = c.np.prod(resolutions)
            lower = upper - c.np.prod(resolutions[:-1])
            ref_raveled = c.np.arange(lower, upper)
            to_test = multi[..., -1]
            assert c.np.array_equal(ref_raveled, to_test)

            # test first and last slice that's orthogonal to the first dim
            full_ids = c.np.arange(c.np.prod(resolutions))

            # first slice
            ref_raveled = full_ids[0 :: resolutions[0]]
            to_test = multi[0, ...]
            assert c.np.array_equal(ref_raveled, to_test)
            # syntex check - this should hold in np.ndarray.__getitem__ syntax
            to_test = multi[0, :]
            assert c.np.array_equal(ref_raveled, to_test)

            # last slice
            ref_raveled = full_ids[int(resolutions[0] - 1) :: resolutions[0]]
            to_test = multi[-1, ...]
            assert c.np.array_equal(ref_raveled, to_test)


if __name__ == "__main__":
    c.unittest.main()
