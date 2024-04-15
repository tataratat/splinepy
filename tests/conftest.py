import contextlib
import os
import re

import numpy as np
import pytest

import splinepy


def error_log(*args):
    """used instead of print within the test"""
    splinepy.utils.log.error(*args)


@pytest.fixture
def np_rng():
    return np.random.default_rng()


# query points
@pytest.fixture
def queries_2D():
    return [
        [0.01, 0.01],
        [0.01, 0.5],
        [0.9, 0.1],
        [0.8, 0.7],
        [0.4, 0.99],
    ]


@pytest.fixture
def queries_3D():
    return [
        [0.1, 0.1, 0.1],
        [0.734, 0.525, 0.143],
        [0.9666, 0.991, 0.003],
        [0.5623, 0.0089, 0.99],
        [0.0431, 0.2, 0.523],
    ]


# initializing a spline should be a test itself, so provide `dict_spline`
# this is "iga-book"'s fig 2.15.
@pytest.fixture
def bspline_2p2d():
    return splinepy.BSpline(
        degrees=[2, 2],
        knot_vectors=[
            [0, 0, 0, 0.5, 1, 1, 1],
            [0, 0, 0, 1, 1, 1],
        ],
        control_points=[
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
        ],
    )


# half-half circle.
@pytest.fixture
def nurbs_2p2d():
    return splinepy.NURBS(
        degrees=[2, 1],
        knot_vectors=[
            [0, 0, 0, 1, 1, 1],
            [0, 0, 1, 1],
        ],
        control_points=[
            [-1.0, 0.0],
            [-1.0, 1.0],
            [0.0, 1.0],
            [-2.0, 0.0],
            [-2.0, 2.0],
            [0.0, 2.0],
        ],
        weights=[
            [1.0],
            [2**-0.5],
            [1.0],
            [1.0],
            [2**-0.5],
            [1.0],
        ],
    )


@pytest.fixture
def bezier_2p2d():
    return splinepy.Bezier(
        degrees=[2, 1],
        control_points=[
            [-1.0, 0.0],
            [-1.0, 1.0],
            [0.0, 1.0],
            [-2.0, 0.0],
            [-2.0, 2.0],
            [0.0, 2.0],
        ],
    )


@pytest.fixture
def rational_bezier_2p2d():
    return splinepy.RationalBezier(
        degrees=[2, 1],
        control_points=[
            [-1.0, 0.0],
            [-1.0, 1.0],
            [0.0, 1.0],
            [-2.0, 0.0],
            [-2.0, 2.0],
            [0.0, 2.0],
        ],
        weights=[
            [1.0],
            [2**-0.5],
            [1.0],
            [1.0],
            [2**-0.5],
            [1.0],
        ],
    )


@pytest.fixture
def bezier_3p3d():
    return splinepy.Bezier(
        degrees=[1, 1, 1],
        control_points=[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, -1.0, 1.0],
            [1.0, 0.0, 1.0],
            [-1.0, 1.0, 2.0],
            [2.0, 2.0, 2.0],
        ],
    )


@pytest.fixture
def rational_bezier_3p3d():
    return splinepy.RationalBezier(
        degrees=[1, 1, 1],
        control_points=[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, -1.0, 1.0],
            [1.0, 0.0, 1.0],
            [-1.0, 1.0, 2.0],
            [2.0, 2.0, 2.0],
        ],
        weights=[1.0] * 8,
    )


@pytest.fixture
def bspline_3p3d():
    return splinepy.BSpline(
        degrees=[1, 1, 1],
        control_points=[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, -1.0, 1.0],
            [1.0, 0.0, 1.0],
            [-1.0, 1.0, 2.0],
            [2.0, 2.0, 2.0],
        ],
        knot_vectors=[
            [0.0, 0.0, 1.0, 1.0],
            [0.0, 0.0, 1.0, 1.0],
            [0.0, 0.0, 1.0, 1.0],
        ],
    )


@pytest.fixture
def nurbs_3p3d():
    return splinepy.NURBS(
        degrees=[1, 1, 1],
        control_points=[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, -1.0, 1.0],
            [1.0, 0.0, 1.0],
            [-1.0, 1.0, 2.0],
            [2.0, 2.0, 2.0],
        ],
        weights=[1.0] * 8,
        knot_vectors=[
            [0.0, 0.0, 1.0, 1.0],
            [0.0, 0.0, 1.0, 1.0],
            [0.0, 0.0, 1.0, 1.0],
        ],
    )


@pytest.fixture
def raster():
    def _raster(bounds, resolutions):
        """prepares raster points using np.meshgrid"""
        l_bounds, u_bounds = bounds[0], bounds[1]
        pts = np.meshgrid(
            *[
                np.linspace(lo, up, re)
                for lo, up, re in zip(l_bounds, u_bounds, resolutions)
            ],
            indexing="ij",
        )
        # return pts
        return np.hstack([p.reshape(-1, 1) for p in pts[::-1]])

    return _raster


@pytest.fixture
def nd_box(raster):
    def _nd_box(dim):
        """creates simple box in n-d"""
        ds = [1 for _ in range(dim)]
        cps = raster(
            [[0 for _ in range(dim)], [1 for _ in range(dim)]],
            [2 for _ in range(dim)],
        )
        kvs = [[0, 0, 1, 1] for _ in range(dim)]
        ws = np.ones((len(cps), 1))
        return {
            "degrees": ds,
            "control_points": cps,
            "knot_vectors": kvs,
            "weights": ws,
        }

    return _nd_box


@pytest.fixture
def to_tmpf():
    def _to_tmpf(tmpd):
        """given tmpd, returns tmpf"""
        return os.path.join(tmpd, "nqv248p90")

    return _to_tmpf


@pytest.fixture
def are_splines_equal():
    def _are_splines_equal(a, b, print_=False):
        """returns True if Splines are equivalent"""
        if a.whatami != b.whatami:
            return False
        for req_prop in a.required_properties:
            if req_prop == "knot_vectors":
                for aa, bb in zip(a.knot_vectors, b.knot_vectors):
                    if not np.allclose(aa.numpy(), bb.numpy()):
                        if print_:
                            error_log("a.kvs", a.kvs)
                            error_log("b.kvs", b.kvs)
                        return False
            elif not np.allclose(getattr(a, req_prop), getattr(b, req_prop)):
                if print_:
                    error_log(f"a.{req_prop}", getattr(a, req_prop))
                    error_log(f"b.{req_prop}", getattr(b, req_prop))
                return False
        return True

    return _are_splines_equal


@pytest.fixture
def are_items_close():
    def _are_items_close(a, b):
        """returns True if items in a and b are close"""
        all_close = True

        for i, (aa, bb) in enumerate(zip(a, b)):
            if not all(np.isclose(aa, bb)):
                # print to inform
                error_log(f"elements in index-{i} are not close")
                error_log(f"  from first: {aa}")
                error_log(f"  from second: {bb}")

                all_close = False

        return all_close

    return _are_items_close


@pytest.fixture
def are_items_same():
    def _are_items_same(a, b):
        """returns True if items in a and b are same"""
        all_same = True

        for i, (aa, bb) in enumerate(zip(a, b)):
            if aa != bb:
                # print to inform
                error_log(f"element in index-{i} are not same")
                error_log(f"  from first: {aa}")
                error_log(f"  from second: {bb}")

                all_same = False

        return all_same

    return _are_items_same


@pytest.fixture
def are_stripped_lines_same():
    def _are_stripped_lines_same(a, b, ignore_order=False):
        """returns True if items in a and b same, preceding and tailing whitespaces
        are ignored and strings are joined"""
        all_same = True

        for i, (line_a, line_b) in enumerate(zip(a, b)):
            # check stripped string
            stripped_a, stripped_b = line_a.strip(), line_b.strip()

            # print general info
            if stripped_a != stripped_b:
                error_log(f"stripped line at index-{i} are not the same")
                error_log(f"  from first: {line_a}")
                error_log(f"  from second: {line_b}")

            # give one more chance if ignore_order
            if stripped_a != stripped_b and ignore_order:
                error_log("  checking again, while ignoring word order:")

                # This is meant for attributes
                delimiters = r" |\>|\<|\t|,"
                splitted_a = list(
                    filter(None, re.split(delimiters, stripped_a))
                )
                splitted_b = list(
                    filter(None, re.split(delimiters, stripped_b))
                )
                # first, len check
                len_a, len_b = len(splitted_a), len(splitted_b)
                if len(splitted_a) != len(splitted_b):
                    error_log(
                        f"    different word counts: a-{len_a}, b-{len_b}"
                    )
                    all_same = False
                else:
                    # word order
                    a_to_b = []
                    nums_b = None
                    for word_a in splitted_a:
                        try:
                            a_to_b.append(splitted_b.index(word_a))
                        except BaseException:
                            try:
                                num_a = float(word_a)
                            except ValueError:
                                pass
                            else:
                                if nums_b is None:
                                    nums_b = []
                                    for idx, num_b in enumerate(splitted_b):
                                        with contextlib.suppress(ValueError):
                                            nums_b.append((idx, float(num_b)))
                                for idx, num_b in nums_b:
                                    if np.isclose(num_a, num_b):
                                        a_to_b.append(idx)
                                        break
                                else:
                                    error_log(
                                        f"    second does not contain ({word_a})"
                                    )
                                    all_same = False

        return all_same

    return _are_stripped_lines_same
