import os
import unittest

import gustaf as gus
import numpy as np

import splinepy

__all__ = [
    "unittest",
    "np",
    "splinepy",
    "gus",
]

# abbreviation
# z: bezier
# r: rational bezier
# b: bspline
# n: nurbs


# initializing a spline should be a test itself, so provide `dict_spline`
# this is "iga-book"'s fig 2.15.
def b2p2d():
    return dict(
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


b2P2D = b2p2d()


# half-half circle.
def n2p2d():
    return dict(
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


n2P2D = n2p2d()


def n2p2d_quarter_circle():
    """explicit function for quarter circle
    incase n2p2d changes in the future..."""
    return n2p2d()


#
def z2p2d():
    return dict(
        degrees=n2P2D["degrees"], control_points=n2P2D["control_points"]
    )


z2P2D = z2p2d()


#
def r2p2d():
    return dict(
        degrees=n2P2D["degrees"],
        control_points=n2P2D["control_points"],
        weights=n2P2D["weights"],
    )


r2P2D = r2p2d()


def all2p2d():
    z = splinepy.Bezier(**z2p2d())
    r = splinepy.RationalBezier(**r2p2d())
    b = splinepy.BSpline(**b2p2d())
    n = splinepy.NURBS(**n2p2d())
    return [z, r, b, n]


# 3D
def z3p3d():
    return dict(
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


z3P3D = z3p3d()


def r3p3d():
    return dict(
        **z3P3D,
        weights=[1.0] * len(z3P3D["control_points"]),
    )


r3P3D = r3p3d()


def b3p3d():
    return dict(
        **z3P3D,
        knot_vectors=[
            [0.0, 0.0, 1.0, 1.0],
            [0.0, 0.0, 1.0, 1.0],
            [0.0, 0.0, 1.0, 1.0],
        ],
    )


b3P3D = b3p3d()


def n3p3d():
    return dict(
        **r3P3D,
        knot_vectors=b3P3D["knot_vectors"],
    )


n3P3D = n3p3d()


def all3p3d():
    z = splinepy.Bezier(**z3p3d())
    r = splinepy.RationalBezier(**r3p3d())
    b = splinepy.BSpline(**b3p3d())
    n = splinepy.NURBS(**n3p3d())
    return [z, r, b, n]


# query points
q2D = [
    [0.01, 0.01],
    [0.01, 0.5],
    [0.9, 0.1],
    [0.8, 0.7],
    [0.4, 0.99],
]

q3D = [
    [0.1, 0.1, 0.1],
    [0.734, 0.525, 0.143],
    [0.9666, 0.991, 0.003],
    [0.5623, 0.0089, 0.99],
    [0.0431, 0.2, 0.523],
]


def raster(bounds, resolutions):
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


def nd_box(dim):
    """creates simple box in nd"""
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


def to_tmpf(tmpd):
    """given tmpd, returns tmpf"""
    return os.path.join(tmpd, "nqv248p90")


def are_splines_equal(a, b):
    """returns True if Splines are equivalent"""
    if not a.whatami == b.whatami:
        return False
    for req_prop in a.required_properties:
        if req_prop == "knot_vectors":
            for aa, bb in zip(a.knot_vectors, b.knot_vectors):
                if not np.allclose(aa.numpy(), bb.numpy()):
                    return False
        else:
            if not np.allclose(getattr(a, req_prop), getattr(b, req_prop)):
                return False
    return True


def are_items_close(a, b):
    """returns True if items in a and b are close"""
    all_close = True

    for i, (aa, bb) in enumerate(zip(a, b)):
        this_is_close = all(np.isclose(aa, bb))
        if not this_is_close:
            # print to inform
            print(f"elements in index-{i} are not close")
            print(f"  from first: {aa}")
            print(f"  from second: {bb}")

            all_close = False

    return all_close


def are_items_same(a, b):
    """returns True if items in a and b are same"""
    all_same = True

    for i, (aa, bb) in enumerate(zip(a, b)):
        this_is_same = aa == bb
        if not this_is_same:
            # print to inform
            print(f"element in index-{i} are not same")
            print(f"  from first: {aa}")
            print(f"  from second: {bb}")

            all_same = False

    return all_same


def are_stripped_lines_same(a, b, ignore_order=False):
    """returns True if items in a and b same, preceding and tailing whitespaces
    are ignored and strings are joined"""
    all_same = True

    for i, (line_a, line_b) in enumerate(zip(a, b)):
        # check stripped string
        stripped_a, stripped_b = line_a.strip(), line_b.strip()
        this_is_same = stripped_a == stripped_b

        # print general info
        if not this_is_same:
            print(f"stripped line at index-{i} are not the same")
            print(f"  from first: {line_a}")
            print(f"  from second: {line_b}")

        # give one more chance if ignore_order
        if not this_is_same and ignore_order:
            print("  checking again, while ignoring word order:")

            splitted_a, splitted_b = stripped_a.split(), stripped_b.split()

            # first, len check
            len_a, len_b = len(splitted_a), len(splitted_b)
            if len(splitted_a) != len(splitted_b):
                print(f"    different word counts: a-{len_a}, b-{len_b}")
                all_same = False
            else:
                # word order
                a_to_b = list()
                for word_a in splitted_a:
                    try:
                        a_to_b.append(splitted_b.index(word_a))
                    except BaseException:
                        print(f"    second does not contain ({word_a})")
                        all_same = False

    return all_same
