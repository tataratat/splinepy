import unittest

import numpy as np

import splinepy

__all__ = [
        "unittest",
        "np",
        "splinepy",
]

# abbreviation
# z: bezier
# r: rational bezier
# b: bspline
# n: nurbs

# initializing a spline should be a test itself, so provide `dict_spline`
# this is "iga-book"'s fig 2.15.
b2P2D = dict(
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
n2P2D = dict(
        degrees=[2, 1],
        knot_vectors=[
                [0, 0, 0, 1, 1, 1],
                [0, 0, 1, 1],
        ],
        control_points=[
                [-1., 0.],
                [-1., 1.],
                [0., 1.],
                [-2., 0.],
                [-2., 2.],
                [0., 2.],
        ],
        weights=[
                [1.],
                [2**-.5],
                [1.],
                [1.],
                [2**-.5],
                [1.],
        ],
)

#
z2P2D = dict(degrees=n2P2D["degrees"], control_points=n2P2D["control_points"])

#
r2P2D = dict(
        degrees=n2P2D["degrees"],
        control_points=n2P2D["control_points"],
        weights=n2P2D["weights"]
)

# 3D
z3P3D = dict(
        degrees=[1, 1, 1],
        control_points=[
                [0., 0., 0.],
                [1., 0., 0.],
                [0., 1., 0.],
                [1., 1., 0.],
                [0., -1., 1.],
                [1., 0., 1.],
                [-1., 1., 2.],
                [2., 2., 2.],
        ],
)

r3P3D = dict(
        **z3P3D,
        weights=[1.] * len(z3P3D["control_points"]),
)

b3P3D = dict(
        **z3P3D,
        knot_vectors=[
                [0., 0., 1., 1.],
                [0., 0., 1., 1.],
                [0., 0., 1., 1.],
        ],
)

n3P3D = dict(
        **r3P3D,
        knot_vectors=b3P3D["knot_vectors"],
)

# query points
q2D = [
        [.01, .01],
        [.01, .5],
        [.9, .1],
        [.8, .7],
        [.4, .99],
]

q3D = [
        [.1, .1, .1],
        [.734, .525, .143],
        [.9666, .991, .003],
        [.5623, .0089, .99],
        [.0431, .2, .523],
]


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
        this_is_same = (aa == bb)
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
        this_is_same = (stripped_a == stripped_b)

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
