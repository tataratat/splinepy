import argparse

import gustaf as gus
import numpy as np
import vedo

import splinepy

# this will cause a deprecated warning for newer version of vedo.
# for backward compatibility, we keep this one.
vedo.settings.screenshot_transparent_background = 1

if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        description=(
            "Creates a fanfare. Default option corresponds to "
            "the tataratat logo."
        )
    )
    arg_parser.add_argument(
        "--r_mouth",
        dest="r_mouth",
        default=0.2,
        help="Mouth radius for mouth piece.",
    )
    arg_parser.add_argument(
        "--r_shaft",
        dest="r_shaft",
        default=0.1,
        help="Shaft radius for mouth piece.",
    )
    arg_parser.add_argument(
        "--r_rounding",
        dest="r_rounding",
        default=0.05,
        help="Radius offset for mouth piece rounding.",
    )
    arg_parser.add_argument(
        "--l_mouth",
        dest="l_mouth",
        default=0.4,
        help="Length of the mouth piece.",
    )
    arg_parser.add_argument(
        "--l_cone",
        dest="l_cone",
        default=1.2,
        help="Length of the cone.",
    )
    arg_parser.add_argument(
        "--r_cone",
        dest="r_cone",
        default=0.8,
        help="Radius of the cone.",
    )
    arg_parser.add_argument(
        "--l_shaft",
        dest="l_shaft",
        default=4.0,
        help="Length of the shaft.",
    )
    arg_parser.add_argument(
        "--r_curve",
        dest="r_curve",
        default=0.3,
        help="Radius of connecting curves.",
    )
    args = arg_parser.parse_args()

    # Let's draw a fanfare
    fanfare = []

    # cast all args to float
    # Mouthpiece
    r_mouth = float(args.r_mouth)
    r_shaft = float(args.r_shaft)
    r_rounding = float(args.r_rounding)
    l_mouth = float(args.l_mouth)
    # Cone
    l_cone = float(args.l_cone)
    r_cone = float(args.r_cone)
    # Shaft
    l_shaft = float(args.l_shaft)
    curved_radius = float(args.r_curve)

    # Create fanfare part by part
    mouth_piece = splinepy.BSpline(
        knot_vectors=[[0, 0, 0, 0.25, 0.5, 0.5, 0.75, 1, 1, 1]],
        degrees=[2],
        control_points=[
            [l_mouth, r_mouth - r_rounding],
            [l_mouth + r_rounding, r_mouth - r_rounding],
            [l_mouth + r_rounding, r_mouth],
            [l_mouth, r_mouth],
            [0.5 * l_mouth, r_mouth],
            [0.5 * l_mouth, r_shaft],
            [0.0, r_shaft],
        ],
    ).nurbs.create.revolved(
        axis=[1, 0, 0],
        center=[0, 0, 0],
        angle=360,
        degree=True,
        n_knot_spans=4,
    )

    isqrt2 = 2 ** (-0.5)
    cross_section = splinepy.NURBS(
        degrees=[2],
        knot_vectors=[[0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1]],
        control_points=(
            np.array(
                [
                    [0.0, 0, 0],
                    [0.0, 0.5, -0.5],
                    [0.0, 1, 0],
                    [0.0, 1 + 0.5, 0.5],
                    [0.0, 1, 1],
                    [0.0, 0.5, 1 + 0.5],
                    [0.0, 0, 1],
                    [0.0, -0.5, 0.5],
                    [0, 0, 0],
                ]
            )
            + np.array([0, -0.5, -0.5])
        )
        * (r_shaft * 2 * isqrt2),
        weights=[1, isqrt2, 1, isqrt2, 1, isqrt2, 1, isqrt2, 1],
    )

    shaft = cross_section.create.extruded(extrusion_vector=[-l_shaft, 0, 0])

    lower_shaft = shaft.copy()
    lower_shaft.control_points += [0, -2 * curved_radius, 0]
    upper_shaft = shaft.copy()
    upper_shaft.control_points += [
        0,
        -2 * curved_radius + 2 * isqrt2 * curved_radius,
        -2 * isqrt2 * curved_radius,
    ]

    cross_section.control_points += [-l_shaft, 0, 0]

    curved_piece = cross_section.nurbs.create.revolved(
        axis=[0.0, 0, 1], center=[-l_shaft, -curved_radius, 0], angle=180
    )
    cross_section.control_points += [l_shaft, -2 * curved_radius, 0]
    second_curved_piece = cross_section.nurbs.create.revolved(
        axis=[0, 1.0, 1.0],
        center=[
            0,
            -2 * curved_radius + isqrt2 * curved_radius,
            -isqrt2 * curved_radius,
        ],
        angle=180,
    )

    cone_part = splinepy.BSpline(
        degrees=[2],
        knot_vectors=[[0, 0, 0, 0.5, 1, 1, 1]],
        control_points=[
            [0.0, r_shaft, 0.0],
            [-0.25 * l_cone, r_shaft, 0.0],
            [-l_cone, 0.5 * (r_cone + r_shaft), 0.0],
            [-l_cone, r_cone, 0.0],
        ],
    ).nurbs.create.revolved(
        axis=[1, 0, 0],
        center=[0, 0, 0],
        angle=360,
        degree=True,
    )
    cone_part.control_points += [
        -l_shaft,
        -2 * curved_radius + 2 * isqrt2 * curved_radius,
        -2 * isqrt2 * curved_radius,
    ]

    fanfare.append(mouth_piece)
    fanfare.append(curved_piece)
    fanfare.append(second_curved_piece)
    fanfare.append(upper_shaft)
    fanfare.append(cone_part)
    fanfare.append(shaft)
    fanfare.append(lower_shaft)

    # show fanfare - prepare only surface and full spline components
    face_showables = []
    spline_showables = []
    for f in fanfare:
        showables = f.showable(resolutions=500, c=(255, 215, 0))
        faces = showables["spline"]
        face_showables.append(faces)
        spline_showables.extend(list(showables.values()))

    cam = {
        "pos": (-11.90193, 5.995478, 3.057389),
        "focalPoint": (-2.817062, -0.6220471, -0.3948362),
        "viewup": (0.3271027, 0.7493888, -0.5756911),
        "distance": 11.75773,
        "clippingRange": (6.030973, 19.21880),
    }

    gus.show(face_showables, spline_showables, cam=cam)
