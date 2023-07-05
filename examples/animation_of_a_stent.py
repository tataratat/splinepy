import gustaf as gus
import imageio
import numpy as np
import vedo

import splinepy as sp

# ------------------
# Animation settings
# ------------------

# Global parameters
fps = 15
vedo.settings.screenshot_transparent_background = 1
export_resolution = [1440, 1088]
sample_resolution = 6
duration_per_scene = 2
n_frames_per_scene = int(duration_per_scene * fps)
increment_per_frame = 1 / (n_frames_per_scene - 1)

# Initialize writer
write_ms = imageio.get_writer(
    "animation_stent.mp4",
    fps=fps,
)

# Camera settings
camera_1 = dict(
    position=(5.0, 0.7, -0.1),
    focal_point=(0.2, 0.5, -0.1),
    viewup=(0, 1.00000, 0),
    distance=5.0,
    clipping_range=(4.74874, 5.48835),
)
camera_2 = dict(
    position=(3, 0.6, 0.5),
    focal_point=(0, 0.5, 0.50),
    viewup=(0.0, 1.0, 0.0),
    distance=2.73205,
    clipping_range=(2.55757, 2.95828),
)
camera_3 = dict(
    position=(8.2, 7.3, -0.3),
    focal_point=(-0.448800, 1.59083, 0.0400746),
    viewup=(-0.545630, 0.837155, -0.0381893),
    distance=10.3750,
    clipping_range=(5.78697, 14.8182),
)


# Helper function to interpolate smoothly between two camera settings
def interpolate_cameras(unravel_factor, start, end):
    cam = dict()
    for key in ["position", "focal_point", "viewup", "clipping_range"]:
        cam[key] = tuple(
            [
                (1 - unravel_factor) * p + unravel_factor * pp
                for p, pp in zip(start[key], end[key])
            ]
        )
    cam["distance"] = (1 - unravel_factor) * start[
        "distance"
    ] + unravel_factor * end["distance"]
    return cam


# -------------------
# Generate Geometries
# -------------------


# Create the basis tile
def Microtile(midline_position):
    return [
        sp.Bezier(
            degrees=[3],
            control_points=[
                [0.0, 0.0],
                [0.0, midline_position],
                [0.5, midline_position],
                [0.5, 0.5],
            ],
        ),
        sp.Bezier(
            degrees=[3],
            control_points=[
                [1.0, 0.0],
                [1.0, midline_position],
                [0.5, midline_position],
                [0.5, 0.5],
            ],
        ),
        sp.Bezier(
            degrees=[3],
            control_points=[
                [0.5, 0.5],
                [0.5, 1 - midline_position],
                [0.0, 1 - midline_position],
                [0.0, 1.0],
            ],
        ),
        sp.Bezier(
            degrees=[3],
            control_points=[
                [0.5, 0.5],
                [0.5, 1 - midline_position],
                [1.0, 1 - midline_position],
                [1.0, 1.0],
            ],
        ),
    ]


# Create Macrospline
height = 4
radius = 1
tiling = [8, 8]
macro_spline = sp.Bezier(
    degrees=[1], control_points=[[radius, 0, 0], [radius, height, 0]]
).nurbs.create.revolved(
    axis=[0, 1, 0],
    center=[0, 0, 0],
    angle=359.9,
    degree=True,
    n_knot_spans=tiling[0],
)
macro_spline.normalize_knot_vectors()
macro_spline.insert_knots(0, [i / tiling[1] for i in range(1, tiling[1])])

# Unraveld spline
mac_spline_unraveled = sp.NURBS(**macro_spline.todict())
mac_spline_unraveled.control_points[:, 2] = np.arccos(
    macro_spline.cps[:, 0]
    / np.linalg.norm(macro_spline.cps[:, [True, False, True]], axis=1)
)
mac_spline_unraveled.control_points[
    (macro_spline.control_points[:, 2] - 1e-12 > 0), 2
] = (
    2 * np.pi
    - mac_spline_unraveled.control_points[
        (macro_spline.control_points[:, 2] - 1e-12 > 0), 2
    ]
)
mac_spline_unraveled.control_points[:, 0] = 0.0
mac_spline_unraveled.cps[:, 2] /= np.max(mac_spline_unraveled.cps[:, 2])
mac_spline_unraveled.cps[:, 1] /= np.max(mac_spline_unraveled.cps[:, 1])


# Unravel function to interpolate splines
def interpolate_between_splines(unravel_factor, target, start):
    """
    Unravel factor must be (0,1) interpolatest between unraveld and real
    geometry
    """
    macro_spline_c = type(target)(**target.todict())
    macro_spline_c.cps = (
        unravel_factor * target.cps + (1 - unravel_factor) * start.cps
    )
    return macro_spline_c


# Unit Square representation for sliding effect
unit_cube_representation = mac_spline_unraveled.extract.beziers()[0]
unit_cube_mapped = unit_cube_representation.copy()
unit_cube_mapped.cps[:, [1, 2]] /= np.max(
    unit_cube_mapped.cps[:, [1, 2]], axis=0
)
unit_cube_mapped.cps += [0, 0, -1.3]

# Scene 1 - Show setup
unit_cube_mapped.show_options["alpha"] = 0.1
unit_cube_mapped.show_options["knots"] = False
unit_cube_mapped.show_options["c"] = "lightgray"
mac_spline_unraveled.show_options["knots"] = True
mac_spline_unraveled.show_options["c"] = "lightgray"


microtile_lines = Microtile(0.25)
microtile = [unit_cube_mapped.compose(line) for line in microtile_lines]
for i in range(n_frames_per_scene):
    plt = gus.show(
        [unit_cube_mapped, mac_spline_unraveled, *microtile],
        knots=True,
        control_points=False,
        # c="lightgray",
        lighting="off",
        offscreen=True,
        resolution=sample_resolution,
        size=export_resolution,
        cam=camera_1,
        azimuth=i * increment_per_frame,
    )
    write_ms.append_data(plt.screenshot(asarray=True))

# Scene 2 - Move Microtile into parametric domain
for i in range(n_frames_per_scene):
    current = interpolate_between_splines(
        increment_per_frame * i,
        unit_cube_representation,
        unit_cube_mapped,
    )
    microtile = [current.compose(line) for line in microtile_lines]
    plt = gus.show(
        [unit_cube_mapped, mac_spline_unraveled, *microtile],
        knots=True,
        control_points=False,
        # c="lightgray",
        lighting="off",
        offscreen=True,
        resolution=sample_resolution,
        size=export_resolution,
        cam=camera_1,
    )
    write_ms.append_data(plt.screenshot(asarray=True))

# Scene 3 - shift perspective
microtile = [
    unit_cube_representation.compose(line) for line in microtile_lines
]
for i in range(n_frames_per_scene):
    unit_cube_mapped.show_options["alpha"] = 0.1 * (n_frames_per_scene - i)
    plt = gus.show(
        [unit_cube_mapped, mac_spline_unraveled, *microtile],
        knots=True,
        control_points=False,
        # c="lightgray",
        lighting="off",
        offscreen=True,
        resolution=sample_resolution,
        size=export_resolution,
        cam=interpolate_cameras(i * increment_per_frame, camera_1, camera_2),
    )
    write_ms.append_data(plt.screenshot(asarray=True))

# Scene 4 - build MS
generator = sp.microstructure.Microstructure()
generator.microtile = microtile_lines
generator.deformation_function = mac_spline_unraveled
generator.tiling = [1, 1]
ms = generator.create()

for i in range(n_frames_per_scene):
    alpha_h = (i * increment_per_frame) ** 2
    for m in ms:
        m.show_options["alpha"] = alpha_h
    plt = gus.show(
        [mac_spline_unraveled, *microtile, *ms],
        knots=True,
        control_points=False,
        # c="lightgray",
        lighting="off",
        offscreen=True,
        resolution=sample_resolution,
        size=export_resolution,
        cam=camera_2,
    )
    write_ms.append_data(plt.screenshot(asarray=True))

# Scene 5 - Show MS in parametric view
plt = gus.show(
    [*ms],
    knots=True,
    control_points=False,
    # c="lightgray",
    lighting="off",
    offscreen=True,
    resolution=sample_resolution,
    size=export_resolution,
    cam=camera_2,
)
for i in range(n_frames_per_scene):
    write_ms.append_data(plt.screenshot(asarray=True))

# Scene 5 - Map Para->Phys with MS
for i in range(n_frames_per_scene):
    generator.deformation_function = interpolate_between_splines(
        i * increment_per_frame, macro_spline, mac_spline_unraveled
    )
    plt = gus.show(
        [*generator.create()],
        knots=True,
        control_points=False,
        # c="lightgray",
        lighting="off",
        offscreen=True,
        resolution=sample_resolution,
        size=export_resolution,
        cam=interpolate_cameras(i * increment_per_frame, camera_2, camera_3),
    )
    write_ms.append_data(plt.screenshot(asarray=True))


# Scene 6 - rotate results
generator.deformation_function = macro_spline
ms = generator.create()

ss, cc = np.sin(0.05), np.cos(0.05)
rot = np.array([[cc, 0, ss], [0, 1, 0], [-ss, 0, cc]])

# Rotate the resulting cylinder
for i in range(n_frames_per_scene * 3):
    camera_3["position"] = tuple(np.matmul(rot, camera_3["position"]).tolist())
    camera_3["viewup"] = tuple(np.matmul(rot, camera_3["viewup"]).tolist())
    camera_3["focal_point"] = tuple(
        np.matmul(rot, camera_3["focal_point"]).tolist()
    )
    plt = gus.show(
        [*ms],
        knots=True,
        control_points=False,
        # c="lightgray",
        lighting="off",
        offscreen=True,
        resolution=sample_resolution,
        size=export_resolution,
        cam=camera_3,
    )
    write_ms.append_data(plt.screenshot(asarray=True))

# Scene 6a - Continue rotation, but modify
increment_per_frame = 0.2 / n_frames_per_scene
for i in range(n_frames_per_scene):
    camera_3["position"] = tuple(np.matmul(rot, camera_3["position"]).tolist())
    camera_3["viewup"] = tuple(np.matmul(rot, camera_3["viewup"]).tolist())
    camera_3["focal_point"] = tuple(
        np.matmul(rot, camera_3["focal_point"]).tolist()
    )

    generator.microtile = Microtile(0.25 + i * increment_per_frame)
    ms = generator.create()
    plt = gus.show(
        [*ms],
        knots=True,
        control_points=False,
        # c="lightgray",
        lighting="off",
        offscreen=True,
        resolution=sample_resolution,
        size=export_resolution,
        cam=camera_3,
    )
    write_ms.append_data(plt.screenshot(asarray=True))

# Scene 6b - Continue rotation, but modify
for i in range(n_frames_per_scene * 2):
    camera_3["position"] = tuple(np.matmul(rot, camera_3["position"]).tolist())
    camera_3["viewup"] = tuple(np.matmul(rot, camera_3["viewup"]).tolist())
    camera_3["focal_point"] = tuple(
        np.matmul(rot, camera_3["focal_point"]).tolist()
    )
    generator.microtile = Microtile(0.45 - i * increment_per_frame)
    ms = generator.create()
    plt = gus.show(
        [*ms],
        knots=True,
        control_points=False,
        # c="lightgray",
        lighting="off",
        offscreen=True,
        resolution=sample_resolution,
        size=export_resolution,
        cam=camera_3,
    )
    write_ms.append_data(plt.screenshot(asarray=True))

# Scene 6c - Continue rotation, but modify
increment_per_frame = 0.2 / n_frames_per_scene
for i in range(n_frames_per_scene):
    camera_3["position"] = tuple(np.matmul(rot, camera_3["position"]).tolist())
    camera_3["viewup"] = tuple(np.matmul(rot, camera_3["viewup"]).tolist())
    camera_3["focal_point"] = tuple(
        np.matmul(rot, camera_3["focal_point"]).tolist()
    )

    generator.microtile = Microtile(0.05 + i * increment_per_frame)
    ms = generator.create()
    plt = gus.show(
        [*ms],
        knots=True,
        control_points=False,
        # c="lightgray",
        lighting="off",
        offscreen=True,
        resolution=sample_resolution,
        size=export_resolution,
        cam=camera_3,
    )
    write_ms.append_data(plt.screenshot(asarray=True))
