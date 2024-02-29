"""Create gus"""

import gustaf as gus

import splinepy


def example(
    with_pupiles: bool = False,
    with_dark_background: bool = False,
    dark_with_white_oval_background: bool = False,
    with_alpha_background: bool = False,
):
    if with_alpha_background:
        import vedo

        vedo.settings.screenshot_transparent_background = 1

    l_eye = splinepy.helpme.create.disk(0.15)
    l_eye.cps += [0.6, 0.2]
    l_eye.show_options["c"] = "black"
    l_eye.show_options["control_points"] = False
    l_eye.show_options["knots"] = False

    r_eye = splinepy.helpme.create.disk(0.15)
    r_eye.cps += [-0.6, 0.2]
    r_eye.show_options["c"] = "black"
    r_eye.show_options["control_points"] = False
    r_eye.show_options["knots"] = False

    x_offset = 0.59
    y_offset = 0.16

    l_pupil = splinepy.helpme.create.disk(0.05)
    l_pupil.cps += [x_offset, y_offset]
    l_pupil.show_options["c"] = "white"
    l_pupil.show_options["control_points"] = False
    l_pupil.show_options["knots"] = False

    r_pupil = splinepy.helpme.create.disk(0.05)
    r_pupil.cps += [-x_offset, y_offset]
    r_pupil.show_options["c"] = "white"
    r_pupil.show_options["control_points"] = False
    r_pupil.show_options["knots"] = False

    upper_lip = splinepy.Bezier(
        [5, 1],
        [
            [-0.8, -0.1],
            [-0.2, -0.4],
            [-0.5, 0.5],
            [0.5, 0.5],
            [0.2, -0.4],
            [0.8, -0.1],
            [-0.75, -0.15],
            [-0.2, -0.55],
            [-0.5, 0.05],
            [0.5, 0.05],
            [0.2, -0.55],
            [0.75, -0.15],
        ],
    )
    upper_lip.show_options["c"] = "orange7"
    upper_lip.show_options["control_points"] = False
    upper_lip.show_options["knots"] = False

    inner_mouth = splinepy.Bezier(
        [5, 1],
        [
            *upper_lip.cps[6:],
            [-0.75, -0.15],
            [-0.6, -0.325],
            [-0.4, -0.525],
            [0.4, -0.525],
            [0.6, -0.325],
            [0.75, -0.15],
        ],
    )
    inner_mouth.show_options["c"] = "orange5"
    inner_mouth.show_options["control_points"] = False
    inner_mouth.show_options["knots"] = False

    lower_lip = splinepy.Bezier(
        [5, 1],
        [
            *inner_mouth.cps[6:],
            [-0.8, -0.1],
            [-0.6, -0.4],
            [-0.4, -0.6],
            [0.4, -0.6],
            [0.6, -0.4],
            [0.8, -0.1],
        ],
    )
    lower_lip.show_options["c"] = "orange7"
    lower_lip.show_options["control_points"] = False
    lower_lip.show_options["knots"] = False

    item_to_show = [l_eye, r_eye, upper_lip, inner_mouth, lower_lip]
    if with_pupiles:
        item_to_show += [l_pupil, r_pupil]

    gus.show(
        item_to_show,
        lighting="off",
        background="white",
        close=False,
    )

    if with_dark_background:
        background = []

        if dark_with_white_oval_background:
            height = 0.6
            width = 0

            x_shift = 0.3
            y_shift = -0.02

            background = [
                splinepy.helpme.create.disk(height),
                splinepy.helpme.create.disk(height),
                splinepy.helpme.create.box(width + 2 * x_shift, 2 * height),
            ]
            for spline in background:
                spline.show_options["control_points"] = False
                spline.show_options["knots"] = False

            background[0].cps += [-x_shift, y_shift]
            background[1].cps += [width + x_shift, y_shift]
            background[2].cps += [-x_shift, -height + y_shift]
            for elem in background:
                elem.show_options["c"] = "white"
                elem.show_options["alpha"] = 0.25
        else:
            l_eye.show_options["c"] = "white"
            r_eye.show_options["c"] = "white"
            l_pupil.show_options["c"] = "black"
            r_pupil.show_options["c"] = "black"

        gus.show(
            [*background, *item_to_show],
            lighting="off",
            background="black",
            close=False,
        )


if __name__ == "__main__":
    example()
