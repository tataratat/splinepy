import gustaf as gus
import numpy as np

import splinepy


def bspline2p3d():
    """
    Creates new bspline of of parametric dimension 2, physical dimension 3
    """
    ds2 = [2, 2]

    # define knot vectors
    kvs2 = [
        [0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
        [0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
    ]

    # define control points
    cps2 = np.array(
        [
            [0, 0, 0],
            [0, 1, 0],
            [1, 1.5, 0],
            [3, 1.5, 0],
            [-1, 0, 0],
            [-1, 2, 0],
            [1, 4, 0],
            [3, 4, 0],
            [-2, 0, 0],
            [-2, 2, 0],
            [1, 5, 0],
            [3, 5, 2],
        ]
    )

    # init bspline
    return splinepy.BSpline(
        degrees=ds2,
        knot_vectors=kvs2,
        control_points=cps2,
    )


if __name__ == "__main__":
    # turn on debug logs
    # gus.utils.log.configure(debug=True)

    # define spline_data
    # 1. see coordinates's norm
    b = bspline2p3d()
    b.spline_data["me"] = b
    b.show_options["data"] = "me"
    gus.show(["1. Show norm of coordinates.", b])

    # 1.1 default scalarbar
    b = bspline2p3d()
    b.spline_data["me"] = b
    b.show_options["data"] = "me"
    b.show_options["scalarbar"] = True
    gus.show(["1.1. Show 1. plus scalarbar", b])

    # 2. see coordinate's norm and as arrow
    b = bspline2p3d()
    b.spline_data["me"] = b
    b.show_options["data"] = "me"
    b.show_options["arrow_data"] = "me"
    gus.show(
        ["2. Show coordinate norm as scalar field and coordinate as arrows", b]
    )

    # 3. see coordinates norm and as arrows only on specified places
    b = bspline2p3d()
    b.spline_data["me"] = b
    b.show_options["data"] = "me"
    b.show_options["arrow_data"] = "me"
    b.show_options["arrow_data_on"] = np.random.random((100, 2))  # para_coords
    gus.show(
        ["3. Show coordinates norm and as arrows on 100 random points.", b]
    )

    # 4. see 3. with parametric_view
    b = bspline2p3d()
    b.spline_data["me"] = b
    b.show_options["data"] = "me"
    b.show_options["arrow_data"] = "me"
    b.show_options["arrow_data_on"] = np.random.random((100, 2))  # para_coords
    b.show_options["scalarbar"] = True
    gus.show(
        ["4. Show 3. and 3. in parametric space view with scalarbar.", b],
        b.create.parametric_view(),
    )
    # b.show(parametric_space=True)  # does the same

    # 5. plot function defines how to evaluate at given parametric location
    def plot_func(data, on):
        """
        callback to evaluate derivatives
        """
        return data.derivative(on, [0, 1])

    # use adaptor to defie a data
    plot_func_data = splinepy.SplineDataAdaptor(b, function=plot_func)
    # the rest is the same
    b = bspline2p3d()
    b.spline_data["der01"] = plot_func_data
    b.show_options["arrow_data"] = "der01"
    gus.show(
        [
            "5. Show partial derivative of seconda parametric dimension.\n"
            "This uses a callback function and SplineDataAdaptor.\n"
            "Also achievable with derivative spline.",
            b,
        ]
    )

    # remove on to sample same way as spline.
    # however, gold
    b = bspline2p3d()
    b.spline_data["der01"] = plot_func_data
    b.show_options["arrow_data"] = "der01"
    b.show_options["arrow_data_on"] = np.random.random((100, 2))  # para_coords
    b.show_options["arrow_data_color"] = "gold"
    gus.show(
        [
            "5.1. Same data as 5. However, on 100 random places with gold "
            "arrows",
            b,
        ]
    )

    # 6. plot on predefined place - this will be only available as arrow data
    locations = np.repeat(np.linspace(0, 1, 15), 2).reshape(-1, 2)
    values = np.repeat(np.linspace(1, 2, 15), 3).reshape(-1, 3) + [0, 0, 2]
    fixed_data = splinepy.SplineDataAdaptor(values, locations=locations)
    b = bspline2p3d()
    b.spline_data["fixed"] = fixed_data
    b.show_options["arrow_data"] = "fixed"
    gus.show(
        [
            "6. Show arbitrary array data on predefined locations using "
            "SplineDataAdaptor.",
            b,
        ]
    )

    # fixed location data can't be shown in other requested locations.
    # following won't work
    # b.show_options["arrow_data_on"] = locations[:5]
    # b.show()

    # 7. plot any data with a function
    # some manually defined deformed spline
    deformed = bspline2p3d()  # minimal copy - properties and cached data only
    deformed.cps[11, -1] -= 4
    deformed.cps *= [5.5, 5.5, 5.5]
    deformed.cps += [-5, 0, 8]
    deformed.cps[0, [0, -1]] += 4
    deformed.show_options["c"] = "hotpink"

    # define callback
    def func(self_and_deformed, on):
        """
        callback to sample displacements.
        """
        # unpack data
        self, deformed = self_and_deformed
        return deformed.evaluate(on) - self.evaluate(on)

    # use adaptor - data is used as the first arg for callback
    deformed_data = splinepy.SplineDataAdaptor(
        data=(b, deformed),
        function=func,
    )
    b.spline_data["deformed"] = deformed_data
    b.show_options["arrow_data"] = "deformed"
    # arrows are always automatically scaled. for this one, let's not
    b.show_options["arrow_data_scale"] = 1
    b.show_options["arrow_data_on"] = locations
    # let's see in parametric space
    p_view = b.create.parametric_view()  # shallow copies data and options
    p_view.show_options.pop("arrow_data_scale")  # we want automatic scaling
    p_view.show_options["arrow_data_color"] = "gold"
    p_view.show_options["data"] = "deformed"
    # plot side by side
    gus.show(
        [
            "7. Original spline, deformed spline\n"
            "and pointwise correspondence on selected locations.",
            b,
            deformed,
        ],
        [
            "Parametric view of displacements.\n"
            "Arrows are directly transferred with automatic rescaling.\n"
            "(no inverse mapping)",
            p_view,
        ],
    )
    # say, deformed has changed again - plots should adapt automatically
    deformed.cps[:, [1, 2]] = deformed.cps[:, [2, 1]]
    gus.show(
        [
            "7.1. Same setup as 7. with inplace changes in deformed spline.",
            b,
            deformed,
        ],
        [
            "Parametric view of displacements.\n"
            "Arrows are directly transferred with automatic rescaling.\n"
            "(no inverse mapping)",
            p_view,
        ],
    )

    # 8. fixed location data that uses callback
    # predefind some locations
    bottom = splinepy.helpme.create.arc(
        radius=0.5, angle=-180
    )  # zero centered
    bottom.cps += [0.5, 0.55]
    circle1 = splinepy.helpme.create.circle(radius=0.1)
    circle2 = circle1.copy()
    circle1.cps += [0.25, 0.75]
    circle2.cps += [0.75, 0.75]
    n = 30
    locations = np.vstack(
        (bottom.sample(n), circle1.sample(n), circle2.sample(n))
    )

    # callback
    def heights(bot_c1_c2_n):
        bot, c1, c2, n_repeat = bot_c1_c2_n

        values = np.repeat(
            [
                [0, 0, bot],
                [0, 0, c1],
                [0, 0, c2],
            ],
            n_repeat,
            axis=0,
        )
        return values

    nice_data = splinepy.SplineDataAdaptor(
        data=(4, 2, 1, n),
        locations=locations,
        function=heights,
    )
    disc = splinepy.helpme.create.disk(2, angle=123)
    disc.normalize_knot_vectors()
    disc.spline_data["nice"] = nice_data
    disc.show_options["arrow_data"] = "nice"
    disc.show_options["arrow_data_color"] = "jet"
    gus.show(
        [
            "8. Showing arbitrary data with callback.\n"
            "Any array-like multi-dim data can be plotted as long as\n"
            "len(data) == len(locations) holds.",
            disc,
        ]
    )

    # 9. plot derivative
    bezier = splinepy.Bezier(
        degrees=[3], control_points=[[0, 0], [1, 0], [1, 2], [3, 2]]
    )
    bezier.spline_data["deriv"] = bezier.derivative_spline([1])
    bezier.show_options["arrow_data"] = "deriv"
    bezier.show_options["arrow_data_on"] = np.linspace(0, 1, 7).reshape(-1, 1)
    bezier.show_options["control_points"] = True
    bezier.show_options["arrow_data_scale"] = 0.3
    bezier.show_options["control_mesh"] = False
    bezier.show_options["control_point_ids"] = False
    splinepy.show(
        [
            "9. Show the derivative of the of the line as a vector plor\n"
            "Close form derivatives can be computed for Bezier types",
            bezier,
        ]
    )
