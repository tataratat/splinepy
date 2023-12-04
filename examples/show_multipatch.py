import gustaf as gus
import numpy as np

import splinepy

if __name__ == "__main__":
    """
    x1
    *---*---*
    | 2 | 3 |
    *---*---*
    | 0 | 1 |
    *---*---* -> x0
    """
    box0 = splinepy.helpme.create.box(1, 2)
    x0_max, x1_max = box0.cps.max(axis=0)

    box1 = box0.copy()
    box1.cps[:, 0] += x0_max

    box2 = box0.copy()
    box2.cps[:, 1] += x1_max

    box3 = box0.copy()
    box3.cps += [x0_max, x1_max]

    # plot scalar data
    m = splinepy.Multipatch([box0, box1, box2, box3])
    m.spline_data["t"] = m
    m.show_options["data_name"] = "t"
    m.show_options["lighting"] = "off"
    m.show()

    # plot arrow data
    m.show_options["arrow_data"] = "t"
    m.show()

    # plot arrow data with callback
    def plot_eval(data, resolutions=None, on=None):
        """
        evaluates at given location
        """
        if resolutions is not None:
            q = gus.create.vertices.raster(
                box0.parametric_bounds, resolutions
            ).vertices
            return data.evaluate(q)
        elif on is not None:
            return data.evaluate(on)

    def plot_jacs(data, resolutions=None, on=None):
        """
        evaluates at given location
        """
        if resolutions is not None:
            q = gus.create.vertices.raster(
                box0.parametric_bounds, resolutions
            ).vertices
            return np.linalg.det(
                np.vstack(d.jacobian(q) for d in data.patches)
            )
        elif on is not None:
            return np.linalg.det(data.jacobians(on))

    # use adaptor to plot values at specific places
    # this works because multipatch's evaluate() evaluates at
    # each spline. you need to make sure that your queries are within
    # each bounds.
    m.spline_data["t2"] = splinepy.SplineDataAdaptor(m, function=plot_eval)
    m.show_options["arrow_data"] = "t2"

    # form query and set evaluation points for arrow_data
    circle = splinepy.helpme.create.circle(radius=0.3)
    circle.cps -= circle.cps.min(axis=0)
    queries = circle.sample(20)[:-1]
    m.show_options["arrow_data_on"] = queries
    m.show()

    # Generate MS
    generator = splinepy.microstructure.Microstructure()
    generator.deformation_function = splinepy.Bezier(
        degrees=[2, 1],
        control_points=[[0, 0], [2, 1], [4, 0], [0, 2], [2, 4], [4, 2]],
    )
    generator.microtile = splinepy.microstructure.tiles.Cross2D()
    generator.tiling = [3, 2]

    m = splinepy.Multipatch(generator.create(center_expansion=1.2))
    m.spline_data["detJ"] = splinepy.SplineDataAdaptor(m, function=plot_jacs)
    m.show_options["data_name"] = "detJ"
    m.show_options["lighting"] = "off"
    m.show()
