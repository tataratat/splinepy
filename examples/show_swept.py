import gustaf as gus
import numpy as np

import splinepy

if __name__ == "__main__":
    # trajectory
    # define degrees
    ds1 = [2]
    # define knot vectors
    kvs1 = [[0.0, 0.0, 0.0, 1.0, 1.0, 1.0]]
    # define control points
    cps1 = np.array(
        [
            [0.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
            [4.0, 0.0, 0.0],
        ]
    )
    # init trajectory as bspline
    trajectory = splinepy.BSpline(
        degrees=ds1,
        knot_vectors=kvs1,
        control_points=cps1,
    )

    # define sections along trajectory
    nsect = len(kvs1[0])-ds1[0]-1

    # compute parameter values for inserting cross sections
    par_value = np.zeros((nsect,1))
    for i in range(nsect):
        par_value[i][0] = (kvs1[0][i+1]+kvs1[0][i+2])/2

    # evaluate trajectory at these parameter values
    evals = trajectory.evaluate(par_value)
 
    # cross section
    # define degrees
    ds_cs = [1]

    # define knot vectors
    kvs_cs = [[0.0, 0.0, 1.0, 1.0]]

    # define control points
    cps_cs = np.array(
        [
            [0.0, 0.0, -1.0],
            [0.0, 0.0, 1.0],
        ]
    )

    # init cross section as bspline
    cross_section = splinepy.BSpline(
        degrees=ds_cs,
        knot_vectors=kvs_cs,
        control_points=cps_cs,
    )

    # set cross section CPs along trajectory
    cps_eval = []
    for eval_point in evals:
        # for-loop representing the nr. of cross section CPs
        # note: only working for this special case
        for i in range(len(cps_cs)): 
            current_cps = []
            current_cps.append(eval_point[0])
            current_cps.append(eval_point[1])
            current_cps.append(cps_cs[i][2])
            cps_eval.append(current_cps)

    fitting_points = np.array(cps_eval)
    
    # fit surface
    interpolated_surface, _ = splinepy.helpme.fit.surface(
        fitting_points=fitting_points,
        size=[2, 3],
        n_control_points=[2, 3],
        degrees=[1, 2],
    )
    
    # testing derivative and evaluation
    # derivative = trajectory.derivative([[0.5], [1]], [1])
    # evals = trajectory.evaluate([[0.5]])

    # show
    def der(data, on):
        return data.derivative(on, [1])

    trajectory.spline_data["der"] = splinepy.SplineDataAdaptor(
        trajectory, function=der
    )
    trajectory.show_options["arrow_data"] = "der"
    trajectory.show_options["arrow_data_on"] = np.linspace(0, 1, 20).reshape(
        -1, 1
    )

    gus.show(
        ["Trajectory", trajectory],
        ["Cross section", cross_section],
        resolution=50,
    )

    gus.show(
        ["Surface", interpolated_surface],
        resolution=50,
    )
