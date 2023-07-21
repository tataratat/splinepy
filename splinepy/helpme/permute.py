import numpy as np

from splinepy import utils


def parametric_axes(spline, permutation_list, inplace=True):
    """
    Permutates the parametric dimensions

    This function can be used, e.g., to  interchange the parametric
    dimensions xi and eta in order to have them in the right orientation
    (applications in boundary condition definition or mfem export)

    Parameters
    ----------
    spline: Spline
    permutation_list : list
        New order of parametric dimensions
    inplace: bool
        Default is True. If True, modifies spline inplace, else, returns
        a modified_spline.

    Returns
    -------
    modified_spline : type(spline)
        spline with reordered parametric dimensions. iff `inplace=True`.
    """
    # Data collector for new spline object
    dict_spline = {}

    # Sanity checks
    if not isinstance(permutation_list, list):
        raise ValueError("Permutation list incomprehensive")
    if not set(range(spline.para_dim)) == set(permutation_list):
        raise ValueError("Permutation list invalid")

    utils.log.debug("Permuting parametric axes...")

    # Update knot_vectors where applicable
    if "knot_vectors" in spline.required_properties:
        dict_spline["knot_vectors"] = [
            spline.knot_vectors[permutation_list[i]]
            for i in range(spline.para_dim)
        ]
    # Update degrees
    dict_spline["degrees"] = [
        spline.degrees[permutation_list[i]] for i in range(spline.para_dim)
    ]

    # Retrieve control mesh resolutions
    ctps_dims = spline.control_mesh_resolutions
    new_ctps_dims = [
        ctps_dims[permutation_list[i]] for i in range(spline.para_dim)
    ]
    n_ctps = spline.control_points.shape[0]

    # Map global to local index
    # i_glob = i + n_i * j  + n_i * n_j * k ...
    local_indices = np.empty([n_ctps, spline.para_dim], dtype=int)
    global_indices = np.arange(n_ctps, dtype=int)
    for i_p in range(spline.para_dim):
        local_indices[:, i_p] = global_indices % ctps_dims[i_p]
        global_indices -= local_indices[:, i_p]
        global_indices = np.floor_divide(global_indices, ctps_dims[i_p])

    # Reorder indices
    local_indices[:] = local_indices[:, permutation_list]

    # Rearange global to local
    global_indices = np.matmul(
        local_indices, np.cumprod([1] + new_ctps_dims)[0:-1]
    )
    # Get inverse mapping
    global_indices = np.argsort(global_indices)
    if "weights" in spline.required_properties:
        dict_spline["weights"] = spline.weights[global_indices]
    dict_spline["control_points"] = spline.control_points[global_indices, :]

    if inplace:
        utils.log.debug("  applying permutation inplace")
        spline._new_core(**dict_spline)

        return None

    else:
        utils.log.debug("  returning permuted spline")
        return type(spline)(**dict_spline)
