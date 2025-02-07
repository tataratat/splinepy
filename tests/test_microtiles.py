from inspect import getfullargspec

import numpy as np
from pytest import mark, raises, skip

import splinepy.microstructure as ms
from splinepy.utils.data import cartesian_product as _cartesian_product

EPS = 1e-8

all_tile_classes = list(ms.tiles.everything().values())
# Tile classes where closure should be tested
tile_classes_with_closure = [
    tile_class
    for tile_class in all_tile_classes
    if tile_class._closure_directions is not None
]
# Tile classes where sensitivities should be tested
tile_classes_with_sensitivities = [
    tile_class
    for tile_class in all_tile_classes
    if tile_class._sensitivities_implemented
]

# Skip certain tile classes for parameters testing
skip_tiles = {
    ms.tiles.EllipsVoid: "control points easily lie outside unitcube",
    ms.tiles.Snappy: "has no 'parameters' implemented",
}


def check_control_points(tile_patches):
    """Helper function. Check if all of tile's control points all lie within unit
    square/cube"""
    # Go through all patches
    for tile_patch in tile_patches:
        cps = tile_patch.control_points
        assert np.all((cps >= 0.0 - EPS) & (cps <= 1.0 + EPS)), (
            "Control points of tile must lie inside the unit square/cube. "
            + "Found points outside bounds: "
            f"{cps[~((cps >= 0.0 - EPS) & (cps <= 1.0 + EPS))]}"
        )


def make_bounds_feasible(bounds):
    """Helper function. Bounds are understood as open set of min. and max. values.
    Therefore, convert bounds to open set of these values.

    Parameters
    ------------
    bounds: list<list<float>>
        Values of bounds

    Returns
    ----------
    feasible_bounds: np.ndarray
        Values of bounds
    """
    feasible_bounds = [
        [min_value + EPS, max_value - EPS] for min_value, max_value in bounds
    ]
    return np.array(feasible_bounds)


@mark.parametrize("tile_class", all_tile_classes)
def test_tile_class(tile_class):
    """Checks if all tile classes have the appropriate members and functions.

    Parameters
    ---------
    tile_class: tile class in splinepy.microstructure.tiles
        Microtile
    """
    # Create instance of class
    tile_instance = tile_class()
    tile_name = tile_class.__name__

    required_class_variables = {
        "_para_dim": int,
        "_dim": int,
        "_evaluation_points": np.ndarray,
        "_n_info_per_eval_point": int,
        "_sensitivities_implemented": bool,
        "_parameter_bounds": list,
        "_parameters_shape": tuple,
    }

    # Get tile class' objects
    members = [
        attr for attr in dir(tile_instance) if not attr.startswith("__")
    ]

    # Class must have function create_tile()
    assert hasattr(
        tile_instance, "create_tile"
    ), f"Tile class {tile_name} must have create_tile() method"

    # Tile must be able to take parameters and sensitivities as input
    create_parameters = getfullargspec(tile_instance.create_tile).args
    for required_param in ["parameters", "parameter_sensitivities"]:
        assert required_param in create_parameters, (
            f"{tile_name}.create_tile() must have '{required_param}' as an "
            "input parameter"
        )

    # Ensure closure can be handled correctly
    if "closure" in create_parameters:
        assert "_closure_directions" in members, (
            f"Tile class {tile_name} has closure ability. The available closure "
            + "directions are missing"
        )
        assert hasattr(
            tile_instance, "_closing_tile"
        ), f"Tile class {tile_name} has closure ability but no _closing_tile() function"

    # Check if tile class has all required variables and they are the correct type
    for required_variable, var_type in required_class_variables.items():
        assert (
            required_variable in members
        ), f"Tile class {tile_name} needs to have member variable '{required_variable}'"
        assert isinstance(
            getattr(tile_instance, required_variable), var_type
        ), (
            f"Variable {required_variable} needs to be of type {var_type} and not"
            f"{required_variable}"
        )

    # Check default parameter value if there is one
    if tile_instance._parameters_shape != ():
        # Assert that there is a default value
        assert hasattr(tile_instance, "_default_parameter_value"), (
            f'{tile_name} must have "_default_parameter_value" as a class '
            " attribute."
        )
        # Check the default value's type
        default_value = tile_instance._default_parameter_value
        if isinstance(default_value, np.ndarray):
            # Check the dimensions
            assert default_value.shape == tile_instance._parameters_shape, (
                f"Default parameter values for tile {tile_name} has the wrong"
                " dimensions"
            )
            # Check if default values are within bounds
            default_value = default_value.ravel()
        elif not isinstance(default_value, float):
            raise ValueError(
                f"Default parameter value for tile {tile_name} must either be "
                "a float or a numpy array"
            )

        # Check if default values are within bounds
        parameter_bounds = np.asarray(tile_instance._parameter_bounds)
        lower_bounds = parameter_bounds[:, 0]
        upper_bounds = parameter_bounds[:, 1]
        assert np.all(
            (default_value > lower_bounds) & (default_value < upper_bounds)
        ), f"Default parameter value of tile {tile_name} is not within bounds"


@mark.parametrize("tile_class", all_tile_classes)
def test_tile_bounds(tile_class):
    """Test if tile is still in unit cube at the bounds. Checks default and also
    non-default parameter values.

    Parameters
    ---------
    tile_class: tile class in splinepy.microstructure.tiles
        Microtile
    """
    tile_creator = tile_class()
    # Create tile with default parameters
    tile_patches, _ = tile_creator.create_tile()
    check_control_points(tile_patches)

    # Skip certain classes for testing
    if tile_class in skip_tiles:
        skip(
            "Known issue: skip bound test for non-default parameter values for tile "
            f"{tile_class.__name__}. Reason: {skip_tiles[tile_class]}"
        )

    # Go through all extremes of parameters and ensure that also they create
    # tiles within unit square/cube
    feasible_parameter_bounds = make_bounds_feasible(
        tile_creator._parameter_bounds
    )
    all_parameter_bounds = _cartesian_product(feasible_parameter_bounds)
    for parameter_extremes in all_parameter_bounds:
        tile_patches, _ = tile_creator.create_tile(
            parameters=parameter_extremes.reshape(
                tile_creator._parameters_shape
            )
        )
        check_control_points(tile_patches)


@mark.parametrize("tile_class", tile_classes_with_closure)
def test_tile_closure(tile_class):
    """Check if closing tiles also lie in unit cube.

    Parameters
    ---------
    tile_class: tile class in splinepy.microstructure.tiles
        Microtile
    """
    tile_name = tile_class.__name__

    # Skip tile if if does not support closure
    if tile_class._closure_directions is None:
        skip(f"Tile {tile_name} does not have a closure implementation. Skip")
    tile_creator = tile_class()
    # Go through all implemented closure directions
    for closure_direction in tile_creator._closure_directions:
        tile_patches, sensitivities = tile_creator.create_tile(
            closure=closure_direction
        )
        assert sensitivities is None, (
            f"Expected sensitivities to be None for closure {closure_direction} "
            + f"when no sensitivities are requested, got {sensitivities}"
        )

        check_control_points(tile_patches)

    # Also check non-default parameters
    # Skip certain classes for testing
    if tile_class in skip_tiles:
        skip(
            "Known issue: skip closure test for non-default parameter for tile "
            f"{tile_name}. Reason: {skip_tiles[tile_class]}"
        )

    # Go through all extremes of parameters and ensure that also they create
    # tiles within unit square/cube
    feasible_parameter_bounds = make_bounds_feasible(
        tile_creator._parameter_bounds
    )
    all_parameter_bounds = _cartesian_product(feasible_parameter_bounds)
    # Go through all implemented closure directions
    for closure_direction in tile_creator._closure_directions:
        for parameter_extremes in all_parameter_bounds:
            tile_patches, _ = tile_creator.create_tile(
                parameters=parameter_extremes.reshape(
                    tile_creator._parameters_shape
                ),
                closure=closure_direction,
            )
            check_control_points(tile_patches)


@mark.parametrize("tile_class", tile_classes_with_sensitivities)
def test_tile_derivatives(
    tile_class,
    np_rng,
    heps,
    n_test_points,
    fd_derivative_stepsizes_and_weights,
):
    """Testing the correctness of the tile derivatives using Finite Differences.
    This includes every closure and no closure, every parameter and every patch
    by evaluating at random points and for random parameters.

    Parameters
    ---------
    tile_class: tile class in splinepy.microstructure.tiles
        Microtile
    np_rng: numpy.random._generator.Generator
        Default random number generator
    heps: float
        Perturbation size for finite difference evaluation. Defined in conftest.py
    n_test_points: int
        Number of testing points in the parametric domain. Defined in conftest.py
    fd_derivative_stepsizes_and_weights: dict<int:int>
        Stepsizes and weights for finite difference scheme. Defined in conftest.py
    """
    tile_creator = tile_class()
    # Skip test if tile class has no implemented sensitivities
    if not tile_creator._sensitivities_implemented:
        skip(f"Tile {tile_class.__name__} has no sensitivities implemented")

    def generate_random_parameters(tile_creator, np_rng):
        """Generate random parameters within bounds"""
        parameter_bounds = np.array(tile_creator._parameter_bounds)
        parameters = parameter_bounds[:, 0] + np_rng.random(
            len(parameter_bounds)
        ) * np.ptp(parameter_bounds, axis=1)
        return parameters.reshape(tile_creator._parameters_shape)

    parameters = generate_random_parameters(tile_creator, np_rng)

    # Test no closure as well as ...
    closure_directions = [None]
    # ... every closure implemented
    if tile_creator._closure_directions is not None:
        closure_directions += tile_creator._closure_directions

    # Retrieve shape values of parameters
    n_eval_points = tile_creator._evaluation_points.shape[0]
    n_info_per_eval_point = tile_creator._n_info_per_eval_point

    def derivative_finite_difference_evaluation(
        tile_creator, parameters, i_parameter, closure, eval_points, n_patches
    ):
        """Compute the derivative using fourth-order accurate centered finite
        differences

        Parameters
        ------------
        tile_creator: instance of microtile
            Microtile class
        parameters: np.ndarray
            Tile parameter values
        i_parameter: int
            Index of parameter, on which the FD derivative calculation should be
            performed on
        closure: None/str
            Closure direction of tile
        eval_points: np.ndarray
            Evaluation points on where to evaluate the derivative on
        n_patches: int
            Number of patches of tile
        """
        # Initialize array for FD evaluation
        fd_sensitivities = np.zeros(
            (n_patches, len(eval_points), tile_creator.dim)
        )

        # Go through the
        for stepsize, weighting in fd_derivative_stepsizes_and_weights.items():
            # Perturb parameter with respective stepsize
            parameters_perturbed = parameters.copy()
            parameters_perturbed[:, i_parameter] += stepsize * heps
            # Create patches with perturbed parameter value
            splines_perturbed, _ = tile_creator.create_tile(
                parameters=parameters_perturbed, closure=closure
            )
            fd_sensitivities += np.array(
                [
                    weighting / heps * spl.evaluate(eval_points)
                    for spl in splines_perturbed
                ]
            )

        return fd_sensitivities

    # Test each closure direction
    for closure in closure_directions:
        # Evaluate tile with given parameter and closure configuration
        splines_orig, _ = tile_creator.create_tile(
            parameters=parameters, closure=closure
        )
        n_patches = len(splines_orig)
        # Set evaluation points as random spots in the parametric space
        eval_points = np_rng.random((n_test_points, splines_orig[0].para_dim))
        # Go through all the parameters individually
        for i_parameter in range(n_info_per_eval_point):
            # Get implemented derivatives w.r.t. one parameter
            parameter_sensitivities = np.zeros(
                (n_eval_points, n_info_per_eval_point, 1)
            )
            parameter_sensitivities[:, i_parameter, :] = 1
            _, derivatives = tile_creator.create_tile(
                parameters=parameters,
                parameter_sensitivities=parameter_sensitivities,
                closure=closure,
            )
            deriv_evaluations = [
                deriv.evaluate(eval_points) for deriv in derivatives[0]
            ]
            # Perform finite difference evaluation
            fd_sensitivities = derivative_finite_difference_evaluation(
                tile_creator,
                parameters,
                i_parameter,
                closure,
                eval_points,
                n_patches,
            )
            # Check every patch
            for i_patch, deriv_orig, deriv_fd in zip(
                range(n_patches), deriv_evaluations, fd_sensitivities
            ):
                message = (
                    f"Implemented derivative calculation for tile class {tile_class} "
                    f"with closure {closure}, parameter {i_parameter+1}/"
                    f"{n_info_per_eval_point} at patch {i_patch+1}/{n_patches} does not"
                    f" match the derivative obtained using Finite Differences at the "
                    f"following evaluation points:\n {eval_points}\nImplemented "
                    f"derivative:\n{deriv_orig}\nFinite Difference derivative:\n"
                    f"{deriv_fd}"
                )
                assert np.allclose(deriv_orig, deriv_fd), message


@mark.parametrize("tile_class", all_tile_classes)
def test_invalid_parameter_values(tile_class, big_perturbation):
    """Testing whether the tile class correctly raises an error if invalid parameters
    are given. Current tests include too low or too high parameter values and wrong
    shapes of the parameter array.

    Parameters
    ----------
    tile_class: tile class in splinepy.microstructure.tiles
        Microtile class
    """
    tile_creator = tile_class()
    # For certain tiles skip tests
    if len(tile_creator._parameter_bounds) == 0:
        skip(
            f"Skip check for invalid parameters for tile {tile_class.__name__} "
            "since there are no parameter bounds implemented for this tile"
        )

    parameter_bounds = np.asarray(tile_creator._parameter_bounds)

    # Check if tile class correctly raises an error if parameter values are too low
    parameters_too_low = (
        parameter_bounds[:, 0].reshape(tile_creator._parameters_shape)
        - big_perturbation
    )
    with raises(ValueError) as exc_info_low:
        tile_creator.create_tile(parameters=parameters_too_low)
    # Check if the exception message calls TileBase.check_params()
    assert "The following parameters are out of bounds: " in str(
        exc_info_low.value
    ), (
        f"Tile class {tile_class.__name__} must call TileBase.check_params() and raise",
        " a ValueError if parameters values are too low",
    )

    # Check the same if parameter values are too high
    parameters_too_high = (
        parameter_bounds[:, 1].reshape(tile_creator._parameters_shape)
        + big_perturbation
    )
    with raises(ValueError) as exc_info_high:
        tile_creator.create_tile(parameters=parameters_too_high)
    # Check if the exception message calls TileBase.check_params()
    assert "The following parameters are out of bounds: " in str(
        exc_info_high.value
    ), (
        f"Tile class {tile_class.__name__} must call TileBase.check_params() and raise",
        " a ValueError if parameters values are too high",
    )

    # Test if error is correctly thrown if the parameter shape is incompatible
    # Take parameters in the middle of the bounds and double the array size
    parameters_middle = np.mean(parameter_bounds, axis=1).reshape(
        tile_creator._parameters_shape
    )
    parameters_middle = np.tile(parameters_middle, [2, 1])
    # Check if the error is correctly raised
    with raises(TypeError) as exc_info_size:
        tile_creator.create_tile(parameters=parameters_middle)
    assert "Mismatch in parameter size, expected" in str(
        exc_info_size.value
    ), (
        f"Tile class {tile_class.__name__} must call TileBase.check_params() and raise",
        " a TypeError if the given parameters have the wrong shape",
    )
