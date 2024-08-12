from inspect import getfullargspec

import numpy as np

import splinepy.microstructure as ms
from splinepy.utils.data import cartesian_product as _cartesian_product

EPS = 1e-8

all_tile_classes = list(ms.tiles.everything().values())


def check_control_points(tile_patches):
    """Check if all of tile's control points all lie within unit square/cube"""
    # Go through all patches
    for tile_patch in tile_patches:
        cps = tile_patch.control_points
        assert np.all(
            (cps >= 0.0 - EPS) & (cps <= 1.0 + EPS)
        ), "Control points of tile must lie inside the unit square/cube"


def make_bounds_feasible(bounds):
    """Bounds are given are understood as open set of min. and max. values. Therefore,
    convert bounds to open set of these values.

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


def test_tile_class():
    """Checks if all tile classes have the appropriate members and functions."""
    required_class_variables = {
        "_para_dim": int,
        "_dim": int,
        "_evaluation_points": np.ndarray,
        "_n_info_per_eval_point": int,
        "_sensitivities_implemented": bool,
    }

    # Go through all available tiles
    for tile_class in all_tile_classes:
        # Get tile class' objects
        members = [
            attr for attr in dir(tile_class) if not attr.startswith("__")
        ]

        # Class must have function create_tile()
        assert hasattr(
            tile_class, "create_tile"
        ), "Tile class must have create_tile()"
        # Tile must be able to take parameters and sensitivities as input
        create_parameters = getfullargspec(tile_class.create_tile).args
        for required_param in ["parameters", "parameter_sensitivities"]:
            assert (
                required_param in create_parameters
            ), f"create_tile() must have '{required_param}' as an input parameter"

        # Ensure closure can be correctly handled
        if "closure" in create_parameters:
            assert "_closure_directions" in members, (
                "Tile class has closure ability. The available closure directions "
                + "are missing"
            )
            assert hasattr(
                tile_class, "_closing_tile"
            ), "Tile class has closure ability but no _closing_tile() function!"

        # Check if tile class has all required variables and they are the correct type
        for required_variable, var_type in required_class_variables.items():
            assert (
                required_variable in members
            ), f"Tile class needs to have member variable '{required_variable}'"
            assert isinstance(
                eval(f"tile_class.{required_variable}"), var_type
            ), f"Variable {required_variable} needs to be of type {var_type}"


def test_tile_bounds():
    """Test if tile is still in unit cube at the bounds. Currently only checks
    default parameter values."""
    # Skip certain tile classes for testing
    skip_tiles = [
        ms.tiles.EllipsVoid,  # Control points easily lie outside unitcube
        ms.tiles.Snappy,  # Has no parameters
    ]

    for tile_class in all_tile_classes:
        # Skip certain classes for testing
        skip_tile_class = False
        for skip_tile in skip_tiles:
            if tile_class is skip_tile:
                skip_tile_class = True
                break
        if skip_tile_class:
            continue

        tile_creator = tile_class()
        # Create tile with default parameters
        tile_patches, _ = tile_creator.create_tile()
        check_control_points(tile_patches)

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


def test_tile_derivatives():
    """Testing the correctness of the tile derivatives using Finite Differences"""
    # TODO
    for tile_class in all_tile_classes:
        tile_creator = tile_class()
        # Skip test if tile class has no implemented sensitivities
        if not tile_creator._sensitivities_implemented:
            continue
        pass


def test_tile_closure():
    """Check if closing tiles also lie in unit cube. Maybe also check derivatives."""
    # TODO: check closure also for extreme values of parameters
    for tile_class in all_tile_classes:
        # Skip tile if if does not support closure
        if "_closure_directions" not in dir(tile_class):
            continue
        tile_creator = tile_class()
        # Go through all implemented closure directions
        for closure_direction in tile_creator._closure_directions:
            tile_patches, sensitivities = tile_creator.create_tile(
                closure=closure_direction
            )
            assert (
                sensitivities is None
            ), "Ensure sensitivities for closure are None if no sensitivities are asked"

            check_control_points(tile_patches)
