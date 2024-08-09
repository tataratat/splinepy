from inspect import getfullargspec

import numpy as np

import splinepy.microstructure as ms

EPS = 1e-8

all_tile_classes = list(ms.tiles.everything().values())


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
    # TODO: also check non-default parameters

    def check_control_points(tile_patches):
        """Check if all of tile's control points all lie within unit square/cube"""
        # Go through all patches
        for tile_patch in tile_patches:
            cps = tile_patch.control_points
            assert np.all(
                (cps >= 0.0) & (cps <= 1.0)
            ), "Control points of tile must lie inside the unit square/cube"

    for tile_class in all_tile_classes:
        tile_creator = tile_class()
        # Create tile with default parameters
        tile_patches, _ = tile_creator.create_tile()
        check_control_points(tile_patches)


def test_tile_derivatives():
    """Testing the correctness of the tile derivatives using Finite Differences"""
    # TODO
    pass


def test_tile_closure():
    """Check if closing tiles also lie in unit cube. Maybe also check derivatives."""
    # TODO
    pass
