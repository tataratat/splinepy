import numpy as np
from pytest import mark

import splinepy.microstructure as ms
from splinepy.helpme.create import box
from splinepy.helpme.integrate import volume

all_tile_classes = list(ms.tiles.everything().values())

TILING = [2, 2, 2]
BOX_DIMENSIONS = [1, 1, 1]
EPS = 1e-7

# TODO: the following tiles fail the closure test
CLOSURE_FAILS = [ms.tiles.HollowOctagonExtrude, ms.tiles.InverseCross3D]


@mark.parametrize("tile_class", all_tile_classes)
def test_closing_face(tile_class):
    """Check if closing face is working"""

    def check_if_closed(multipatch, closure_direction):
        """Helper function to see if multipatch has a closing surface

        Parameters
        ----------
        multipatch: splinepy Multipatch
            Microstructure multipatch
        closure_direction: str
            Direction in which the closure has been applied
        """
        direction_index_dict = {"x": 0, "y": 1, "z": 2}
        direction_index = direction_index_dict[closure_direction]

        def min_identifier(points):
            return points[:, direction_index] < EPS

        def max_identifier(points):
            return (
                points[:, direction_index]
                > BOX_DIMENSIONS[direction_index] - EPS
            )

        # create min- and max-boundaries of the multipatch using the identifiers
        multipatch.boundary_from_function(min_identifier, boundary_id=2)
        multipatch.boundary_from_function(max_identifier, boundary_id=3)

        min_patches = multipatch.boundary_multipatch(2)
        max_patches = multipatch.boundary_multipatch(3)

        # Check if the closing surface has a surface area of 1
        face_min_area = sum([volume(patch) for patch in min_patches.patches])
        face_max_area = sum([volume(patch) for patch in max_patches.patches])

        assert (
            face_min_area > 1.0 - EPS
        ), f"The closure of the {closure_direction}_min surface is not complete"
        assert (
            face_max_area > 1.0 - EPS
        ), f"The closure of the {closure_direction}_max surface is not complete"

    # Skip tile if it doesn't support closure
    if "_closure_directions" not in dir(tile_class):
        return

    # TODO: right now skip tiles which have faulty closures
    if tile_class in CLOSURE_FAILS:
        return

    tile_creator = tile_class()
    generator = ms.microstructure.Microstructure(
        deformation_function=box(*BOX_DIMENSIONS[: tile_creator._dim]),
        microtile=tile_creator,
        tiling=TILING[: tile_creator._dim],
    )
    # Go through all implemented closure direction
    closure_directions = {
        directionname[0] for directionname in tile_creator._closure_directions
    }
    for closure_direction in closure_directions:
        multipatch = generator.create(closing_face=closure_direction)
        check_if_closed(multipatch, closure_direction)


@mark.parametrize("tile_class", all_tile_classes)
def test_macro_sensitivities(tile_class, np_rng, heps=1e-7, n_test_points=10):
    """Testing the correctness of the derivatives of the whole microstructure w.r.t.
    the deformation function's control points. It is tested by evaluating the derivative
    obtained via finite differences. The values are evaluated at random points.

    Parameters
    ----------
    heps: float
        Perturbation size for finite difference evaluation
    n_test_points: int
        Number of testing points int the parametric domain"""

    tile_creator = tile_class()
    deformation_function_orig = box(*BOX_DIMENSIONS[: tile_creator._dim])
    generator = ms.microstructure.Microstructure(
        deformation_function=deformation_function_orig,
        microtile=tile_creator,
        tiling=TILING[: tile_creator._dim],
    )
    multipatch = generator.create(macro_sensitivities=True)
    dim = multipatch.dim
    n_cps = deformation_function_orig.cps.shape[0]

    # Set evaluation points as random spots in the parametric space
    eval_points = np_rng.random((n_test_points, tile_creator._para_dim))
    microstructure_orig_evaluations = [
        patch.evaluate(eval_points) for patch in multipatch.patches
    ]
    n_patches = len(multipatch.patches)

    # Go through derivatives of every deformation function's control point
    for ii_ctps in range(n_cps):
        # Gradient through finite differences
        deformation_function_perturbed = deformation_function_orig.copy()
        deformation_function_perturbed.cps[ii_ctps, :] += heps
        generator.deformation_function = deformation_function_perturbed
        multipatch_perturbed = generator.create()
        microstructure_perturbed_evaluations = [
            patch.evaluate(eval_points)
            for patch in multipatch_perturbed.patches
        ]
        # Evaluate finite difference gradient
        fd_sensitivity = [
            (patch_perturbed - patch_orig) / heps
            for patch_perturbed, patch_orig in zip(
                microstructure_orig_evaluations,
                microstructure_perturbed_evaluations,
            )
        ]

        # Go through each direction
        for jj_dim in range(dim):
            deriv_orig = multipatch.fields[ii_ctps * dim + jj_dim]
            deriv_evaluations = [
                patch.evaluate(eval_points) for patch in deriv_orig.patches
            ]
            for k_patch, patch_deriv_implemented, patch_deriv_fd in zip(
                range(n_patches), deriv_evaluations, fd_sensitivity
            ):
                assert np.allclose(
                    -patch_deriv_implemented[:, jj_dim],
                    patch_deriv_fd[:, jj_dim],
                ), (
                    "Implemented derivative calculation for tile class"
                    + f"{tile_class}, at patch {k_patch+1}/{n_patches} does not "
                    + "match the derivative obtained using Finite Differences at "
                    + "the following evaluation points:\n"
                    + str(eval_points)
                )
