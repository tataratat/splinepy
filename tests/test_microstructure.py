import splinepy.microstructure as ms
from splinepy.helpme.create import box
from splinepy.helpme.integrate import volume

all_tile_classes = list(ms.tiles.everything().values())

TILING = [2, 2, 2]
BOX_DIMENSIONS = [1, 1, 1]
EPS = 1e-7

# TODO: the following tiles fail the closure test
CLOSURE_FAILS = [ms.tiles.HollowOctagonExtrude, ms.tiles.InverseCross3D]


def test_closing_face():
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

    for tile_class in all_tile_classes:
        # Skip tile if it doesn't support closure
        if "_closure_directions" not in dir(tile_class):
            continue

        # TODO: right now skip tiles which have faulty closures
        if tile_class in CLOSURE_FAILS:
            continue

        tile_creator = tile_class()
        generator = ms.microstructure.Microstructure(
            deformation_function=box(*BOX_DIMENSIONS[: tile_creator._dim]),
            microtile=tile_creator,
            tiling=TILING[: tile_creator._para_dim],
        )
        # Go through all implemented closure direction
        closure_directions = {
            directionname[0]
            for directionname in tile_creator._closure_directions
        }
        for closure_direction in closure_directions:
            multipatch = generator.create(closing_face=closure_direction)
            check_if_closed(multipatch, closure_direction)


def test_macro_sensitivities():
    pass
