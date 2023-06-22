"""splinepy/helpme/ffd.py.

Freeform Deformation!
"""
import gustaf as gus

from splinepy import settings
from splinepy._base import SplinepyBase
from splinepy.helpme.create import from_bounds
from splinepy.spline import Spline


def _rescaled_vertices(vertices, scale, offset):
    return (vertices * (1 / scale)) + offset


def _spline_para_dim_and_mesh_dim_matches(spline, mesh):
    return spline.para_dim == mesh.vertices.shape[1]


class FFD(SplinepyBase):
    """
    Free-form deformation is a method used to deform an object by a
    deformation function. In our case the object is given via a mesh, the
    currently supported mesh-types are gustaf.Vertices and its subclasses,
    and the deformations function by a spline, supported splines are
    splinepy.Spline and its subclasses.
    The splines parametric dimension will be scaled in to a unit-hypercube
    as well as the original meshes vertices.

    The FFD class provides functions to modify the spline by completely
    overwriting the spline whole spline or parts of it. To obtain the
    deformed mesh mapped into the latest spline, retrieve the mesh
    attribute.

    Parameters
    ----------
    mesh: gustaf.Vertices
      Mesh used in the FFD. Defaults to None.
    spline: splinepy.Spline
      Spline used in the FFD. Defaults to None.
    padding: float
      Padding factor to scale mesh into spline's parametric space.
      Ideally, this should be 0, but to avoid floating point issues during
      evaluation, we apply default settings.TOLERANCE

    Returns
    -------
    None
    """

    def __init__(
        self,
        mesh=None,
        spline=None,
        padding=None,
    ):
        # init variables
        self._mesh = None
        self._spline = None
        self._q_vertices = None
        self._q_scale = None
        self._q_offset = None

        # use setters for attr
        if padding is None:
            self.padding = settings.TOLERANCE
        if mesh is not None:
            self.mesh = mesh
        if spline is not None:
            self.spline = spline

    @property
    def mesh(self):
        """Returns copy of current mesh. Before copying, it applies
        deformation.

        Returns
        -------
        current_mesh: gustaf.Vertices or derived
            Current Mesh with the deformation according to the current spline.
        """
        # is there mesh?
        if self._mesh is None:
            raise ValueError("Please set mesh first.")

        # evaluate new vertices
        current_mesh = type(self._mesh)(
            vertices=self._spline.evaluate(self._q_vertices), copy=False
        )

        # apply connectivity if applicable
        if hasattr(self._mesh, "elements"):
            current_mesh.elements = self._mesh.elements.copy()

        # setback to default
        current_mesh.setter_copies = True

        return current_mesh

    @mesh.setter
    def mesh(self, mesh):
        """Sets mesh. If it is first time, the copy of it will be saved as
        original mesh. If spline is already defined and in transformed status,
        it applies transformation directly.

        Parameters
        -----------
        mesh: gustaf.Vertices
            Mesh used for the FFD

        Returns
        --------
        None
        """
        if mesh is None:
            self._mesh = None
            self._q_vertices = None
            self._q_offset = None
            self._q_scale = None
            return None

        # type check
        if not isinstance(mesh, gus.Vertices):
            raise TypeError("mesh should be an instance of gus.Vertices")

        # warn and erase existing spline if dimension doesn't match
        if (
            self._spline is not None
            and not _spline_para_dim_and_mesh_dim_matches(self._spline, mesh)
        ):
            self._logw(
                "Mismatch between mesh vertices' dim and set spline's",
                "para_dim.",
                "Deleting existing spline",
            )

            self._spline = None

        # create a spline if there isn't one
        if self._spline is None:
            # Define a default spline if mesh is given but no spline
            par_dim = mesh.vertices.shape[1]
            # try to expand bounds same ratio as
            bounds = mesh.bounds().copy()
            self.spline = from_bounds([[0] * par_dim, [1] * par_dim], bounds)
            self._logd("created default spline")

        self._logd("Setting mesh.")
        self._logd("Mesh Info:")
        self._logd(f"  Vertices: {mesh.vertices.shape}.")
        self._logd(f"  Bounds: {mesh.bounds()}.")

        # keep original copy
        self._mesh = mesh.copy()

        self._logd("Fitting mesh into spline's parametric space.")

        self._q_vertices = self._mesh.vertices.copy()

        # get mesh bound to scale
        scaling_bounds = self._mesh.bounds().copy()

        # apply padding
        scaling_bounds[0] -= self.padding
        scaling_bounds[1] += self.padding

        # save mesh offset and scale
        self._q_offset = scaling_bounds[0]
        self._q_scale = 1 / (scaling_bounds[1] - scaling_bounds[0])

        # scale and offset vertices coordinates
        self._q_vertices -= self._q_offset
        self._q_vertices *= self._q_scale

        self._logd("Successfully scaled mesh vertices!")

    @property
    def spline(self):
        """Returns a copy of the spline. Please use the setter to explicitly
        make changes to the spline.

        Parameters
        -----------
        None

        Returns
        --------
        self._spline: Spline
        """
        return self._spline

    @spline.setter
    def spline(self, spline):
        """Sets spline. The spline parametric range bounds will be converted
        into the bounds [0,1]^para_dim.

        Parameters
        -----------
        spline: Spline
            New Spline for the next deformation

        Returns
        --------
        None
        """
        if spline is None:
            self._spline = None
            return None

        if not isinstance(spline, Spline):
            raise TypeError(
                "spline should be an instance of splinepy.Spline class"
            )

        if (
            self._mesh is not None
            and not _spline_para_dim_and_mesh_dim_matches(spline, self._mesh)
        ):
            raise ValueError(
                "Spline's para_dim should match mesh vertices' dim."
            )

        # copy input
        self._spline = spline.copy()

        # normalize if needed
        if spline.has_knot_vectors:
            self._spline.normalize_knot_vectors()

    @property
    def padding(self):
        """
        Padding ratio for fitting mesh inside spline's parametric dimension.
        It should be a positive number.

        Parameters
        ----------
        None

        Returns
        -------
        padding: float
        """
        return self._padding

    @padding.setter
    def padding(self, padding):
        """
        Padding Setter.
        """
        if padding < 0:
            raise ValueError("Padding should be a positive number.")

        if padding < settings.TOLERANCE:
            raise ValueError(
                "Padding is too small - "
                "should be bigger than the tolerance ({settings.TOLERANCE})."
            )

        self._padding = padding

        # if there's q_vertices, reset mesh, so that new q_vertices are formed
        if self._q_vertices is not None:
            self.mesh = self._mesh

    def show(self, **kwargs):
        """Visualize. Shows the deformed mesh and the current spline. Currently
        visualization is limited to vedo.

        Parameters
        ----------
        title: str
            Title of the vedo window. Defaults to "gustaf - FFD".
        return_showable: bool
            If true returns a dict of the showable items. Defaults to False.
        return_gustaf: bool
            Return dict of gustaf discrete objects, for example,
            {Vertices, Edges, Faces}, instead of opening a window.
            Defaults to False.
        kwargs: Any
            Arbitrary keyword arguments. These are passed onto the vedo
            functions. Please be aware, that no checking of these are performed
            in this function.

        Returns
        -------
        Any:
            Returns, if applicable, the vedo plotter. 'close=False' as argument
            to get the plotter.
        """
        if self._spline is None and self._mesh is None:
            raise ValueError("Please set a mesh before calling show()")
        return_showable = kwargs.pop("return_showable", False)
        return_gustaf = kwargs.pop("return_gustaf", False)
        title = kwargs.pop("title", "gustaf - FFD")

        if return_gustaf and return_showable:
            raise ValueError(
                "Either one of following params can be True: "
                "{return_gustaf, return_showable} "
                "You've set both True."
            )

        things_to_show = self.showable(**kwargs)

        if return_gustaf or return_showable:
            return things_to_show

        return gus.show(
            ["Original Mesh", things_to_show["original_mesh"]],
            [
                "Spline and Mesh",
                *things_to_show["spline"].values(),
                things_to_show["mesh"],
            ],
            ["Deformed Mesh", things_to_show["mesh"]],
            title=title,
        )

    def showable(self, **kwargs):
        """Returns a dictionary of showable items to describe the FFD at the
        current state.

        See show() for more information. This function redirects to it
        directly with the return_showable keyword set to True.
        """
        if self._spline is None and self._mesh is None:
            raise ValueError("Please set a mesh before calling show()")

        return_gustaf = kwargs.get("return_gustaf", False)

        things_to_show = dict()

        # let's show faces at most, since volumes can take awhile
        mesh = self._mesh
        if mesh.kind == "volume":
            mesh = gus.Faces(
                mesh.const_vertices,
                mesh.faces()[mesh.single_faces()],
                copy=False,
            )

        # original mesh
        things_to_show["original_mesh"] = mesh

        # current spline
        self._spline.show_options["alpha"] = 0.3
        things_to_show["spline"] = self._spline

        # deformed mesh
        deformed = self.mesh
        things_to_show["mesh"] = deformed

        if return_gustaf:
            return things_to_show

        # let's turn everything into showable and return
        for k, v in things_to_show.items():
            things_to_show[k] = v.showable(**kwargs)

        return things_to_show
