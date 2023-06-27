import numpy as np

from splinepy.bezier import Bezier
from splinepy.microstructure.tiles.tilebase import TileBase


class Ellipsvoid(TileBase):
    """Void in form of an ellipse set into a unit cell.

    Parametrization acts on the elipsoid's orientation as well as on its
    expansion (radii in x and y/z direction).

    Used in ADAMM test case for Hutchinson.

    See create_tile for more information
    """

    def __init__(self):
        self._dim = 3
        self._evaluation_points = np.array(
            [
                [0.5, 0.5, 0.5],
            ]
        )
        self._n_info_per_eval_point = 4

        # Aux values
        c0 = 0.5 / 3**0.5
        aux_1 = 3**0.5 / 2 * 0.5
        aux_0 = aux_1 * 2
        self._sphere_ctps = np.array(
            [
                [-c0, -c0, -c0],
                [0, -aux_1, -aux_1],
                [c0, -c0, -c0],
                [-aux_1, 0, -aux_1],
                [0, 0, -aux_0],
                [aux_1, 0, -aux_1],
                [-c0, c0, -c0],
                [0, aux_1, -aux_1],
                [c0, c0, -c0],
                # Second layer
                [-aux_1, -aux_1, 0],
                [0, -aux_0, 0],
                [aux_1, -aux_1, 0],
                [-aux_0, 0, 0],
                [0, 0, 0],
                [aux_0, 0, 0],
                [-aux_1, aux_1, 0],
                [0, aux_0, 0],
                [aux_1, aux_1, 0],
                # Third layer
                [-c0, -c0, c0],
                [0, -aux_1, aux_1],
                [c0, -c0, c0],
                [-aux_1, 0, aux_1],
                [0, 0, aux_0],
                [aux_1, 0, aux_1],
                [-c0, c0, c0],
                [0, aux_1, aux_1],
                [c0, c0, c0],
            ]
        )

    def _rotation_matrix_x(self, angle):
        cc, ss = np.cos(angle), np.sin(angle)
        return np.array([[1, 0, 0], [0, cc, ss], [0, -ss, cc]])

    def _rotation_matrix_x_deriv(self, angle):
        cc, ss = np.cos(angle), np.sin(angle)
        return np.array([[0, 0, 0], [0, -ss, cc], [0, -cc, -ss]])

    def _rotation_matrix_y(self, angle):
        cc, ss = np.cos(angle), np.sin(angle)
        return np.array([[cc, 0, -ss], [0, 1, 0], [ss, 0, cc]])

    def _rotation_matrix_y_deriv(self, angle):
        cc, ss = np.cos(angle), np.sin(angle)
        return np.array([[-ss, 0, -cc], [0, 1, 0], [cc, 0, -ss]])

    def create_tile(
        self,
        parameters=None,
        parameter_sensitivities=None,
        **kwargs,
    ):
        """Create a full cuboid with an ellipsoid void in the center, that has
        a rotation applied in the form

        T = RotY * RotX * Tz

        Parameters are stored as a 4 dimensional field. The entries represent:
        p_0 -> radius along x-axis
        p_1 -> radius along y/z axis
        p_2 -> rotation around x-axis
        p_3 -> rotation around y-axis

        Parameters
        ----------
        parameters : np.array
          Defines 4 values necessary for tile definition
          1. expansion in ellipsoidal void in x-axis
          2. expansion of ellipsoidal void in y- and z-axis
          3. rotation around x-axis
          4. rotation around y- and z-axis
        parameter_sensitivities: np.array
          Sensitivities of all 4 values with respect to some number of external
          design variables

        Returns
        -------
        microtile_list : list(splines)
        """  # set to default if nothing is given
        if parameters is None:
            parameters = np.array([0.5, 0.5, 0, 0]).reshape(1, 1, 4)

        # Create center ellipsoid
        # RotY * RotX * DIAG(r_x, r_yz) * T_base
        [r_x, r_yz, rot_x, rot_y] = parameters.flatten()

        # Check if user requests derivative splines
        if self.check_param_derivatives(parameter_sensitivities):
            n_derivatives = parameter_sensitivities.shape[2]
        else:
            n_derivatives = 0

        # Prepare auxiliary matrices and values
        diag = np.diag([r_x, r_yz, r_yz])
        rotMx = self._rotation_matrix_x(rot_x)
        rotMy = self._rotation_matrix_y(rot_y)

        # Start assembly
        derivatives = []
        splines = []
        for i_derivative in range(n_derivatives + 1):
            spline_list = []
            # Constant auxiliary values
            if i_derivative == 0:
                v_one_half = 0.5
                v_one = 1.0
                v_zero = 0.0
                ctps = np.einsum(
                    "ij,jk,kl,ml->mi",
                    rotMy,
                    rotMx,
                    diag,
                    self._sphere_ctps.copy(),
                    optimize=True,
                )
            else:
                [dr_x, dr_yz, drot_x, drot_y] = parameter_sensitivities[
                    0, :, i_derivative - 1
                ].flatten()
                v_one_half = 0.0
                v_one = 0.0
                v_zero = 0.0
                ddiag = np.diag([dr_x, dr_yz, dr_yz])
                drotMx = self._rotation_matrix_x_deriv(rot_x)
                drotMy = self._rotation_matrix_y_deriv(rot_y)
                ctps_mat = (
                    drot_y
                    * np.einsum(
                        "ij,jk,kl->il",
                        drotMy,
                        rotMx,
                        diag,
                    )
                    + drot_x * np.einsum("ij,jk,kl->il", rotMy, drotMx, diag)
                    + np.einsum("ij,jk,kl->il", rotMy, rotMx, ddiag)
                )
                ctps = np.einsum("ij,kj->ki", ctps_mat, self._sphere_ctps)

            ctps += [v_one_half, v_one_half, v_one_half]
            # Start the assembly
            spline_list.append(
                Bezier(
                    degrees=[2, 2, 1],
                    control_points=[
                        [v_zero, v_zero, v_zero],
                        [v_one_half, v_zero, v_zero],
                        [v_one, v_zero, v_zero],
                        [v_zero, v_one_half, v_zero],
                        [v_one_half, v_one_half, v_zero],
                        [v_one, v_one_half, v_zero],
                        [v_zero, v_one, v_zero],
                        [v_one_half, v_one, v_zero],
                        [v_one, v_one, v_zero],
                        ctps[0, :],
                        ctps[1, :],
                        ctps[2, :],
                        ctps[3, :],
                        ctps[4, :],
                        ctps[5, :],
                        ctps[6, :],
                        ctps[7, :],
                        ctps[8, :],
                    ],
                )
            )
            spline_list.append(
                Bezier(
                    degrees=[2, 2, 1],
                    control_points=[
                        ctps[18, :],
                        ctps[19, :],
                        ctps[20, :],
                        ctps[21, :],
                        ctps[22, :],
                        ctps[23, :],
                        ctps[24, :],
                        ctps[25, :],
                        ctps[26, :],
                        [v_zero, v_zero, v_one],
                        [v_one_half, v_zero, v_one],
                        [v_one, v_zero, v_one],
                        [v_zero, v_one_half, v_one],
                        [v_one_half, v_one_half, v_one],
                        [v_one, v_one_half, v_one],
                        [v_zero, v_one, v_one],
                        [v_one_half, v_one, v_one],
                        [v_one, v_one, v_one],
                    ],
                )
            )
            # Y-Axis
            spline_list.append(
                Bezier(
                    degrees=[2, 1, 2],
                    control_points=[
                        [v_zero, v_zero, v_zero],
                        [v_one_half, v_zero, v_zero],
                        [v_one, v_zero, v_zero],
                        ctps[0, :],
                        ctps[1, :],
                        ctps[2, :],
                        [v_zero, v_zero, v_one_half],
                        [v_one_half, v_zero, v_one_half],
                        [v_one, v_zero, v_one_half],
                        ctps[9, :],
                        ctps[10, :],
                        ctps[11, :],
                        [v_zero, v_zero, v_one],
                        [v_one_half, v_zero, v_one],
                        [v_one, v_zero, v_one],
                        ctps[18, :],
                        ctps[19, :],
                        ctps[20, :],
                    ],
                )
            )
            spline_list.append(
                Bezier(
                    degrees=[2, 1, 2],
                    control_points=[
                        ctps[6, :],
                        ctps[7, :],
                        ctps[8, :],
                        [v_zero, v_one, v_zero],
                        [v_one_half, v_one, v_zero],
                        [v_one, v_one, v_zero],
                        ctps[15, :],
                        ctps[16, :],
                        ctps[17, :],
                        [v_zero, v_one, v_one_half],
                        [v_one_half, v_one, v_one_half],
                        [v_one, v_one, v_one_half],
                        ctps[24, :],
                        ctps[25, :],
                        ctps[26, :],
                        [v_zero, v_one, v_one],
                        [v_one_half, v_one, v_one],
                        [v_one, v_one, v_one],
                    ],
                )
            )
            # Z-Axis
            spline_list.append(
                Bezier(
                    degrees=[1, 2, 2],
                    control_points=[
                        [v_zero, v_zero, v_zero],
                        ctps[0, :],
                        [v_zero, v_one_half, v_zero],
                        ctps[3, :],
                        [v_zero, v_one, v_zero],
                        ctps[6, :],
                        [v_zero, v_zero, v_one_half],
                        ctps[9, :],
                        [v_zero, v_one_half, v_one_half],
                        ctps[12, :],
                        [v_zero, v_one, v_one_half],
                        ctps[15, :],
                        [v_zero, v_zero, v_one],
                        ctps[18, :],
                        [v_zero, v_one_half, v_one],
                        ctps[21, :],
                        [v_zero, v_one, v_one],
                        ctps[24, :],
                    ],
                )
            )
            spline_list.append(
                Bezier(
                    degrees=[1, 2, 2],
                    control_points=[
                        ctps[2, :],
                        [v_one, v_zero, v_zero],
                        ctps[5, :],
                        [v_one, v_one_half, v_zero],
                        ctps[8, :],
                        [v_one, v_one, v_zero],
                        ctps[11, :],
                        [v_one, v_zero, v_one_half],
                        ctps[14, :],
                        [v_one, v_one_half, v_one_half],
                        ctps[17, :],
                        [v_one, v_one, v_one_half],
                        ctps[20, :],
                        [v_one, v_zero, v_one],
                        ctps[23, :],
                        [v_one, v_one_half, v_one],
                        ctps[26, :],
                        [v_one, v_one, v_one],
                    ],
                )
            )

            # Pass to output
            if i_derivative == 0:
                splines = spline_list.copy()
            else:
                derivatives.append(spline_list)

        # Return results
        if i_derivative == 0:
            return splines
        else:
            return (splines, derivatives)
