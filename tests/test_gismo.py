import tempfile
try:
    from . import common as c
except BaseException:
    import common as c

_gismo_export_ref_3d = [
        '<xml>\n',
        '  <MultiPatch parDim="3" id="0">\n',
        '    <patches type="id_range">100 103</patches>\n',
        '    <interfaces>101 3 100 4 0 1 2 1 1 1\n',
        '103 1 100 2 0 1 2 1 1 1\n',
        '102 1 101 2 0 1 2 1 1 1\n',
        '103 4 102 3 0 1 2 1 1 1</interfaces>\n',
        '    <boundary>100 1\n',
        '100 3\n',
        '100 5\n',
        '100 6\n',
        '101 1\n',
        '101 4\n',
        '101 5\n',
        '101 6\n',
        '102 2\n',
        '102 4\n',
        '102 5\n',
        '102 6\n',
        '103 2\n',
        '103 3\n',
        '103 5\n',
        '103 6</boundary>\n',
        '  </MultiPatch>\n',
        '  <boundaryConditions multipatch="0" id="1">\n',
        '    <!--Please fill boundary conditions here-->\n',
        '    <bc type="Dirichlet" unknown="0">0 3\n',
        '0 5\n',
        '0 6\n',
        '1 4\n',
        '1 5\n',
        '1 6\n',
        '2 4\n',
        '2 5\n',
        '2 6\n',
        '3 3\n',
        '3 5\n',
        '3 6</bc>\n',
        '    <bc type="Dirichlet" unknown="0">0 1\n',
        '1 1</bc>\n',
        '  </boundaryConditions>\n',
        '  <Geometry type="TensorBSpline3" id="100">\n',
        '    <Basis type="TensorBSplineBasis3">\n',
        '      <Basis type="BSplineBasis" index="0">\n',
        (
                '        <KnotVector degree="2">0.0 0.0 0.0'
                ' 1.0 1.0 1.0</KnotVector>\n'
        ),
        '      </Basis>\n',
        '      <Basis type="BSplineBasis" index="1">\n',
        '        <KnotVector degree="1">0.0 0.0 1.0 1.0</KnotVector>\n',
        '      </Basis>\n',
        '      <Basis type="BSplineBasis" index="2">\n',
        (
                '        <KnotVector degree="2">0.0 0.0 '
                '0.0 1.0 1.0 1.0</KnotVector>\n'
        ),
        '      </Basis>\n',
        '    </Basis>\n',
        '    <coefs geoDim="3">0.0 0.0 0.0\n',
        '0.5 0.0 0.0\n',
        '1.0 0.0 0.0\n',
        '0.0 1.0 0.0\n',
        '0.5 1.0 0.0\n',
        '1.0 1.0 0.0\n',
        '0.0 0.0 0.5\n',
        '0.5 0.0 0.5\n',
        '1.0 0.0 0.5\n',
        '0.0 1.0 0.5\n',
        '0.5 1.0 0.5\n',
        '1.0 1.0 0.5\n',
        '0.0 0.0 1.0\n',
        '0.5 0.0 1.0\n',
        '1.0 0.0 1.0\n',
        '0.0 1.0 1.0\n',
        '0.5 1.0 1.0\n',
        '1.0 1.0 1.0</coefs>\n',
        '  </Geometry>\n',
        '  <Geometry type="TensorBSpline3" id="101">\n',
        '    <Basis type="TensorBSplineBasis3">\n',
        '      <Basis type="BSplineBasis" index="0">\n',
        (
                '        <KnotVector degree="2">0.0 0.0 '
                '0.0 1.0 1.0 1.0</KnotVector>\n'
        ),
        '      </Basis>\n',
        '      <Basis type="BSplineBasis" index="1">\n',
        '        <KnotVector degree="1">0.0 0.0 0.5 1.0 1.0</KnotVector>\n',
        '      </Basis>\n',
        '      <Basis type="BSplineBasis" index="2">\n',
        (
                '        <KnotVector degree="2'
                '">0.0 0.0 0.0 1.0 1.0 1.0</KnotVector>\n'
        ),
        '      </Basis>\n',
        '    </Basis>\n',
        '    <coefs geoDim="3">0.0 1.0 0.0\n',
        '0.5 1.0 0.0\n',
        '1.0 1.0 0.0\n',
        '0.0 1.5 0.0\n',
        '0.5 1.5 0.0\n',
        '1.0 1.5 0.0\n',
        '0.0 2.0 0.0\n',
        '0.5 2.0 0.0\n',
        '1.0 2.0 0.0\n',
        '0.0 1.0 0.5\n',
        '0.5 1.0 0.5\n',
        '1.0 1.0 0.5\n',
        '0.0 1.5 0.5\n',
        '0.5 1.5 0.5\n',
        '1.0 1.5 0.5\n',
        '0.0 2.0 0.5\n',
        '0.5 2.0 0.5\n',
        '1.0 2.0 0.5\n',
        '0.0 1.0 1.0\n',
        '0.5 1.0 1.0\n',
        '1.0 1.0 1.0\n',
        '0.0 1.5 1.0\n',
        '0.5 1.5 1.0\n',
        '1.0 1.5 1.0\n',
        '0.0 2.0 1.0\n',
        '0.5 2.0 1.0\n',
        '1.0 2.0 1.0</coefs>\n',
        '  </Geometry>\n',
        '  <Geometry type="TensorNurbs3" id="102">\n',
        '    <Basis type="TensorNurbsBasis3">\n',
        '      <Basis type="TensorBSplineBasis3">\n',
        '        <Basis type="BSplineBasis" index="0">\n',
        '          <KnotVector degree="1">0.0 0.0 1.0 1.0</KnotVector>\n',
        '        </Basis>\n',
        '        <Basis type="BSplineBasis" index="1">\n',
        '          <KnotVector degree="1">0.0 0.0 0.5 1.0 1.0</KnotVector>\n',
        '        </Basis>\n',
        '        <Basis type="BSplineBasis" index="2">\n',
        (
                '          <KnotVector degree="2">0.0 0.'
                '0 0.0 1.0 1.0 1.0</KnotVector>\n'
        ),
        '        </Basis>\n',
        '      </Basis>\n',
        '      <weights>1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0</weights>\n',
        '    </Basis>\n',
        '    <coefs geoDim="3">1.0 1.0 0.0\n',
        '2.0 1.0 0.0\n',
        '1.0 1.5 0.0\n',
        '2.0 1.5 0.0\n',
        '1.0 2.0 0.0\n',
        '2.0 2.0 0.0\n',
        '1.0 1.0 0.5\n',
        '2.0 1.0 0.5\n',
        '1.0 1.5 0.5\n',
        '2.0 1.5 0.5\n',
        '1.0 2.0 0.5\n',
        '2.0 2.0 0.5\n',
        '1.0 1.0 1.0\n',
        '2.0 1.0 1.0\n',
        '1.0 1.5 1.0\n',
        '2.0 1.5 1.0\n',
        '1.0 2.0 1.0\n',
        '2.0 2.0 1.0</coefs>\n',
        '  </Geometry>\n',
        '  <Geometry type="TensorNurbs3" id="103">\n',
        '    <Basis type="TensorNurbsBasis3">\n',
        '      <Basis type="TensorBSplineBasis3">\n',
        '        <Basis type="BSplineBasis" index="0">\n',
        '          <KnotVector degree="1">0.0 0.0 1.0 1.0</KnotVector>\n',
        '        </Basis>\n',
        '        <Basis type="BSplineBasis" index="1">\n',
        '          <KnotVector degree="1">0.0 0.0 1.0 1.0</KnotVector>\n',
        '        </Basis>\n',
        '        <Basis type="BSplineBasis" index="2">\n',
        (
                '          <KnotVector degree="2">0.0 0.0 0.0'
                ' 1.0 1.0 1.0</KnotVector>\n'
        ),
        '        </Basis>\n',
        '      </Basis>\n',
        '      <weights>1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0</weights>\n',
        '    </Basis>\n',
        '    <coefs geoDim="3">1.0 0.0 0.0\n',
        '2.0 0.0 0.0\n',
        '1.0 1.0 0.0\n',
        '2.0 1.0 0.0\n',
        '1.0 0.0 0.5\n',
        '2.0 0.0 0.5\n',
        '1.0 1.0 0.5\n',
        '2.0 1.0 0.5\n',
        '1.0 0.0 1.0\n',
        '2.0 0.0 1.0\n',
        '1.0 1.0 1.0\n',
        '2.0 1.0 1.0</coefs>\n',
        '  </Geometry>\n',
        '</xml>'
]

_gismo_export_ref_2d = [
        '<xml>\n',
        '  <MultiPatch parDim="2" id="0">\n',
        '    <patches type="id_range">100 103</patches>\n',
        '    <interfaces>101 1 100 2 0 1 1 1\n',
        '102 3 100 4 0 1 1 1\n',
        '103 3 101 4 0 1 1 1\n',
        '103 1 102 2 0 1 1 1</interfaces>\n',
        '    <boundary>100 1\n',
        '100 3\n',
        '101 2\n',
        '101 3\n',
        '102 1\n',
        '102 4\n',
        '103 2\n',
        '103 4</boundary>\n',
        '  </MultiPatch>\n',
        '  <boundaryConditions multipatch="0" id="1">\n',
        '    <!--Please fill boundary conditions here-->\n',
        '    <bc type="Dirichlet" unknown="0">0 3\n',
        '1 3\n',
        '2 4\n',
        '3 4</bc>\n',
        '    <bc type="Dirichlet" unknown="0">0 1\n',
        '2 1</bc>\n',
        '  </boundaryConditions>\n',
        '  <Geometry type="TensorBSpline2" id="100">\n',
        '    <Basis type="TensorBSplineBasis2">\n',
        '      <Basis type="BSplineBasis" index="0">\n',
        (
                '        <KnotVector degree="2">0.0 0.0 '
                '0.0 1.0 1.0 1.0</KnotVector>\n'
        ),
        '      </Basis>\n',
        '      <Basis type="BSplineBasis" index="1">\n',
        '        <KnotVector degree="1">0.0 0.0 1.0 1.0</KnotVector>\n',
        '      </Basis>\n',
        '    </Basis>\n',
        '    <coefs geoDim="2">0.0 0.0\n',
        '0.5 0.0\n',
        '1.0 0.0\n',
        '0.0 1.0\n',
        '0.5 1.0\n',
        '1.0 1.0</coefs>\n',
        '  </Geometry>\n',
        '  <Geometry type="TensorNurbs2" id="101">\n',
        '    <Basis type="TensorNurbsBasis2">\n',
        '      <Basis type="TensorBSplineBasis2">\n',
        '        <Basis type="BSplineBasis" index="0">\n',
        '          <KnotVector degree="1">0.0 0.0 1.0 1.0</KnotVector>\n',
        '        </Basis>\n',
        '        <Basis type="BSplineBasis" index="1">\n',
        '          <KnotVector degree="1">0.0 0.0 1.0 1.0</KnotVector>\n',
        '        </Basis>\n',
        '      </Basis>\n',
        '      <weights>1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0</weights>\n',
        '    </Basis>\n',
        '    <coefs geoDim="2">1.0 0.0\n',
        '2.0 0.0\n',
        '1.0 1.0\n',
        '2.0 1.0</coefs>\n',
        '  </Geometry>\n',
        '  <Geometry type="TensorBSpline2" id="102">\n',
        '    <Basis type="TensorBSplineBasis2">\n',
        '      <Basis type="BSplineBasis" index="0">\n',
        (
                '        <KnotVector degree="2">0.0 0.0 0.0'
                ' 1.0 1.0 1.0</KnotVector>\n'
        ),
        '      </Basis>\n',
        '      <Basis type="BSplineBasis" index="1">\n',
        '        <KnotVector degree="1">0.0 0.0 0.5 1.0 1.0</KnotVector>\n',
        '      </Basis>\n',
        '    </Basis>\n',
        '    <coefs geoDim="2">0.0 1.0\n',
        '0.5 1.0\n',
        '1.0 1.0\n',
        '0.0 1.5\n',
        '0.5 1.5\n',
        '1.0 1.5\n',
        '0.0 2.0\n',
        '0.5 2.0\n',
        '1.0 2.0</coefs>\n',
        '  </Geometry>\n',
        '  <Geometry type="TensorNurbs2" id="103">\n',
        '    <Basis type="TensorNurbsBasis2">\n',
        '      <Basis type="TensorBSplineBasis2">\n',
        '        <Basis type="BSplineBasis" index="0">\n',
        '          <KnotVector degree="1">0.0 0.0 1.0 1.0</KnotVector>\n',
        '        </Basis>\n',
        '        <Basis type="BSplineBasis" index="1">\n',
        '          <KnotVector degree="1">0.0 0.0 0.5 1.0 1.0</KnotVector>\n',
        '        </Basis>\n',
        '      </Basis>\n',
        '      <weights>1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0\n',
        '1.0</weights>\n',
        '    </Basis>\n',
        '    <coefs geoDim="2">1.0 1.0\n',
        '2.0 1.0\n',
        '1.0 1.5\n',
        '2.0 1.5\n',
        '1.0 2.0\n',
        '2.0 2.0</coefs>\n',
        '  </Geometry>\n',
        '</xml>'
]


class gismoExportTest(c.unittest.TestCase):

    def test_gismo_export(self):
        """
        Test gismo export routine
        """
        from sys import version as python_version

        # Define some splines
        bez_el0 = c.splinepy.Bezier(
                degrees=[1, 1],
                control_points=[[0, 0], [1, 0], [0, 1], [1, 1]]
        )
        rbz_el1 = c.splinepy.RationalBezier(
                degrees=[1, 1],
                control_points=[[1, 0], [2, 0], [1, 1], [2, 1]],
                weights=[1, 1, 1, 1]
        )
        bsp_el2 = c.splinepy.BSpline(
                degrees=[1, 1],
                control_points=[[0, 1], [1, 1], [0, 2], [1, 2]],
                knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1]]
        )
        nur_el3 = c.splinepy.NURBS(
                degrees=[1, 1],
                control_points=[[1, 1], [2, 1], [1, 2], [2, 2]],
                weights=[1, 1, 1, 1],
                knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1]]
        )

        # Make it more tricky
        bez_el0.elevate_degree(0)
        bsp_el2.elevate_degree(0)

        bsp_el2.insert_knots(1, [0.5])
        nur_el3.insert_knots(1, [0.5])

        # Init multipatch
        multipatch = c.splinepy.Multipatch(
                splines=[bez_el0, rbz_el1, bsp_el2, nur_el3]
        )

        # Define some functions for boundary identification
        def is_bottom(x):
            return x[:, 0] < 0.01

        def is_top(x):
            return x[:, 0] > 1.99

        # Add boundary
        multipatch.add_boundary_with_function(is_bottom)
        multipatch.add_boundary_with_function(is_top)

        # Test Output
        with tempfile.NamedTemporaryFile() as tmpf:
            c.splinepy.io.gismo.export(tmpf.name, multipatch)

            with open(tmpf.name, "r") as tmp_read:
                if int(python_version.split('.')[1]) >= 9:
                    self.assertTrue(
                            c.are_items_same(
                                    _gismo_export_ref_2d, tmp_read.readlines()
                            )
                    )
                else:
                    self.assertTrue(
                            c.are_string_items_same(
                                    _gismo_export_ref_2d, tmp_read.readlines()
                            )
                    )

        # Test Also 3D Meshes
        bez_el0 = c.splinepy.Bezier(
                degrees=[1, 1, 1],
                control_points=[
                        [0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1],
                        [1, 0, 1], [0, 1, 1], [1, 1, 1]
                ]
        )
        rbz_el1 = c.splinepy.RationalBezier(
                degrees=[1, 1, 1],
                control_points=[
                        [1, 0, 0], [2, 0, 0], [1, 1, 0], [2, 1, 0], [1, 0, 1],
                        [2, 0, 1], [1, 1, 1], [2, 1, 1]
                ],
                weights=[1] * 8
        )
        bsp_el2 = c.splinepy.BSpline(
                degrees=[1, 1, 1],
                control_points=[
                        [0, 1, 0], [1, 1, 0], [0, 2, 0], [1, 2, 0], [0, 1, 1],
                        [1, 1, 1], [0, 2, 1], [1, 2, 1]
                ],
                knot_vectors=[[0, 0, 1, 1]] * 3
        )
        nur_el3 = c.splinepy.NURBS(
                degrees=[1, 1, 1],
                control_points=[
                        [1, 1, 0], [2, 1, 0], [1, 2, 0], [2, 2, 0], [1, 1, 1],
                        [2, 1, 1], [1, 2, 1], [2, 2, 1]
                ],
                weights=[1] * 8,
                knot_vectors=[[0, 0, 1, 1]] * 3
        )

        # Make it more tricky
        bez_el0.elevate_degree(0)
        bsp_el2.elevate_degree(0)
        bsp_el2.insert_knots(1, [0.5])
        nur_el3.insert_knots(1, [0.5])
        for s in [bez_el0, bsp_el2, nur_el3, rbz_el1]:
            s.elevate_degree(2)

        # Define some functions for boundary identification
        def is_bottom(x):
            return x[:, 0] < 0.01

        def is_top(x):
            return x[:, 0] > 1.99

        # Add boundary
        multipatch = c.splinepy.Multipatch(
                splines=[bez_el0, bsp_el2, nur_el3, rbz_el1]
        )
        multipatch.add_boundary_with_function(is_bottom)
        multipatch.add_boundary_with_function(is_top)

        # Test output
        with tempfile.NamedTemporaryFile() as tmpf:
            c.splinepy.io.gismo.export(tmpf.name, multipatch)
            with open(tmpf.name, "r") as tmp_read:
                if int(python_version.split('.')[1]) >= 9:
                    self.assertTrue(
                            c.are_items_same(
                                    _gismo_export_ref_3d, tmp_read.readlines()
                            )
                    )
                else:
                    self.assertTrue(
                            c.are_string_items_same(
                                    _gismo_export_ref_3d, tmp_read.readlines()
                            )
                    )


if __name__ == "__main__":
    c.unittest.main()
