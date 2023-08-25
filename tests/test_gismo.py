import tempfile
from sys import version as python_version

try:
    from . import common as c
except BaseException:
    import common as c

_gismo_export_ref_2d = [
    (
        '<xml><MultiPatch parDim="2" id="0"><patches type="id_range">100 '
        "103</patches><interfaces>100 2 101 1 0 1 1 1\n"
    ),
    "100 4 102 3 0 1 1 1\n",
    "101 4 103 3 0 1 1 1\n",
    "102 2 103 1 0 1 1 1</interfaces><boundary>100 1\n",
    "100 3\n",
    "101 2\n",
    "101 3\n",
    "102 1\n",
    "102 4\n",
    "103 2\n",
    (
        '103 4</boundary></MultiPatch><boundaryConditions multipatch="4" id="1'
        '"><!--Please fill boundary conditions here--><bc type="Dirichlet" unk'
        'nown="0">0 3\n'
    ),
    "1 3\n",
    "2 4\n",
    '3 4</bc><bc type="Dirichlet" unknown="0">0 1\n',
    '2 1</bc><bc type="Dirichlet" unknown="0">1 2\n',
    (
        '3 2</bc></boundaryConditions><Geometry type="TensorBSpline2" id="100"'
        '><Basis type="TensorBSplineBasis2"><Basis type="BSplineBasis" index="'
        '0"><KnotVector degree="2" format="ASCII">0.0 0.0 0.0 1.0 1.0 1.0</Kno'
        'tVector></Basis><Basis type="BSplineBasis" index="1"><KnotVector degr'
        'ee="1" format="ASCII">0.0 0.0 '
        "1.0 1.0</KnotVector></Basis></Basis><"
        'coefs geoDim="2" format="ASCII">0.0 0.0\n'
    ),
    "0.5 0.0\n",
    "1.0 0.0\n",
    "0.0 1.0\n",
    "0.5 1.0\n",
    (
        '1.0 1.0</coefs></Geometry><Geometry type="TensorNurbs2" id="101"><Bas'
        'is type="TensorNurbsBasis2"><Basis type="TensorBSplineBasis2"><Basis '
        'type="BSplineBasis" index="0"><KnotVector degree="1" format="ASCII">0'
        ".0 0.0 1.0 1.0<"
        '/KnotVector></Basis><Basis type="BSplineBasis" index="1"><KnotVector '
        'degree="1" format="ASCII">0.0 0.0 1.0 1.0</KnotVector></Basis></Basis'
        '><weights format="ASCII">1.0\n'
    ),
    "1.0\n",
    "1.0\n",
    '1.0</weights></Basis><coefs geoDim="2" format="ASCII">1.0 0.0\n',
    "2.0 0.0\n",
    "1.0 1.0\n",
    (
        '2.0 1.0</coefs></Geometry><Geometry type="TensorBSpline2" id="102"><B'
        'asis type="TensorBSplineBasis2"><Basis type="BSplineBasis" index="0">'
        '<KnotVector degree="2" format="ASCII">0.0 0.0 0.0 1.0 1.0 1.0</KnotVe'
        "ctor></Basis><B"
        'asis type="BSplineBasis" index="1"><KnotVector degree="1" format="ASC'
        'II">0.0 0.0 0.5'
        ' 1.0 1.0</KnotVector></Basis></Basis><coefs geoDim="2" format="ASCII"'
        ">0.0 1.0\n"
    ),
    "0.5 1.0\n",
    "1.0 1.0\n",
    "0.0 1.5\n",
    "0.5 1.5\n",
    "1.0 1.5\n",
    "0.0 2.0\n",
    "0.5 2.0\n",
    (
        '1.0 2.0</coefs></Geometry><Geometry type="TensorNurbs2" id="103"><Bas'
        'is type="TensorNurbsBasis2"><Basis type="TensorBSplineBasis2"><Basis '
        'type="BSplineBasis" index="0"><KnotVector degree="1" format="ASCII">0'
        ".0 0.0 1.0 1.0<"
        '/KnotVector></Basis><Basis type="BSplineBasis" index="1"><KnotVector '
        'degree="1" format="ASCII">0.0 0.0 0.5 1.0 1.0</KnotVector></Basis></B'
        "asis><weights f"
        'ormat="ASCII">1'
        ".0\n"
    ),
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    '1.0</weights></Basis><coefs geoDim="2" format="ASCII">1.0 1.0\n',
    "2.0 1.0\n",
    "1.0 1.5\n",
    "2.0 1.5\n",
    "1.0 2.0\n",
    "2.0 2.0</coefs></Geometry></xml>\n",
]

_gismo_export_ref_2d_indent = [
    "<xml>\n",
    '  <MultiPatch parDim="2" id="0">\n',
    '    <patches type="id_range">100 103</patches>\n',
    "    <interfaces>100 2 101 1 0 1 1 1\n",
    "100 4 102 3 0 1 1 1\n",
    "101 4 103 3 0 1 1 1\n",
    "102 2 103 1 0 1 1 1</interfaces>\n",
    "    <boundary>100 1\n",
    "100 3\n",
    "101 2\n",
    "101 3\n",
    "102 1\n",
    "102 4\n",
    "103 2\n",
    "103 4</boundary>\n",
    "  </MultiPatch>\n",
    '  <boundaryConditions multipatch="4" id="1">\n',
    "    <!--Please fill boundary conditions here-->\n",
    '    <bc type="Dirichlet" unknown="0">0 3\n',
    "1 3\n",
    "2 4\n",
    "3 4</bc>\n",
    '    <bc type="Dirichlet" unknown="0">0 1\n',
    "2 1</bc>\n",
    '    <bc type="Dirichlet" unknown="0">1 2\n',
    "3 2</bc>\n",
    "  </boundaryConditions>\n",
    '  <Geometry type="TensorBSpline2" id="100">\n',
    '    <Basis type="TensorBSplineBasis2">\n',
    '      <Basis type="BSplineBasis" index="0">\n',
    '        <KnotVector degree="2" format="ASCII">0.0 0.0 0.0 1.0 1.0 1.0</Kn'
    "otVector>\n",
    "      </Basis>\n",
    '      <Basis type="BSplineBasis" index="1">\n',
    (
        '        <KnotVector degree="1" format="ASCII">0.0 0.0 1.0 1.0</KnotVe'
        "ctor>\n"
    ),
    "      </Basis>\n",
    "    </Basis>\n",
    '    <coefs geoDim="2" format="ASCII">0.0 0.0\n',
    "0.5 0.0\n",
    "1.0 0.0\n",
    "0.0 1.0\n",
    "0.5 1.0\n",
    "1.0 1.0</coefs>\n",
    "  </Geometry>\n",
    '  <Geometry type="TensorNurbs2" id="101">\n',
    '    <Basis type="TensorNurbsBasis2">\n',
    '      <Basis type="TensorBSplineBasis2">\n',
    '        <Basis type="BSplineBasis" index="0">\n',
    (
        '          <KnotVector degree="1" format="ASCII">0.0 0.0 1.0 1.0</KnotVector>\n'
    ),
    "        </Basis>\n",
    '        <Basis type="BSplineBasis" index="1">\n',
    (
        '          <KnotVector degree="1" format="ASCII">0.0 0.0 1.0 1.0</KnotVector>\n'
    ),
    "        </Basis>\n",
    "      </Basis>\n",
    '      <weights format="ASCII">1.0\n',
    "1.0\n",
    "1.0\n",
    "1.0</weights>\n",
    "    </Basis>\n",
    '    <coefs geoDim="2" format="ASCII">1.0 0.0\n',
    "2.0 0.0\n",
    "1.0 1.0\n",
    "2.0 1.0</coefs>\n",
    "  </Geometry>\n",
    '  <Geometry type="TensorBSpline2" id="102">\n',
    '    <Basis type="TensorBSplineBasis2">\n',
    '      <Basis type="BSplineBasis" index="0">\n',
    (
        '        <KnotVector degree="2" format="ASCII">0.0 0.0 0.0 1.0 1.0 1.0'
        "</KnotVector>\n"
    ),
    "      </Basis>\n",
    '      <Basis type="BSplineBasis" index="1">\n',
    (
        '        <KnotVector degree="1" format="ASCII">0.0 0.0 0.5 1.0 1.0</Kn'
        "otVector>\n"
    ),
    "      </Basis>\n",
    "    </Basis>\n",
    '    <coefs geoDim="2" format="ASCII">0.0 1.0\n',
    "0.5 1.0\n",
    "1.0 1.0\n",
    "0.0 1.5\n",
    "0.5 1.5\n",
    "1.0 1.5\n",
    "0.0 2.0\n",
    "0.5 2.0\n",
    "1.0 2.0</coefs>\n",
    "  </Geometry>\n",
    '  <Geometry type="TensorNurbs2" id="103">\n',
    '    <Basis type="TensorNurbsBasis2">\n',
    '      <Basis type="TensorBSplineBasis2">\n',
    '        <Basis type="BSplineBasis" index="0">\n',
    (
        '          <KnotVector degree="1" format="ASCII">0.0 0.0 1.0 1.0</Knot'
        "Vector>\n"
    ),
    "        </Basis>\n",
    '        <Basis type="BSplineBasis" index="1">\n',
    (
        '          <KnotVector degree="1" format="ASCII">0.0 0.0 0.5 1.0 1.0</'
        "KnotVector>\n"
    ),
    "        </Basis>\n",
    "      </Basis>\n",
    '      <weights format="ASCII">1.0\n',
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0</weights>\n",
    "    </Basis>\n",
    '    <coefs geoDim="2" format="ASCII">1.0 1.0\n',
    "2.0 1.0\n",
    "1.0 1.5\n",
    "2.0 1.5\n",
    "1.0 2.0\n",
    "2.0 2.0</coefs>\n",
    "  </Geometry>\n",
    "</xml>",
]

_gismo_export_ref_2d_labeled = [
    (
        '<xml><MultiPatch parDim="2" id="0"><patches type="id_range">100 103</'
        "patches><interfaces>100 2 101 1 0 1 1 1\n"
    ),
    "100 4 102 3 0 1 1 1\n",
    "101 4 103 3 0 1 1 1\n",
    '102 2 103 1 0 1 1 1</interfaces><boundary name="BID1">100 3\n',
    "101 3\n",
    "102 4\n",
    '103 4</boundary><boundary name="BID2">100 1\n',
    '102 1</boundary><boundary name="BID3">101 2\n',
    (
        '103 2</boundary></MultiPatch><Geometry type="TensorBSpline2" id="100"'
        '><Basis type="TensorBSplineBasis2"><Basis type="BSplineBasis" index="'
        '0"><KnotVector degree="2" format="ASCII">0.0 0.0 0.0 1.0 1.0 1.0</Kno'
        'tVector></Basis><Basis type="BSplineBasis" index="1"><KnotVector degr'
        'ee="1" format="ASCII">0.0 0.0 '
        '1.0 1.0</KnotVector></Basis></Basis><coefs geoDim="2" format="ASCII">'
        "0.0 0.0\n"
    ),
    "0.5 0.0\n",
    "1.0 0.0\n",
    "0.0 1.0\n",
    "0.5 1.0\n",
    (
        '1.0 1.0</coefs></Geometry><Geometry type="TensorNurbs2" id="101"><Bas'
        'is type="TensorNurbsBasis2"><Basis type="TensorBSplineBasis2"><Basis '
        'type="BSplineBasis" index="0"><KnotVector degree="1" format="ASCII">0'
        '.0 0.0 1.0 1.0</KnotVector></Basis><Basis type="BSplineBasis" index="'
        '1"><KnotVector degree="1" format="ASCII">0.0 0.0 1.0 1.0</KnotVector>'
        '</Basis></Basis><weights format="ASCII">1.0\n'
    ),
    "1.0\n",
    "1.0\n",
    '1.0</weights></Basis><coefs geoDim="2" format="ASCII">1.0 0.0\n',
    "2.0 0.0\n",
    "1.0 1.0\n",
    (
        '2.0 1.0</coefs></Geometry><Geometry type="TensorBSpline2" id="102"><B'
        'asis type="TensorBSplineBasis2"><Basis type="BSplineBasis" index="0">'
        '<KnotVector degree="2" format="ASCII">0.0 0.0 0.0 1.0 1.0 1.0</KnotVe'
        'ctor></Basis><Basis type="BSplineBasis" index="1"><KnotVector degree='
        '"1" format="ASCII">0.0 0.0 0.5 1.0 1.0</KnotVector></Basis></Basis><c'
        'oefs geoDim="2" format="ASCII">0.0 1.0\n'
    ),
    "0.5 1.0\n",
    "1.0 1.0\n",
    "0.0 1.5\n",
    "0.5 1.5\n",
    "1.0 1.5\n",
    "0.0 2.0\n",
    "0.5 2.0\n",
    (
        '1.0 2.0</coefs></Geometry><Geometry type="TensorNurbs2" id="103"><Bas'
        'is type="TensorNurbsBasis2"><Basis type="TensorBSplineBasis2"><Basis '
        'type="BSplineBasis" index="0"><KnotVector degree="1" format="ASCII">0'
        '.0 0.0 1.0 1.0</KnotVector></Basis><Basis type="BSplineBasis" index="'
        '1"><KnotVector degree="1" format="ASCII">0.0 0.0 0.5 1.0 1.0</KnotVec'
        'tor></Basis></Basis><weights format="ASCII">1.0\n'
    ),
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    '1.0</weights></Basis><coefs geoDim="2" format="ASCII">1.0 1.0\n',
    "2.0 1.0\n",
    "1.0 1.5\n",
    "2.0 1.5\n",
    "1.0 2.0\n",
    "2.0 2.0</coefs></Geometry></xml>",
]

_gismo_export_ref_2d_indent_labeled = [
    "<xml>\n",
    '  <MultiPatch parDim="2" id="0">\n',
    '    <patches type="id_range">100 103</patches>\n',
    "    <interfaces>100 2 101 1 0 1 1 1\n",
    "100 4 102 3 0 1 1 1\n",
    "101 4 103 3 0 1 1 1\n",
    "102 2 103 1 0 1 1 1</interfaces>\n",
    '    <boundary name="BID1">100 3\n',
    "101 3\n",
    "102 4\n",
    "103 4</boundary>\n",
    '    <boundary name="BID2">100 1\n',
    "102 1</boundary>\n",
    '    <boundary name="BID3">101 2\n',
    "103 2</boundary>\n",
    "  </MultiPatch>\n",
    '  <Geometry type="TensorBSpline2" id="100">\n',
    '    <Basis type="TensorBSplineBasis2">\n',
    '      <Basis type="BSplineBasis" index="0">\n',
    (
        '        <KnotVector degree="2" format="ASCII">0.0 0.0 0.0 1.0 1.0 1.0'
        "</KnotVector>\n"
    ),
    "      </Basis>\n",
    '      <Basis type="BSplineBasis" index="1">\n',
    (
        '        <KnotVector degree="1" format="ASCII">0.0 0.0 1.0 1.0</KnotVe'
        "ctor>\n"
    ),
    "      </Basis>\n",
    "    </Basis>\n",
    '    <coefs geoDim="2" format="ASCII">0.0 0.0\n',
    "0.5 0.0\n",
    "1.0 0.0\n",
    "0.0 1.0\n",
    "0.5 1.0\n",
    "1.0 1.0</coefs>\n",
    "  </Geometry>\n",
    '  <Geometry type="TensorNurbs2" id="101">\n',
    '    <Basis type="TensorNurbsBasis2">\n',
    '      <Basis type="TensorBSplineBasis2">\n',
    '        <Basis type="BSplineBasis" index="0">\n',
    (
        '          <KnotVector degree="1" format="ASCII">0.0 0.0 1.0 1.0</KnotVector>\n'
    ),
    "        </Basis>\n",
    '        <Basis type="BSplineBasis" index="1">\n',
    (
        '          <KnotVector degree="1" format="ASCII">0.0 0.0 1.0 1.0</KnotVector>\n'
    ),
    "        </Basis>\n",
    "      </Basis>\n",
    '      <weights format="ASCII">1.0\n',
    "1.0\n",
    "1.0\n",
    "1.0</weights>\n",
    "    </Basis>\n",
    '    <coefs geoDim="2" format="ASCII">1.0 0.0\n',
    "2.0 0.0\n",
    "1.0 1.0\n",
    "2.0 1.0</coefs>\n",
    "  </Geometry>\n",
    '  <Geometry type="TensorBSpline2" id="102">\n',
    '    <Basis type="TensorBSplineBasis2">\n',
    '      <Basis type="BSplineBasis" index="0">\n',
    (
        '        <KnotVector degree="2" format="ASCII">0.0 0.0 0.0 1.0 1.0 1.0'
        "</KnotVector>\n"
    ),
    "      </Basis>\n",
    '      <Basis type="BSplineBasis" index="1">\n',
    (
        '        <KnotVector degree="1" format="ASCII">0.0 0.0 0.5 1.0 1.0</Kn'
        "otVector>\n"
    ),
    "      </Basis>\n",
    "    </Basis>\n",
    '    <coefs geoDim="2" format="ASCII">0.0 1.0\n',
    "0.5 1.0\n",
    "1.0 1.0\n",
    "0.0 1.5\n",
    "0.5 1.5\n",
    "1.0 1.5\n",
    "0.0 2.0\n",
    "0.5 2.0\n",
    "1.0 2.0</coefs>\n",
    "  </Geometry>\n",
    '  <Geometry type="TensorNurbs2" id="103">\n',
    '    <Basis type="TensorNurbsBasis2">\n',
    '      <Basis type="TensorBSplineBasis2">\n',
    '        <Basis type="BSplineBasis" index="0">\n',
    (
        '          <KnotVector degree="1" format="ASCII">0.0 0.0 1.0 1.0</Knot'
        "Vector>\n"
    ),
    "        </Basis>\n",
    '        <Basis type="BSplineBasis" index="1">\n',
    (
        '          <KnotVector degree="1" format="ASCII">0.0 0.0 0.5 1.0 1.0</'
        "KnotVector>\n"
    ),
    "        </Basis>\n",
    "      </Basis>\n",
    '      <weights format="ASCII">1.0\n',
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0</weights>\n",
    "    </Basis>\n",
    '    <coefs geoDim="2" format="ASCII">1.0 1.0\n',
    "2.0 1.0\n",
    "1.0 1.5\n",
    "2.0 1.5\n",
    "1.0 2.0\n",
    "2.0 2.0</coefs>\n",
    "  </Geometry>\n",
    "</xml>",
]

_gismo_export_ref_3d = [
    (
        '<xml><MultiPatch parDim="3" id="0"><patches type="id_range">100 103</'
        "patches><interfaces>100 2 103 1 0 1 2 1 1 1\n"
    ),
    "100 4 101 3 0 1 2 1 1 1\n",
    "101 2 102 1 0 1 2 1 1 1\n",
    "102 3 103 4 0 1 2 1 1 1</interfaces><boundary>100 1\n",
    "100 3\n",
    "100 5\n",
    "100 6\n",
    "101 1\n",
    "101 4\n",
    "101 5\n",
    "101 6\n",
    "102 2\n",
    "102 4\n",
    "102 5\n",
    "102 6\n",
    "103 2\n",
    "103 3\n",
    "103 5\n",
    (
        '103 6</boundary></MultiPatch><boundaryConditions multipatch="4" id="1'
        '"><!--Please fill boundary conditions here--><bc type="Dirichlet" unk'
        'nown="0">0 3\n'
    ),
    "0 5\n",
    "0 6\n",
    "1 4\n",
    "1 5\n",
    "1 6\n",
    "2 4\n",
    "2 5\n",
    "2 6\n",
    "3 3\n",
    "3 5\n",
    '3 6</bc><bc type="Dirichlet" unknown="0">0 1\n',
    '1 1</bc><bc type="Dirichlet" unknown="0">2 2\n',
    (
        '3 2</bc></boundaryConditions><Geometry type="TensorBSpline3" id="100"'
        '><Basis type="TensorBSplineBasis3"><Basis type="BSplineBasis" index="'
        '0"><KnotVector degree="2" format="ASCII">0.0 0.0 0.0 1.0 1.0 1.0</Kno'
        'tVector></Basis><Basis type="BSplineBasis" index="1"><KnotVector degr'
        'ee="1" format="ASCII">0.0 0.0 1.0 1.0</KnotVector></Basis><Basis type'
        '="BSplineBasis" index="2"><KnotVector degree="2" format="ASCII">0.0 0'
        '.0 0.0 1.0 1.0 1.0</KnotVector></Basis></Basis><coefs geoDim="3" form'
        'at="ASCII">0.0 0.0 0.0\n'
    ),
    "0.5 0.0 0.0\n",
    "1.0 0.0 0.0\n",
    "0.0 1.0 0.0\n",
    "0.5 1.0 0.0\n",
    "1.0 1.0 0.0\n",
    "0.0 0.0 0.5\n",
    "0.5 0.0 0.5\n",
    "1.0 0.0 0.5\n",
    "0.0 1.0 0.5\n",
    "0.5 1.0 0.5\n",
    "1.0 1.0 0.5\n",
    "0.0 0.0 1.0\n",
    "0.5 0.0 1.0\n",
    "1.0 0.0 1.0\n",
    "0.0 1.0 1.0\n",
    "0.5 1.0 1.0\n",
    (
        '1.0 1.0 1.0</coefs></Geometry><Geometry type="TensorBSpline3" id="101'
        '"><Basis type="TensorBSplineBasis3"><Basis type="BSplineBasis" index='
        '"0"><KnotVector degree="2" format="ASCII">0.0 0.0 0.0 1.0 1.0 1.0</Kn'
        'otVector></Basis><Basis type="BSplineBasis" index="1"><KnotVector deg'
        'ree="1" format="ASCII">0.0 0.0 0.5 1.0 1.0</KnotVector></Basis><Basis'
        ' type="BSplineBasis" index="2"><KnotVector degree="2" format="ASCII">'
        '0.0 0.0 0.0 1.0 1.0 1.0</KnotVector></Basis></Basis><coefs geoDim="3"'
        ' format="ASCII">0.0 1.0 0.0\n'
    ),
    "0.5 1.0 0.0\n",
    "1.0 1.0 0.0\n",
    "0.0 1.5 0.0\n",
    "0.5 1.5 0.0\n",
    "1.0 1.5 0.0\n",
    "0.0 2.0 0.0\n",
    "0.5 2.0 0.0\n",
    "1.0 2.0 0.0\n",
    "0.0 1.0 0.5\n",
    "0.5 1.0 0.5\n",
    "1.0 1.0 0.5\n",
    "0.0 1.5 0.5\n",
    "0.5 1.5 0.5\n",
    "1.0 1.5 0.5\n",
    "0.0 2.0 0.5\n",
    "0.5 2.0 0.5\n",
    "1.0 2.0 0.5\n",
    "0.0 1.0 1.0\n",
    "0.5 1.0 1.0\n",
    "1.0 1.0 1.0\n",
    "0.0 1.5 1.0\n",
    "0.5 1.5 1.0\n",
    "1.0 1.5 1.0\n",
    "0.0 2.0 1.0\n",
    "0.5 2.0 1.0\n",
    (
        '1.0 2.0 1.0</coefs></Geometry><Geometry type="TensorNurbs3" id="102">'
        '<Basis type="TensorNurbsBasis3"><Basis type="TensorBSplineBasis3"><Ba'
        'sis type="BSplineBasis" index="0"><KnotVector degree="1" format="ASCI'
        'I">0.0 0.0 1.0 1.0</KnotVector></Basis><Basis type="BSplineBasis" ind'
        'ex="1"><KnotVector degree="1" format="ASCII">0.0 0.0 0.5 1.0 1.0</Kno'
        'tVector></Basis><Basis type="BSplineBasis" index="2"><KnotVector degr'
        'ee="2" format="ASCII">0.0 0.0 0.0 1.0 1.0 1.0</KnotVector></Basis></B'
        'asis><weights format="ASCII">1.0\n'
    ),
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    '1.0</weights></Basis><coefs geoDim="3" format="ASCII">1.0 1.0 0.0\n',
    "2.0 1.0 0.0\n",
    "1.0 1.5 0.0\n",
    "2.0 1.5 0.0\n",
    "1.0 2.0 0.0\n",
    "2.0 2.0 0.0\n",
    "1.0 1.0 0.5\n",
    "2.0 1.0 0.5\n",
    "1.0 1.5 0.5\n",
    "2.0 1.5 0.5\n",
    "1.0 2.0 0.5\n",
    "2.0 2.0 0.5\n",
    "1.0 1.0 1.0\n",
    "2.0 1.0 1.0\n",
    "1.0 1.5 1.0\n",
    "2.0 1.5 1.0\n",
    "1.0 2.0 1.0\n",
    (
        '2.0 2.0 1.0</coefs></Geometry><Geometry type="TensorNurbs3" id="103">'
        '<Basis type="TensorNurbsBasis3"><Basis type="TensorBSplineBasis3"><Ba'
        'sis type="BSplineBasis" index="0"><KnotVector degree="1" format="ASCI'
        'I">0.0 0.0 1.0 1.0</KnotVector></Basis><Basis type="BSplineBasis" ind'
        'ex="1"><KnotVector degree="1" format="ASCII">0.0 0.0 1.0 1.0</KnotVec'
        'tor></Basis><Basis type="BSplineBasis" index="2"><KnotVector degree="'
        '2" format="ASCII">0.0 0.0 0.0 1.0 1.0 1.0</KnotVector></Basis></Basis'
        '><weights format="ASCII">1.0\n'
    ),
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    '1.0</weights></Basis><coefs geoDim="3" format="ASCII">1.0 0.0 0.0\n',
    "2.0 0.0 0.0\n",
    "1.0 1.0 0.0\n",
    "2.0 1.0 0.0\n",
    "1.0 0.0 0.5\n",
    "2.0 0.0 0.5\n",
    "1.0 1.0 0.5\n",
    "2.0 1.0 0.5\n",
    "1.0 0.0 1.0\n",
    "2.0 0.0 1.0\n",
    "1.0 1.0 1.0\n",
    "2.0 1.0 1.0</coefs></Geometry></xml>",
]

_gismo_export_ref_3d_indent = [
    "<xml>\n",
    '  <MultiPatch parDim="3" id="0">\n',
    '    <patches type="id_range">100 103</patches>\n',
    "    <interfaces>100 2 103 1 0 1 2 1 1 1\n",
    "100 4 101 3 0 1 2 1 1 1\n",
    "101 2 102 1 0 1 2 1 1 1\n",
    "102 3 103 4 0 1 2 1 1 1</interfaces>\n",
    "    <boundary>100 1\n",
    "100 3\n",
    "100 5\n",
    "100 6\n",
    "101 1\n",
    "101 4\n",
    "101 5\n",
    "101 6\n",
    "102 2\n",
    "102 4\n",
    "102 5\n",
    "102 6\n",
    "103 2\n",
    "103 3\n",
    "103 5\n",
    "103 6</boundary>\n",
    "  </MultiPatch>\n",
    '  <boundaryConditions multipatch="4" id="1">\n',
    "    <!--Please fill boundary conditions here-->\n",
    '    <bc type="Dirichlet" unknown="0">0 3\n',
    "0 5\n",
    "0 6\n",
    "1 4\n",
    "1 5\n",
    "1 6\n",
    "2 4\n",
    "2 5\n",
    "2 6\n",
    "3 3\n",
    "3 5\n",
    "3 6</bc>\n",
    '    <bc type="Dirichlet" unknown="0">0 1\n',
    "1 1</bc>\n",
    '    <bc type="Dirichlet" unknown="0">2 2\n',
    "3 2</bc>\n",
    "  </boundaryConditions>\n",
    '  <Geometry type="TensorBSpline3" id="100">\n',
    '    <Basis type="TensorBSplineBasis3">\n',
    '      <Basis type="BSplineBasis" index="0">\n',
    (
        '        <KnotVector degree="2" format="ASCII">0.0 0.0 0.0 1.0 1.0 1.0'
        "</KnotVector>\n"
    ),
    "      </Basis>\n",
    '      <Basis type="BSplineBasis" index="1">\n',
    (
        '        <KnotVector degree="1" format="ASCII">0.0 0.0 1.0 1.0</KnotVector>\n'
    ),
    "      </Basis>\n",
    '      <Basis type="BSplineBasis" index="2">\n',
    (
        '        <KnotVector degree="2" format="ASCII">0.0 0.0 0.0 1.0 1.0 1.0'
        "</KnotVector>\n"
    ),
    "      </Basis>\n",
    "    </Basis>\n",
    '    <coefs geoDim="3" format="ASCII">0.0 0.0 0.0\n',
    "0.5 0.0 0.0\n",
    "1.0 0.0 0.0\n",
    "0.0 1.0 0.0\n",
    "0.5 1.0 0.0\n",
    "1.0 1.0 0.0\n",
    "0.0 0.0 0.5\n",
    "0.5 0.0 0.5\n",
    "1.0 0.0 0.5\n",
    "0.0 1.0 0.5\n",
    "0.5 1.0 0.5\n",
    "1.0 1.0 0.5\n",
    "0.0 0.0 1.0\n",
    "0.5 0.0 1.0\n",
    "1.0 0.0 1.0\n",
    "0.0 1.0 1.0\n",
    "0.5 1.0 1.0\n",
    "1.0 1.0 1.0</coefs>\n",
    "  </Geometry>\n",
    '  <Geometry type="TensorBSpline3" id="101">\n',
    '    <Basis type="TensorBSplineBasis3">\n',
    '      <Basis type="BSplineBasis" index="0">\n',
    (
        '        <KnotVector degree="2" format="ASCII">0.0 0.0 0.0 1.0 1.0 1.0'
        "</KnotVector>\n"
    ),
    "      </Basis>\n",
    '      <Basis type="BSplineBasis" index="1">\n',
    (
        '        <KnotVector degree="1" format="ASCII">0.0 0.0 0.5 1.0 1.0</Kn'
        "otVector>\n"
    ),
    "      </Basis>\n",
    '      <Basis type="BSplineBasis" index="2">\n',
    (
        '        <KnotVector degree="2" format="ASCII">0.0 0.0 0.0 1.0 1.0 1.0'
        "</KnotVector>\n"
    ),
    "      </Basis>\n",
    "    </Basis>\n",
    '    <coefs geoDim="3" format="ASCII">0.0 1.0 0.0\n',
    "0.5 1.0 0.0\n",
    "1.0 1.0 0.0\n",
    "0.0 1.5 0.0\n",
    "0.5 1.5 0.0\n",
    "1.0 1.5 0.0\n",
    "0.0 2.0 0.0\n",
    "0.5 2.0 0.0\n",
    "1.0 2.0 0.0\n",
    "0.0 1.0 0.5\n",
    "0.5 1.0 0.5\n",
    "1.0 1.0 0.5\n",
    "0.0 1.5 0.5\n",
    "0.5 1.5 0.5\n",
    "1.0 1.5 0.5\n",
    "0.0 2.0 0.5\n",
    "0.5 2.0 0.5\n",
    "1.0 2.0 0.5\n",
    "0.0 1.0 1.0\n",
    "0.5 1.0 1.0\n",
    "1.0 1.0 1.0\n",
    "0.0 1.5 1.0\n",
    "0.5 1.5 1.0\n",
    "1.0 1.5 1.0\n",
    "0.0 2.0 1.0\n",
    "0.5 2.0 1.0\n",
    "1.0 2.0 1.0</coefs>\n",
    "  </Geometry>\n",
    '  <Geometry type="TensorNurbs3" id="102">\n',
    '    <Basis type="TensorNurbsBasis3">\n',
    '      <Basis type="TensorBSplineBasis3">\n',
    '        <Basis type="BSplineBasis" index="0">\n',
    (
        '          <KnotVector degree="1" format="ASCII">0.0 0.0 1.0 1.0</Knot'
        "Vector>\n"
    ),
    "        </Basis>\n",
    '        <Basis type="BSplineBasis" index="1">\n',
    (
        '          <KnotVector degree="1" format="ASCII">0.0 0.0 0.5 1.0 1.0</'
        "KnotVector>\n"
    ),
    "        </Basis>\n",
    '        <Basis type="BSplineBasis" index="2">\n',
    (
        '          <KnotVector degree="2" format="ASCII">0.0 0.0 0.0 1.0 1.0 1'
        ".0</KnotVector>\n"
    ),
    "        </Basis>\n",
    "      </Basis>\n",
    '      <weights format="ASCII">1.0\n',
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0</weights>\n",
    "    </Basis>\n",
    '    <coefs geoDim="3" format="ASCII">1.0 1.0 0.0\n',
    "2.0 1.0 0.0\n",
    "1.0 1.5 0.0\n",
    "2.0 1.5 0.0\n",
    "1.0 2.0 0.0\n",
    "2.0 2.0 0.0\n",
    "1.0 1.0 0.5\n",
    "2.0 1.0 0.5\n",
    "1.0 1.5 0.5\n",
    "2.0 1.5 0.5\n",
    "1.0 2.0 0.5\n",
    "2.0 2.0 0.5\n",
    "1.0 1.0 1.0\n",
    "2.0 1.0 1.0\n",
    "1.0 1.5 1.0\n",
    "2.0 1.5 1.0\n",
    "1.0 2.0 1.0\n",
    "2.0 2.0 1.0</coefs>\n",
    "  </Geometry>\n",
    '  <Geometry type="TensorNurbs3" id="103">\n',
    '    <Basis type="TensorNurbsBasis3">\n',
    '      <Basis type="TensorBSplineBasis3">\n',
    '        <Basis type="BSplineBasis" index="0">\n',
    (
        '          <KnotVector degree="1" format="ASCII">0.0 0.0 1.0 1.0</Knot'
        "Vector>\n"
    ),
    "        </Basis>\n",
    '        <Basis type="BSplineBasis" index="1">\n',
    (
        '          <KnotVector degree="1" format="ASCII">0.0 0.0 1.0 1.0</Knot'
        "Vector>\n"
    ),
    "        </Basis>\n",
    '        <Basis type="BSplineBasis" index="2">\n',
    (
        '          <KnotVector degree="2" format="ASCII">0.0 0.0 0.0 1.0 1.0 1'
        ".0</KnotVector>\n"
    ),
    "        </Basis>\n",
    "      </Basis>\n",
    '      <weights format="ASCII">1.0\n',
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0</weights>\n",
    "    </Basis>\n",
    '    <coefs geoDim="3" format="ASCII">1.0 0.0 0.0\n',
    "2.0 0.0 0.0\n",
    "1.0 1.0 0.0\n",
    "2.0 1.0 0.0\n",
    "1.0 0.0 0.5\n",
    "2.0 0.0 0.5\n",
    "1.0 1.0 0.5\n",
    "2.0 1.0 0.5\n",
    "1.0 0.0 1.0\n",
    "2.0 0.0 1.0\n",
    "1.0 1.0 1.0\n",
    "2.0 1.0 1.0</coefs>\n",
    "  </Geometry>\n",
    "</xml>",
]

_gismo_export_ref_binary = [
    (
        '<xml><MultiPatch parDim="3" id="0"><patches type="id_range">100 100</'
        "patches><interfaces /><boundary>100 1\n"
    ),
    "100 2\n",
    "100 3\n",
    "100 4\n",
    "100 5\n",
    (
        '100 6</boundary></MultiPatch><Geometry type="TensorNurbs3" id="100"><'
        'Basis type="TensorNurbsBasis3"><Basis type="TensorBSplineBasis3"><Bas'
        'is type="BSplineBasis" index="0"><KnotVector degree="2" format="B64Fl'
        'oat64">AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA6IsZ1lGpzD/IhucXrjnhP36cNePwj+'
        "c/wUPr0iaA7j8AAAAAAADwPwAAAAAAAPA/AAAAAAAA8D8=</KnotVector></Basis><B"
        'asis type="BSplineBasis" index="1"><KnotVector degree="2" format="B64'
        'Float64">AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAASA5HnqbExj+085Xz6HPOP5jLCkb1'
        "49o/gpPXhXDp5z8AAAAAAADwPwAAAAAAAPA/AAAAAAAA8D8=</KnotVector></Basis>"
        '<Basis type="BSplineBasis" index="2"><KnotVector degree="3" format="B'
        '64Float64">AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADg+v5vc+iiP7TFaR'
        "6eW98/vuFXRKUt6D9h7nNQ533sPwAAAAAAAPA/AAAAAAAA8D8AAAAAAADwPwAAAAAAAPA"
        '/</KnotVector></Basis></Basis><weights format="B64Float64">1Yyh/QFO5T'
        "8yU5fTK1bpP4bnlmKAp+8/NIJStq9DwT/Lqxsx/fbqP4C8zkCYbMI/Ov0iA+Kg3T/YnbD"
        "0stPlP5DQNw2vtt4/MZevQNo86D/4yoqRipHbP+oL2df3vd0/CLQkTid61D/YoUlIb3DT"
        "P7jlz2VovLY/04dbjgAI5z8MVTShYYfSPxiUid9Pndo/gO9V0yzZ0D8T1FxF8bHhP/Uep"
        "mzHuO4/Ur7PFpdM0T8i9ky03vXhP/QcNwd8h98/WeSrBh4Y4z8fkRDQATzhP5AQgqeNrN"
        "I/XbnGoUdL7T/AP1+2TFHNPzIUPd/HfNA/UJXAolG+xz+mkLHPP0LVP0Siu02uZuA/vGq"
        "b2+A22T9KedEoOZ7tP1JbcRC9p90/F5klE8vv5z/1lvr3nmrgP3E+9NW44uM/vYQVBcQ6"
        "6j+eEt+rjO/SP7a/o+zL0NM/NFLFxgEG6z9AXpWJDAPMP4R70jAF9cM/YK1g+J7npD9Wc"
        "Fjb5XHgP/ZgzC7Kzdg/2X/Iy8f56z9EF6pF+9jJP1t4P69KS+8/qF2UUIdHuD9ep3FTZR"
        "rjPxWkvFle0+Y/DE/04/ukyj82uHDpnaDrP3Ivrx8GttE/XFwWfDCI7z+oTk+gQ2LiP6x"
        "+OXAJ3tg/EnscgLal1D8AiRKttH1pPxhQBqW73MI/+BuhWvKR0z+1Ey+t+1TqP9jQMVla"
        "7Lo/9L/6+1Ur3z/I8SqEusjfP3fW3pi9POE/WvMUhm0S0D+wxVyRxOHGP1MknxsVN+o/x"
        "InAbnqVyz/MeknlRLTHPxOMovpdFuQ/fkASFU7L1j8eAp367FzsPyAIVrF9/Jk/JDDfrw"
        "0B1z+QEbIRUCqsP/MVF4VLFeQ/C+Jb2rHq5z9WQj2S/H/cPw+Kcyl8qOk/KQRpTepF6D8"
        "GRznc9vniP47OCWkrv9w/CKGQwasWsz/415XUm7u9P0pRgCxYBuc/LGtoPkKP2z9Ytbou"
        "0izGP5r32fjR4dQ/7wtvuQuX7T/wCX37zyPtP/IBQlpF3eU/MBckvmGuzD8gCim9dqejP"
        "1jZJx8bV7M/mtRYJS5R5D+rUzxNcyziPxXolQMGO+A/Y3X7bMf94z/mK0O6xVzrPy1kF+"
        "I08+o/kJm+PdBHxz9HUmFdNHPhP6ErR+K2E+M/8LuMKpkBsz/gyDIVamC4PyJWVVDtm+0"
        "/HCmtG9Vo7j+KiQP7XOTrPxxpeQDqwuY/eCG4QXqTsT8UU0V+6j7QP0xcqp8bcso/aL7Q"
        "q975tD/qRc4VSEvXP3xkZrqheNM/yWPKjAY26z9t1F2GfPHqP0w6cvCHCNs/2eZyZNLB4"
        "j9KdT//QeHsP55NTBo8h+A/48gRX1yv5z8/6NNzoeLhPyZx6u1ALdk/OQ7GwyjY5T+4il"
        "jLSUjDP4srMstI8+Q/8C0tuF0W2D8ssryZ8BbjP8CttkUsiuA/Jp8TDWML7z+ylqzoz4H"
        "sP6Q2GF6z9uM/dXuNPm5Z4j+6/r03jRfWPyq7QmNbvOU/oPFU1WDe5T+P1XqLIhbtPwj2"
        "QySOmM8/MOXJZPPauD+imEkPPDzjP6TyE98+3t4/iI+HlKskzT+GC74R1R3dPxq2vhMk/"
        "es/sP/iUj2yqD/M+6uL6cbYP559VWolG9c/cPg3lvVy4j/s5sxRAjHvP/haj3TocOk/z8"
        "5fYmCw6D/QdH/CorTdPzq/pEIYdeE/oeu6m72z5D/bbviJcXLgPw5rzX8dPuY//BcOV9f"
        "/2j+kxbhHGLzKP+jrDqO7Mbo/QkEWkvED5T8ZHdUOFHvhP3L8ov0fqd4/uEDGARM3vj9s"
        "RjNTt6zXP8BQ4XhMNYg/9EEesAL27z+A9M57fcx5P+Oozq0MiuI/R3XPhbKD5D9mC+xuB"
        "9rcP5jdIXf/F8E/QGRE561kxD9srctX+oTDPwgc+qABcrA/eWMeB+Bf7T+y2rFaw3vqP8"
        "QysNODC8M/8O2uCmHhoD+LGwlLqCniP6DMJWdHeJ0/WxqGkxoi5z9QLUaAMsCqP57BYW5"
        "p2OU/ctag+lk12D+NH/m1LyfmP8S6BhUO8OA/Jt7lgBYI6z8HzKUBMyLiP/pKKUnJftw/"
        "9B+kJfce5z9YCwHa+MTYP819+OhCg+Q/QWn/eqz47z+UGjSzqkvmP2ApJO3mzM4/06lBd"
        "N7k7z9s+fQ3rxfhPwrO8Pt2gd0/3Y0tB1LL4D+Vyce0kH3uP/e1ywVlB+Y/4r5O2ygO6z"
        "8CNV2fSljfPxDkMSvfK6o/vkirOGJ01T92KSyeVujXP2x9HufMHtM/Zcm51ol66D/0CSz"
        "vgg7KP/Q3obFUwuI/XZwQNeCg4T9AuDEWWhnrPxP4t4wGSuc/XtmtEO676T+HBSkeT0no"
        "P+zJ+21Dkes/A3dP/BMo7T8hNutfhpLhP2gY6eQ0JOA/+NAtIDVB3z8O1AVBpsvmPzhS4"
        "jkA4eA/WHajr5nN0z/enUvNxKHjP4jL65517+w/tZvhKJ716z/dPjl9SPnsP7gKLSgNlr"
        "I/JIskwkhd5D+3wU+wW5/lPyjhIG2t8tI/4poL9A/S0j/o6K9P8Aq/P+5r+kDQDeM/BFw"
        "SOo487z+zWMHaRMHuP+ekjk285Oc/AjfrHrjC0j+08vIXA2nDP1qxOROyrtg/PjlCcw2q"
        "2D9wDir4bvS0P06jtYxys+8/PLtuOUk+0j+SGkeu0C/uP9ahXGHd6O8/gCuIC+n2oT9Ec"
        "Q4dxtbGP1LGKwQ8W9c/dscRmBFm7D/AOi0ikd2jP3TgSo3zqek/zUHGdezy4T+v/CXy2B"
        "foP8BDwSJRstI/5/yCxNmc5j+5S6b0GJLoP2Bw//xmus8/mA79LMdfwz9fBCIr9WnuP8K"
        "88Km72eM/lrQaF+B35T8w7LUAZ8PuP1nGV0yGgOM/HE2rB5AyyT+p7IGbECThPy6hF428"
        "edU/2AzzGOkJwj/FRrX4P+vqP3B6Kbyegs8/LAlPLO01wz9kjmYKz0LMP+CjdSG57a8/U"
        "CI83m/U3T8IXw871EvvP3yfOhDJQs4/wEY/S+eanz/Q+DI2BC/pPxhX9RSQac4/CHof3L"
        "rxuD9YfhkD6d+wP1L50H+aUO8/sFQCl9dk1j88Z5t+uUzEP1JKKS1Sqeo/wvwnalf00j+"
        "cDHiiNTHhP/xm0CoMN9M/HGQ90u/J2z83MO8ji8njP8SynW14hdo/EAG4eYH/1z843uqN"
        "IWjmP/wUYFfXkOM/rWTSkxTw5T/0/Jryu2HnP8xe+6i4btw/kKgfTnKqzj91/9tQNXrtP"
        "2L03DSGftQ/4gcqN85i7D9qE8j+uOHuP1gIIfAKILg/4PQqqjJMpD/WL2ttBzPSPyCitL"
        "MGq70//wTKgu8d6D+LXOsf6q3uP8AutQvFj8s/KGmBGLhmtj9HLDhIa7TvP1N4LrrJW+s"
        "/sWG2cgqN4j/8GMBjp3DCP1ys+L68XuY/58JBQa/D6j8So1hEEfnsP/QKEKsQFOE/BCnS"
        "NqTjzj/tLZkVdLjjP6SLYz8dssA/8PSfVdndoD95cQjTI97qPxooDcILatM/SJr/oQ4s0"
        "D8hvFzfOXfoP5gsAH61rcQ/gLprKFn3qj8bFcSAbIzmP14njVt/j+Y/JQggzc2x4D+0+U"
        "VYCkThP57FXKv7Gtk/jiCu0i3Z3T8gpRb9mSizP9AMT8lw2uY/qI594fe96T8nUIaTp2X"
        "pP0RPZ+AkIdI/O2pTr0XZ7T/hZ9NofJrkP57g5OZb19Y/1ecRQW1X5T8i+zfTjUThPyHs"
        "8CmL9+A/IHPgW0Vp6D+bG1+PXjvrPzYFMWWg+t0/FktELrvx1T+qljQsHDPuPzmuJPsX2"
        "OE/mDVsFzju4j8WpV/VogTrP9gLVX0DgrQ/Zzp0HyTh7D/gOUCGEv6zPyaO4Bisc+k/nf"
        "8NWiCk7j9aUDsv0T3ZP/ECJNd/GO4/LSdNfGsv7D8WNLPWhsHaP4VZfdmWBOQ/QE0huG/"
        "u6T9esPnxxfXRP2DzBS4Djsk/KjUuEvPo3j/Yo9JXUN+3P1Cd28aMb7Y/kKpo8aBttz/U"
        "3EuHoFzhP6AnuSzWR+E/AC6hg7MC7T+4f12OOaixP6e7XtO60+w/hxQKcgwo4T+uBNp6N"
        "hLUP+tnsutNquM/RqMvMUPq4j/Ki7SK8vTnP/iyxyVQ2L4/bdQ3ZhUJ5T9yUtLDWFLtP5"
        'wlqP9SHdQ/h8ES0Oxp6j+D/xEoSPzgPw==</weights></Basis><coefs geoDim="3"'
        ' format="B64Float64">gL3Tcl5Cez8QXdqroIyuP2vsN+mADeI/6j3BghQi6D+XxVD4'
        "+jvhP1ASop6UFNc/QEVvgeCszT8qr7m9r6fXP3DSI/R7m+E/dfMj7iQE7D9ASZ73bHe3P"
        "+yn4K56JNo/8OBwuglJ2z/A+/ZCJ9qjP7OZ6lESnuw/QPYQeuOv0z/AtpZkMhWAP5Cas1"
        "nI8sM/SE/rfzVi5T83tqfKW87oP8QTDEoZJdQ/pgSnuLfz4z8VGHFm44LtP0Xdtjg/Z+s"
        "/3VxGMZos4j+yylDaNkjvP+/yp2aPVO4/BKsyi0fR7z8UuNEhFTHaP/upqq8m1eM/XGvg"
        "d2Jw1j/kfnk7P97mP/RwzPAHsso/NieQlms36z9sUWN3KQLlPxzi1uTyhNo/0IXKavrLy"
        "z+QD71GLxLYP4AHTbMwOso/CMdxY0Lm7z+pJlVBxnzpPyIdeh6CrdM/ETEQvg/95j8H06"
        "dY1c3mP/S2KTzj1sM/nAQwUST17z/accrAqBTjP43Tb+iwB+E/fqsRuljw2j8Saj3FZcX"
        "oPwAL4j3V9NA/AIfN184UYz+Cdhm+Yb7nP754y/g7PtQ/zbsqA3qh6T+NVE5JLCvpPyCR"
        "SN7V1ro/eGC5Npcs2j/kWe+oja3cP4z43xaR48k/ANpTysNWXT8o5I5wxxPVP4JERofHx"
        "eE/pmt7nGK55D+cp5UusqPDP5CnBTKaZ8A/eG7A9Mv1xD964OQIMt3WP0qBSNEeSdw/QG"
        "enSJhFmD+g4xsk6bXVP0SmfHWee84/oDQqqjtQqD/uNefd5/jkPzRIPA0KC8k/1BAbloB"
        "f6z/qzvT5AEzTPwc7oYQCku0/XPPQSAJO6j+zrVWpcZniPws13k+Zzeo/7t6PAeOT1j+h"
        "R4oa3aHrPzxlrrwBM94/aEhRAwwY0j9fhBorhILnP66hRhl9+eI/fHW5NTmpwj8gkMwDR"
        "g+1P9ApohoJ58U/hqxm99yr6T9y7CDvhM/aP4BE0Lw5csU/OA6Eb14y0j9wzdMiR3DNP3"
        "7Bhk6t6Nc/YDc69k1Q0z8EVIDnDQTIP0xOnWSVks4/WE20r6zX6D/I73FS8cPfP4EBniG"
        "lj+M/4E26JKYlvz/YssHS6CzDPzO7XZFFoeI/gPAmBIM/qD/5N6mW4j3nP+sOxQVsteY/"
        "8TMp5Jxk7z8zHr1zO3nmP1jWHtV9meg/8Ecv7QH7qz/QNGUuJUXZPx6ZGlggud4/FaUae"
        "xDP7j/u+/EG4nfYP16Hn65nNuU/StvguWf16T+g/mieUYDeP6ieXmYQcNU/onGVHlsJ4z"
        "/SBuVb9BHtP6dPPCd/2Ok/9K9D0/szzT/aMHU9kIndP4aipwJXDOA/VjDhOaP27z8cyhw"
        "YGKbYP+B13nYFN80/sBQmWtWn2D9D3DAkmR/gP9jWplXLScc/4FfnB13Vlj9TJJ39+i/g"
        "P863wUtl4dA/PNYy+RkPwD825JSulantP2F1+OARuO4/b2LM9KJN6z8yhvVjRXnsP5Jw2"
        "LdlHuM/WivZ2dr70D8ZO/iqia/tPxv1A17hh+g/KHlGUzkm1j9Mwj9IbWfgP9PAx4xzYe"
        "o/3nuq/x3P1T/46rsG62vIP/ijtNTVlbw/KqrjqiGz7D9RyxXaJxzmP8VgacZjH+4/uSP"
        "kVJli6z/U6ONaBi/MPz/ODOLJOOU/fLiKjQk42z/jMXk+pKntP0em2CBvxu8/KLza0ApY"
        "zT8Ia6M5hivjPzqkcW354OI/gKMR2sfytj+g9aZ53u/kP0ArsIbJEK0/VM8jOzXr1T/q1"
        "Gh/tEXWP6DOtUKWBc0/Avg6jx3H4z9QOWJoccbBPyS2cy9A7ck/UKUxrwrrvz/Y0huJ13"
        "jTP2dD/vYwhe4/zloVMVgT1j+4b/BHE/HLP593JcI3ZeA/m0YQw7GL4D902VX6W/zgP1w"
        "tKlTJk8w/TkSS6l2Y4j+aMhyy1TzWP2TE1AG6M+Y/MIMvU0Hswz+UyC5871vEP4hZBJ1D"
        "bL0/f+YCtFGu6z9gVh/4taTbP2BJ5ug0gck/Z/IcOsWm4D/O72Vvp+fmPzK0SmhNudA/y"
        "h8g4ECE5j/zrIiYNkjkP5yunfqkeOA/v/X4laYU4T8iiYecqXfaP5wTD9ENutQ/xK0pVi"
        "Hl3D+AuruxLHTAP/wrSEr+we8/SUDeZOps7z+0vj3yiernP3hahu8C8MI/wbuZaGK97j9"
        "AHqj3+C3CP7WV9yQy2+4/WThad+2u6j/vYP+fcrTiP1KFbNihzd8/gw9AczMq6D9QZex9"
        "JVW7P9AvtLSH29U/kfn/g9YV7T83ZnaxdfHhP8QLX5uEyuo/KJ7FZOrK6j86+yMLEfjfP"
        "/Y5KphvitQ/tkaJYq/o1T8ApCxOTSbUP2T+HT8N08s/tCeyTJXv3j9guEH2cXfGP5oiS5"
        "VKAeQ/yF8nvoUIxT+EH9h/bwXiP5MA3UQie+g/JlUqMd3u2T+UBWfVy4TlPygz4TMGzNY"
        "/uTpEk7mA7j9xEmOiOxfuP0X2qbB4ZuI/MFa+m9r80j+p19xsgYXvPyhnDMvLz9Q/AHe9"
        "gykLvz/2npoN7TbqP1QMj2hLeOE//A0PgVJXwD90KEmubbvfP9F/tvaVkeA/BRn3FubH4"
        "z9ws8MuH6DXP6jeR9Wkm8Q/Y77BAIwG7T8kRC3iCgzXP6BxyLYays0/SBUeVSp2zz+w5/"
        "MpsTy+P8kmqYOOfOg/tLntbeRm6T8wFdlZ6/G8P3DeAS9p7ug/JiSJm8Mo5z9gvlQmfkL"
        "WP70dNoIJmeo/bPdDGS1jxz9iVrlEhuvmPwQYn4f0/ME/dFj1jpApxz/ggjpFHF7eP7BO"
        "KU/8adA/bvEmX2qm2j9g6Kodb1rqP1TIi3s/pOk/vtYZ1v3W3z9QuhUCTsnFP0zVz0HhC"
        "ck/WEuG0mk0vT+nWfCMZIHsPzyzfsdfpcw/ahF6bC2v3T+2nHHx0MHaP4mG59zscOo/kL"
        "4y6Whnxj/oNh3gk2XhP0gjSkaI4e8/EPgPg/bfzj9y81QDqqXSP56vlRsnpdI//PIDpjp"
        "Mzj8O2UPa3h7qP8oxd/cUrtM/usdaJPLV7j8QzZtivRvLP2xhkwUV78Y/M9DHSOpo5z+v"
        "QAXJgw7nPwjwQfmVH88/FEfdz7Ov4D/0ydAjnRXGP3YhGToxpNA/AGJ6fO9U6D+Hrv00E"
        "zHpP2O2N4Tn4eE/MN9tBlyR4D8RZgGYts/sP6AEL/HvjaI/mHl6JZbK2D/ipTDVzPbkPw"
        "Mcyr07NeU/PJLyDoO45z/ESyIsf0bXP8WxUmPS1+g//pWX3R4N5z+Q0yfSNMTRP8SW915"
        "S8s0/UeQS0Uxd6j+QzNATXqvWP5T88vC1Gc4/DKOgOvFp7j909/67pTjNPzxZnCr9t9s/"
        "GNcD4uLV3z+HDnLg4aHgP+Dtskz+N94/2INitaT/wz9EOmFjtmTWP7nmPLx8Z+o/6JuEg"
        "ADByz+GuRGcL0zhP3C7F6EfXr0/oJ4KpNnqnT/YvpRvKdTWPwQ2YDT+oc8/iiEEgT5t2D"
        "+aLfQXlzLZP5AhT4uAw7o/ZQJDrYnb6j/MgMvrHK/DP/C1xSWbEME/YlVKVSM37D/ot0v"
        "Vkmu5PykK/89sEO4/xLvU48UI4j99iweDErvjPzHeUUL25OQ/EMmuXrqYvT9ctab2PEfG"
        "P31oT8h7oew/2KhWKAQ5uD/ExvnJ1XvSPzDijxxb+Ow/KyAvTWmj4T/IJqXEPr/rPyJAR"
        "XFdi+E/oA+6ukEUkT9NaDGIDL/uPwghsxXS2No/qvX/kD2P7D+CGZN+tmjTP8NxEspRd+"
        "w/FAs9Ueg52T+kvo3ECcbGP3hdnxwGvsE/glZuqAcB2j9IOIyvlF7XPzQQuIHZeOQ/dNE"
        "Bz685wT9X6ZG5MPbjP9oA467mr+U/HADj0A2Wwj+vPXNo6anoP8BCFtMhc98/gWoeV490"
        "7j9T8ZMbE7TkP0in8U61AOo/BiOjTYed4z8iEzp7IvjkP/B8ZwyCS9k/DttQI8SE0D+2x"
        "Ygu7PPXP17ZXJ0ERdM/eEzUhTJMtj/RFuqagOvsP5C3jEpsKOQ/cDO+B9gYqj/cgdAPYN"
        "3tPywmd9XaztY/ygnqU+fR1z84X8SREhHnP6DzKhwF/9M/JMjyDIPt6T9comMUbuTdP0R"
        "rTfiIKMk/gfaUlN0y7j866US4ngnfP40E/VLK0+0/1kxjYxu22D+c+47V3ObHPwhfuSrt"
        "Q8U/SCmDZwePvj8lUvKZVlHuPwD+6XtUaO4/KF40gAOJvz+ilGXfhyXhP4zYsNerZ9Y/6"
        "jvHwTmY3D8QeM/Hr8rrPxSu7k8Jccg/CKxvFCRDvj+ojD+9VG/bP3jKCsZPis0/ZiYLBm"
        "jO5j8g3NFRa+aiP8LBn+adG9A/8A5/tZdY3j/17ElI4Y3rP9WFt3Ci+eE/ZP03I/C+2T8"
        "ayUFfraHtP7vDawrrC+M/SLbG6xt16D9/ZZhGZ4jvPyZP8Jm/geg/4FXebaq31z8wfct0"
        "pBnTPxrjolAwm+w/Uk4NyttK4T/Qjxb7rcmpP/4X3dAjReE/gOVAV1Wx4T+sElf+zobpP"
        "+VzbFxjQ+Y/emNPlepy3D+3N4RD3u7pP3dgZVpBkec/NmSrF1Ua0T80QNbIphbpP7Qw6d"
        "wnjuA/aAZkCCbJ0D8wAIO7FnjrP58h6spqiu8/KuAjuJ514D/cNT96VmvnP+Q+8lzgIck"
        "/ws6hNFcz5T+wJ/8fwMPSP6gD4SMr/8I/Sg2ihUMH3D9knr9ikF7QPzgg2nftWL8/hC35"
        "6o1R3T98aIqb963PP0ghECH1gu8/p9PzCdEr5T/AscIwQ6jbP3Oo8+lRru4/QMYSr5CPo"
        "T+a9bFj11PWP+9u/MhebOk/IdFVehcS6D9YQV3mT66zPzJACBxhMtE/apG4K5rn3T/mdE"
        "+uTCrpP5RuEHgrFOk/8Nn0721Fuj8Ubxyf3JnNP6cW6kqLxeQ/UN4/or5Kwj9+YYBPeCP"
        "XP0dyerRUN+U/ABj+KQ5/hj8lZ1PrLpbqP3yyZU3wfcQ/0LVou4D/vj/TYmrlCKTqPzE5"
        "WBrk8us/3ITKPTRo6z+QnIBRVc/OP2izek5vL8A/G8jcw7WJ6j9EbSUlPm7RP2MY6/Ll1"
        "eM/WWUk7Dfg4D8JJLy4v8XrPxgKg5KEk8A/6hmWVIvj3D+dnyijE+znPxS8jRIbqcE/Qg"
        "mrKtKZ0D8/mP6dbQLrP7ppB/Uh0+4/RPBbehoU3j8RQKDIWq3vPzDoGYqJMbk/oJstQf/"
        "xyT+GNBKtK2jcP0Z4Jqg5UO8/gNp48OiZ6j+KomDTjUTdP8wzYgKlP9g/fJn4ScZqyz+E"
        "XgrCbZHJP11ZhLs8z+M/dKgdodZc3j8MeYyi9MrHP2PLlj4RPuU/k2ww9z475T/KslNGK"
        "jjWP7JNxMaVMuo/cAX4Ftkt2D/k5GiwnXfuP2qEd1iGKtw/6wCQGfvC6j/fNAs/IX7lPy"
        "zIsMyPvtg/snHNZmzk4z9ybDQzeoXgP3yQ+IBdquw/TML9+tgy5T/h4IUCWTXnP3HZB7U"
        "Qo+M/YgVCPIj74D8T6TkUyxbqPwBBr78VD8o/Bl3DWUTr2D+QcbbGolvFP+gH42PBoO8/"
        "iJkh+Ps93T/gKvcaHEjDP7zVqkB4etA/iBb7HW9/xD9Y27x2ai+yP2RyqT9h5do/MOJCw"
        "FTDpj+DAxbTQubtPynHwu25s+0/1iWHO0Uy2j9mQIt5mpPYP5e1jtN/xOU/8S3TI0/a6j"
        "9UHKQWVHDNP4Ak9bE8cpA/3glYyXRN7D/xNbWN0rniPwDM/BqoVVw/qkLPZe0R6j/4IiA"
        "q7WKxP+K/s1hNI98/iNuoL3LG1j9FYuFWrjPtPyxuL/BqT9A/DGVx34Yv7j8AlQHN5QPV"
        "P6zLR28mtM4/oIRQQOBSqD8kkpYhGPzSP1jRo4z5Pek//FHD4Evb5D8GMJU/H8PTPzAex"
        "mJrXdc/UNRQbRAG7D9A5PAZYsiEP0DgSxJtf7M/kpifl5Dz0z/MkJ3y8afsP2nrXA0ooO"
        "0/FdIJJysM6T9UtCbHS6HOP5/rmI1CyO8/cntmr1rH3D9MiiJBG+/oP+AjYtef8dE/wns"
        "30CKS1j961PT8fyzjP4BSZTUrR8c/0UcdID9e6j8IsG312sXLP6A8eKwwxMc/YPRucC0o"
        "kD/SHqLh5cLVP0B3n+pFZ+c/iDT3KNFiuD+6Hb1RvA7aP5DSX24fQs8/DNBunWCb5j/Wl"
        "+bU2bXcP7B9ynyLYdc/yM4qb9w37j8ANHBOrsBBPyGTqWYbJuY/JFYJ3R0k6T8aPpscSZ"
        "TsPxvohSUe5+Q/R8w+Tk8d5T96CR3B2pHnPwCtKevPWu8/odl8KmIh6z9YJKet75zCP+s"
        "g6rgjSew/GWTVLcKf7T84U4d/SqbgP46+h/g5Yuo/srJy6hP01j8tqVtNu1vlP+cIGlz0"
        "DeA/TGeIlQj+1D+qeDXakDbYP7fGwxdrXeM/bvUcajO06D/m5gcPRgzoP0jo1fgWi8I/t"
        "GkctxP/1j/aNpxX96nlPwQEM8wDndQ/ONn8mAz06D/M3hApAK/fP1SnQ6GYkMw/lzYSuT"
        "Qf7z+O27MUUE7rP/+Rzv7NjOI/dIO2+XFV7D9kAihVwLTNPxQDLf7CT88/AKZdr2d75z8"
        "7zfA8VoHlP/D8aJ6fyKY/YCYGj8fT6j8efuEMparqP0rog1FG8OU/JNNm4czM6z/Y9tjz"
        "uEjPP+sf0DmepeY/QBlnWNZQyj9mpppZMm3fP/sE6Obfku0/nHunTjDhyz/EvEl7LqbWP"
        "wKn8uxd/N4/GEwp/jwEyT/YefrrQT7dP7VUnEEtcuM/KhT1gaP07z/IChaf9AnQPwryni"
        "FsSNo/cTPMoaNc5T9sS7orxv7UPzfCzajnGu4/RkAYc5/y2j8Q9r8ey3jGP2zsW7BlMtk"
        "/k5qko/0M7j9/kTvlrDrmP6P3VhV3Y+E/0KebQ/I6xT9YX1cpewDsPyiTe1jWPro/jItz"
        "yO7rwz9KJylnCYThP9NfbsY+quU/NFB4v7jD7j8SKD5/SHrdP+COBTkIr6k/GA9iN1a90"
        "D9YAYCuvmvOP/CE+HYpoLg/vLjJe7FP4T9j92Hd09HnP1SueoZkS9o/8qijSDy37z9Szo"
        "onJPLlP2TocLxpoM4/IPCu17sHsD+kxPd/lXnFP2pR972szuk/ZAU8AYK8yz/IPa1S6zj"
        "NP/KKFT4Mzdw/syJoSyHE6j9ukb0MMNjdP/TaFlQ23NI/q8ezL+5c5j/XEqlKp0bqP8J5"
        "T+/Zr+g/w4ej2M3Z7z8grO8pyl+pPyQgb+0E5MA/YKbCCFRFnz8GQN4IJGPkP/KdX3guc"
        "to//H5h8WEIzj9nLKw3kVviP8zUN6+ZQcs/UN+1EybJvD8UtXjhFY/JP2BVc3R+C5M/Hh"
        "cDLres0T/4xBjSTibsP3y271QfLt8/NQslO5Do5z95AXHMgtzvPzyVPljnMMo/BibRLVi"
        "o2j/ajmZ5hSjeP3C04U2zv9Q/eB8ffwDf5z8u05syFvPQPz8cz0i2t+c/+2T5enXa4j/M"
        "ZTpZSR/dP6BkZ7TsI9E/KQ5a/vHl7z/UJpwAsSTUP8fRUxu4juE/1m1GroAv4T/st1Rbd"
        "1nnPxj5byYHB+I/7Q8xIbfd7j8SrMfgbdnVP+ijxdOu7Ow/aiNPUB0y3D8ZHk/5n0LqP/"
        "SNgRKT1+8/gGDKP+GRez+NGhASTpLoPwCFXnPUZME/CM3OxWpf5z9wUAOPkyapP1W/TSJ"
        "mEOM/UgC8kDl+3D/aSSiFXHjiP2Ga1H3al+U/hDYsfSy2zD9aFuvaK0rUPzK/VFKxRuk/"
        "hDCoCpI65z+x/3DIfnHoP7QaQ62Vst0/gtZZhvEW5z9evxTe6DHqPwgwo/P/tOM/wDOrB"
        "4bCzD/mFTcdfDPtP5SCDYA9ft8/D72PRtSA7T/eMwiDNEbbPwx6FOwREOo/OG1RGyxw1T"
        "/gd/CTjkbWP7Y4WihgJt0/YjLen3ZH3j8UDpTYrTnnP5RSPdgfZdA/aMS/WwpJuT/o8kH"
        "YQUO/P6mwQ8FRaew/9guQz0vf1j8EET7b7GjWP2yOW0JmGuA/Z5Lco1+66T/FtZHYIpDh"
        "PywdU0q6++A/xIIsIaIs2z+wooaVN1fvP/iZt4+P0cI/mwaik74r4D+fZF7bUBDkP4of+"
        "zJmndE/NqdWEErx0D82CqCBDBroP8sH/5nWLu0/+/862jgy5D+oPFQ4MqvTPxoDJJhL99"
        "s/NBVMaxa6zz9busnSZ9TlP0wlICFr/8Y/uiC8wtJM6T/A5rgqBozJP8g4SWLgadE/bKs"
        "DY6egxD9EgeuF1mPjP/XWkvUs9u4/jqXmQLpT1z9mcY+Z61vcP0ZWiDw5nNM/0o9nAdvu"
        "0j92ytc1CtTUPwBEda35L9s/+LaxQFFI2z//fcMxLT7qP4YlY0AeHdg/dNvCzKbw5j/d1"
        "QNdyWTpP7T9cm6thek/6FD0NliI3z/A8CoR6XytPwDAJb0XiDc/DBTobmlwyD/uLNPr39"
        "biPzP3kSOWbO0/ALyzg53UgT/eK5AfGoHmP4h3qbBI1M0/EqLIFbkK2z+QmiSIjKjLPwi"
        "ODHzovto/MC0ih/pA2T8MVpnQu5DtP/T1f8rK0No/KZFFR0Ou4T8Tv4IJJUblP5pMVApf"
        "TOQ/fQuId73E7j8I3peptyq7P9jHf9pNutI/ANx/qMpX5D/xVFvaWJ/pPyQ1ropVb+g/i"
        "MRCp8tL1z+aDuJew2zRPwQcwCwPZ+s/OMASZtVNsz94ULf57OXiP8DQDFgJv4Y/UMHm9+"
        "fusD/9+bRh4pbpP0BZw6u9Fuo/dGRMyuuv3j+uwGGicLXeP8wsT7CjgNo/oHgMlWw81D9"
        "OQMVR0+TZP7Z7yxk8k9U/gF/MHor+fD9p4XbQX+/vP3xHECoxt+0/OPlQJpHKzj8QSziv"
        "piPYP8cyQGxLbe4/VtkFd/cJ3j+Afn3LV3zvPzBj7Rwettg/Eo0tww533z8EONA6WxrbP"
        "8tKxwdmmeI/jYD6KLMy4D//OY/3+xbjP7R4AlhEwe4/kGWJ9iE42z+Ofove9G7QP7aWAH"
        "4vj+o/+mfhyaxW3j9odZFNrcrhPzDdZtibHrg/iA7OVoBa0j84jlFNLr3tPwRleCB27d4"
        "/qtcwOjqu3T+A/S+ofaGUP6BAcA5/15A/Bf6KE2Al4T8qcLWCmwHqP5ks9LVir+U/i4G+"
        "bfld4z/WALQRmR3mP0DjkCthybA/dNn2XRwF7T8oFJkOjL25P+CV2ohCFJ4/KZsYM7kc5"
        "z+soFMcZ+DjPzCjA9o+t88/gIri1l3Zez+wliLu8ozWPzRzB9QXzcQ/iBDtnM7OtT/gg6"
        "gfT5/gP5X3yh8svOk/kLg8hnAV7T+P32oOBXLtP21xX3cTOu8/JPOjOnZv2j/qsaKPd77"
        "lP5nUZFosROo/JwIwFB/26D8G7reFjk7tP0A3bK/mEp8/9mgNSrht4T9YQlJRoB2zPz/J"
        "+BYeCOI/esii0wCV5z8O4aLQsEPXP5DiPKIRV9Q/hEznyI3q2z94Hldj3FfMP4xt5CX06"
        "Mg/dsnr2PvE6T+gomy09snOPwQ9GkPbuu8/iqL67CKE4z/qbwhSOWDjP3jYOxh7ku0/wl"
        "AccM6L4j9XhxC12OLiP1jvT6YOfug/AHKPaeQNgz+JuC/6RJnpP9qi7T2AVeM/pOq1bLB"
        "9yj9oiRxrdcLgP9VousUm/ew/SIYL9Ejsvj9LBx0vj4noP+g3CwGr6eI/JJzQGz2v4z94"
        "68le+MG3P8OQzYyvieE/eLfedleiyT+C/7b4PE/QP5qfyRGs1+o/jaijKFkL4T+lv3Lnm"
        "03qP/6O4iZkM9o/QHCBOjU7tT+4zqyDSgPvPxZ7y4x3JtU/KKA0F5posz/ghVRgu4WyPy"
        "AGjo2Jbts/Je70wDy/5T9F+lbrK8frP7aKslN5M98/wGGK6xOZzz/ZHcmfGETkP+ifnia"
        "uoNQ/ENRm0y0f6D8Kq90yhK7YP+7r3Btu5Nw/YPzCYs+Bvz84fXcJvqXVP9yZQljPYc4/"
        "cGk7zow1uT8yYCh28CjiPxD4JIE8frM/4vAh3Cts5j+ArK5JvOW5P47NDEn5bOU/CCKsf"
        "LcRzz9TYjc7d/TuP8I+sEyk1uk/Y2B+tffN5T84mFge9164P4lkyk49ROg/nJcUIPSWyD"
        "8uZfqfkjvlP157AhTfZ+E/43CbDn+k4z/J/KNS9nXkPw6n9aARZNg/yaoMYwKs6T8EWr0"
        "Gq6PbP+LsyS7xYO8/EM6sSHY2yD8UKhRZsMDfP+RSfjHa7+8/rHCbCSHY1D9jfJk8ApPr"
        "P6+ourmiwuY/MMrzwOQRzT98kuclb2TvP95UekRN5+I/zO4/e+q1yT9IM4ebN1XdPyEhe"
        "+orGOA/EbbYCa2c6j9o4zM+LfrWP8Cx2RRik5A/hEC0ScFb3T9WJVN3zSrZPwtJj0Vziu"
        "s/7hhZtMuH5D+QaW257+XpP+zKCTp3X90/4BhLaHZkpj9vWKZENP7lP2S52RktHuQ/gIM"
        "nsL9eyz8JZ45ancLtPxD/hXZbz78/3uIA/gun3z8TGX0zJfTrP5+9JWa2MOA/59mAsggT"
        "5z/ieLUBomfpP8608gHYftw/XtxOsp332j/273ETK3jiP07KgmMotuI/oHBBlg6wkz+s2"
        "7dN1RDMP+BU6Zyu+NQ/IJrGUA9f2D/gWVc2CyjEPwxXPsk5sOk/Ti1LWjwE5j/UnA4EN2"
        "TrP0D/a0Es85I/e1n09BSi4D+6Q2M+xaLgP8gr+hAo9M0/BgFoSdSl6T/McqcWqM7XP1R"
        "pV3/TG8Y/0LGQIWdr4T/O28TvuZvsP7G+noT7mOA/eIQ8Vrdl4j/Mx8LjnzTfP0i5UosH"
        "+Oc/d+6K9PaW4T/gvi39qUbhP+UiC5yK/O8/gHfNA/RC0z8FdbwPQjHiP0sHfrVdwOQ/w"
        "1n1F8S17z/WRgwmS27cP05zfyD5cuM/ixyoV+PH7T9Mua04FGHsP+vgT3CIN+A/mGdd3e"
        "fB4D+7DQf+JzXnPwBeTuGPR2M/y2JHn2G07z/QcbctadntPyLyRvgJR+I/uESGMd6yuz9"
        "YQOAlR4rvP2SgbKDAu+I/qI5O4FlszD8GshsseHfgP/CkY+gkh+o/a1lnzIJ14T8wSYfs"
        "ih3QP3zmGLEcQeM/3tAd9zHr5D+BzNb9oKfuP0cVv7aJK+s/0hvivx3f1z/Ca7Z5HU7ZP"
        "+M23iXXxu0/sNb6gJsq3T/IabGjC3O6P5spwcs02OE/1qJyvHNX6T/NNHXizsTiP8DVO1"
        "awx8w/HdrzlDXv4D9nymN2si3kP/Kklgpyk9Q/3oT073kB6j912q5uMATlP6gjuKbLRrY"
        "/J/KrnMMW7D+Qu1R2ZTHXP9Sh5SP7yd0/DoS6rr2j5D+ZAoXGQl/nP1ndPnh7vec/3Aa9"
        "PZBd7T91z3zt8tLnP2zem0Pv4tg/FJKH7kQE0D/Ec1GGn6vAPwwhjLrdX84/O2LjQ75i6"
        "z+MCao1G8XJPwg2Ky6km9Q/KeT2Ri1N7j/uhRUEO+vrPyB/NS7N2Zs/ilmQD/wo4T+Ohj"
        "vNuFTuP8gKHCH0gcs/eHfw06LT1D/+XXICdjzQP+wnGu7zW+4/Ix22NR/06T9YZq8/ZUe"
        "4P7Z9jdVc3tY/x7HLILb46D/8JPaYBdzjP3GYJfDKxO0/qq1H5g9v4z8EYer4WcnIP4C+"
        "lKf8IXM/qrTyiHYe1D9Sh/H0iWXhP2j5ZZN2Ibs/MK8aplnq0z/Y7+XYvFvqP4+Mt5w5Q"
        "Oo/t+9sVRIn4j/M0vlrJlnAP4kRKaN+P+o/8+3LAitZ6j8U6nCvGz7AP7CsAylJmsE/pn"
        "RYlPNi7j+EAki0YwvWP/ToILOygtc/6ivlLQGx3T9AnQ/SA+/iP+hO+sCNqtA/TEzEpZH"
        "V1z94N6GwfE7OP4q0VjjkudA/cB1+9/GPtz8oF1kfuKDVPxCT/hARAMU/xoHLhUke4z/4"
        "4ZLPAcPYP20RnX61yOg/ACU7IHaK7j+uxkuBDQHgPzaWEVJjw+c/tpfJQXck4D/4pRWxg"
        "Q3VPyJjWIqoN9s/zFytIW1Uzj9+7VdYTQXVPwAF655x8Lc/Z89IbXrp6D/RNPrb4CrnPy"
        "nH1vRe8+w/NKTyjZJf3z++rUcK1oTRPyfC8FzkNuI/bKjEA2Nx4D+x2pEV5wboP1Tem6C"
        "paNM//MmBwFan0D+ofWJQmzrmP/TeiQSYHds/lP9UcZpk6T82F/PJRkbfP/Av11y7nME/"
        "7w/bt4jx7z+qu7sKDzbpP+RUw5jC/cU/fNMZIFzD1z+Y4k2zXDXOP4QXnoY298Q/5Cf5p"
        "Sxx2D+Qb3YZhcraP559bwHrEeQ/KAvJphv91D8AfnomqkhpP0BIDZk5qNo/tO0AlSwS1z"
        "/4f0N9XUq4P9CSk+ZzcNc/LMqybRA2yT+/umPNjeTjP7QTZWywnu0/2BkzkVPXsD/vhlI"
        "3n4LoP96qwtURqdQ/+PXLuzBM3z+g8RoerZG8P94btgB9h+E/ghyz68Iu7z+wpCV06LTF"
        "P5pSRYGaMdE/mMiUym4d0z8wiZl8L5GuP+A8I+OF9NE/7CvTAKr/0z+KdjN5+lbXP6ZN5"
        "+PMZtw/</coefs></Geometry></xml>"
    ),
]


class gismoExportTest(c.unittest.TestCase):
    def test_gismo_export(self):
        """
        Test gismo export routine
        """
        if int(python_version.split(".")[1]) < 8:
            c.splinepy.utils.log.info(
                "gismo export is only tested here from python3.8+. "
                "Skipping test, because current version is: "
                f"{python_version}"
            )
            return True

        ###########
        # 2D Mesh #
        ###########
        # Define some splines
        bez_el0 = c.splinepy.Bezier(
            degrees=[1, 1], control_points=[[0, 0], [1, 0], [0, 1], [1, 1]]
        )
        rbz_el1 = c.splinepy.RationalBezier(
            degrees=[1, 1],
            control_points=[[1, 0], [2, 0], [1, 1], [2, 1]],
            weights=[1, 1, 1, 1],
        )
        bsp_el2 = c.splinepy.BSpline(
            degrees=[1, 1],
            control_points=[[0, 1], [1, 1], [0, 2], [1, 2]],
            knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1]],
        )
        nur_el3 = c.splinepy.NURBS(
            degrees=[1, 1],
            control_points=[[1, 1], [2, 1], [1, 2], [2, 2]],
            weights=[1, 1, 1, 1],
            knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1]],
        )

        # Make it more tricky
        bez_el0.elevate_degrees(0)
        bsp_el2.elevate_degrees(0)

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
        multipatch.boundary_from_function(is_bottom)
        multipatch.boundary_from_function(is_top)

        # Test Output
        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = c.to_tmpf(tmpd)
            c.splinepy.io.gismo.export(
                tmpf,
                multipatch=multipatch,
                indent=False,
                labeled_boundaries=False,
            )

            with open(tmpf) as tmp_read:
                self.assertTrue(
                    c.are_stripped_lines_same(
                        _gismo_export_ref_2d, tmp_read.readlines(), True
                    )
                )

        # for python version > 3.9, test indented version
        if int(python_version.split(".")[1]) >= 9:
            with tempfile.TemporaryDirectory() as tmpd:
                tmpf = c.to_tmpf(tmpd)
                c.splinepy.io.gismo.export(
                    tmpf,
                    multipatch=multipatch,
                    indent=True,
                    labeled_boundaries=False,
                )

                with open(tmpf) as tmp_read:
                    self.assertTrue(
                        c.are_items_same(
                            _gismo_export_ref_2d_indent, tmp_read.readlines()
                        )
                    )

        ########################
        # 2D Mesh - new format #
        ########################
        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = c.to_tmpf(tmpd)
            c.splinepy.io.gismo.export(
                tmpf,
                multipatch=multipatch,
                indent=False,
                labeled_boundaries=True,
            )

            with open(tmpf) as tmp_read:
                self.assertTrue(
                    c.are_stripped_lines_same(
                        _gismo_export_ref_2d_labeled,
                        tmp_read.readlines(),
                        True,
                    )
                )

        # for python version > 3.9, test indented version
        if int(python_version.split(".")[1]) >= 9:
            with tempfile.TemporaryDirectory() as tmpd:
                tmpf = c.to_tmpf(tmpd)
                c.splinepy.io.gismo.export(
                    tmpf,
                    multipatch=multipatch,
                    indent=True,
                    labeled_boundaries=True,
                )

                with open(tmpf) as tmp_read:
                    self.assertTrue(
                        c.are_items_same(
                            _gismo_export_ref_2d_indent_labeled,
                            tmp_read.readlines(),
                        )
                    )

        ###########
        # 3D Mesh #
        ###########

        # Test Also 3D Meshes
        bez_el0 = c.splinepy.Bezier(
            degrees=[1, 1, 1],
            control_points=[
                [0, 0, 0],
                [1, 0, 0],
                [0, 1, 0],
                [1, 1, 0],
                [0, 0, 1],
                [1, 0, 1],
                [0, 1, 1],
                [1, 1, 1],
            ],
        )
        rbz_el1 = c.splinepy.RationalBezier(
            degrees=[1, 1, 1],
            control_points=[
                [1, 0, 0],
                [2, 0, 0],
                [1, 1, 0],
                [2, 1, 0],
                [1, 0, 1],
                [2, 0, 1],
                [1, 1, 1],
                [2, 1, 1],
            ],
            weights=[1] * 8,
        )
        bsp_el2 = c.splinepy.BSpline(
            degrees=[1, 1, 1],
            control_points=[
                [0, 1, 0],
                [1, 1, 0],
                [0, 2, 0],
                [1, 2, 0],
                [0, 1, 1],
                [1, 1, 1],
                [0, 2, 1],
                [1, 2, 1],
            ],
            knot_vectors=[[0, 0, 1, 1]] * 3,
        )
        nur_el3 = c.splinepy.NURBS(
            degrees=[1, 1, 1],
            control_points=[
                [1, 1, 0],
                [2, 1, 0],
                [1, 2, 0],
                [2, 2, 0],
                [1, 1, 1],
                [2, 1, 1],
                [1, 2, 1],
                [2, 2, 1],
            ],
            weights=[1] * 8,
            knot_vectors=[[0, 0, 1, 1]] * 3,
        )

        # Make it more tricky
        bez_el0.elevate_degrees(0)
        bsp_el2.elevate_degrees(0)
        bsp_el2.insert_knots(1, [0.5])
        nur_el3.insert_knots(1, [0.5])
        for s in [bez_el0, bsp_el2, nur_el3, rbz_el1]:
            s.elevate_degrees(2)

        # Define some functions for boundary identification
        def is_bottom(x):
            return x[:, 0] < 0.01

        def is_top(x):
            return x[:, 0] > 1.99

        # Add boundary
        multipatch = c.splinepy.Multipatch(
            splines=[bez_el0, bsp_el2, nur_el3, rbz_el1]
        )
        multipatch.boundary_from_function(is_bottom)
        multipatch.boundary_from_function(is_top)

        # Test output
        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = c.to_tmpf(tmpd)
            c.splinepy.io.gismo.export(
                tmpf,
                multipatch=multipatch,
                indent=False,
                labeled_boundaries=False,
            )
            with open(tmpf) as tmp_read:
                self.assertTrue(
                    c.are_stripped_lines_same(
                        _gismo_export_ref_3d, tmp_read.readlines(), True
                    )
                )

        # for python version > 3.9, test indented version
        if int(python_version.split(".")[1]) >= 9:
            with tempfile.TemporaryDirectory() as tmpd:
                tmpf = c.to_tmpf(tmpd)
                c.splinepy.io.gismo.export(
                    tmpf,
                    multipatch=multipatch,
                    indent=True,
                    labeled_boundaries=False,
                )

                with open(tmpf) as tmp_read:
                    self.assertTrue(
                        c.are_items_same(
                            _gismo_export_ref_3d_indent, tmp_read.readlines()
                        )
                    )

    def test_gismo_import(self):
        """
        Test gismo export routine
        """
        if int(python_version.split(".")[1]) < 8:
            c.splinepy.utils.log.info(
                "gismo export is only tested here from python3.8+. "
                "Skipping test, because current version is: "
                f"{python_version}"
            )
            return True

        # Define some splines
        bsp_el2 = c.splinepy.BSpline(
            degrees=[1, 1],
            control_points=[[0, 0], [1, 0], [0, 1], [1, 1]],
            knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1]],
        )
        nur_el3 = c.splinepy.NURBS(
            degrees=[1, 1],
            control_points=[[1, 0], [2, 0], [1, 1], [2, 1]],
            weights=[1, 1, 1, 1],
            knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1]],
        )

        # Make it more tricky
        bsp_el2.elevate_degrees(0)
        bsp_el2.elevate_degrees(1)
        nur_el3.elevate_degrees(0)
        nur_el3.elevate_degrees(1)
        bsp_el2.insert_knots(1, [0.5])
        nur_el3.insert_knots(1, [0.5])

        # Test Output against input
        multipatch_geometry = c.splinepy.Multipatch([bsp_el2, nur_el3])
        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = c.to_tmpf(tmpd)
            c.splinepy.io.gismo.export(
                tmpf,
                multipatch=multipatch_geometry,
                indent=False,
                labeled_boundaries=False,
            )
            multipatch_geometry_loaded = c.splinepy.io.gismo.load(
                tmpf, load_options=False
            )
            self.assertTrue(
                all(
                    c.are_splines_equal(a, b)
                    for a, b in zip(
                        multipatch_geometry.patches,
                        multipatch_geometry_loaded.patches,
                    )
                )
            )
            self.assertTrue(
                c.np.allclose(
                    multipatch_geometry.interfaces,
                    multipatch_geometry_loaded.interfaces,
                )
            )

            # Now with modified boundaries
            multipatch_geometry.boundaries_from_continuity()
        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = c.to_tmpf(tmpd)
            c.splinepy.io.gismo.export(
                tmpf,
                multipatch=multipatch_geometry,
                indent=False,
                labeled_boundaries=True,
            )
            multipatch_geometry_loaded = c.splinepy.io.gismo.load(
                tmpf, load_options=False
            )
            self.assertTrue(
                all(
                    c.are_splines_equal(a, b)
                    for a, b in zip(
                        multipatch_geometry.patches,
                        multipatch_geometry_loaded.patches,
                    )
                )
            )
            self.assertTrue(
                c.np.allclose(
                    multipatch_geometry.interfaces,
                    multipatch_geometry_loaded.interfaces,
                )
            )

    def test_gismo_import_with_options(self):
        """
        Test gismo export routine
        """
        if int(python_version.split(".")[1]) < 8:
            c.splinepy.utils.log.info(
                "gismo export is only tested here from python3.8+. "
                "Skipping test, because current version is: "
                f"{python_version}"
            )
            return True

        # Define some splines
        bsp_el2 = c.splinepy.BSpline(
            degrees=[1, 1],
            control_points=[[0, 0], [1, 0], [0, 1], [1, 1]],
            knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1]],
        )
        nur_el3 = c.splinepy.NURBS(
            degrees=[1, 1],
            control_points=[[1, 0], [2, 0], [1, 1], [2, 1]],
            weights=[1, 1, 1, 1],
            knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1]],
        )

        # Make it more tricky
        bsp_el2.elevate_degrees(0)
        bsp_el2.elevate_degrees(1)
        nur_el3.elevate_degrees(0)
        nur_el3.elevate_degrees(1)
        bsp_el2.insert_knots(1, [0.5])
        nur_el3.insert_knots(1, [0.5])

        # Test Output against input
        multipatch_geometry = c.splinepy.Multipatch([bsp_el2, nur_el3])

        # Set some options
        gismo_options = [
            {
                "tag": "OptionNumber1",
                "attributes": {"Mambo": "No. 5", "id": "5"},
                "text": "One, two, three, four, five\nEverybody in the car, "
                "so come on, let's ride",
                "children": [
                    {
                        "tag": "AnotherOne",
                        "text": "0 ,0",
                        "attributes": {},
                        "children": [],
                    }
                ],
            }
        ]
        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = c.to_tmpf(tmpd)
            c.splinepy.io.gismo.export(
                tmpf,
                multipatch=multipatch_geometry,
                indent=False,
                labeled_boundaries=False,
                options=gismo_options,
            )
            (
                multipatch_geometry_loaded,
                gismo_options_loaded,
            ) = c.splinepy.io.gismo.load(tmpf, load_options=True)
            self.assertTrue(
                all(
                    c.are_splines_equal(a, b)
                    for a, b in zip(
                        multipatch_geometry.patches,
                        multipatch_geometry_loaded.patches,
                    )
                )
            )
            self.assertTrue(
                c.np.allclose(
                    multipatch_geometry.interfaces,
                    multipatch_geometry_loaded.interfaces,
                )
            )
            self.assertEqual(gismo_options_loaded, gismo_options)

    def test_gismo_io_binary(self):
        """Test the base64 io-routines"""
        # We test this with just one (big, 3D) spline
        nurbs = c.splinepy.NURBS(
            degrees=[1, 1, 1],
            knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1], [0, 0, 1, 1]],
            control_points=c.np.ones((8, 3)),
            weights=c.np.ones((8, 1)),
        )
        nurbs.elevate_degrees([0, 1, 2, 2])
        c.np.random.seed(19284918)
        for i in range(3):
            nurbs.insert_knots(i, c.np.random.random(4))

        # Randomize points
        nurbs.cps = c.np.random.random(nurbs.cps.shape)
        nurbs.weights = c.np.random.random(nurbs.weights.shape)

        # Create a multipatch geometry
        multipatch_geometry = c.splinepy.Multipatch([nurbs])

        # Export
        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = c.to_tmpf(tmpd)
            c.splinepy.io.gismo.export(
                tmpf,
                multipatch=multipatch_geometry,
                indent=False,
                labeled_boundaries=False,
                as_base64=True,
            )
            (
                multipatch_geometry_loaded,
                gismo_options_loaded,
            ) = c.splinepy.io.gismo.load(tmpf, load_options=True)
            with open(tmpf) as tmp_read:
                self.assertTrue(
                    c.are_items_same(
                        _gismo_export_ref_binary, tmp_read.readlines()
                    )
                )
            self.assertTrue(
                all(
                    c.are_splines_equal(a, b)
                    for a, b in zip(
                        multipatch_geometry.patches,
                        multipatch_geometry_loaded.patches,
                    )
                )
            )
            self.assertTrue(
                c.np.allclose(
                    multipatch_geometry.interfaces,
                    multipatch_geometry_loaded.interfaces,
                )
            )


if __name__ == "__main__":
    c.unittest.main()
