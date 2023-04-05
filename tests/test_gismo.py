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
        '0"><KnotVector degree="2">0.0 0.0 0.0 1.0 1.0 1.0</KnotVector></Basis'
        '><Basis type="BSplineBasis" index="1"><KnotVector degree="1">0.0 0.0 '
        '1.0 1.0</KnotVector></Basis></Basis><coefs geoDim="2">0.0 0.0\n'
    ),
    "0.5 0.0\n",
    "1.0 0.0\n",
    "0.0 1.0\n",
    "0.5 1.0\n",
    (
        '1.0 1.0</coefs></Geometry><Geometry type="TensorNurbs2" id="101"><Bas'
        'is type="TensorNurbsBasis2"><Basis type="TensorBSplineBasis2"><Basis '
        'type="BSplineBasis" index="0"><KnotVector degree="1">0.0 0.0 1.0 1.0<'
        '/KnotVector></Basis><Basis type="BSplineBasis" index="1"><KnotVector '
        'degree="1">0.0 0.0 1.0 1.0</KnotVector></Basis></Basis><weights>1.0\n'
    ),
    "1.0\n",
    "1.0\n",
    '1.0</weights></Basis><coefs geoDim="2">1.0 0.0\n',
    "2.0 0.0\n",
    "1.0 1.0\n",
    (
        '2.0 1.0</coefs></Geometry><Geometry type="TensorBSpline2" id="102"><B'
        'asis type="TensorBSplineBasis2"><Basis type="BSplineBasis" index="0">'
        '<KnotVector degree="2">0.0 0.0 0.0 1.0 1.0 1.0</KnotVector></Basis><B'
        'asis type="BSplineBasis" index="1"><KnotVector degree="1">0.0 0.0 0.5'
        ' 1.0 1.0</KnotVector></Basis></Basis><coefs geoDim="2">0.0 1.0\n'
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
        'type="BSplineBasis" index="0"><KnotVector degree="1">0.0 0.0 1.0 1.0<'
        '/KnotVector></Basis><Basis type="BSplineBasis" index="1"><KnotVector '
        'degree="1">0.0 0.0 0.5 1.0 1.0</KnotVector></Basis></Basis><weights>1'
        ".0\n"
    ),
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    '1.0</weights></Basis><coefs geoDim="2">1.0 1.0\n',
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
    '        <KnotVector degree="2">0.0 0.0 0.0 1.0 1.0 1.0</KnotVector>\n',
    "      </Basis>\n",
    '      <Basis type="BSplineBasis" index="1">\n',
    '        <KnotVector degree="1">0.0 0.0 1.0 1.0</KnotVector>\n',
    "      </Basis>\n",
    "    </Basis>\n",
    '    <coefs geoDim="2">0.0 0.0\n',
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
    '          <KnotVector degree="1">0.0 0.0 1.0 1.0</KnotVector>\n',
    "        </Basis>\n",
    '        <Basis type="BSplineBasis" index="1">\n',
    '          <KnotVector degree="1">0.0 0.0 1.0 1.0</KnotVector>\n',
    "        </Basis>\n",
    "      </Basis>\n",
    "      <weights>1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0</weights>\n",
    "    </Basis>\n",
    '    <coefs geoDim="2">1.0 0.0\n',
    "2.0 0.0\n",
    "1.0 1.0\n",
    "2.0 1.0</coefs>\n",
    "  </Geometry>\n",
    '  <Geometry type="TensorBSpline2" id="102">\n',
    '    <Basis type="TensorBSplineBasis2">\n',
    '      <Basis type="BSplineBasis" index="0">\n',
    '        <KnotVector degree="2">0.0 0.0 0.0 1.0 1.0 1.0</KnotVector>\n',
    "      </Basis>\n",
    '      <Basis type="BSplineBasis" index="1">\n',
    '        <KnotVector degree="1">0.0 0.0 0.5 1.0 1.0</KnotVector>\n',
    "      </Basis>\n",
    "    </Basis>\n",
    '    <coefs geoDim="2">0.0 1.0\n',
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
    '          <KnotVector degree="1">0.0 0.0 1.0 1.0</KnotVector>\n',
    "        </Basis>\n",
    '        <Basis type="BSplineBasis" index="1">\n',
    '          <KnotVector degree="1">0.0 0.0 0.5 1.0 1.0</KnotVector>\n',
    "        </Basis>\n",
    "      </Basis>\n",
    "      <weights>1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0</weights>\n",
    "    </Basis>\n",
    '    <coefs geoDim="2">1.0 1.0\n',
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
        '0"><KnotVector degree="2">0.0 0.0 0.0 1.0 1.0 1.0</KnotVector></Basis'
        '><Basis type="BSplineBasis" index="1"><KnotVector degree="1">0.0 0.0 '
        '1.0 1.0</KnotVector></Basis></Basis><coefs geoDim="2">0.0 0.0\n'
    ),
    "0.5 0.0\n",
    "1.0 0.0\n",
    "0.0 1.0\n",
    "0.5 1.0\n",
    (
        '1.0 1.0</coefs></Geometry><Geometry type="TensorNurbs2" id="101"><Bas'
        'is type="TensorNurbsBasis2"><Basis type="TensorBSplineBasis2"><Basis '
        'type="BSplineBasis" index="0"><KnotVector degree="1">0.0 0.0 1.0 1.0<'
        '/KnotVector></Basis><Basis type="BSplineBasis" index="1"><KnotVector '
        'degree="1">0.0 0.0 1.0 1.0</KnotVector></Basis></Basis><weights>1.0\n'
    ),
    "1.0\n",
    "1.0\n",
    '1.0</weights></Basis><coefs geoDim="2">1.0 0.0\n',
    "2.0 0.0\n",
    "1.0 1.0\n",
    (
        '2.0 1.0</coefs></Geometry><Geometry type="TensorBSpline2" id="102"><B'
        'asis type="TensorBSplineBasis2"><Basis type="BSplineBasis" index="0">'
        '<KnotVector degree="2">0.0 0.0 0.0 1.0 1.0 1.0</KnotVector></Basis><B'
        'asis type="BSplineBasis" index="1"><KnotVector degree="1">0.0 0.0 0.5'
        ' 1.0 1.0</KnotVector></Basis></Basis><coefs geoDim="2">0.0 1.0\n'
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
        'type="BSplineBasis" index="0"><KnotVector degree="1">0.0 0.0 1.0 1.0<'
        '/KnotVector></Basis><Basis type="BSplineBasis" index="1"><KnotVector '
        'degree="1">0.0 0.0 0.5 1.0 1.0</KnotVector></Basis></Basis><weights>1'
        ".0\n"
    ),
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    '1.0</weights></Basis><coefs geoDim="2">1.0 1.0\n',
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
    '        <KnotVector degree="2">0.0 0.0 0.0 1.0 1.0 1.0</KnotVector>\n',
    "      </Basis>\n",
    '      <Basis type="BSplineBasis" index="1">\n',
    '        <KnotVector degree="1">0.0 0.0 1.0 1.0</KnotVector>\n',
    "      </Basis>\n",
    "    </Basis>\n",
    '    <coefs geoDim="2">0.0 0.0\n',
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
    '          <KnotVector degree="1">0.0 0.0 1.0 1.0</KnotVector>\n',
    "        </Basis>\n",
    '        <Basis type="BSplineBasis" index="1">\n',
    '          <KnotVector degree="1">0.0 0.0 1.0 1.0</KnotVector>\n',
    "        </Basis>\n",
    "      </Basis>\n",
    "      <weights>1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0</weights>\n",
    "    </Basis>\n",
    '    <coefs geoDim="2">1.0 0.0\n',
    "2.0 0.0\n",
    "1.0 1.0\n",
    "2.0 1.0</coefs>\n",
    "  </Geometry>\n",
    '  <Geometry type="TensorBSpline2" id="102">\n',
    '    <Basis type="TensorBSplineBasis2">\n',
    '      <Basis type="BSplineBasis" index="0">\n',
    '        <KnotVector degree="2">0.0 0.0 0.0 1.0 1.0 1.0</KnotVector>\n',
    "      </Basis>\n",
    '      <Basis type="BSplineBasis" index="1">\n',
    '        <KnotVector degree="1">0.0 0.0 0.5 1.0 1.0</KnotVector>\n',
    "      </Basis>\n",
    "    </Basis>\n",
    '    <coefs geoDim="2">0.0 1.0\n',
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
    '          <KnotVector degree="1">0.0 0.0 1.0 1.0</KnotVector>\n',
    "        </Basis>\n",
    '        <Basis type="BSplineBasis" index="1">\n',
    '          <KnotVector degree="1">0.0 0.0 0.5 1.0 1.0</KnotVector>\n',
    "        </Basis>\n",
    "      </Basis>\n",
    "      <weights>1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0</weights>\n",
    "    </Basis>\n",
    '    <coefs geoDim="2">1.0 1.0\n',
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
        "patches><interfaces>100 4 101 3 0 1 2 1 1 1\n"
    ),
    "101 2 102 1 0 1 2 1 1 1\n",
    "100 2 103 1 0 1 2 1 1 1\n",
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
        '0"><KnotVector degree="2">0.0 0.0 0.0 1.0 1.0 1.0</KnotVector></Basis'
        '><Basis type="BSplineBasis" index="1"><KnotVector degree="1">0.0 0.0 '
        '1.0 1.0</KnotVector></Basis><Basis type="BSplineBasis" index="2"><Kno'
        'tVector degree="2">0.0 0.0 0.0 1.0 1.0 1.0</KnotVector></Basis></Basi'
        's><coefs geoDim="3">0.0 0.0 0.0\n'
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
        '"0"><KnotVector degree="2">0.0 0.0 0.0 1.0 1.0 1.0</KnotVector></Basi'
        's><Basis type="BSplineBasis" index="1"><KnotVector degree="1">0.0 0.0'
        ' 0.5 1.0 1.0</KnotVector></Basis><Basis type="BSplineBasis" index="2"'
        '><KnotVector degree="2">0.0 0.0 0.0 1.0 1.0 1.0</KnotVector></Basis><'
        '/Basis><coefs geoDim="3">0.0 1.0 0.0\n'
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
        'sis type="BSplineBasis" index="0"><KnotVector degree="1">0.0 0.0 1.0 '
        '1.0</KnotVector></Basis><Basis type="BSplineBasis" index="1"><KnotVec'
        'tor degree="1">0.0 0.0 0.5 1.0 1.0</KnotVector></Basis><Basis type="B'
        'SplineBasis" index="2"><KnotVector degree="2">0.0 0.0 0.0 1.0 1.0 1.0'
        "</KnotVector></Basis></Basis><weights>1.0\n"
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
    '1.0</weights></Basis><coefs geoDim="3">1.0 1.0 0.0\n',
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
        'sis type="BSplineBasis" index="0"><KnotVector degree="1">0.0 0.0 1.0 '
        '1.0</KnotVector></Basis><Basis type="BSplineBasis" index="1"><KnotVec'
        'tor degree="1">0.0 0.0 1.0 1.0</KnotVector></Basis><Basis type="BSpli'
        'neBasis" index="2"><KnotVector degree="2">0.0 0.0 0.0 1.0 1.0 1.0</Kn'
        "otVector></Basis></Basis><weights>1.0\n"
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
    '1.0</weights></Basis><coefs geoDim="3">1.0 0.0 0.0\n',
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
    "    <interfaces>100 4 101 3 0 1 2 1 1 1\n",
    "101 2 102 1 0 1 2 1 1 1\n",
    "100 2 103 1 0 1 2 1 1 1\n",
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
    '        <KnotVector degree="2">0.0 0.0 0.0 1.0 1.0 1.0</KnotVector>\n',
    "      </Basis>\n",
    '      <Basis type="BSplineBasis" index="1">\n',
    '        <KnotVector degree="1">0.0 0.0 1.0 1.0</KnotVector>\n',
    "      </Basis>\n",
    '      <Basis type="BSplineBasis" index="2">\n',
    '        <KnotVector degree="2">0.0 0.0 0.0 1.0 1.0 1.0</KnotVector>\n',
    "      </Basis>\n",
    "    </Basis>\n",
    '    <coefs geoDim="3">0.0 0.0 0.0\n',
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
    '        <KnotVector degree="2">0.0 0.0 0.0 1.0 1.0 1.0</KnotVector>\n',
    "      </Basis>\n",
    '      <Basis type="BSplineBasis" index="1">\n',
    '        <KnotVector degree="1">0.0 0.0 0.5 1.0 1.0</KnotVector>\n',
    "      </Basis>\n",
    '      <Basis type="BSplineBasis" index="2">\n',
    '        <KnotVector degree="2">0.0 0.0 0.0 1.0 1.0 1.0</KnotVector>\n',
    "      </Basis>\n",
    "    </Basis>\n",
    '    <coefs geoDim="3">0.0 1.0 0.0\n',
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
    '          <KnotVector degree="1">0.0 0.0 1.0 1.0</KnotVector>\n',
    "        </Basis>\n",
    '        <Basis type="BSplineBasis" index="1">\n',
    '          <KnotVector degree="1">0.0 0.0 0.5 1.0 1.0</KnotVector>\n',
    "        </Basis>\n",
    '        <Basis type="BSplineBasis" index="2">\n',
    '          <KnotVector degree="2">0.0 0.0 0.0 1.0 1.0 1.0</KnotVector>\n',
    "        </Basis>\n",
    "      </Basis>\n",
    "      <weights>1.0\n",
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
    '    <coefs geoDim="3">1.0 1.0 0.0\n',
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
    '          <KnotVector degree="1">0.0 0.0 1.0 1.0</KnotVector>\n',
    "        </Basis>\n",
    '        <Basis type="BSplineBasis" index="1">\n',
    '          <KnotVector degree="1">0.0 0.0 1.0 1.0</KnotVector>\n',
    "        </Basis>\n",
    '        <Basis type="BSplineBasis" index="2">\n',
    '          <KnotVector degree="2">0.0 0.0 0.0 1.0 1.0 1.0</KnotVector>\n',
    "        </Basis>\n",
    "      </Basis>\n",
    "      <weights>1.0\n",
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
    '    <coefs geoDim="3">1.0 0.0 0.0\n',
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
        bsp_el2.elevate_degree(0)
        bsp_el2.elevate_degree(1)
        nur_el3.elevate_degree(0)
        nur_el3.elevate_degree(1)
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
                    [
                        c.are_splines_equal(a, b)
                        for a, b in zip(
                            multipatch_geometry.splines,
                            multipatch_geometry_loaded.splines,
                        )
                    ]
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
                    [
                        c.are_splines_equal(a, b)
                        for a, b in zip(
                            multipatch_geometry.splines,
                            multipatch_geometry_loaded.splines,
                        )
                    ]
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
        bsp_el2.elevate_degree(0)
        bsp_el2.elevate_degree(1)
        nur_el3.elevate_degree(0)
        nur_el3.elevate_degree(1)
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
                    [
                        c.are_splines_equal(a, b)
                        for a, b in zip(
                            multipatch_geometry.splines,
                            multipatch_geometry_loaded.splines,
                        )
                    ]
                )
            )
            self.assertTrue(
                c.np.allclose(
                    multipatch_geometry.interfaces,
                    multipatch_geometry_loaded.interfaces,
                )
            )
            self.assertEqual(gismo_options_loaded, gismo_options)


if __name__ == "__main__":
    c.unittest.main()
