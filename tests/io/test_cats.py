import tempfile
from sys import version as python_version

import splinepy

cats_no_indent_export = [
    '<SplineList SplineType="1" NumberOfSplines="4"><!--generated by splinepy'
    ' https://github.com/tataratat/splinepy--><SplineEntry splDim="'
    '2" spaceDim="2" numOfCntrlPntVars="2" numCntrlPnts="6" numOfEleVars="0" c'
    'losed="0 0"><cntrlPntVarNames>x y</cntrlPntVarNames><cntrlPntVars>0.0 0.0'
    " 0.5 0.0 1.0 0.0 0.0 1.0 0.5 1.0 1.0 1.0</cntrlPntVars><deg>2 1</deg><knt"
    "Vecs><kntVec>0.0 0.0 0.0 1.0 1.0 1.0</kntVec><kntVec>0.0 0.0 1.0 1.0</"
    'kntVec></kntVecs></SplineEntry><SplineEntry splDim="2" spaceDim="2" numO'
    'fCntrlPntVars="2" numCntrlPnts="4" numOfEleVars="0" closed="0 0"><cntrlPn'
    "tVarNames>x y</cntrlPntVarNames><cntrlPntVars>1.0 0.0 2.0 0.0 1.0 1.0 2.0"
    " 1.0</cntrlPntVars><deg>1 1</deg><kntVecs><kntVec>0.0 0.0 1.0 1.0</knotV"
    "ec><kntVec>0.0 0.0 1.0 1.0</kntVec></kntVecs><wght>1.0 1.0 1.0 1.0</wgh"
    't></SplineEntry><SplineEntry splDim="2" spaceDim="2" numOfCntrlPntVars="2'
    '" numCntrlPnts="9" numOfEleVars="0" closed="0 0"><cntrlPntVarNames>x y</c'
    "ntrlPntVarNames><cntrlPntVars>0.0 1.0 0.5 1.0 1.0 1.0 0.0 1.5 0.5 1.5 1.0"
    " 1.5 0.0 2.0 0.5 2.0 1.0 2.0</cntrlPntVars><deg>2 1</deg><kntVecs><knotVe"
    "c>0.0 0.0 0.0 1.0 1.0 1.0</kntVec><kntVec>0.0 0.0 0.5 1.0 1.0</kntVec>"
    '</kntVecs></SplineEntry><SplineEntry splDim="2" spaceDim="2" numOfCntrlPn'
    'tVars="2" numCntrlPnts="6" numOfEleVars="0" closed="0 0"><cntrlPntVarName'
    "s>x y</cntrlPntVarNames><cntrlPntVars>1.0 1.0 2.0 1.0 1.0 1.5 2.0 1.5 1.0"
    " 2.0 2.0 2.0</cntrlPntVars><deg>1 1</deg><kntVecs><kntVec>0.0 0.0 1.0 1."
    "0</kntVec><kntVec>0.0 0.0 0.5 1.0 1.0</kntVec></kntVecs><wght>1.0 1.0 "
    "1.0 1.0 1.0 1.0</wght></SplineEntry></SplineList>"
]

cats_indent_export = [
    '<SplineList SplineType="1" NumberOfSplines="4">\n',
    "  <!--generated by splinepy https://github.com/tataratat/splinepy-->\n",
    '  <SplineEntry splDim="2" spaceDim="2" numOfCntrlPntVars="2" '
    'numCntrlPnts="6" numOfEleVars="0" closed="0 0">\n',
    "    <cntrlPntVarNames>x y</cntrlPntVarNames>\n",
    "    <cntrlPntVars>0.0 0.0\n",
    "0.5 0.0\n",
    "1.0 0.0\n",
    "0.0 1.0\n",
    "0.5 1.0\n",
    "1.0 1.0</cntrlPntVars>\n",
    "    <deg>2\n",
    "1</deg>\n",
    "    <kntVecs>\n",
    "      <kntVec>0.0\n",
    "0.0\n",
    "0.0\n",
    "1.0\n",
    "1.0\n",
    "1.0</kntVec>\n",
    "      <kntVec>0.0\n",
    "0.0\n",
    "1.0\n",
    "1.0</kntVec>\n",
    "    </kntVecs>\n",
    "  </SplineEntry>\n",
    '  <SplineEntry splDim="2" spaceDim="2" numOfCntrlPntVars="2" '
    'numCntrlPnts="4" numOfEleVars="0" closed="0 0">\n',
    "    <cntrlPntVarNames>x y</cntrlPntVarNames>\n",
    "    <cntrlPntVars>1.0 0.0\n",
    "2.0 0.0\n",
    "1.0 1.0\n",
    "2.0 1.0</cntrlPntVars>\n",
    "    <deg>1\n",
    "1</deg>\n",
    "    <kntVecs>\n",
    "      <kntVec>0.0\n",
    "0.0\n",
    "1.0\n",
    "1.0</kntVec>\n",
    "      <kntVec>0.0\n",
    "0.0\n",
    "1.0\n",
    "1.0</kntVec>\n",
    "    </kntVecs>\n",
    "    <wght>1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0</wght>\n",
    "  </SplineEntry>\n",
    '  <SplineEntry splDim="2" spaceDim="2" numOfCntrlPntVars="2" '
    'numCntrlPnts="9" numOfEleVars="0" closed="0 0">\n',
    "    <cntrlPntVarNames>x y</cntrlPntVarNames>\n",
    "    <cntrlPntVars>0.0 1.0\n",
    "0.5 1.0\n",
    "1.0 1.0\n",
    "0.0 1.5\n",
    "0.5 1.5\n",
    "1.0 1.5\n",
    "0.0 2.0\n",
    "0.5 2.0\n",
    "1.0 2.0</cntrlPntVars>\n",
    "    <deg>2\n",
    "1</deg>\n",
    "    <kntVecs>\n",
    "      <kntVec>0.0\n",
    "0.0\n",
    "0.0\n",
    "1.0\n",
    "1.0\n",
    "1.0</kntVec>\n",
    "      <kntVec>0.0\n",
    "0.0\n",
    "0.5\n",
    "1.0\n",
    "1.0</kntVec>\n",
    "    </kntVecs>\n",
    "  </SplineEntry>\n",
    '  <SplineEntry splDim="2" spaceDim="2" numOfCntrlPntVars="2" '
    'numCntrlPnts="6" numOfEleVars="0" closed="0 0">\n',
    "    <cntrlPntVarNames>x y</cntrlPntVarNames>\n",
    "    <cntrlPntVars>1.0 1.0\n",
    "2.0 1.0\n",
    "1.0 1.5\n",
    "2.0 1.5\n",
    "1.0 2.0\n",
    "2.0 2.0</cntrlPntVars>\n",
    "    <deg>1\n",
    "1</deg>\n",
    "    <kntVecs>\n",
    "      <kntVec>0.0\n",
    "0.0\n",
    "1.0\n",
    "1.0</kntVec>\n",
    "      <kntVec>0.0\n",
    "0.0\n",
    "0.5\n",
    "1.0\n",
    "1.0</kntVec>\n",
    "    </kntVecs>\n",
    "    <wght>1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0\n",
    "1.0</wght>\n",
    "  </SplineEntry>\n",
    "</SplineList>",
]


def test_cats_export(to_tmpf, are_items_same, are_stripped_lines_same):
    """Test cats export routine"""
    ###########
    # 2D Mesh #
    ###########
    # Define some splines
    bez_el0 = splinepy.Bezier(
        degrees=[1, 1], control_points=[[0, 0], [1, 0], [0, 1], [1, 1]]
    )
    rbz_el1 = splinepy.RationalBezier(
        degrees=[1, 1],
        control_points=[[1, 0], [2, 0], [1, 1], [2, 1]],
        weights=[1, 1, 1, 1],
    )
    bsp_el2 = splinepy.BSpline(
        degrees=[1, 1],
        control_points=[[0, 1], [1, 1], [0, 2], [1, 2]],
        knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1]],
    )
    nur_el3 = splinepy.NURBS(
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
    multipatch = splinepy.Multipatch(
        splines=[bez_el0, rbz_el1, bsp_el2, nur_el3]
    )

    # Test Output
    with tempfile.TemporaryDirectory() as tmpd:
        tmpf = to_tmpf(tmpd)
        splinepy.io.cats.export(
            tmpf, spline_list=multipatch, indent=False, make_rational=False
        )

        with open(tmpf) as tmp_read:
            assert are_stripped_lines_same(
                cats_no_indent_export, tmp_read.readlines(), True
            )

    # for python version > 3.9, test indented version
    if int(python_version.split(".")[1]) >= 9:
        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = to_tmpf(tmpd)
            splinepy.io.cats.export(
                tmpf, spline_list=multipatch, indent=True, make_rational=False
            )

            with open(tmpf) as tmp_read:
                assert are_items_same(cats_indent_export, tmp_read.readlines())


def test_cats_import(to_tmpf, are_splines_equal):
    """Test cats import routine"""
    # Define some splines
    bsp_el2 = splinepy.BSpline(
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
    nur_el3 = splinepy.NURBS(
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
    bsp_el2.elevate_degrees(0)
    bsp_el2.elevate_degrees(1)
    nur_el3.elevate_degrees(0)
    nur_el3.elevate_degrees(1)
    bsp_el2.insert_knots(1, [0.5])
    nur_el3.insert_knots(1, [0.5])

    # Test Output against input
    multipatch_geometry = [bsp_el2, nur_el3]
    with tempfile.TemporaryDirectory() as tmpd:
        tmpf = to_tmpf(tmpd)
        splinepy.io.cats.export(
            tmpf,
            spline_list=multipatch_geometry,
            indent=False,
            make_rational=False,
        )
        multipatch_geometry_loaded = splinepy.io.cats.load(tmpf)
        assert all(
            are_splines_equal(a, b)
            for a, b in zip(
                multipatch_geometry,
                multipatch_geometry_loaded,
            )
        )

    # Make all splines rational
    with tempfile.TemporaryDirectory() as tmpd:
        tmpf = to_tmpf(tmpd)
        splinepy.io.cats.export(
            tmpf,
            spline_list=multipatch_geometry,
            indent=False,
        )
        multipatch_geometry_loaded = splinepy.io.cats.load(tmpf)
        assert all(
            are_splines_equal(a.nurbs, b)
            for a, b in zip(
                multipatch_geometry,
                multipatch_geometry_loaded,
            )
        )
