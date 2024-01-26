import tempfile
from sys import version as python_version

try:
    from ..examples import common as c
except BaseException:
    import common as c

cats_no_indent_export = [
    '<SplineList SplineType="1" NumberOfSplines="4"><SplineEntry splDim="'
    '2" spaceDim="2" numOfCntrlPntVars="2" numCntrlPnts="6" numOfEleVars="0" c'
    'losed="0"><cntrlPntVarNames>x y</cntrlPntVarNames><cntrlPntVars>0.0 0.0'
    " 0.5 0.0 1.0 0.0 0.0 1.0 0.5 1.0 1.0 1.0</cntrlPntVars><deg>2 1</deg><knt"
    "Vecs><knotVec>0.0 0.0 0.0 1.0 1.0 1.0</knotVec><knotVec>0.0 0.0 1.0 1.0</"
    'knotVec></kntVecs></SplineEntry><SplineEntry splDim="2" spaceDim="2" numO'
    'fCntrlPntVars="2" numCntrlPnts="4" numOfEleVars="0" closed="0"><cntrlPntV'
    "arNames>x y</cntrlPntVarNames><cntrlPntVars>1.0 0.0 2.0 0.0 1.0 1.0 2.0"
    " 1.0</cntrlPntVars><deg>1 1</deg><kntVecs><knotVec>0.0 0.0 1.0 1.0</knotV"
    "ec><knotVec>0.0 0.0 1.0 1.0</knotVec></kntVecs><wght>1.0 1.0 1.0 1.0</wgh"
    't></SplineEntry><SplineEntry splDim="2" spaceDim="2" numOfCntrlPntVars="2'
    '" numCntrlPnts="9" numOfEleVars="0" closed="0"><cntrlPntVarNames>x y</c'
    "ntrlPntVarNames><cntrlPntVars>0.0 1.0 0.5 1.0 1.0 1.0 0.0 1.5 0.5 1.5 1.0"
    " 1.5 0.0 2.0 0.5 2.0 1.0 2.0</cntrlPntVars><deg>2 1</deg><kntVecs><knotVe"
    "c>0.0 0.0 0.0 1.0 1.0 1.0</knotVec><knotVec>0.0 0.0 0.5 1.0 1.0</knotVec>"
    '</kntVecs></SplineEntry><SplineEntry splDim="2" spaceDim="2" numOfCntrlPn'
    'tVars="2" numCntrlPnts="6" numOfEleVars="0" closed="0"><cntrlPntVarNames>'
    "x y</cntrlPntVarNames><cntrlPntVars>1.0 1.0 2.0 1.0 1.0 1.5 2.0 1.5 1.0"
    " 2.0 2.0 2.0</cntrlPntVars><deg>1 1</deg><kntVecs><knotVec>0.0 0.0 1.0 1."
    "0</knotVec><knotVec>0.0 0.0 0.5 1.0 1.0</knotVec></kntVecs><wght>1.0 1.0 "
    "1.0 1.0 1.0 1.0</wght></SplineEntry></SplineList>"
]

cats_indent_export = [
    '<SplineList SplineType="1" NumberOfSplines="4">\n',
    (
        '  <SplineEntry splDim="2" spaceDim="2" numOfCntrlPntVars="2" numCnt'
        'rlPnts="6" numOfEleVars="0" closed="0">\n'
    ),
    "    <cntrlPntVarNames>x y</cntrlPntVarNames>\n",
    "    <cntrlPntVars>0.0 0.0\n",
    "0.5 0.0\n",
    "1.0 0.0\n",
    "0.0 1.0\n",
    "0.5 1.0\n",
    "1.0 1.0</cntrlPntVars>\n",
    "    <deg>2 1</deg>\n",
    "    <kntVecs>\n",
    "      <knotVec>0.0 0.0 0.0 1.0 1.0 1.0</knotVec>\n",
    "      <knotVec>0.0 0.0 1.0 1.0</knotVec>\n",
    "    </kntVecs>\n",
    "  </SplineEntry>\n",
    (
        '  <SplineEntry splDim="2" spaceDim="2" numOfCntrlPntVars="2" numCnt'
        'rlPnts="4" numOfEleVars="0" closed="0">\n'
    ),
    "    <cntrlPntVarNames>x y</cntrlPntVarNames>\n",
    "    <cntrlPntVars>1.0 0.0\n",
    "2.0 0.0\n",
    "1.0 1.0\n",
    "2.0 1.0</cntrlPntVars>\n",
    "    <deg>1 1</deg>\n",
    "    <kntVecs>\n",
    "      <knotVec>0.0 0.0 1.0 1.0</knotVec>\n",
    "      <knotVec>0.0 0.0 1.0 1.0</knotVec>\n",
    "    </kntVecs>\n",
    "    <wght>1.0 1.0 1.0 1.0</wght>\n",
    "  </SplineEntry>\n",
    (
        '  <SplineEntry splDim="2" spaceDim="2" numOfCntrlPntVars="2" numCnt'
        'rlPnts="9" numOfEleVars="0" closed="0">\n'
    ),
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
    "    <deg>2 1</deg>\n",
    "    <kntVecs>\n",
    "      <knotVec>0.0 0.0 0.0 1.0 1.0 1.0</knotVec>\n",
    "      <knotVec>0.0 0.0 0.5 1.0 1.0</knotVec>\n",
    "    </kntVecs>\n",
    "  </SplineEntry>\n",
    (
        '  <SplineEntry splDim="2" spaceDim="2" numOfCntrlPntVars="2" numCnt'
        'rlPnts="6" numOfEleVars="0" closed="0">\n'
    ),
    "    <cntrlPntVarNames>x y</cntrlPntVarNames>\n",
    "    <cntrlPntVars>1.0 1.0\n",
    "2.0 1.0\n",
    "1.0 1.5\n",
    "2.0 1.5\n",
    "1.0 2.0\n",
    "2.0 2.0</cntrlPntVars>\n",
    "    <deg>1 1</deg>\n",
    "    <kntVecs>\n",
    "      <knotVec>0.0 0.0 1.0 1.0</knotVec>\n",
    "      <knotVec>0.0 0.0 0.5 1.0 1.0</knotVec>\n",
    "    </kntVecs>\n",
    "    <wght>1.0 1.0 1.0 1.0 1.0 1.0</wght>\n",
    "  </SplineEntry>\n",
    "</SplineList>",
]


class CATSioTest(c.unittest.TestCase):
    def test_cats_export(self):
        """
        Test cats export routine
        """
        if int(python_version.split(".")[1]) < 8:
            c.splinepy.utils.log.info(
                "cats export is only tested here from python3.8+. "
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

        # Test Output
        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = c.to_tmpf(tmpd)
            c.splinepy.io.cats.export(
                tmpf,
                spline_list=multipatch,
                indent=False,
            )

            with open(tmpf) as tmp_read:
                self.assertTrue(
                    c.are_stripped_lines_same(
                        cats_no_indent_export, tmp_read.readlines(), True
                    )
                )

        # for python version > 3.9, test indented version
        if int(python_version.split(".")[1]) >= 9:
            with tempfile.TemporaryDirectory() as tmpd:
                tmpf = c.to_tmpf(tmpd)
                c.splinepy.io.cats.export(
                    tmpf,
                    spline_list=multipatch,
                    indent=True,
                )

                with open(tmpf) as tmp_read:
                    self.assertTrue(
                        c.are_items_same(
                            cats_indent_export, tmp_read.readlines()
                        )
                    )

    def test_cats_import(self):
        """
        Test cats export routine
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
        bsp_el2.elevate_degrees(0)
        bsp_el2.elevate_degrees(1)
        nur_el3.elevate_degrees(0)
        nur_el3.elevate_degrees(1)
        bsp_el2.insert_knots(1, [0.5])
        nur_el3.insert_knots(1, [0.5])

        # Test Output against input
        multipatch_geometry = [bsp_el2, nur_el3]
        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = c.to_tmpf(tmpd)
            c.splinepy.io.cats.export(
                tmpf,
                spline_list=multipatch_geometry,
                indent=False,
            )
            multipatch_geometry_loaded = c.splinepy.io.cats.load(tmpf)
            self.assertTrue(
                all(
                    c.are_splines_equal(a, b)
                    for a, b in zip(
                        multipatch_geometry,
                        multipatch_geometry_loaded,
                    )
                )
            )


if __name__ == "__main__":
    c.unittest.main()
