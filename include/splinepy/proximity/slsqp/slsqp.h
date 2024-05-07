#pragma once
/*      ALGORITHM 733, COLLECTED ALGORITHMS FROM ACM. */
/*      TRANSACTIONS ON MATHEMATICAL SOFTWARE, */
/*      VOL. 20, NO. 3, SEPTEMBER, 1994, PP. 262-281. */
/*      https://doi.org/10.1145/192115.192124 */

/*      https://web.archive.org/web/20170106155705/http://permalink.gmane.org/gmane.comp.python.scientific.dev
el/6725 */
/*      ------ */
/*      From: Deborah Cotton <cotton@hq.acm.org> */
/*      Date: Fri, 14 Sep 2007 12:35:55 -0500 */
/*      Subject: RE: Algorithm License requested */
/*      To: Alan Isaac */

/*      Prof. Issac, */

/*      In that case, then because the author consents to [the ACM] releasing */
/*      the code currently archived at http://www.netlib.org/toms/733 under the
 */
/*      BSD license, the ACM hereby releases this code under the BSD license. */

/*      Regards, */

/*      Deborah Cotton, Copyright & Permissions */
/*      ACM Publications */
/*      2 Penn Plaza, Suite 701** */
/*      New York, NY 10121-0701 */
/*      permissions@acm.org */
/*      212.869.7440 ext. 652 */
/*      Fax. 212.869.0481 */
/*      ------ */

#ifdef __cplusplus
extern "C" {
#endif

/* *********************************************************************** */
/*                              optimizer                               * */
/* *********************************************************************** */
/*   SLSQP       S EQUENTIAL  L EAST  SQ UARES  P ROGRAMMING */
/*            TO SOLVE GENERAL NONLINEAR OPTIMIZATION PROBLEMS */
/* *********************************************************************** */
/* *                                                                     * */
/* *                                                                     * */
/* *            A NONLINEAR PROGRAMMING METHOD WITH                      * */
/* *            QUADRATIC  PROGRAMMING  SUBPROBLEMS                      * */
/* *                                                                     * */
/* *                                                                     * */
/* *  THIS SUBROUTINE SOLVES THE GENERAL NONLINEAR PROGRAMMING PROBLEM   * */
/* *                                                                     * */
/* *            MINIMIZE    F(X)                                         * */
/* *                                                                     * */
/* *            SUBJECT TO  C (X) .EQ. 0  ,  J = 1,...,MEQ               * */
/* *                         J                                           * */
/* *                                                                     * */
/* *                        C (X) .GE. 0  ,  J = MEQ+1,...,M             * */
/* *                         J                                           * */
/* *                                                                     * */
/* *                        XL .LE. X .LE. XU , I = 1,...,N.             * */
/* *                          I      I       I                           * */
/* *                                                                     * */
/* *  THE ALGORITHM IMPLEMENTS THE METHOD OF HAN AND POWELL              * */
/* *  WITH BFGS-UPDATE OF THE B-MATRIX AND L1-TEST FUNCTION              * */
/* *  WITHIN THE STEPLENGTH ALGORITHM.                                   * */
/* *                                                                     * */
/* *    PARAMETER DESCRIPTION:                                           * */
/* *    ( * MEANS THIS PARAMETER WILL BE CHANGED DURING CALCULATION )    * */
/* *                                                                     * */
/* *    M              IS THE TOTAL NUMBER OF CONSTRAINTS, M .GE. 0      * */
/* *    MEQ            IS THE NUMBER OF EQUALITY CONSTRAINTS, MEQ .GE. 0 * */
/* *    LA             SEE A, LA .GE. MAX(M,1)                           * */
/* *    N              IS THE NUMBER OF VARIABLES, N .GE. 1               * */
/* *  * X()            X() STORES THE CURRENT ITERATE OF THE N VECTOR X  * */
/* *                   ON ENTRY X() MUST BE INITIALIZED. ON EXIT X()     * */
/* *                   STORES THE SOLUTION VECTOR X IF MODE = 0.         * */
/* *    XL()           XL() STORES AN N VECTOR OF LOWER BOUNDS XL TO X.  * */
/* *                   ELEMENTS MAY BE NAN TO INDICATE NO LOWER BOUND.   * */
/* *    XU()           XU() STORES AN N VECTOR OF UPPER BOUNDS XU TO X.  * */
/* *                   ELEMENTS MAY BE NAN TO INDICATE NO UPPER BOUND.   * */
/* *    F              IS THE VALUE OF THE OBJECTIVE FUNCTION.           * */
/* *    C()            C() STORES THE M VECTOR C OF CONSTRAINTS,         * */
/* *                   EQUALITY CONSTRAINTS (IF ANY) FIRST.              * */
/* *                   DIMENSION OF C MUST BE GREATER OR EQUAL LA,       * */
/* *                   which must be GREATER OR EQUAL MAX(1,M).          * */
/* *    G()            G() STORES THE N VECTOR G OF PARTIALS OF THE      * */
/* *                   OBJECTIVE FUNCTION; DIMENSION OF G MUST BE        * */
/* *                   GREATER OR EQUAL N+1.                             * */
/* *    A(),LA,M,N     THE LA BY N + 1 ARRAY A() STORES                  * */
/* *                   THE M BY N MATRIX A OF CONSTRAINT NORMALS.        * */
/* *                   A() HAS FIRST DIMENSIONING PARAMETER LA,          * */
/* *                   WHICH MUST BE GREATER OR EQUAL MAX(1,M).          * */
/* *    F,C,G,A        MUST ALL BE SET BY THE USER BEFORE EACH CALL.     * */
/* *  * ACC            ABS(ACC) CONTROLS THE FINAL ACCURACY.             * */
/* *                   IF ACC .LT. ZERO AN EXACT LINESEARCH IS PERFORMED,* */
/* *                   OTHERWISE AN ARMIJO-TYPE LINESEARCH IS USED.      * */
/* *  * ITER           PRESCRIBES THE MAXIMUM NUMBER OF ITERATIONS.      * */
/* *                   ON EXIT ITER INDICATES THE NUMBER OF ITERATIONS.  * */
/* *  * MODE           MODE CONTROLS CALCULATION:                        * */
/* *                   REVERSE COMMUNICATION IS USED IN THE SENSE THAT   * */
/* *                   THE PROGRAM IS INITIALIZED BY MODE = 0; THEN IT IS* */
/* *                   TO BE CALLED REPEATEDLY BY THE USER UNTIL A RETURN* */
/* *                   WITH MODE .NE. IABS(1) TAKES PLACE.               * */
/* *                   IF MODE = -1 GRADIENTS HAVE TO BE CALCULATED,     * */
/* *                   WHILE WITH MODE = 1 FUNCTIONS HAVE TO BE CALCULATED */
/* *                   MODE MUST NOT BE CHANGED BETWEEN SUBSEQUENT CALLS * */
/* *                   OF SQP.                                           * */
/* *                   EVALUATION MODES:                                 * */
/* *        MODE = -1: GRADIENT EVALUATION, (G&A)                        * */
/* *                0: ON ENTRY: INITIALIZATION, (F,G,C&A)               * */
/* *                   ON EXIT : REQUIRED ACCURACY FOR SOLUTION OBTAINED * */
/* *                1: FUNCTION EVALUATION, (F&C)                        * */
/* *                                                                     * */
/* *                   FAILURE MODES:                                    * */
/* *                2: NUMBER OF EQUALITY CONSTRAINTS LARGER THAN N      * */
/* *                3: MORE THAN 3*N ITERATIONS IN LSQ SUBPROBLEM        * */
/* *                4: INEQUALITY CONSTRAINTS INCOMPATIBLE               * */
/* *                5: SINGULAR MATRIX E IN LSQ SUBPROBLEM               * */
/* *                6: SINGULAR MATRIX C IN LSQ SUBPROBLEM               * */
/* *                7: RANK-DEFICIENT EQUALITY CONSTRAINT SUBPROBLEM HFTI* */
/* *                8: POSITIVE DIRECTIONAL DERIVATIVE FOR LINESEARCH    * */
/* *                9: MORE THAN ITER ITERATIONS IN SQP                  * */
/* *             >=10: WORKING SPACE W OR JW TOO SMALL,                  * */
/* *                   W SHOULD BE ENLARGED TO L_W=MODE/1000             * */
/* *                   JW SHOULD BE ENLARGED TO L_JW=MODE-1000*L_W       * */
/* *  * W(), L_W       W() IS A ONE DIMENSIONAL WORKING SPACE,           * */
/* *                   THE LENGTH L_W OF WHICH SHOULD BE AT LEAST        * */
/* *                   (3*N1+M)*(N1+1)                        for LSQ    * */
/* *                  +(N1-MEQ+1)*(MINEQ+2) + 2*MINEQ         for LSI    * */
/* *                  +(N1+MINEQ)*(N1-MEQ) + 2*MEQ + N1       for LSEI   * */
/* *                  + N1*N/2 + 2*M + 3*N + 3*N1 + 1         for SLSQPB * */
/* *                   with MINEQ = M - MEQ + 2*N1  &  N1 = N+1          * */
/* *        NOTICE:    FOR PROPER DIMENSIONING OF W IT IS RECOMMENDED TO * */
/* *                   COPY THE FOLLOWING STATEMENTS INTO THE HEAD OF    * */
/* *                   THE CALLING PROGRAM (AND REMOVE THE COMMENT C)    * */
/* ####################################################################### */
/*     INTEGER LEN_W, LEN_JW, M, N, N1, MEQ, MINEQ */
/*     PARAMETER (M=... , MEQ=... , N=...  ) */
/*     PARAMETER (N1= N+1, MINEQ= M-MEQ+N1+N1) */
/*     PARAMETER (LEN_W= */
/*    $           (3*N1+M)*(N1+1) */
/*    $          +(N1-MEQ+1)*(MINEQ+2) + 2*MINEQ */
/*    $          +(N1+MINEQ)*(N1-MEQ) + 2*MEQ + N1 */
/*    $          +(N+1)*N/2 + 2*M + 3*N + 3*N1 + 1, */
/*    $           LEN_JW=MINEQ) */
/*     DOUBLE PRECISION W(LEN_W) */
/*     INTEGER          JW(LEN_JW) */
/* ####################################################################### */
/* *                   THE FIRST M+N+N*N1/2 ELEMENTS OF W MUST NOT BE    * */
/* *                   CHANGED BETWEEN SUBSEQUENT CALLS OF SLSQP.        * */
/* *                   ON RETURN W(1) ... W(M) CONTAIN THE MULTIPLIERS   * */
/* *                   ASSOCIATED WITH THE GENERAL CONSTRAINTS, WHILE    * */
/* *                   W(M+1) ... W(M+N(N+1)/2) STORE THE CHOLESKY FACTOR* */
/* *                   L*D*L(T) OF THE APPROXIMATE HESSIAN OF THE        * */
/* *                   LAGRANGIAN COLUMNWISE DENSE AS LOWER TRIANGULAR   * */
/* *                   UNIT MATRIX L WITH D IN ITS 'DIAGONAL' and        * */
/* *                   W(M+N(N+1)/2+N+2 ... W(M+N(N+1)/2+N+2+M+2N)       * */
/* *                   CONTAIN THE MULTIPLIERS ASSOCIATED WITH ALL       * */
/* *                   ALL CONSTRAINTS OF THE QUADRATIC PROGRAM FINDING  * */
/* *                   THE SEARCH DIRECTION TO THE SOLUTION X*           * */
/* *  * JW(), L_JW     JW() IS A ONE DIMENSIONAL INTEGER WORKING SPACE   * */
/* *                   THE LENGTH L_JW OF WHICH SHOULD BE AT LEAST       * */
/* *                   MINEQ                                             * */
/* *                   with MINEQ = M - MEQ + 2*N1  &  N1 = N+1          * */
/* *                                                                     * */
/* *  THE USER HAS TO PROVIDE THE FOLLOWING SUBROUTINES:                 * */
/* *     LDL(N,A,Z,SIG,W) :   UPDATE OF THE LDL'-FACTORIZATION.          * */
/* *     LINMIN(A,B,F,TOL) :  LINESEARCH ALGORITHM IF EXACT = 1          * */
/* *     LSQ(M,MEQ,LA,N,NC,C,D,A,B,XL,XU,X,LAMBDA,W,....) :              * */
/* *                                                                     * */
/* *        SOLUTION OF THE QUADRATIC PROGRAM                            * */
/* *                QPSOL IS RECOMMENDED:                                * */
/* *     PE GILL, W MURRAY, MA SAUNDERS, MH WRIGHT:                      * */
/* *     USER'S GUIDE FOR SOL/QPSOL:                                     * */
/* *     A FORTRAN PACKAGE FOR QUADRATIC PROGRAMMING,                    * */
/* *     TECHNICAL REPORT SOL 83-7, JULY 1983                            * */
/* *     DEPARTMENT OF OPERATIONS RESEARCH, STANFORD UNIVERSITY          * */
/* *     STANFORD, CA 94305                                              * */
/* *     QPSOL IS THE MOST ROBUST AND EFFICIENT QP-SOLVER                * */
/* *     AS IT ALLOWS WARM STARTS WITH PROPER WORKING SETS               * */
/* *                                                                     * */
/* *     IF IT IS NOT AVAILABLE USE LSEI, A CONSTRAINT LINEAR LEAST      * */
/* *     SQUARES SOLVER IMPLEMENTED USING THE SOFTWARE HFTI, LDP, NNLS   * */
/* *     FROM C.L. LAWSON, R.J.HANSON: SOLVING LEAST SQUARES PROBLEMS,   * */
/* *     PRENTICE HALL, ENGLEWOOD CLIFFS, 1974.                          * */
/* *     LSEI COMES WITH THIS PACKAGE, together with all necessary SR's. * */
/* *                                                                     * */
/* *     TOGETHER WITH A COUPLE OF SUBROUTINES FROM BLAS LEVEL 1         * */
/* *                                                                     * */
/* *     SQP IS HEAD SUBROUTINE FOR BODY SUBROUTINE SQPBDY               * */
/* *     IN WHICH THE ALGORITHM HAS BEEN IMPLEMENTED.                    * */
/* *                                                                     * */
/* *  IMPLEMENTED BY: DIETER KRAFT, DFVLR OBERPFAFFENHOFEN               * */
/* *  as described in Dieter Kraft: A Software Package for               * */
/* *                                Sequential Quadratic Programming     * */
/* *                                DFVLR-FB 88-28, 1988                 * */
/* *  which should be referenced if the user publishes results of SLSQP  * */
/* *                                                                     * */
/* *  DATE:           APRIL - OCTOBER, 1981.                             * */
/* *  STATUS:         DECEMBER, 31-ST, 1984.                             * */
/* *  STATUS:         MARCH   , 21-ST, 1987, REVISED TO FORTRAN 77       * */
/* *  STATUS:         MARCH   , 20-th, 1989, REVISED TO MS-FORTRAN       * */
/* *  STATUS:         APRIL   , 14-th, 1989, HESSE   in-line coded       * */
/* *  STATUS:         FEBRUARY, 28-th, 1991, FORTRAN/2 Version 1.04      * */
/* *                                         accepts Statement Functions * */
/* *  STATUS:         MARCH   ,  1-st, 1991, tested with SALFORD         * */
/* *                                         FTN77/386 COMPILER VERS 2.40* */
/* *                                         in protected mode           * */
/* *                                                                     * */
/* *********************************************************************** */
/* *                                                                     * */
/* *  Copyright 1991: Dieter Kraft, FHM                                  * */
/* *                                                                     * */
/* *********************************************************************** */
/*     dim(W) =         N1*(N1+1) + MEQ*(N1+1) + MINEQ*(N1+1)  for LSQ */
/*                    +(N1-MEQ+1)*(MINEQ+2) + 2*MINEQ          for LSI */
/*                    +(N1+MINEQ)*(N1-MEQ) + 2*MEQ + N1        for LSEI */
/*                    + N1*N/2 + 2*M + 3*N +3*N1 + 1           for SLSQPB */
/*                      with MINEQ = M - MEQ + 2*N1  &  N1 = N+1 */
/*   CHECK LENGTH OF WORKING ARRAYS */
int slsqp(int* m,
          int* meg,
          int* la,
          int* n,
          double* x,
          double* xl,
          double* xu,
          double* f,
          double* c,
          double* g,
          double* a,
          double* acc,
          int* iter,
          int* mode,
          double* w,
          int* l_w,
          int* jw,
          int* l_jw,
          double* alpha,
          double* f0,
          double* gs,
          double* h1,
          double* h2,
          double* h3,
          double* h4,
          double* t,
          double* t0,
          double* tol,
          int* iexact,
          int* incons,
          int* ireset,
          int* itermx,
          int* line,
          int* n1,
          int* n2,
          int* n3);

/*   NONLINEAR PROGRAMMING BY SOLVING SEQUENTIALLY QUADRATIC PROGRAMS */
/*        -  L1 - LINE SEARCH,  POSITIVE DEFINITE  BFGS UPDATE  - */
/*                      BODY SUBROUTINE FOR SLSQP */
/*     dim(W) =         N1*(N1+1) + MEQ*(N1+1) + MINEQ*(N1+1)  for LSQ */
/*                     +(N1-MEQ+1)*(MINEQ+2) + 2*MINEQ */
/*                     +(N1+MINEQ)*(N1-MEQ) + 2*MEQ + N1       for LSEI */
/*                      with MINEQ = M - MEQ + 2*N1  &  N1 = N+1 */
int slsqpb(int* m,
           int* meq,
           int* la,
           int* n,
           double* x,
           double* xl,
           double* xu,
           double* f,
           double* c__,
           double* g,
           double* a,
           double* acc,
           int* iter,
           int* mode,
           double* r__,
           double* l,
           double* x0,
           double* mu,
           double* s,
           double* u,
           double* v,
           double* w,
           int* iw,
           double* alpha,
           double* f0,
           double* gs,
           double* h1,
           double* h2,
           double* h3,
           double* h4,
           double* t,
           double* t0,
           double* tol,
           int* iexact,
           int* incons,
           int* ireset,
           int* itermx,
           int* line,
           int* n1,
           int* n2,
           int* n3);

/*   MINIMIZE with respect to X */
/*             ||E*X - F|| */
/*                                      1/2  T */
/*   WITH UPPER TRIANGULAR MATRIX E = +D   *L , */
/*                                      -1/2  -1 */
/*                     AND VECTOR F = -D    *L  *G, */
/*  WHERE THE UNIT LOWER TRIDIANGULAR MATRIX L IS STORED COLUMNWISE */
/*  DENSE IN THE N*(N+1)/2 ARRAY L WITH VECTOR D STORED IN ITS */
/* 'DIAGONAL' THUS SUBSTITUTING THE ONE-ELEMENTS OF L */
/*   SUBJECT TO */
/*             A(J)*X - B(J) = 0 ,         J=1,...,MEQ, */
/*             A(J)*X - B(J) >=0,          J=MEQ+1,...,M, */
/*             XL(I) <= X(I) <= XU(I),     I=1,...,N, */
/*     ON ENTRY, THE USER HAS TO PROVIDE THE ARRAYS L, G, A, B, XL, XU. */
/*     WITH DIMENSIONS: L(N*(N+1)/2), G(N), A(LA,N), B(M), XL(N), XU(N) */
/*     THE WORKING ARRAY W MUST HAVE AT LEAST THE FOLLOWING DIMENSION: */
/*     DIM(W) =        (3*N+M)*(N+1)                        for LSQ */
/*                    +(N-MEQ+1)*(MINEQ+2) + 2*MINEQ        for LSI */
/*                    +(N+MINEQ)*(N-MEQ) + 2*MEQ + N        for LSEI */
/*                      with MINEQ = M - MEQ + 2*N */
/*     ON RETURN, NO ARRAY WILL BE CHANGED BY THE SUBROUTINE. */
/*     X     STORES THE N-DIMENSIONAL SOLUTION VECTOR */
/*     Y     STORES THE VECTOR OF LAGRANGE MULTIPLIERS OF DIMENSION */
/*           M+N+N (CONSTRAINTS+LOWER+UPPER BOUNDS) */
/*     MODE  IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANINGS: */
/*          MODE=1: SUCCESSFUL COMPUTATION */
/*               2: ERROR RETURN BECAUSE OF WRONG DIMENSIONS (N<1) */
/*               3: ITERATION COUNT EXCEEDED BY NNLS */
/*               4: INEQUALITY CONSTRAINTS INCOMPATIBLE */
/*               5: MATRIX E IS NOT OF FULL RANK */
/*               6: MATRIX C IS NOT OF FULL RANK */
/*               7: RANK DEFECT IN HFTI */
/*     coded            Dieter Kraft, april 1987 */
/*     revised                        march 1989 */
int lsq(const int* m,
        const int* meq,
        const int* n,
        const int* nl,
        const int* la,
        const double* l,
        const double* g,
        const double* a,
        const double* b,
        const double* xl,
        const double* xu,
        double* x,
        double* y,
        double* w,
        int* jw,
        int* mode);

/*     FOR MODE=1, THE SUBROUTINE RETURNS THE SOLUTION X OF */
/*     EQUALITY & INEQUALITY CONSTRAINED LEAST SQUARES PROBLEM LSEI : */
/*                MIN ||E*X - F|| */
/*                 X */
/*                S.T.  C*X  = D, */
/*                      G*X >= H. */
/*     USING QR DECOMPOSITION & ORTHOGONAL BASIS OF NULLSPACE OF C */
/*     CHAPTER 23.6 OF LAWSON & HANSON: SOLVING LEAST SQUARES PROBLEMS. */
/*     THE FOLLOWING DIMENSIONS OF THE ARRAYS DEFINING THE PROBLEM */
/*     ARE NECESSARY */
/*     DIM(E) :   FORMAL (LE,N),    ACTUAL (ME,N) */
/*     DIM(F) :   FORMAL (LE  ),    ACTUAL (ME  ) */
/*     DIM(C) :   FORMAL (LC,N),    ACTUAL (MC,N) */
/*     DIM(D) :   FORMAL (LC  ),    ACTUAL (MC  ) */
/*     DIM(G) :   FORMAL (LG,N),    ACTUAL (MG,N) */
/*     DIM(H) :   FORMAL (LG  ),    ACTUAL (MG  ) */
/*     DIM(X) :   FORMAL (N   ),    ACTUAL (N   ) */
/*     DIM(W) :   2*MC+ME+(ME+MG)*(N-MC)  for LSEI */
/*              +(N-MC+1)*(MG+2)+2*MG     for LSI */
/*     DIM(JW):   MAX(MG,L) */
/*     ON ENTRY, THE USER HAS TO PROVIDE THE ARRAYS C, D, E, F, G, AND H. */
/*     ON RETURN, ALL ARRAYS WILL BE CHANGED BY THE SUBROUTINE. */
/*     X     STORES THE SOLUTION VECTOR */
/*     XNORM STORES THE RESIDUUM OF THE SOLUTION IN EUCLIDEAN NORM */
/*     W     STORES THE VECTOR OF LAGRANGE MULTIPLIERS IN ITS FIRST */
/*           MC+MG ELEMENTS */
/*     MODE  IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANINGS: */
/*          MODE=1: SUCCESSFUL COMPUTATION */
/*               2: ERROR RETURN BECAUSE OF WRONG DIMENSIONS (N<1) */
/*               3: ITERATION COUNT EXCEEDED BY NNLS */
/*               4: INEQUALITY CONSTRAINTS INCOMPATIBLE */
/*               5: MATRIX E IS NOT OF FULL RANK */
/*               6: MATRIX C IS NOT OF FULL RANK */
/*               7: RANK DEFECT IN HFTI */
/*     18.5.1981, DIETER KRAFT, DFVLR OBERPFAFFENHOFEN */
/*     20.3.1987, DIETER KRAFT, DFVLR OBERPFAFFENHOFEN */
int lsei(double* c__,
         double* d__,
         double* e,
         double* f,
         double* g,
         double* h__,
         const int* lc,
         const int* mc,
         const int* le,
         const int* me,
         const int* lg,
         const int* mg,
         const int* n,
         double* x,
         double* xnrm,
         double* w,
         int* jw,
         int* mode);

/*     FOR MODE=1, THE SUBROUTINE RETURNS THE SOLUTION X OF */
/*     INEQUALITY CONSTRAINED LINEAR LEAST SQUARES PROBLEM: */
/*                    MIN ||E*X-F|| */
/*                     X */
/*                    S.T.  G*X >= H */
/*     THE ALGORITHM IS BASED ON QR DECOMPOSITION AS DESCRIBED IN */
/*     CHAPTER 23.5 OF LAWSON & HANSON: SOLVING LEAST SQUARES PROBLEMS */
/*     THE FOLLOWING DIMENSIONS OF THE ARRAYS DEFINING THE PROBLEM */
/*     ARE NECESSARY */
/*     DIM(E) :   FORMAL (LE,N),    ACTUAL (ME,N) */
/*     DIM(F) :   FORMAL (LE  ),    ACTUAL (ME  ) */
/*     DIM(G) :   FORMAL (LG,N),    ACTUAL (MG,N) */
/*     DIM(H) :   FORMAL (LG  ),    ACTUAL (MG  ) */
/*     DIM(X) :   N */
/*     DIM(W) :   (N+1)*(MG+2) + 2*MG */
/*     DIM(JW):   LG */
/*     ON ENTRY, THE USER HAS TO PROVIDE THE ARRAYS E, F, G, AND H. */
/*     ON RETURN, ALL ARRAYS WILL BE CHANGED BY THE SUBROUTINE. */
/*     X     STORES THE SOLUTION VECTOR */
/*     XNORM STORES THE RESIDUUM OF THE SOLUTION IN EUCLIDEAN NORM */
/*     W     STORES THE VECTOR OF LAGRANGE MULTIPLIERS IN ITS FIRST */
/*           MG ELEMENTS */
/*     MODE  IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANINGS: */
/*          MODE=1: SUCCESSFUL COMPUTATION */
/*               2: ERROR RETURN BECAUSE OF WRONG DIMENSIONS (N<1) */
/*               3: ITERATION COUNT EXCEEDED BY NNLS */
/*               4: INEQUALITY CONSTRAINTS INCOMPATIBLE */
/*               5: MATRIX E IS NOT OF FULL RANK */
/*     03.01.1980, DIETER KRAFT: CODED */
/*     20.03.1987, DIETER KRAFT: REVISED TO FORTRAN 77 */
int lsi(double* e,
        double* f,
        double* g,
        double* h__,
        const int* le,
        const int* me,
        const int* lg,
        const int* mg,
        const int* n,
        double* x,
        double* xnorm,
        double* w,
        int* jw,
        int* mode);

/*                     T */
/*     MINIMIZE   1/2 X X    SUBJECT TO   G * X >= H. */
/*       C.L. LAWSON, R.J. HANSON: 'SOLVING LEAST SQUARES PROBLEMS' */
/*       PRENTICE HALL, ENGLEWOOD CLIFFS, NEW JERSEY, 1974. */
/*     PARAMETER DESCRIPTION: */
/*     G(),MG,M,N   ON ENTRY G() STORES THE M BY N MATRIX OF */
/*                  LINEAR INEQUALITY CONSTRAINTS. G() HAS FIRST */
/*                  DIMENSIONING PARAMETER MG */
/*     H()          ON ENTRY H() STORES THE M VECTOR H REPRESENTING */
/*                  THE RIGHT SIDE OF THE INEQUALITY SYSTEM */
/*     REMARK: G(),H() WILL NOT BE CHANGED DURING CALCULATIONS BY LDP */
/*     X()          ON ENTRY X() NEED NOT BE INITIALIZED. */
/*                  ON EXIT X() STORES THE SOLUTION VECTOR X IF MODE=1. */
/*     XNORM        ON EXIT XNORM STORES THE EUCLIDEAN NORM OF THE */
/*                  SOLUTION VECTOR IF COMPUTATION IS SUCCESSFUL */
/*     W()          W IS A ONE DIMENSIONAL WORKING SPACE, THE LENGTH */
/*                  OF WHICH SHOULD BE AT LEAST (M+2)*(N+1) + 2*M */
/*                  ON EXIT W() STORES THE LAGRANGE MULTIPLIERS */
/*                  ASSOCIATED WITH THE CONSTRAINTS */
/*                  AT THE SOLUTION OF PROBLEM LDP */
/*     INDEX()      INDEX() IS A ONE DIMENSIONAL INTEGER WORKING SPACE */
/*                  OF LENGTH AT LEAST M */
/*     MODE         MODE IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING */
/*                  MEANINGS: */
/*          MODE=1: SUCCESSFUL COMPUTATION */
/*               2: ERROR RETURN BECAUSE OF WRONG DIMENSIONS (N.LE.0) */
/*               3: ITERATION COUNT EXCEEDED BY NNLS */
/*               4: INEQUALITY CONSTRAINTS INCOMPATIBLE */
int ldp(const double* g,
        const int* mg,
        const int* m,
        const int* n,
        const double* h__,
        double* x,
        double* xnorm,
        double* w,
        int* index,
        int* mode);

/*     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY: */
/*     'SOLVING LEAST SQUARES PROBLEMS'. PRENTICE-HALL.1974 */
/*      **********   NONNEGATIVE LEAST SQUARES   ********** */
/*     GIVEN AN M BY N MATRIX, A, AND AN M-VECTOR, B, COMPUTE AN */
/*     N-VECTOR, X, WHICH SOLVES THE LEAST SQUARES PROBLEM */
/*                  A*X = B  SUBJECT TO  X >= 0 */
/*     A(),MDA,M,N */
/*            MDA IS THE FIRST DIMENSIONING PARAMETER FOR THE ARRAY,A(). */
/*            ON ENTRY A()  CONTAINS THE M BY N MATRIX,A. */
/*            ON EXIT A() CONTAINS THE PRODUCT Q*A, */
/*            WHERE Q IS AN M BY M ORTHOGONAL MATRIX GENERATED */
/*            IMPLICITLY BY THIS SUBROUTINE. */
/*            EITHER M>=N OR M<N IS PERMISSIBLE. */
/*            THERE IS NO RESTRICTION ON THE RANK OF A. */
/*     B()    ON ENTRY B() CONTAINS THE M-VECTOR, B. */
/*            ON EXIT B() CONTAINS Q*B. */
/*     X()    ON ENTRY X() NEED NOT BE INITIALIZED. */
/*            ON EXIT X() WILL CONTAIN THE SOLUTION VECTOR. */
/*     RNORM  ON EXIT RNORM CONTAINS THE EUCLIDEAN NORM OF THE */
/*            RESIDUAL VECTOR. */
/*     W()    AN N-ARRAY OF WORKING SPACE. */
/*            ON EXIT W() WILL CONTAIN THE DUAL SOLUTION VECTOR. */
/*            W WILL SATISFY W(I)=0 FOR ALL I IN SET P */
/*            AND W(I)<=0 FOR ALL I IN SET Z */
/*     Z()    AN M-ARRAY OF WORKING SPACE. */
/*     INDEX()AN INTEGER WORKING ARRAY OF LENGTH AT LEAST N. */
/*            ON EXIT THE CONTENTS OF THIS ARRAY DEFINE THE SETS */
/*            P AND Z AS FOLLOWS: */
/*            INDEX(1)    THRU INDEX(NSETP) = SET P. */
/*            INDEX(IZ1)  THRU INDEX (IZ2)  = SET Z. */
/*            IZ1=NSETP + 1 = NPP1, IZ2=N. */
/*     MODE   THIS IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANING: */
/*            1    THE SOLUTION HAS BEEN COMPUTED SUCCESSFULLY. */
/*            2    THE DIMENSIONS OF THE PROBLEM ARE WRONG, */
/*                 EITHER M <= 0 OR N <= 0. */
/*            3    ITERATION COUNT EXCEEDED, MORE THAN 3*N ITERATIONS. */
int nnls(double* a,
         const int* mda,
         const int* m,
         const int* n,
         double* b,
         double* x,
         double* rnorm,
         double* w,
         double* z__,
         int* index,
         int* mode);

/*     RANK-DEFICIENT LEAST SQUARES ALGORITHM AS DESCRIBED IN: */
/*     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY, 1973 JUN 12 */
/*     TO APPEAR IN 'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974 */
/*     A(*,*),MDA,M,N   THE ARRAY A INITIALLY CONTAINS THE M x N MATRIX A */
/*                      OF THE LEAST SQUARES PROBLEM AX = B. */
/*                      THE FIRST DIMENSIONING PARAMETER MDA MUST SATISFY */
/*                      MDA >= M. EITHER M >= N OR M < N IS PERMITTED. */
/*                      THERE IS NO RESTRICTION ON THE RANK OF A. */
/*                      THE MATRIX A WILL BE MODIFIED BY THE SUBROUTINE. */
/*     B(*,*),MDB,NB    IF NB = 0 THE SUBROUTINE WILL MAKE NO REFERENCE */
/*                      TO THE ARRAY B. IF NB > 0 THE ARRAY B() MUST */
/*                      INITIALLY CONTAIN THE M x NB MATRIX B  OF THE */
/*                      THE LEAST SQUARES PROBLEM AX = B AND ON RETURN */
/*                      THE ARRAY B() WILL CONTAIN THE N x NB SOLUTION X. */
/*                      IF NB>1 THE ARRAY B() MUST BE DOUBLE SUBSCRIPTED */
/*                      WITH FIRST DIMENSIONING PARAMETER MDB>=MAX(M,N), */
/*                      IF NB=1 THE ARRAY B() MAY BE EITHER SINGLE OR */
/*                      DOUBLE SUBSCRIPTED. */
/*     TAU              ABSOLUTE TOLERANCE PARAMETER FOR PSEUDORANK */
/*                      DETERMINATION, PROVIDED BY THE USER. */
/*     KRANK            PSEUDORANK OF A, SET BY THE SUBROUTINE. */
/*     RNORM            ON EXIT, RNORM(J) WILL CONTAIN THE EUCLIDEAN */
/*                      NORM OF THE RESIDUAL VECTOR FOR THE PROBLEM */
/*                      DEFINED BY THE J-TH COLUMN VECTOR OF THE ARRAY B. */
/*     H(), G()         ARRAYS OF WORKING SPACE OF LENGTH >= N. */
/*     IP()             INTEGER ARRAY OF WORKING SPACE OF LENGTH >= N */
/*                      RECORDING PERMUTATION INDICES OF COLUMN VECTORS */
int hfti(double* a,
         const int* mda,
         const int* m,
         const int* n,
         double* b,
         const int* mdb,
         const int* nb,
         const double* tau,
         int* krank,
         double* rnorm,
         double* h__,
         double* g,
         int* ip);

/*     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY, 1973 JUN 12 */
/*     TO APPEAR IN 'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974 */
/*     CONSTRUCTION AND/OR APPLICATION OF A SINGLE */
/*     HOUSEHOLDER TRANSFORMATION  Q = I + U*(U**T)/B */
/*     MODE    = 1 OR 2   TO SELECT ALGORITHM  H1  OR  H2 . */
/*     LPIVOT IS THE INDEX OF THE PIVOT ELEMENT. */
/*     L1,M   IF L1 <= M   THE TRANSFORMATION WILL BE CONSTRUCTED TO */
/*            ZERO ELEMENTS INDEXED FROM L1 THROUGH M. */
/*            IF L1 > M THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION. */
/*     U(),IUE,UP */
/*            ON ENTRY TO H1 U() STORES THE PIVOT VECTOR. */
/*            IUE IS THE STORAGE INCREMENT BETWEEN ELEMENTS. */
/*            ON EXIT FROM H1 U() AND UP STORE QUANTITIES DEFINING */
/*            THE VECTOR U OF THE HOUSEHOLDER TRANSFORMATION. */
/*            ON ENTRY TO H2 U() AND UP */
/*            SHOULD STORE QUANTITIES PREVIOUSLY COMPUTED BY H1. */
/*            THESE WILL NOT BE MODIFIED BY H2. */
/*     C()    ON ENTRY TO H1 OR H2 C() STORES A MATRIX WHICH WILL BE */
/*            REGARDED AS A SET OF VECTORS TO WHICH THE HOUSEHOLDER */
/*            TRANSFORMATION IS TO BE APPLIED. */
/*            ON EXIT C() STORES THE SET OF TRANSFORMED VECTORS. */
/*     ICE    STORAGE INCREMENT BETWEEN ELEMENTS OF VECTORS IN C(). */
/*     ICV    STORAGE INCREMENT BETWEEN VECTORS IN C(). */
/*     NCV    NUMBER OF VECTORS IN C() TO BE TRANSFORMED. */
/*            IF NCV <= 0 NO OPERATIONS WILL BE DONE ON C(). */
int h12(const int* mode,
        const int* lpivot,
        const int* l1,
        const int* m,
        double* u,
        const int* iue,
        double* up,
        double* c__,
        const int* ice,
        const int* icv,
        const int* ncv);

/*   LDL     LDL' - RANK-ONE - UPDATE */
/*   PURPOSE: */
/*           UPDATES THE LDL' FACTORS OF MATRIX A BY RANK-ONE MATRIX */
/*           SIGMA*Z*Z' */
/*   INPUT ARGUMENTS: (* MEANS PARAMETERS ARE CHANGED DURING EXECUTION) */
/*     N     : ORDER OF THE COEFFICIENT MATRIX A */
/*   * A     : POSITIVE DEFINITE MATRIX OF DIMENSION N; */
/*             ONLY THE LOWER TRIANGLE IS USED AND IS STORED COLUMN BY */
/*             COLUMN AS ONE DIMENSIONAL ARRAY OF DIMENSION N*(N+1)/2. */
/*   * Z     : VECTOR OF DIMENSION N OF UPDATING ELEMENTS */
/*     SIGMA : SCALAR FACTOR BY WHICH THE MODIFYING DYADE Z*Z' IS */
/*             MULTIPLIED */
/*   OUTPUT ARGUMENTS: */
/*     A     : UPDATED LDL' FACTORS */
/*   WORKING ARRAY: */
/*     W     : VECTOR OP DIMENSION N (USED ONLY IF SIGMA .LT. ZERO) */
/*   METHOD: */
/*     THAT OF FLETCHER AND POWELL AS DESCRIBED IN : */
/*     FLETCHER,R.,(1974) ON THE MODIFICATION OF LDL' FACTORIZATION. */
/*     POWELL,M.J.D.      MATH.COMPUTATION 28, 1067-1078. */
/*   IMPLEMENTED BY: */
/*     KRAFT,D., DFVLR - INSTITUT FUER DYNAMIK DER FLUGSYSTEME */
/*               D-8031  OBERPFAFFENHOFEN */
/*   STATUS: 15. JANUARY 1980 */
/*   SUBROUTINES REQUIRED: NONE */
int ldl(const int* n, double* a, double* z__, const double* sigma, double* w);

/*   LINMIN  LINESEARCH WITHOUT DERIVATIVES */
/*   PURPOSE: */
/*  TO FIND THE ARGUMENT LINMIN WHERE THE FUNCTION F TAKES IT'S MINIMUM */
/*  ON THE INTERVAL AX, BX. */
/*  COMBINATION OF GOLDEN SECTION AND SUCCESSIVE QUADRATIC INTERPOLATION. */
/*   INPUT ARGUMENTS: (* MEANS PARAMETERS ARE CHANGED DURING EXECUTION) */
/* *MODE   SEE OUTPUT ARGUMENTS */
/*  AX     LEFT ENDPOINT OF INITIAL INTERVAL */
/*  BX     RIGHT ENDPOINT OF INITIAL INTERVAL */
/*  F      FUNCTION VALUE AT LINMIN WHICH IS TO BE BROUGHT IN BY */
/*         REVERSE COMMUNICATION CONTROLLED BY MODE */
/*  TOL    DESIRED LENGTH OF INTERVAL OF UNCERTAINTY OF FINAL RESULT */
/*   OUTPUT ARGUMENTS: */
/*  LINMIN ABSCISSA APPROXIMATING THE POINT WHERE F ATTAINS A MINIMUM */
/*  MODE   CONTROLS REVERSE COMMUNICATION */
/*         MUST BE SET TO 0 INITIALLY, RETURNS WITH INTERMEDIATE */
/*         VALUES 1 AND 2 WHICH MUST NOT BE CHANGED BY THE USER, */
/*         ENDS WITH CONVERGENCE WITH VALUE 3. */
/*   WORKING ARRAY: */
/*  NONE */
/*   METHOD: */
/*  THIS FUNCTION SUBPROGRAM IS A SLIGHTLY MODIFIED VERSION OF THE */
/*  ALGOL 60 PROCEDURE LOCALMIN GIVEN IN */
/*  R.P. BRENT: ALGORITHMS FOR MINIMIZATION WITHOUT DERIVATIVES, */
/*              PRENTICE-HALL (1973). */
/*   IMPLEMENTED BY: */
/*     KRAFT, D., DFVLR - INSTITUT FUER DYNAMIK DER FLUGSYSTEME */
/*                D-8031  OBERPFAFFENHOFEN */
/*   STATUS: 31. AUGUST  1984 */
/*   SUBROUTINES REQUIRED: NONE */
/*  EPS = SQUARE - ROOT OF MACHINE PRECISION */
/*  C = GOLDEN SECTION RATIO = (3-SQRT(5))/2 */
double linmin(int* mode,
              const double* ax,
              const double* bx,
              const double* f,
              const double* tol);

/*     CONSTANT TIMES A VECTOR PLUS A VECTOR. */
/*     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */
int daxpy_sl(const int* n,
             const double* da,
             const double* dx,
             const int* incx,
             double* dy,
             const int* incy);

/*     COPIES A VECTOR, X, TO A VECTOR, Y. */
/*     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */
int dcopy(const int* n,
          const double* dx,
          const int* incx,
          double* dy,
          const int* incy);

/*     FORMS THE DOT PRODUCT OF TWO VECTORS. */
/*     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */
double ddot_sl(const int* n,
               const double* dx,
               const int* incx,
               const double* dy,
               const int* incy);

/*      DNRM1 - COMPUTES THE I-NORM OF A VECTOR */
/*              BETWEEN THE ITH AND THE JTH ELEMENTS */
/*      INPUT - */
/*      N       LENGTH OF VECTOR */
/*      X       VECTOR OF LENGTH N */
/*      I       INITIAL ELEMENT OF VECTOR TO BE USED */
/*      J       FINAL ELEMENT TO USE */
/*      OUTPUT - */
/*      DNRM1   NORM */
double dnrm1(const int* n, const double* x, const int* i__, const int* j);

/* Function Body */
/*     EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX() WITH STORAGE */
/*     INCREMENT INCX . */
/*     IF    N .LE. 0 RETURN WITH RESULT = 0. */
/*     IF N .GE. 1 THEN INCX MUST BE .GE. 1 */
/*           C.L.LAWSON, 1978 JAN 08 */
/*     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE */
/*     HOPEFULLY APPLICABLE TO ALL MACHINES. */
/*         CUTLO = MAXIMUM OF  SQRT(U/EPS)   OVER ALL KNOWN MACHINES. */
/*         CUTHI = MINIMUM OF  SQRT(V)       OVER ALL KNOWN MACHINES. */
/*     WHERE */
/*         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1. */
/*         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT) */
/*         V   = LARGEST  NO.            (OVERFLOW  LIMIT) */
/*     BRIEF OUTLINE OF ALGORITHM.. */
/*     PHASE 1    SCANS ZERO COMPONENTS. */
/*     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO */
/*     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO */
/*     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M */
/*     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX. */
/*     VALUES FOR CUTLO AND CUTHI.. */
/*     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER */
/*     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS.. */
/*     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE */
/*                   UNIVAC AND DEC AT 2**(-103) */
/*                   THUS CUTLO = 2**(-51) = 4.44089E-16 */
/*     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC. */
/*                   THUS CUTHI = 2**(63.5) = 1.30438E19 */
/*     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC. */
/*                   THUS CUTLO = 2**(-33.5) = 8.23181D-11 */
/*     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19 */
/*     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 / */
/*     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 / */
double dnrm2(const int* n, const double* dx, const int* incx);

/*     APPLIES A PLANE ROTATION. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */
int dsrot(const int* n,
          double* dx,
          const int* incx,
          double* dy,
          const int* incy,
          const double* c__,
          const double* s);

/*     CONSTRUCT GIVENS PLANE ROTATION. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */
/*                    MODIFIED 9/27/86. */
int dsrotg(double* da, double* db, double* c__, double* s);

/*     SCALES A VECTOR BY A CONSTANT. */
/*     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */
int dscal_sl(const int* n, const double* da, double* dx, const int* incx);

int bound(const int* n, double* x, const double* xl, const double* xu);

#ifdef __cplusplus
} /* extern "C */
#endif
