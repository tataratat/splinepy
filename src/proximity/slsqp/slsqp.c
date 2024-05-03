#include "splinepy/proximity/slsqp/slsqp.h"

#ifdef __cplusplus
#include <cmath>
using std::fabs, std::sqrt, std::copysign;
#else

/* fabs, sqrt, copysign */
#include <math.h>

#endif

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

/* Table of constant values */

const static int c__0 = 0;
const static int c__1 = 1;
const static int c__2 = 2;

/* Subroutine */ int slsqp(int* m,
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
                           double* w,
                           int* l_w__,
                           int* jw,
                           int* l_jw__,
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
                           int* n3) {
  /* System generated locals */
  int a_dim1, a_offset, i__1, i__2;

  /* Local variables */
  int il, im, ir, is, iu, iv, iw, ix, mineq;

  /* Parameter adjustments */
  --c__;
  a_dim1 = *la;
  a_offset = 1 + a_dim1;
  a -= a_offset;
  --g;
  --xu;
  --xl;
  --x;
  --w;
  --jw;

  /* Function Body */
  *n1 = *n + 1;
  mineq = *m - *meq + *n1 + *n1;
  il = (*n1 * 3 + *m) * (*n1 + 1) + (*n1 - *meq + 1) * (mineq + 2)
       + (mineq << 1) + (*n1 + mineq) * (*n1 - *meq) + (*meq << 1)
       + *n1 * *n / 2 + (*m << 1) + *n * 3 + (*n1 << 2) + 1;
  /* Computing MAX */
  i__1 = mineq, i__2 = *n1 - *meq;
  im = max(i__1, i__2);
  if (*l_w__ < il || *l_jw__ < im) {
    *mode = max(10, il) * 1000;
    *mode += max(10, im);
    return 0;
  }
  /*   PREPARE DATA FOR CALLING SQPBDY  -  INITIAL ADDRESSES IN W */
  im = 1;
  il = im + max(1, *m);
  il = im + *la;
  ix = il + *n1 * *n / 2 + 1;
  ir = ix + *n;
  is = ir + *n + *n + max(1, *m);
  is = ir + *n + *n + *la;
  iu = is + *n1;
  iv = iu + *n1;
  iw = iv + *n1;
  slsqpb(m,
         meq,
         la,
         n,
         &x[1],
         &xl[1],
         &xu[1],
         f,
         &c__[1],
         &g[1],
         &a[a_offset],
         acc,
         iter,
         mode,
         &w[ir],
         &w[il],
         &w[ix],
         &w[im],
         &w[is],
         &w[iu],
         &w[iv],
         &w[iw],
         &jw[1],
         alpha,
         f0,
         gs,
         h1,
         h2,
         h3,
         h4,
         t,
         t0,
         tol,
         iexact,
         incons,
         ireset,
         itermx,
         line,
         n1,
         n2,
         n3);
  return 0;
} /* slsqp */

/* Subroutine */ int slsqpb(int* m,
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
                            int* n3) {
  /* Initialized data */

  const double zero = 0.;
  const double one = 1.;
  const double alfmin = .1;
  const double hun = 100.;
  const double ten = 10.;
  const double two = 2.;

  /* System generated locals */
  int a_dim1, a_offset, i__1, i__2;
  double d__1, d__2;

  /* Local variables */
  int i__, j, k;
  int badlin;

  /* Parameter adjustments */
  --mu;
  --c__;
  --v;
  --u;
  --s;
  --x0;
  --l;
  --r__;
  a_dim1 = *la;
  a_offset = 1 + a_dim1;
  a -= a_offset;
  --g;
  --xu;
  --xl;
  --x;
  --w;
  --iw;

  /* Function Body */
  /*     The badlin flag keeps track whether the SQP problem on the current */
  /*     iteration was inconsistent or not. */
  badlin = 0;
  if (*mode < 0) {
    goto L260;
  } else if (*mode == 0) {
    goto L100;
  } else {
    goto L220;
  }
L100:
  *itermx = *iter;
  if (*acc >= zero) {
    *iexact = 0;
  } else {
    *iexact = 1;
  }
  *acc = fabs(*acc);
  *tol = ten * *acc;
  *iter = 0;
  *ireset = 0;
  *n1 = *n + 1;
  *n2 = *n1 * *n / 2;
  *n3 = *n2 + 1;
  s[1] = zero;
  mu[1] = zero;
  dcopy(n, &s[1], &c__0, &s[1], &c__1);
  dcopy(m, &mu[1], &c__0, &mu[1], &c__1);
/*   RESET BFGS MATRIX */
L110:
  ++(*ireset);
  if (*ireset > 5) {
    goto L255;
  }
  l[1] = zero;
  dcopy(n2, &l[1], &c__0, &l[1], &c__1);
  j = 1;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    l[j] = one;
    j = j + *n1 - i__;
    /* L120: */
  }
/*   MAIN ITERATION : SEARCH DIRECTION, STEPLENGTH, LDL'-UPDATE */
L130:
  ++(*iter);
  *mode = 9;
  if (*iter > *itermx) {
    goto L330;
  }
  /*   SEARCH DIRECTION AS SOLUTION OF QP - SUBPROBLEM */
  dcopy(n, &xl[1], &c__1, &u[1], &c__1);
  dcopy(n, &xu[1], &c__1, &v[1], &c__1);
  d__1 = -one;
  daxpy_sl(n, &d__1, &x[1], &c__1, &u[1], &c__1);
  d__1 = -one;
  daxpy_sl(n, &d__1, &x[1], &c__1, &v[1], &c__1);
  *h4 = one;
  lsq(m,
      meq,
      n,
      n3,
      la,
      &l[1],
      &g[1],
      &a[a_offset],
      &c__[1],
      &u[1],
      &v[1],
      &s[1],
      &r__[1],
      &w[1],
      &iw[1],
      mode);
  /*   AUGMENTED PROBLEM FOR INCONSISTENT LINEARIZATION */

  /*   If it turns out that the original SQP problem is inconsistent, */
  /*   disallow termination with convergence on this iteration, */
  /*   even if the augmented problem was solved. */
  badlin = 0;
  if (*mode == 6) {
    if (*n == *meq) {
      *mode = 4;
    }
  }
  if (*mode == 4) {
    badlin = 1;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
      if (j <= *meq) {
        a[j + *n1 * a_dim1] = -c__[j];
      } else {
        /* Computing MAX */
        d__1 = -c__[j];
        a[j + *n1 * a_dim1] = max(d__1, zero);
      }
      /* L140: */
    }
    s[1] = zero;
    dcopy(n, &s[1], &c__0, &s[1], &c__1);
    *h3 = zero;
    g[*n1] = zero;
    l[*n3] = hun;
    s[*n1] = one;
    u[*n1] = zero;
    v[*n1] = one;
    *incons = 0;
  L150:
    lsq(m,
        meq,
        n1,
        n3,
        la,
        &l[1],
        &g[1],
        &a[a_offset],
        &c__[1],
        &u[1],
        &v[1],
        &s[1],
        &r__[1],
        &w[1],
        &iw[1],
        mode);
    *h4 = one - s[*n1];
    if (*mode == 4) {
      l[*n3] = ten * l[*n3];
      ++(*incons);
      if (*incons > 5) {
        goto L330;
      }
      goto L150;
    } else if (*mode != 1) {
      goto L330;
    }
  } else if (*mode != 1) {
    goto L330;
  }
  /*   UPDATE MULTIPLIERS FOR L1-TEST */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    v[i__] = g[i__] - ddot_sl(m, &a[i__ * a_dim1 + 1], &c__1, &r__[1], &c__1);
    /* L160: */
  }
  *f0 = *f;
  dcopy(n, &x[1], &c__1, &x0[1], &c__1);
  *gs = ddot_sl(n, &g[1], &c__1, &s[1], &c__1);
  *h1 = fabs(*gs);
  *h2 = zero;
  i__1 = *m;
  for (j = 1; j <= i__1; ++j) {
    if (j <= *meq) {
      *h3 = c__[j];
    } else {
      *h3 = zero;
    }
    /* Computing MAX */
    d__1 = -c__[j];
    *h2 += max(d__1, *h3);
    *h3 = (d__1 = r__[j], fabs(d__1));
    /* Computing MAX */
    d__1 = *h3, d__2 = (mu[j] + *h3) / two;
    mu[j] = max(d__1, d__2);
    *h1 += *h3 * (d__1 = c__[j], fabs(d__1));
    /* L170: */
  }
  /*   CHECK CONVERGENCE */
  *mode = 0;
  if (*h1 < *acc && *h2 < *acc && !badlin && *f == *f) {
    goto L330;
  }
  *h1 = zero;
  i__1 = *m;
  for (j = 1; j <= i__1; ++j) {
    if (j <= *meq) {
      *h3 = c__[j];
    } else {
      *h3 = zero;
    }
    /* Computing MAX */
    d__1 = -c__[j];
    *h1 += mu[j] * max(d__1, *h3);
    /* L180: */
  }
  *t0 = *f + *h1;
  *h3 = *gs - *h1 * *h4;
  *mode = 8;
  if (*h3 >= zero) {
    goto L110;
  }
  /*   LINE SEARCH WITH AN L1-TESTFUNCTION */
  *line = 0;
  *alpha = one;
  if (*iexact == 1) {
    goto L210;
  }
/*   INEXACT LINESEARCH */
L190:
  ++(*line);
  *h3 = *alpha * *h3;
  dscal_sl(n, alpha, &s[1], &c__1);
  dcopy(n, &x0[1], &c__1, &x[1], &c__1);
  daxpy_sl(n, &one, &s[1], &c__1, &x[1], &c__1);
  *mode = 1;
  goto L330;
L200:
  if (*h1 <= *h3 / ten || *line > 10) {
    goto L240;
  }
  /* Computing MAX */
  d__1 = *h3 / (two * (*h3 - *h1));
  *alpha = max(d__1, alfmin);
  goto L190;
/*   EXACT LINESEARCH */
L210:
  if (*line != 3) {
    *alpha = linmin(line, &alfmin, &one, t, tol);
    dcopy(n, &x0[1], &c__1, &x[1], &c__1);
    daxpy_sl(n, alpha, &s[1], &c__1, &x[1], &c__1);
    *mode = 1;
    goto L330;
  }
  dscal_sl(n, alpha, &s[1], &c__1);
  goto L240;
/*   CALL FUNCTIONS AT CURRENT X */
L220:
  *t = *f;
  i__1 = *m;
  for (j = 1; j <= i__1; ++j) {
    if (j <= *meq) {
      *h1 = c__[j];
    } else {
      *h1 = zero;
    }
    /* Computing MAX */
    d__1 = -c__[j];
    *t += mu[j] * max(d__1, *h1);
    /* L230: */
  }
  *h1 = *t - *t0;
  switch (*iexact + 1) {
  case 1:
    goto L200;
  case 2:
    goto L210;
  }
/*   CHECK CONVERGENCE */
L240:
  *h3 = zero;
  i__1 = *m;
  for (j = 1; j <= i__1; ++j) {
    if (j <= *meq) {
      *h1 = c__[j];
    } else {
      *h1 = zero;
    }
    /* Computing MAX */
    d__1 = -c__[j];
    *h3 += max(d__1, *h1);
    /* L250: */
  }
  if (((d__1 = *f - *f0, fabs(d__1)) < *acc || dnrm2(n, &s[1], &c__1) < *acc)
      && *h3 < *acc && !badlin && *f == *f) {
    *mode = 0;
  } else {
    *mode = -1;
  }
  goto L330;
/*   CHECK relaxed CONVERGENCE in case of positive directional derivative */
L255:
  *h3 = zero;
  i__1 = *m;
  for (j = 1; j <= i__1; ++j) {
    if (j <= *meq) {
      *h1 = c__[j];
    } else {
      *h1 = zero;
    }
    /* Computing MAX */
    d__1 = -c__[j];
    *h3 += max(d__1, *h1);
    /* L256: */
  }
  if (((d__1 = *f - *f0, fabs(d__1)) < *tol || dnrm2(n, &s[1], &c__1) < *tol)
      && *h3 < *tol && !badlin && *f == *f) {
    *mode = 0;
  } else {
    *mode = 8;
  }
  goto L330;
/*   CALL JACOBIAN AT CURRENT X */
/*   UPDATE CHOLESKY-FACTORS OF HESSIAN MATRIX BY MODIFIED BFGS FORMULA */
L260:
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    u[i__] = g[i__] - ddot_sl(m, &a[i__ * a_dim1 + 1], &c__1, &r__[1], &c__1)
             - v[i__];
    /* L270: */
  }
  /*   L'*S */
  k = 0;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    *h1 = zero;
    ++k;
    i__2 = *n;
    for (j = i__ + 1; j <= i__2; ++j) {
      ++k;
      *h1 += l[k] * s[j];
      /* L280: */
    }
    v[i__] = s[i__] + *h1;
    /* L290: */
  }
  /*   D*L'*S */
  k = 1;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    v[i__] = l[k] * v[i__];
    k = k + *n1 - i__;
    /* L300: */
  }
  /*   L*D*L'*S */
  for (i__ = *n; i__ >= 1; --i__) {
    *h1 = zero;
    k = i__;
    i__1 = i__ - 1;
    for (j = 1; j <= i__1; ++j) {
      *h1 += l[k] * v[j];
      k = k + *n - j;
      /* L310: */
    }
    v[i__] += *h1;
    /* L320: */
  }
  *h1 = ddot_sl(n, &s[1], &c__1, &u[1], &c__1);
  *h2 = ddot_sl(n, &s[1], &c__1, &v[1], &c__1);
  *h3 = *h2 * .2;
  if (*h1 < *h3) {
    *h4 = (*h2 - *h3) / (*h2 - *h1);
    *h1 = *h3;
    dscal_sl(n, h4, &u[1], &c__1);
    d__1 = one - *h4;
    daxpy_sl(n, &d__1, &v[1], &c__1, &u[1], &c__1);
  }
  if (*h1 == 0. || *h2 == 0.) {
    /*         Singular update: reset hessian. */
    goto L110;
  }
  d__1 = one / *h1;
  ldl(n, &l[1], &u[1], &d__1, &v[1]);
  d__1 = -one / *h2;
  ldl(n, &l[1], &v[1], &d__1, &u[1]);
  /*   END OF MAIN ITERATION */
  goto L130;
/*   END OF SLSQPB */
L330:
  return 0;
} /* slsqpb */

/* Subroutine */ int lsq(const int* m,
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
                         int* mode) {
  /* Initialized data */

  const double zero = 0.;
  const double one = 1.;

  /* System generated locals */
  int a_dim1, a_offset, i__1, i__2;
  double d__1;

  /* Local variables */
  int i__, j, i1, i2, i3, i4, m1, n1, n2, n3, ic, id, ie, if__, ig, ih, il, ip,
      iw;
  double diag;
  int mineq;
  double xnorm;
  int nancnt;

  /* Parameter adjustments */
  --y;
  --x;
  --xu;
  --xl;
  --g;
  --l;
  --b;
  a_dim1 = *la;
  a_offset = 1 + a_dim1;
  a -= a_offset;
  --w;
  --jw;

  /* Function Body */
  n1 = *n + 1;
  mineq = *m - *meq;
  m1 = mineq + *n + *n;
  /*  determine whether to solve problem */
  /*  with inconsistent linerarization (n2=1) */
  /*  or not (n2=0) */
  n2 = n1 * *n / 2 + 1;
  if (n2 == *nl) {
    n2 = 0;
  } else {
    n2 = 1;
  }
  n3 = *n - n2;
  /*  RECOVER MATRIX E AND VECTOR F FROM L AND G */
  i2 = 1;
  i3 = 1;
  i4 = 1;
  ie = 1;
  if__ = *n * *n + 1;
  i__1 = n3;
  for (i__ = 1; i__ <= i__1; ++i__) {
    i1 = n1 - i__;
    diag = sqrt(l[i2]);
    w[i3] = zero;
    dcopy(&i1, &w[i3], &c__0, &w[i3], &c__1);
    i__2 = i1 - n2;
    dcopy(&i__2, &l[i2], &c__1, &w[i3], n);
    i__2 = i1 - n2;
    dscal_sl(&i__2, &diag, &w[i3], n);
    w[i3] = diag;
    i__2 = i__ - 1;
    w[if__ - 1 + i__] =
        (g[i__] - ddot_sl(&i__2, &w[i4], &c__1, &w[if__], &c__1)) / diag;
    i2 = i2 + i1 - n2;
    i3 += n1;
    i4 += *n;
    /* L10: */
  }
  if (n2 == 1) {
    w[i3] = l[*nl];
    w[i4] = zero;
    dcopy(&n3, &w[i4], &c__0, &w[i4], &c__1);
    w[if__ - 1 + *n] = zero;
  }
  d__1 = -one;
  dscal_sl(n, &d__1, &w[if__], &c__1);
  ic = if__ + *n;
  id = ic + *meq * *n;
  if (*meq > 0) {
    /*  RECOVER MATRIX C FROM UPPER PART OF A */
    i__1 = *meq;
    for (i__ = 1; i__ <= i__1; ++i__) {
      dcopy(n, &a[i__ + a_dim1], la, &w[ic - 1 + i__], meq);
      /* L20: */
    }
    /*  RECOVER VECTOR D FROM UPPER PART OF B */
    dcopy(meq, &b[1], &c__1, &w[id], &c__1);
    d__1 = -one;
    dscal_sl(meq, &d__1, &w[id], &c__1);
  }
  ig = id + *meq;
  /*  RECOVER MATRIX G FROM LOWER PART OF A */
  /*  The matrix G(mineq+2*n,m1) is stored at w(ig) */
  /*  Not all rows will be filled if some of the upper/lower */
  /*  bounds are unbounded. */
  if (mineq > 0) {
    i__1 = mineq;
    for (i__ = 1; i__ <= i__1; ++i__) {
      dcopy(n, &a[*meq + i__ + a_dim1], la, &w[ig - 1 + i__], &m1);
      /* L30: */
    }
  }
  ih = ig + m1 * *n;
  iw = ih + mineq + (*n << 1);
  if (mineq > 0) {
    /*  RECOVER H FROM LOWER PART OF B */
    /*  The vector H(mineq+2*n) is stored at w(ih) */
    dcopy(&mineq, &b[*meq + 1], &c__1, &w[ih], &c__1);
    d__1 = -one;
    dscal_sl(&mineq, &d__1, &w[ih], &c__1);
  }
  /*  AUGMENT MATRIX G BY +I AND -I, AND, */
  /*  AUGMENT VECTOR H BY XL AND XU */
  /*  NaN value indicates no bound */
  ip = ig + mineq;
  il = ih + mineq;
  nancnt = 0;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    if (xl[i__] == xl[i__]) {
      w[il] = xl[i__];
      i__2 = *n;
      for (j = 1; j <= i__2; ++j) {
        w[ip + m1 * (j - 1)] = 0.;
        /* L41: */
      }
      w[ip + m1 * (i__ - 1)] = 1.;
      ++ip;
      ++il;
    } else {
      ++nancnt;
    }
    /* L40: */
  }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    if (xu[i__] == xu[i__]) {
      w[il] = -xu[i__];
      i__2 = *n;
      for (j = 1; j <= i__2; ++j) {
        w[ip + m1 * (j - 1)] = 0.;
        /* L51: */
      }
      w[ip + m1 * (i__ - 1)] = -1.;
      ++ip;
      ++il;
    } else {
      ++nancnt;
    }
    /* L50: */
  }
  i__1 = max(1, *meq);
  i__2 = m1 - nancnt;
  lsei(&w[ic],
       &w[id],
       &w[ie],
       &w[if__],
       &w[ig],
       &w[ih],
       &i__1,
       meq,
       n,
       n,
       &m1,
       &i__2,
       n,
       &x[1],
       &xnorm,
       &w[iw],
       &jw[1],
       mode);
  if (*mode == 1) {
    /*   restore Lagrange multipliers (only for user-defined variables) */
    dcopy(m, &w[iw], &c__1, &y[1], &c__1);
    /*   set rest of the multipliers to nan (they are not used) */
    if (n3 > 0) {
      y[*m + 1] = 0.;
      y[*m + 1] = 0 / y[*m + 1];
      i__1 = *m + n3 + n3;
      for (i__ = *m + 2; i__ <= i__1; ++i__) {
        y[i__] = y[*m + 1];
        /* L60: */
      }
    }
  }
  bound(n, &x[1], &xl[1], &xu[1]);
  /*   END OF SUBROUTINE LSQ */
  return 0;
} /* lsq */

/* Subroutine */ int lsei(double* c__,
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
                          int* mode) {
  /* Initialized data */

  const double epmach = 2.22e-16;
  const double zero = 0.;

  /* System generated locals */
  int c_dim1, c_offset, e_dim1, e_offset, g_dim1, g_offset, i__1, i__2, i__3;
  double d__1;

  /* Local variables */
  int i__, j, k, l;
  double t;
  int ie, if__, ig, iw, mc1;
  int krank;
  double rnorm[1];

  /* Parameter adjustments */
  --d__;
  --f;
  --h__;
  --x;
  g_dim1 = *lg;
  g_offset = 1 + g_dim1;
  g -= g_offset;
  e_dim1 = *le;
  e_offset = 1 + e_dim1;
  e -= e_offset;
  c_dim1 = *lc;
  c_offset = 1 + c_dim1;
  c__ -= c_offset;
  --w;
  --jw;

  /* Function Body */
  *mode = 2;
  if (*mc > *n) {
    goto L75;
  }
  l = *n - *mc;
  mc1 = *mc + 1;
  iw = (l + 1) * (*mg + 2) + (*mg << 1) + *mc;
  ie = iw + *mc + 1;
  if__ = ie + *me * l;
  ig = if__ + *me;
  /*  TRIANGULARIZE C AND APPLY FACTORS TO E AND G */
  i__1 = *mc;
  for (i__ = 1; i__ <= i__1; ++i__) {
    /* Computing MIN */
    i__2 = i__ + 1;
    j = min(i__2, *lc);
    i__2 = i__ + 1;
    i__3 = *mc - i__;
    h12(&c__1,
        &i__,
        &i__2,
        n,
        &c__[i__ + c_dim1],
        lc,
        &w[iw + i__],
        &c__[j + c_dim1],
        lc,
        &c__1,
        &i__3);
    i__2 = i__ + 1;
    h12(&c__2,
        &i__,
        &i__2,
        n,
        &c__[i__ + c_dim1],
        lc,
        &w[iw + i__],
        &e[e_offset],
        le,
        &c__1,
        me);
    /* L10: */
    i__2 = i__ + 1;
    h12(&c__2,
        &i__,
        &i__2,
        n,
        &c__[i__ + c_dim1],
        lc,
        &w[iw + i__],
        &g[g_offset],
        lg,
        &c__1,
        mg);
  }
  /*  SOLVE C*X=D AND MODIFY F */
  *mode = 6;
  i__2 = *mc;
  for (i__ = 1; i__ <= i__2; ++i__) {
    if ((d__1 = c__[i__ + i__ * c_dim1], fabs(d__1)) < epmach) {
      goto L75;
    }
    i__1 = i__ - 1;
    x[i__] = (d__[i__] - ddot_sl(&i__1, &c__[i__ + c_dim1], lc, &x[1], &c__1))
             / c__[i__ + i__ * c_dim1];
    /* L15: */
  }
  *mode = 1;
  w[mc1] = zero;
  i__2 = *mg - *mc;
  dcopy(&i__2, &w[mc1], &c__0, &w[mc1], &c__1);
  if (*mc == *n) {
    goto L50;
  }
  i__2 = *me;
  for (i__ = 1; i__ <= i__2; ++i__) {
    /* L20: */
    w[if__ - 1 + i__] =
        f[i__] - ddot_sl(mc, &e[i__ + e_dim1], le, &x[1], &c__1);
  }
  /*  STORE TRANSFORMED E & G */
  i__2 = *me;
  for (i__ = 1; i__ <= i__2; ++i__) {
    /* L25: */
    dcopy(&l, &e[i__ + mc1 * e_dim1], le, &w[ie - 1 + i__], me);
  }
  i__2 = *mg;
  for (i__ = 1; i__ <= i__2; ++i__) {
    /* L30: */
    dcopy(&l, &g[i__ + mc1 * g_dim1], lg, &w[ig - 1 + i__], mg);
  }
  if (*mg > 0) {
    goto L40;
  }
  /*  SOLVE LS WITHOUT INEQUALITY CONSTRAINTS */
  *mode = 7;
  k = max(*le, *n);
  t = sqrt(epmach);
  hfti(&w[ie],
       me,
       me,
       &l,
       &w[if__],
       &k,
       &c__1,
       &t,
       &krank,
       rnorm,
       &w[1],
       &w[l + 1],
       &jw[1]);
  /*  HFTI IS MORE GENERIC, BUT WE ONLY CALL IT WITH NB=1, SO RETRIEVE THE */
  /*  SINGLE VALUE WE NEED FROM RNORM HERE */
  *xnrm = rnorm[0];
  dcopy(&l, &w[if__], &c__1, &x[mc1], &c__1);
  if (krank != l) {
    goto L75;
  }
  *mode = 1;
  goto L50;
/*  MODIFY H AND SOLVE INEQUALITY CONSTRAINED LS PROBLEM */
L40:
  i__2 = *mg;
  for (i__ = 1; i__ <= i__2; ++i__) {
    /* L45: */
    h__[i__] -= ddot_sl(mc, &g[i__ + g_dim1], lg, &x[1], &c__1);
  }
  lsi(&w[ie],
      &w[if__],
      &w[ig],
      &h__[1],
      me,
      me,
      mg,
      mg,
      &l,
      &x[mc1],
      xnrm,
      &w[mc1],
      &jw[1],
      mode);
  if (*mc == 0) {
    goto L75;
  }
  t = dnrm2(mc, &x[1], &c__1);
  *xnrm = sqrt(*xnrm * *xnrm + t * t);
  if (*mode != 1) {
    goto L75;
  }
/*  SOLUTION OF ORIGINAL PROBLEM AND LAGRANGE MULTIPLIERS */
L50:
  i__2 = *me;
  for (i__ = 1; i__ <= i__2; ++i__) {
    /* L55: */
    f[i__] = ddot_sl(n, &e[i__ + e_dim1], le, &x[1], &c__1) - f[i__];
  }
  i__2 = *mc;
  for (i__ = 1; i__ <= i__2; ++i__) {
    /* L60: */
    d__[i__] = ddot_sl(me, &e[i__ * e_dim1 + 1], &c__1, &f[1], &c__1)
               - ddot_sl(mg, &g[i__ * g_dim1 + 1], &c__1, &w[mc1], &c__1);
  }
  for (i__ = *mc; i__ >= 1; --i__) {
    /* L65: */
    i__2 = i__ + 1;
    h12(&c__2,
        &i__,
        &i__2,
        n,
        &c__[i__ + c_dim1],
        lc,
        &w[iw + i__],
        &x[1],
        &c__1,
        &c__1,
        &c__1);
  }
  for (i__ = *mc; i__ >= 1; --i__) {
    /* Computing MIN */
    i__2 = i__ + 1;
    j = min(i__2, *lc);
    i__2 = *mc - i__;
    w[i__] =
        (d__[i__] - ddot_sl(&i__2, &c__[j + i__ * c_dim1], &c__1, &w[j], &c__1))
        / c__[i__ + i__ * c_dim1];
    /* L70: */
  }
/*  END OF SUBROUTINE LSEI */
L75:
  return 0;
} /* lsei */

/* Subroutine */ int lsi(double* e,
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
                         int* mode) {
  /* Initialized data */

  const double epmach = 2.22e-16;
  const double one = 1.;

  /* System generated locals */
  int e_dim1, e_offset, g_dim1, g_offset, i__1, i__2, i__3;
  double d__1;

  /* Local variables */
  int i__, j;
  double t;

  /* Parameter adjustments */
  --f;
  --jw;
  --h__;
  --x;
  g_dim1 = *lg;
  g_offset = 1 + g_dim1;
  g -= g_offset;
  e_dim1 = *le;
  e_offset = 1 + e_dim1;
  e -= e_offset;
  --w;

  /* Function Body */
  /*  QR-FACTORS OF E AND APPLICATION TO F */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    /* Computing MIN */
    i__2 = i__ + 1;
    j = min(i__2, *n);
    i__2 = i__ + 1;
    i__3 = *n - i__;
    h12(&c__1,
        &i__,
        &i__2,
        me,
        &e[i__ * e_dim1 + 1],
        &c__1,
        &t,
        &e[j * e_dim1 + 1],
        &c__1,
        le,
        &i__3);
    /* L10: */
    i__2 = i__ + 1;
    h12(&c__2,
        &i__,
        &i__2,
        me,
        &e[i__ * e_dim1 + 1],
        &c__1,
        &t,
        &f[1],
        &c__1,
        &c__1,
        &c__1);
  }
  /*  TRANSFORM G AND H TO GET LEAST DISTANCE PROBLEM */
  *mode = 5;
  i__2 = *mg;
  for (i__ = 1; i__ <= i__2; ++i__) {
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      if (!((d__1 = e[j + j * e_dim1], fabs(d__1)) >= epmach)) {
        goto L50;
      }
      /* L20: */
      i__3 = j - 1;
      g[i__ + j * g_dim1] =
          (g[i__ + j * g_dim1]
           - ddot_sl(&i__3, &g[i__ + g_dim1], lg, &e[j * e_dim1 + 1], &c__1))
          / e[j + j * e_dim1];
    }
    /* L30: */
    h__[i__] -= ddot_sl(n, &g[i__ + g_dim1], lg, &f[1], &c__1);
  }
  /*  SOLVE LEAST DISTANCE PROBLEM */
  ldp(&g[g_offset], lg, mg, n, &h__[1], &x[1], xnorm, &w[1], &jw[1], mode);
  if (*mode != 1) {
    goto L50;
  }
  /*  SOLUTION OF ORIGINAL PROBLEM */
  daxpy_sl(n, &one, &f[1], &c__1, &x[1], &c__1);
  for (i__ = *n; i__ >= 1; --i__) {
    /* Computing MIN */
    i__2 = i__ + 1;
    j = min(i__2, *n);
    /* L40: */
    i__2 = *n - i__;
    x[i__] = (x[i__] - ddot_sl(&i__2, &e[i__ + j * e_dim1], le, &x[j], &c__1))
             / e[i__ + i__ * e_dim1];
  }
  /* Computing MIN */
  i__2 = *n + 1;
  j = min(i__2, *me);
  i__2 = *me - *n;
  t = dnrm2(&i__2, &f[j], &c__1);
  *xnorm = sqrt(*xnorm * *xnorm + t * t);
/*  END OF SUBROUTINE LSI */
L50:
  return 0;
} /* lsi */

/* Subroutine */ int ldp(const double* g,
                         const int* mg,
                         const int* m,
                         const int* n,
                         const double* h__,
                         double* x,
                         double* xnorm,
                         double* w,
                         int* index,
                         int* mode) {
  /* Initialized data */

  const double zero = 0.;
  const double one = 1.;

  /* System generated locals */
  int g_dim1, g_offset, i__1, i__2;
  double d__1;

  /* Local variables */
  int i__, j, n1, if__, iw, iy, iz;
  double fac;
  double rnorm;
  int iwdual;

  /* Parameter adjustments */
  --index;
  --h__;
  --x;
  g_dim1 = *mg;
  g_offset = 1 + g_dim1;
  g -= g_offset;
  --w;

  /* Function Body */
  *mode = 2;
  if (*n <= 0) {
    goto L50;
  }
  /*  STATE DUAL PROBLEM */
  *mode = 1;
  x[1] = zero;
  dcopy(n, &x[1], &c__0, &x[1], &c__1);
  *xnorm = zero;
  if (*m == 0) {
    goto L50;
  }
  iw = 0;
  i__1 = *m;
  for (j = 1; j <= i__1; ++j) {
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
      ++iw;
      /* L10: */
      w[iw] = g[j + i__ * g_dim1];
    }
    ++iw;
    /* L20: */
    w[iw] = h__[j];
  }
  if__ = iw + 1;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    ++iw;
    /* L30: */
    w[iw] = zero;
  }
  w[iw + 1] = one;
  n1 = *n + 1;
  iz = iw + 2;
  iy = iz + n1;
  iwdual = iy + *m;
  /*  SOLVE DUAL PROBLEM */
  nnls(&w[1],
       &n1,
       &n1,
       m,
       &w[if__],
       &w[iy],
       &rnorm,
       &w[iwdual],
       &w[iz],
       &index[1],
       mode);
  if (*mode != 1) {
    goto L50;
  }
  *mode = 4;
  if (rnorm <= zero) {
    goto L50;
  }
  /*  COMPUTE SOLUTION OF PRIMAL PROBLEM */
  fac = one - ddot_sl(m, &h__[1], &c__1, &w[iy], &c__1);
  d__1 = one + fac;
  if (!(d__1 - one > zero)) {
    goto L50;
  }
  *mode = 1;
  fac = one / fac;
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    /* L40: */
    x[j] = fac * ddot_sl(m, &g[j * g_dim1 + 1], &c__1, &w[iy], &c__1);
  }
  *xnorm = dnrm2(n, &x[1], &c__1);
  /*  COMPUTE LAGRANGE MULTIPLIERS FOR PRIMAL PROBLEM */
  w[1] = zero;
  dcopy(m, &w[1], &c__0, &w[1], &c__1);
  daxpy_sl(m, &fac, &w[iy], &c__1, &w[1], &c__1);
/*  END OF SUBROUTINE LDP */
L50:
  return 0;
} /* ldp */

/* Subroutine */ int nnls(double* a,
                          const int* mda,
                          const int* m,
                          const int* n,
                          double* b,
                          double* x,
                          double* rnorm,
                          double* w,
                          double* z__,
                          int* index,
                          int* mode) {
  /* Initialized data */

  const double zero = 0.;
  const double one = 1.;
  const double factor = .01;

  /* System generated locals */
  int a_dim1, a_offset, i__1, i__2;
  double d__1;

  /* Local variables */
  double c__;
  int i__, j, k, l;
  double s, t;
  int ii, jj, ip, iz, jz;
  double up;
  int iz1, iz2, npp1, iter;
  double wmax, alpha, asave;
  int itmax, izmax, nsetp;
  double unorm;

  /* Parameter adjustments */
  --z__;
  --b;
  --index;
  --w;
  --x;
  a_dim1 = *mda;
  a_offset = 1 + a_dim1;
  a -= a_offset;

  /* Function Body */
  /*     revised          Dieter Kraft, March 1983 */
  *mode = 2;
  if (*m <= 0 || *n <= 0) {
    goto L290;
  }
  *mode = 1;
  iter = 0;
  itmax = *n * 3;
  /* STEP ONE (INITIALIZE) */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    /* L100: */
    index[i__] = i__;
  }
  iz1 = 1;
  iz2 = *n;
  nsetp = 0;
  npp1 = 1;
  x[1] = zero;
  dcopy(n, &x[1], &c__0, &x[1], &c__1);
/* STEP TWO (COMPUTE DUAL VARIABLES) */
/* .....ENTRY LOOP A */
L110:
  if (iz1 > iz2 || nsetp >= *m) {
    goto L280;
  }
  i__1 = iz2;
  for (iz = iz1; iz <= i__1; ++iz) {
    j = index[iz];
    /* L120: */
    i__2 = *m - nsetp;
    w[j] = ddot_sl(&i__2, &a[npp1 + j * a_dim1], &c__1, &b[npp1], &c__1);
  }
/* STEP THREE (TEST DUAL VARIABLES) */
L130:
  wmax = zero;
  i__2 = iz2;
  for (iz = iz1; iz <= i__2; ++iz) {
    j = index[iz];
    if (w[j] <= wmax) {
      goto L140;
    }
    wmax = w[j];
    izmax = iz;
  L140:;
  }
  /* .....EXIT LOOP A */
  if (wmax <= zero) {
    goto L280;
  }
  iz = izmax;
  j = index[iz];
  /* STEP FOUR (TEST INDEX J FOR LINEAR DEPENDENCY) */
  asave = a[npp1 + j * a_dim1];
  i__2 = npp1 + 1;
  h12(&c__1,
      &npp1,
      &i__2,
      m,
      &a[j * a_dim1 + 1],
      &c__1,
      &up,
      &z__[1],
      &c__1,
      &c__1,
      &c__0);
  unorm = dnrm2(&nsetp, &a[j * a_dim1 + 1], &c__1);
  t = factor * (d__1 = a[npp1 + j * a_dim1], fabs(d__1));
  d__1 = unorm + t;
  if (d__1 - unorm <= zero) {
    goto L150;
  }
  dcopy(m, &b[1], &c__1, &z__[1], &c__1);
  i__2 = npp1 + 1;
  h12(&c__2,
      &npp1,
      &i__2,
      m,
      &a[j * a_dim1 + 1],
      &c__1,
      &up,
      &z__[1],
      &c__1,
      &c__1,
      &c__1);
  if (z__[npp1] / a[npp1 + j * a_dim1] > zero) {
    goto L160;
  }
L150:
  a[npp1 + j * a_dim1] = asave;
  w[j] = zero;
  goto L130;
/* STEP FIVE (ADD COLUMN) */
L160:
  dcopy(m, &z__[1], &c__1, &b[1], &c__1);
  index[iz] = index[iz1];
  index[iz1] = j;
  ++iz1;
  nsetp = npp1;
  ++npp1;
  i__2 = iz2;
  for (jz = iz1; jz <= i__2; ++jz) {
    jj = index[jz];
    /* L170: */
    h12(&c__2,
        &nsetp,
        &npp1,
        m,
        &a[j * a_dim1 + 1],
        &c__1,
        &up,
        &a[jj * a_dim1 + 1],
        &c__1,
        mda,
        &c__1);
  }
  k = min(npp1, *mda);
  w[j] = zero;
  i__2 = *m - nsetp;
  dcopy(&i__2, &w[j], &c__0, &a[k + j * a_dim1], &c__1);
/* STEP SIX (SOLVE LEAST SQUARES SUB-PROBLEM) */
/* .....ENTRY LOOP B */
L180:
  for (ip = nsetp; ip >= 1; --ip) {
    if (ip == nsetp) {
      goto L190;
    }
    d__1 = -z__[ip + 1];
    daxpy_sl(&ip, &d__1, &a[jj * a_dim1 + 1], &c__1, &z__[1], &c__1);
  L190:
    jj = index[ip];
    /* L200: */
    z__[ip] /= a[ip + jj * a_dim1];
  }
  ++iter;
  if (iter <= itmax) {
    goto L220;
  }
L210:
  *mode = 3;
  goto L280;
/* STEP SEVEN TO TEN (STEP LENGTH ALGORITHM) */
L220:
  alpha = one;
  jj = 0;
  i__2 = nsetp;
  for (ip = 1; ip <= i__2; ++ip) {
    if (z__[ip] > zero) {
      goto L230;
    }
    l = index[ip];
    t = -x[l] / (z__[ip] - x[l]);
    if (alpha < t) {
      goto L230;
    }
    alpha = t;
    jj = ip;
  L230:;
  }
  i__2 = nsetp;
  for (ip = 1; ip <= i__2; ++ip) {
    l = index[ip];
    /* L240: */
    x[l] = (one - alpha) * x[l] + alpha * z__[ip];
  }
  /* .....EXIT LOOP B */
  if (jj == 0) {
    goto L110;
  }
  /* STEP ELEVEN (DELETE COLUMN) */
  i__ = index[jj];
L250:
  x[i__] = zero;
  ++jj;
  i__2 = nsetp;
  for (j = jj; j <= i__2; ++j) {
    ii = index[j];
    index[j - 1] = ii;
    dsrotg(&a[j - 1 + ii * a_dim1], &a[j + ii * a_dim1], &c__, &s);
    t = a[j - 1 + ii * a_dim1];
    dsrot(n, &a[j - 1 + a_dim1], mda, &a[j + a_dim1], mda, &c__, &s);
    a[j - 1 + ii * a_dim1] = t;
    a[j + ii * a_dim1] = zero;
    /* L260: */
    dsrot(&c__1, &b[j - 1], &c__1, &b[j], &c__1, &c__, &s);
  }
  npp1 = nsetp;
  --nsetp;
  --iz1;
  index[iz1] = i__;
  if (nsetp <= 0) {
    goto L210;
  }
  i__2 = nsetp;
  for (jj = 1; jj <= i__2; ++jj) {
    i__ = index[jj];
    if (x[i__] <= zero) {
      goto L250;
    }
    /* L270: */
  }
  dcopy(m, &b[1], &c__1, &z__[1], &c__1);
  goto L180;
/* STEP TWELVE (SOLUTION) */
L280:
  k = min(npp1, *m);
  i__2 = *m - nsetp;
  *rnorm = dnrm2(&i__2, &b[k], &c__1);
  if (npp1 > *m) {
    w[1] = zero;
    dcopy(n, &w[1], &c__0, &w[1], &c__1);
  }
/* END OF SUBROUTINE NNLS */
L290:
  return 0;
} /* nnls */

/* Subroutine */ int hfti(double* a,
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
                          int* ip) {
  /* Initialized data */

  const double zero = 0.;
  const double factor = .001;

  /* System generated locals */
  int a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
  double d__1;

  /* Local variables */
  int i__, j, k, l;
  int jb, kp1;
  double tmp, hmax;
  int lmax, ldiag;

  /* Parameter adjustments */
  --ip;
  --g;
  --h__;
  a_dim1 = *mda;
  a_offset = 1 + a_dim1;
  a -= a_offset;
  --rnorm;
  b_dim1 = *mdb;
  b_offset = 1 + b_dim1;
  b -= b_offset;

  /* Function Body */
  k = 0;
  ldiag = min(*m, *n);
  if (ldiag <= 0) {
    goto L270;
  }
  /*   COMPUTE LMAX */
  i__1 = ldiag;
  for (j = 1; j <= i__1; ++j) {
    if (j == 1) {
      goto L20;
    }
    lmax = j;
    i__2 = *n;
    for (l = j; l <= i__2; ++l) {
      /* Computing 2nd power */
      d__1 = a[j - 1 + l * a_dim1];
      h__[l] -= d__1 * d__1;
      /* L10: */
      if (h__[l] > h__[lmax]) {
        lmax = l;
      }
    }
    d__1 = hmax + factor * h__[lmax];
    if (d__1 - hmax > zero) {
      goto L50;
    }
  L20:
    lmax = j;
    i__2 = *n;
    for (l = j; l <= i__2; ++l) {
      h__[l] = zero;
      i__3 = *m;
      for (i__ = j; i__ <= i__3; ++i__) {
        /* L30: */
        /* Computing 2nd power */
        d__1 = a[i__ + l * a_dim1];
        h__[l] += d__1 * d__1;
      }
      /* L40: */
      if (h__[l] > h__[lmax]) {
        lmax = l;
      }
    }
    hmax = h__[lmax];
  /*   COLUMN INTERCHANGES IF NEEDED */
  L50:
    ip[j] = lmax;
    if (ip[j] == j) {
      goto L70;
    }
    i__2 = *m;
    for (i__ = 1; i__ <= i__2; ++i__) {
      tmp = a[i__ + j * a_dim1];
      a[i__ + j * a_dim1] = a[i__ + lmax * a_dim1];
      /* L60: */
      a[i__ + lmax * a_dim1] = tmp;
    }
    h__[lmax] = h__[j];
  /*   J-TH TRANSFORMATION AND APPLICATION TO A AND B */
  L70:
    /* Computing MIN */
    i__2 = j + 1;
    i__ = min(i__2, *n);
    i__2 = j + 1;
    i__3 = *n - j;
    h12(&c__1,
        &j,
        &i__2,
        m,
        &a[j * a_dim1 + 1],
        &c__1,
        &h__[j],
        &a[i__ * a_dim1 + 1],
        &c__1,
        mda,
        &i__3);
    /* L80: */
    i__2 = j + 1;
    h12(&c__2,
        &j,
        &i__2,
        m,
        &a[j * a_dim1 + 1],
        &c__1,
        &h__[j],
        &b[b_offset],
        &c__1,
        mdb,
        nb);
  }
  /*   DETERMINE PSEUDORANK */
  i__2 = ldiag;
  for (j = 1; j <= i__2; ++j) {
    /* L90: */
    if ((d__1 = a[j + j * a_dim1], fabs(d__1)) <= *tau) {
      goto L100;
    }
  }
  k = ldiag;
  goto L110;
L100:
  k = j - 1;
L110:
  kp1 = k + 1;
  /*   NORM OF RESIDUALS */
  i__2 = *nb;
  for (jb = 1; jb <= i__2; ++jb) {
    /* L130: */
    i__1 = *m - k;
    rnorm[jb] = dnrm2(&i__1, &b[kp1 + jb * b_dim1], &c__1);
  }
  if (k > 0) {
    goto L160;
  }
  i__1 = *nb;
  for (jb = 1; jb <= i__1; ++jb) {
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
      /* L150: */
      b[i__ + jb * b_dim1] = zero;
    }
  }
  goto L270;
L160:
  if (k == *n) {
    goto L180;
  }
  /*   HOUSEHOLDER DECOMPOSITION OF FIRST K ROWS */
  for (i__ = k; i__ >= 1; --i__) {
    /* L170: */
    i__2 = i__ - 1;
    h12(&c__1,
        &i__,
        &kp1,
        n,
        &a[i__ + a_dim1],
        mda,
        &g[i__],
        &a[a_offset],
        mda,
        &c__1,
        &i__2);
  }
L180:
  i__2 = *nb;
  for (jb = 1; jb <= i__2; ++jb) {
    /*   SOLVE K*K TRIANGULAR SYSTEM */
    for (i__ = k; i__ >= 1; --i__) {
      /* Computing MIN */
      i__1 = i__ + 1;
      j = min(i__1, *n);
      /* L210: */
      i__1 = k - i__;
      b[i__ + jb * b_dim1] = (b[i__ + jb * b_dim1]
                              - ddot_sl(&i__1,
                                        &a[i__ + j * a_dim1],
                                        mda,
                                        &b[j + jb * b_dim1],
                                        &c__1))
                             / a[i__ + i__ * a_dim1];
    }
    /*   COMPLETE SOLUTION VECTOR */
    if (k == *n) {
      goto L240;
    }
    i__1 = *n;
    for (j = kp1; j <= i__1; ++j) {
      /* L220: */
      b[j + jb * b_dim1] = zero;
    }
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
      /* L230: */
      h12(&c__2,
          &i__,
          &kp1,
          n,
          &a[i__ + a_dim1],
          mda,
          &g[i__],
          &b[jb * b_dim1 + 1],
          &c__1,
          mdb,
          &c__1);
    }
  /*   REORDER SOLUTION ACCORDING TO PREVIOUS COLUMN INTERCHANGES */
  L240:
    for (j = ldiag; j >= 1; --j) {
      if (ip[j] == j) {
        goto L250;
      }
      l = ip[j];
      tmp = b[l + jb * b_dim1];
      b[l + jb * b_dim1] = b[j + jb * b_dim1];
      b[j + jb * b_dim1] = tmp;
    L250:;
    }
  }
L270:
  *krank = k;
  return 0;
} /* hfti */

/* Subroutine */ int h12(const int* mode,
                         const int* lpivot,
                         const int* l1,
                         const int* m,
                         double* u,
                         const int* iue,
                         double* up,
                         double* c__,
                         const int* ice,
                         const int* icv,
                         const int* ncv) {
  /* Initialized data */

  const double one = 1.;
  const double zero = 0.;

  /* System generated locals */
  int u_dim1, u_offset, i__1, i__2;
  double d__1;

  /* Local variables */
  double b;
  int i__, j, i2, i3, i4;
  double cl, sm;
  int incr;
  double clinv;

  /* Parameter adjustments */
  u_dim1 = *iue;
  u_offset = 1 + u_dim1;
  u -= u_offset;
  --c__;

  /* Function Body */
  if (0 >= *lpivot || *lpivot >= *l1 || *l1 > *m) {
    goto L80;
  }
  cl = (d__1 = u[*lpivot * u_dim1 + 1], fabs(d__1));
  if (*mode == 2) {
    goto L30;
  }
  /*     ****** CONSTRUCT THE TRANSFORMATION ****** */
  i__1 = *m;
  for (j = *l1; j <= i__1; ++j) {
    sm = (d__1 = u[j * u_dim1 + 1], fabs(d__1));
    /* L10: */
    cl = max(sm, cl);
  }
  if (cl <= zero) {
    goto L80;
  }
  clinv = one / cl;
  /* Computing 2nd power */
  d__1 = u[*lpivot * u_dim1 + 1] * clinv;
  sm = d__1 * d__1;
  i__1 = *m;
  for (j = *l1; j <= i__1; ++j) {
    /* L20: */
    /* Computing 2nd power */
    d__1 = u[j * u_dim1 + 1] * clinv;
    sm += d__1 * d__1;
  }
  cl *= sqrt(sm);
  if (u[*lpivot * u_dim1 + 1] > zero) {
    cl = -cl;
  }
  *up = u[*lpivot * u_dim1 + 1] - cl;
  u[*lpivot * u_dim1 + 1] = cl;
  goto L40;
/*     ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C ****** */
L30:
  if (cl <= zero) {
    goto L80;
  }
L40:
  if (*ncv <= 0) {
    goto L80;
  }
  b = *up * u[*lpivot * u_dim1 + 1];
  if (b >= zero) {
    goto L80;
  }
  b = one / b;
  i2 = 1 - *icv + *ice * (*lpivot - 1);
  incr = *ice * (*l1 - *lpivot);
  i__1 = *ncv;
  for (j = 1; j <= i__1; ++j) {
    i2 += *icv;
    i3 = i2 + incr;
    i4 = i3;
    sm = c__[i2] * *up;
    i__2 = *m;
    for (i__ = *l1; i__ <= i__2; ++i__) {
      sm += c__[i3] * u[i__ * u_dim1 + 1];
      /* L50: */
      i3 += *ice;
    }
    if (sm == zero) {
      goto L70;
    }
    sm *= b;
    c__[i2] += sm * *up;
    i__2 = *m;
    for (i__ = *l1; i__ <= i__2; ++i__) {
      c__[i4] += sm * u[i__ * u_dim1 + 1];
      /* L60: */
      i4 += *ice;
    }
  L70:;
  }
L80:
  return 0;
} /* h12 */

/* Subroutine */ int
ldl(const int* n, double* a, double* z__, const double* sigma, double* w) {
  /* Initialized data */

  const double zero = 0.;
  const double one = 1.;
  const double four = 4.;
  const double epmach = 2.22e-16;

  /* System generated locals */
  int i__1, i__2;

  /* Local variables */
  int i__, j;
  double t, u, v;
  int ij;
  double tp, beta, gamma, alpha, delta;

  /* Parameter adjustments */
  --w;
  --z__;
  --a;

  /* Function Body */
  if (*sigma == zero) {
    goto L280;
  }
  ij = 1;
  t = one / *sigma;
  if (*sigma > zero) {
    goto L220;
  }
  /* PREPARE NEGATIVE UPDATE */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    /* L150: */
    w[i__] = z__[i__];
  }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    v = w[i__];
    t += v * v / a[ij];
    i__2 = *n;
    for (j = i__ + 1; j <= i__2; ++j) {
      ++ij;
      /* L160: */
      w[j] -= v * a[ij];
    }
    /* L170: */
    ++ij;
  }
  if (t >= zero) {
    t = epmach / *sigma;
  }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    j = *n + 1 - i__;
    ij -= i__;
    u = w[j];
    w[j] = t;
    /* L210: */
    t -= u * u / a[ij];
  }
L220:
  /* HERE UPDATING BEGINS */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    v = z__[i__];
    delta = v / a[ij];
    if (*sigma < zero) {
      tp = w[i__];
    }
    if (*sigma > zero) {
      tp = t + delta * v;
    }
    alpha = tp / t;
    a[ij] = alpha * a[ij];
    if (i__ == *n) {
      goto L280;
    }
    beta = delta / tp;
    if (alpha > four) {
      goto L240;
    }
    i__2 = *n;
    for (j = i__ + 1; j <= i__2; ++j) {
      ++ij;
      z__[j] -= v * a[ij];
      /* L230: */
      a[ij] += beta * z__[j];
    }
    goto L260;
  L240:
    gamma = t / tp;
    i__2 = *n;
    for (j = i__ + 1; j <= i__2; ++j) {
      ++ij;
      u = a[ij];
      a[ij] = gamma * u + beta * z__[j];
      /* L250: */
      z__[j] -= v * u;
    }
  L260:
    ++ij;
    /* L270: */
    t = tp;
  }
L280:
  return 0;
  /* END OF LDL */
} /* ldl */

double linmin(int* mode,
              const double* ax,
              const double* bx,
              const double* f,
              const double* tol) {
  /* Initialized data */

  const double c__ = .381966011;
  const double eps = 1.5e-8;
  const double zero = 0.;

  /* System generated locals */
  double ret_val, d__1;

  /* Local variables */
  double a, b, d__, e, m, p, q, r__, u, v, w, x, fu, fv, fw, fx, tol1, tol2;

  switch (*mode) {
  case 1:
    goto L10;
  case 2:
    goto L55;
  }
  /*  INITIALIZATION */
  a = *ax;
  b = *bx;
  e = zero;
  v = a + c__ * (b - a);
  w = v;
  x = w;
  ret_val = x;
  *mode = 1;
  goto L100;
/*  MAIN LOOP STARTS HERE */
L10:
  fx = *f;
  fv = fx;
  fw = fv;
L20:
  m = (a + b) * .5;
  tol1 = eps * fabs(x) + *tol;
  tol2 = tol1 + tol1;
  /*  TEST CONVERGENCE */
  if ((d__1 = x - m, fabs(d__1)) <= tol2 - (b - a) * .5) {
    goto L90;
  }
  r__ = zero;
  q = r__;
  p = q;
  if (fabs(e) <= tol1) {
    goto L30;
  }
  /*  FIT PARABOLA */
  r__ = (x - w) * (fx - fv);
  q = (x - v) * (fx - fw);
  p = (x - v) * q - (x - w) * r__;
  q -= r__;
  q += q;
  if (q > zero) {
    p = -p;
  }
  if (q < zero) {
    q = -q;
  }
  r__ = e;
  e = d__;
/*  IS PARABOLA ACCEPTABLE */
L30:
  if (fabs(p) >= (d__1 = q * r__, fabs(d__1)) * .5 || p <= q * (a - x)
      || p >= q * (b - x)) {
    goto L40;
  }
  /*  PARABOLIC INTERPOLATION STEP */
  d__ = p / q;
  /*  F MUST NOT BE EVALUATED TOO CLOSE TO A OR B */
  if (u - a < tol2) {
    d__1 = m - x;
    d__ = copysign(tol1, d__1);
  }
  if (b - u < tol2) {
    d__1 = m - x;
    d__ = copysign(tol1, d__1);
  }
  goto L50;
/*  GOLDEN SECTION STEP */
L40:
  if (x >= m) {
    e = a - x;
  }
  if (x < m) {
    e = b - x;
  }
  d__ = c__ * e;
/*  F MUST NOT BE EVALUATED TOO CLOSE TO X */
L50:
  if (fabs(d__) < tol1) {
    d__ = copysign(tol1, d__);
  }
  u = x + d__;
  ret_val = u;
  *mode = 2;
  goto L100;
L55:
  fu = *f;
  /*  UPDATE A, B, V, W, AND X */
  if (fu > fx) {
    goto L60;
  }
  if (u >= x) {
    a = x;
  }
  if (u < x) {
    b = x;
  }
  v = w;
  fv = fw;
  w = x;
  fw = fx;
  x = u;
  fx = fu;
  goto L85;
L60:
  if (u < x) {
    a = u;
  }
  if (u >= x) {
    b = u;
  }
  if (fu <= fw || w == x) {
    goto L70;
  }
  if (fu <= fv || v == x || v == w) {
    goto L80;
  }
  goto L85;
L70:
  v = w;
  fv = fw;
  w = u;
  fw = fu;
  goto L85;
L80:
  v = u;
  fv = fu;
L85:
  goto L20;
/*  END OF MAIN LOOP */
L90:
  ret_val = x;
  *mode = 3;
L100:
  return ret_val;
  /*  END OF LINMIN */
} /* linmin */

/* ## Following a selection from BLAS Level 1 */
/* Subroutine */ int daxpy_sl(const int* n,
                              const double* da,
                              const double* dx,
                              const int* incx,
                              double* dy,
                              const int* incy) {
  /* System generated locals */
  int i__1;

  /* Local variables */
  int i__, m, ix, iy, mp1;

  /* Parameter adjustments */
  --dy;
  --dx;

  /* Function Body */
  if (*n <= 0) {
    return 0;
  }
  if (*da == 0.) {
    return 0;
  }
  if (*incx == 1 && *incy == 1) {
    goto L20;
  }
  /*        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS */
  /*        NOT EQUAL TO 1 */
  ix = 1;
  iy = 1;
  if (*incx < 0) {
    ix = (-(*n) + 1) * *incx + 1;
  }
  if (*incy < 0) {
    iy = (-(*n) + 1) * *incy + 1;
  }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    dy[iy] += *da * dx[ix];
    ix += *incx;
    iy += *incy;
    /* L10: */
  }
  return 0;
/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */
/*        CLEAN-UP LOOP */
L20:
  m = *n % 4;
  if (m == 0) {
    goto L40;
  }
  i__1 = m;
  for (i__ = 1; i__ <= i__1; ++i__) {
    dy[i__] += *da * dx[i__];
    /* L30: */
  }
  if (*n < 4) {
    return 0;
  }
L40:
  mp1 = m + 1;
  i__1 = *n;
  for (i__ = mp1; i__ <= i__1; i__ += 4) {
    dy[i__] += *da * dx[i__];
    dy[i__ + 1] += *da * dx[i__ + 1];
    dy[i__ + 2] += *da * dx[i__ + 2];
    dy[i__ + 3] += *da * dx[i__ + 3];
    /* L50: */
  }
  return 0;
} /* daxpy_sl */

/* Subroutine */ int dcopy(const int* n,
                           const double* dx,
                           const int* incx,
                           double* dy,
                           const int* incy) {
  /* System generated locals */
  int i__1;

  /* Local variables */
  int i__, m, ix, iy, mp1;

  /* Parameter adjustments */
  --dy;
  --dx;

  /* Function Body */
  if (*n <= 0) {
    return 0;
  }
  if (*incx == 1 && *incy == 1) {
    goto L20;
  }
  /*        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS */
  /*        NOT EQUAL TO 1 */
  ix = 1;
  iy = 1;
  if (*incx < 0) {
    ix = (-(*n) + 1) * *incx + 1;
  }
  if (*incy < 0) {
    iy = (-(*n) + 1) * *incy + 1;
  }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    dy[iy] = dx[ix];
    ix += *incx;
    iy += *incy;
    /* L10: */
  }
  return 0;
/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */
/*        CLEAN-UP LOOP */
L20:
  m = *n % 7;
  if (m == 0) {
    goto L40;
  }
  i__1 = m;
  for (i__ = 1; i__ <= i__1; ++i__) {
    dy[i__] = dx[i__];
    /* L30: */
  }
  if (*n < 7) {
    return 0;
  }
L40:
  mp1 = m + 1;
  i__1 = *n;
  for (i__ = mp1; i__ <= i__1; i__ += 7) {
    dy[i__] = dx[i__];
    dy[i__ + 1] = dx[i__ + 1];
    dy[i__ + 2] = dx[i__ + 2];
    dy[i__ + 3] = dx[i__ + 3];
    dy[i__ + 4] = dx[i__ + 4];
    dy[i__ + 5] = dx[i__ + 5];
    dy[i__ + 6] = dx[i__ + 6];
    /* L50: */
  }
  return 0;
} /* dcopy */

double ddot_sl(const int* n,
               const double* dx,
               const int* incx,
               const double* dy,
               const int* incy) {
  /* System generated locals */
  int i__1;
  double ret_val;

  /* Local variables */
  int i__, m, ix, iy, mp1;
  double dtemp;

  /* Parameter adjustments */
  --dy;
  --dx;

  /* Function Body */
  ret_val = 0.;
  dtemp = 0.;
  if (*n <= 0) {
    return ret_val;
  }
  if (*incx == 1 && *incy == 1) {
    goto L20;
  }
  /*        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS */
  /*          NOT EQUAL TO 1 */
  ix = 1;
  iy = 1;
  if (*incx < 0) {
    ix = (-(*n) + 1) * *incx + 1;
  }
  if (*incy < 0) {
    iy = (-(*n) + 1) * *incy + 1;
  }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    dtemp += dx[ix] * dy[iy];
    ix += *incx;
    iy += *incy;
    /* L10: */
  }
  ret_val = dtemp;
  return ret_val;
/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */
/*        CLEAN-UP LOOP */
L20:
  m = *n % 5;
  if (m == 0) {
    goto L40;
  }
  i__1 = m;
  for (i__ = 1; i__ <= i__1; ++i__) {
    dtemp += dx[i__] * dy[i__];
    /* L30: */
  }
  if (*n < 5) {
    goto L60;
  }
L40:
  mp1 = m + 1;
  i__1 = *n;
  for (i__ = mp1; i__ <= i__1; i__ += 5) {
    dtemp = dtemp + dx[i__] * dy[i__] + dx[i__ + 1] * dy[i__ + 1]
            + dx[i__ + 2] * dy[i__ + 2] + dx[i__ + 3] * dy[i__ + 3]
            + dx[i__ + 4] * dy[i__ + 4];
    /* L50: */
  }
L60:
  ret_val = dtemp;
  return ret_val;
} /* ddot_sl */

double dnrm1(const int* n, const double* x, const int* i__, const int* j) {
  /* Initialized data */

  const double zero = 0.;
  const double one = 1.;

  /* System generated locals */
  int i__1;
  double ret_val, d__1, d__2, d__3;

  /* Local variables */
  int k;
  double sum, temp, scale, snormx;

  /* Parameter adjustments */
  --x;

  /* Function Body */

  snormx = zero;
  i__1 = *j;
  for (k = *i__; k <= i__1; ++k) {
    /* L10: */
    /* Computing MAX */
    d__2 = snormx, d__3 = (d__1 = x[k], fabs(d__1));
    snormx = max(d__2, d__3);
  }
  ret_val = snormx;
  if (snormx == zero) {
    return ret_val;
  }
  scale = snormx;
  if (snormx >= one) {
    scale = sqrt(snormx);
  }
  sum = zero;
  i__1 = *j;
  for (k = *i__; k <= i__1; ++k) {
    temp = zero;
    if ((d__1 = x[k], fabs(d__1)) + scale != scale) {
      temp = x[k] / snormx;
    }
    if (one + temp != one) {
      sum += temp * temp;
    }
    /* L20: */
  }
  sum = sqrt(sum);
  ret_val = snormx * sum;
  return ret_val;
} /* dnrm1 */

double dnrm2(const int* n, const double* dx, const int* incx) {
  /* Initialized data */

  const double zero = 0.;
  const double one = 1.;
  const double cutlo = 8.232e-11;
  const double cuthi = 1.304e19;

  /* System generated locals */
  int i__1, i__2;
  double ret_val, d__1;

  /* Local variables */
  int i__, j, nn;
  double sum, xmax;
  int next;
  double hitest;

  /* Parameter adjustments */
  --dx;

  if (*n > 0) {
    goto L10;
  }
  ret_val = zero;
  goto L300;
L10:
  next = 0;
  sum = zero;
  nn = *n * *incx;
  /*                       BEGIN MAIN LOOP */
  i__ = 1;
L20:
  switch (next) {
  case 0:
    goto L30;
  case 1:
    goto L50;
  case 2:
    goto L70;
  case 3:
    goto L110;
  }
L30:
  if ((d__1 = dx[i__], fabs(d__1)) > cutlo) {
    goto L85;
  }
  next = 1;
  xmax = zero;
/*                        PHASE 1.  SUM IS ZERO */
L50:
  if (dx[i__] == zero) {
    goto L200;
  }
  if ((d__1 = dx[i__], fabs(d__1)) > cutlo) {
    goto L85;
  }
  /*                        PREPARE FOR PHASE 2. */
  next = 2;
  goto L105;
/*                        PREPARE FOR PHASE 4. */
L100:
  i__ = j;
  next = 3;
  sum = sum / dx[i__] / dx[i__];
L105:
  xmax = (d__1 = dx[i__], fabs(d__1));
  goto L115;
/*                   PHASE 2.  SUM IS SMALL. */
/*                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW. */
L70:
  if ((d__1 = dx[i__], fabs(d__1)) > cutlo) {
    goto L75;
  }
/*                   COMMON CODE FOR PHASES 2 AND 4. */
/*                   IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW. */
L110:
  if ((d__1 = dx[i__], fabs(d__1)) <= xmax) {
    goto L115;
  }
  /* Computing 2nd power */
  d__1 = xmax / dx[i__];
  sum = one + sum * (d__1 * d__1);
  xmax = (d__1 = dx[i__], fabs(d__1));
  goto L200;
L115:
  /* Computing 2nd power */
  d__1 = dx[i__] / xmax;
  sum += d__1 * d__1;
  goto L200;
/*                  PREPARE FOR PHASE 3. */
L75:
  sum = sum * xmax * xmax;
/*     FOR REAL OR D.P. SET HITEST = CUTHI/N */
/*     FOR COMPLEX      SET HITEST = CUTHI/(2*N) */
L85:
  hitest = cuthi / (float) (*n);
  /*                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING. */
  i__1 = nn;
  i__2 = *incx;
  for (j = i__; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
    if ((d__1 = dx[j], fabs(d__1)) >= hitest) {
      goto L100;
    }
    /* L95: */
    /* Computing 2nd power */
    d__1 = dx[j];
    sum += d__1 * d__1;
  }
  ret_val = sqrt(sum);
  goto L300;
L200:
  i__ += *incx;
  if (i__ <= nn) {
    goto L20;
  }
  /*              END OF MAIN LOOP. */
  /*              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING. */
  ret_val = xmax * sqrt(sum);
L300:
  return ret_val;
} /* dnrm2 */

/* Subroutine */ int dsrot(const int* n,
                           double* dx,
                           const int* incx,
                           double* dy,
                           const int* incy,
                           const double* c__,
                           const double* s) {
  /* System generated locals */
  int i__1;

  /* Local variables */
  int i__, ix, iy;
  double dtemp;

  /* Parameter adjustments */
  --dy;
  --dx;

  /* Function Body */
  if (*n <= 0) {
    return 0;
  }
  if (*incx == 1 && *incy == 1) {
    goto L20;
  }
  /*       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL */
  /*         TO 1 */
  ix = 1;
  iy = 1;
  if (*incx < 0) {
    ix = (-(*n) + 1) * *incx + 1;
  }
  if (*incy < 0) {
    iy = (-(*n) + 1) * *incy + 1;
  }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    dtemp = *c__ * dx[ix] + *s * dy[iy];
    dy[iy] = *c__ * dy[iy] - *s * dx[ix];
    dx[ix] = dtemp;
    ix += *incx;
    iy += *incy;
    /* L10: */
  }
  return 0;
/*       CODE FOR BOTH INCREMENTS EQUAL TO 1 */
L20:
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    dtemp = *c__ * dx[i__] + *s * dy[i__];
    dy[i__] = *c__ * dy[i__] - *s * dx[i__];
    dx[i__] = dtemp;
    /* L30: */
  }
  return 0;
} /* dsrot */

/* Subroutine */ int dsrotg(double* da, double* db, double* c__, double* s) {
  /* Initialized data */

  const double one = 1.;
  const double zero = 0.;

  /* System generated locals */
  double d__1, d__2;

  /* Local variables */
  double r__, z__, roe, scale;

  roe = *db;
  if (fabs(*da) > fabs(*db)) {
    roe = *da;
  }
  scale = fabs(*da) + fabs(*db);
  if (scale != zero) {
    goto L10;
  }
  *c__ = one;
  *s = zero;
  r__ = zero;
  goto L20;
L10:
  /* Computing 2nd power */
  d__1 = *da / scale;
  /* Computing 2nd power */
  d__2 = *db / scale;
  r__ = scale * sqrt(d__1 * d__1 + d__2 * d__2);
  r__ = copysign(one, roe) * r__;
  *c__ = *da / r__;
  *s = *db / r__;
L20:
  z__ = *s;
  if (fabs(*c__) > zero && fabs(*c__) <= *s) {
    z__ = one / *c__;
  }
  *da = r__;
  *db = z__;
  return 0;
} /* dsrotg */

/* Subroutine */ int
dscal_sl(const int* n, const double* da, double* dx, const int* incx) {
  /* System generated locals */
  int i__1, i__2;

  /* Local variables */
  int i__, m, mp1, nincx;

  /* Parameter adjustments */
  --dx;

  /* Function Body */
  if (*n <= 0) {
    return 0;
  }
  if (*incx == 1) {
    goto L20;
  }
  /*        CODE FOR INCREMENT NOT EQUAL TO 1 */
  nincx = *n * *incx;
  i__1 = nincx;
  i__2 = *incx;
  for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
    dx[i__] = *da * dx[i__];
    /* L10: */
  }
  return 0;
/*        CODE FOR INCREMENT EQUAL TO 1 */
/*        CLEAN-UP LOOP */
L20:
  m = *n % 5;
  if (m == 0) {
    goto L40;
  }
  i__2 = m;
  for (i__ = 1; i__ <= i__2; ++i__) {
    dx[i__] = *da * dx[i__];
    /* L30: */
  }
  if (*n < 5) {
    return 0;
  }
L40:
  mp1 = m + 1;
  i__2 = *n;
  for (i__ = mp1; i__ <= i__2; i__ += 5) {
    dx[i__] = *da * dx[i__];
    dx[i__ + 1] = *da * dx[i__ + 1];
    dx[i__ + 2] = *da * dx[i__ + 2];
    dx[i__ + 3] = *da * dx[i__ + 3];
    dx[i__ + 4] = *da * dx[i__ + 4];
    /* L50: */
  }
  return 0;
} /* dscal_sl */

/* Subroutine */ int
bound(const int* n, double* x, const double* xl, const double* xu) {
  /* System generated locals */
  int i__1;

  /* Local variables */
  int i__;

  /* Parameter adjustments */
  --xu;
  --xl;
  --x;

  /* Function Body */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    /*        Note that xl(i) and xu(i) may be NaN to indicate no bound */
    if (xl[i__] == xl[i__] && x[i__] < xl[i__]) {
      x[i__] = xl[i__];
    } else if (xu[i__] == xu[i__] && x[i__] > xu[i__]) {
      x[i__] = xu[i__];
    }
  }
  return 0;
} /* bound */
