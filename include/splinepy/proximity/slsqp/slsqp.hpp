#ifdef __cplusplus
extern "C" {
#endif
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

int ldl(const int* n, double* a, double* z__, const double* sigma, double* w);

double linmin(int* mode,
              const double* ax,
              const double* bx,
              const double* f,
              const double* tol);

int daxpy_sl(const int* n,
             const double* da,
             const double* dx,
             const int* incx,
             double* dy,
             const int* incy);

int dcopy(const int* n,
          const double* dx,
          const int* incx,
          double* dy,
          const int* incy);

double ddot_sl(const int* n,
               const double* dx,
               const int* incx,
               const double* dy,
               const int* incy);

double dnrm1(const int* n, const double* x, const int* i__, const int* j);

double dnrm2(const int* n, const double* dx, const int* incx);

int dsrot(const int* n,
          double* dx,
          const int* incx,
          double* dy,
          const int* incy,
          const double* c__,
          const double* s);

int dsrotg(double* da, double* db, double* c__, double* s);

int dscal_sl(const int* n, const double* da, double* dx, const int* incx);

int bound(const int* n, double* x, const double* xl, const double* xu);

#ifdef __cplusplus
} /* extern "C */
#endif
