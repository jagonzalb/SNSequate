#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#ifdef HAVE_F77_UNDERSCORE
# define F77_SYMBOL(x)	x ## _
# define F77_QSYMBOL(x)	#x "_"
#else
# define F77_SYMBOL(x)	x
# define F77_QSYMBOL(x) #x
#endif


void F77_SYMBOL(tldbppdensity)(double *y, double *x, int *nrec, int *p, int *npred,
             int *ngrid, double *grid, double* xpred, int *maxn, int *nu, 
             double *alpha, double *lambda, double *tau1, double *tau2, double *psiinv,
             double *s0invm, double *s0inv, int *kk, double *gp, double *beta, 
             double *theta, double *mub, double *sb, int *mcmc, int *nsave,
             double *slice, double *acrate, double *thetasave, double *randsave,
             double *fmean, double *flow, double *fupp, double *meanfpm,
             double *meanfpl, double *meanfph, double *cpo, int *seed, int *iflag,
             double *sbinv, double *workm1, double *workv1, double *workv2, 
             double *workmh1, double *workmh2, double *betal, double *betar,
             double *beta1, double *thetal, double *thetar, double *theta1,
             double *workdpw, double *weight, double *fw, double *fw2, double *fs,
             double *fm, double *worksam, double *v);


static const R_FortranMethodDef FortEntries[] = {
  {"tldbppdensity", (DL_FUNC) &F77_SYMBOL(tldbppdensity),  58},
  {NULL, NULL, 0}
};


void R_init_SNSequate(DllInfo *dll)
{
  //    R_registerRoutines(dll, NULL, NULL, callMethods, NULL);
  R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
