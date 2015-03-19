#ifndef _ALL_R_ROUTINE_H
#define _ALL_R_ROUTINE_H

#include <R.h>
#include <R_ext/Utils.h>
#include <Rinternals.h>
#include <Rdefines.h>

#ifdef __cplusplus
 extern "C" {
#endif

  SEXP phfit_estep_gen_wtime(SEXP Rph, SEXP Rdata, SEXP eps, SEXP ufact);
  SEXP phfit_estep_gen_group(SEXP Rph, SEXP Rbaralpha, SEXP Rdata, SEXP gdatlast, SEXP eps, SEXP ufact);
  SEXP phfit_herlang_estep_wtime(SEXP Rerl, SEXP Rdata);
  SEXP phfit_herlang_estep_group(SEXP Rerl, SEXP Rdata, SEXP gdatlast);
  SEXP mapfit_estep_gen_group(SEXP Rmap, SEXP Rdata, SEXP eps, SEXP ufact);
  SEXP mapfit_hmm_erlang_estep(SEXP Rerl, SEXP Rdata);
  SEXP mapfit_estep_gmmpp(SEXP Rmap, SEXP Rdata, SEXP eps, SEXP divide);

  SEXP phfit_mstep_gen(SEXP Rph, SEXP eres, SEXP Rdata);
  SEXP phfit_mstep_cf1(SEXP Rph, SEXP eres, SEXP Rdata);
  SEXP phfit_herlang_mstep(SEXP Rerl, SEXP eres, SEXP Rdata);
  SEXP mapfit_mstep_gen(SEXP Rmap, SEXP eres, SEXP Rdata);
  SEXP mapfit_hmm_erlang_mstep(SEXP Rerl, SEXP eres, SEXP Rdata);

  SEXP ctmc_st(SEXP RQ, SEXP Rite, SEXP Reps);
  SEXP marsolve(SEXP trans, SEXP Ralpha, SEXP RQ, SEXP Rx, SEXP Rite, SEXP Reps);
  SEXP mpow_mat(SEXP RMA, SEXP Rm);
  SEXP mexp_pade(SEXP Rn, SEXP RMA, SEXP eps);
  SEXP mexp_unif(SEXP Rn, SEXP RMA, SEXP Rt, SEXP eps, SEXP ufact, SEXP atol);
  SEXP mexp_unifvec(SEXP trans, SEXP Rn, SEXP RMA, SEXP Rx, SEXP Rt, SEXP eps, SEXP ufact, SEXP atol);
  SEXP mexp_unifseq(SEXP trans, SEXP Rn, SEXP RMA, SEXP Rx, SEXP Rt, SEXP eps, SEXP ufact, SEXP atol);
  SEXP mexp_kryvec(SEXP trans, SEXP Rn, SEXP RMA, SEXP Rx, SEXP Rt, SEXP Rksub, SEXP ite, SEXP tol, SEXP eps);

#ifdef __cplusplus
}
#endif

#endif
