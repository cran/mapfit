#include "allheader.h"
#include <R_ext/Rdynload.h>

#define CALLDEF(name, n) {#name, (DL_FUNC) &name, n}
#define EXTDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef callMethods[] = {
    CALLDEF(phfit_estep_gen_wtime, 4),
    CALLDEF(phfit_estep_gen_group, 6),
    CALLDEF(phfit_herlang_estep_wtime, 2),
    CALLDEF(phfit_herlang_estep_group, 3),
    CALLDEF(mapfit_estep_gen_group, 4),
    CALLDEF(mapfit_hmm_erlang_estep, 2),
    CALLDEF(mapfit_estep_gmmpp, 4),
    CALLDEF(phfit_mstep_gen, 3),
    CALLDEF(phfit_mstep_cf1, 3),
    CALLDEF(phfit_herlang_mstep, 3),
    CALLDEF(mapfit_mstep_gen, 3),
    CALLDEF(mapfit_hmm_erlang_mstep, 3),
    CALLDEF(ctmc_st, 3),
    CALLDEF(marsolve, 6),
    CALLDEF(mpow_mat, 2),
    CALLDEF(mexp_pade, 3),
    CALLDEF(mexp_unif, 6),
    CALLDEF(mexp_unifvec, 8),
    CALLDEF(mexp_unifseq, 8),
    CALLDEF(mexp_kryvec, 9),
    {NULL, NULL, 0}
};

void R_init_mapfit(DllInfo *dll) {
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
//    R_useDynamicSymbols(dll, FALSE);
}

