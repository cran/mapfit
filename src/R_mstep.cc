
#include "R_common.h"

extern "C" {

  SEXP phfit_mstep_gen(SEXP Rph, SEXP eres, SEXP Rdata) {
    int n = asInteger(getSlot(Rph, "size"));

    SEXP ans_alpha = PROTECT(allocVector(REALSXP, n));
    SEXP ans_xi = PROTECT(allocVector(REALSXP, n));

    sci::vector<double> alpha(n, REAL(AS_NUMERIC(ans_alpha)));
    alpha = REAL(AS_NUMERIC(getSlot(Rph, "alpha")));
    sci::vector<double> xi(n, REAL(AS_NUMERIC(ans_xi)));
    xi = REAL(AS_NUMERIC(getSlot(Rph, "xi")));
    sci::matrix<double>* Qptr = createMatrix(getSlot(Rph, "Q"));
    SEXP ans_Q = PROTECT(allocVector(REALSXP, Qptr->size));
    Qptr->ptr = REAL(AS_NUMERIC(ans_Q));

    double etotal = asReal(getListElement(eres, "etotal"));
    sci::vector<double> eb(n, REAL(AS_NUMERIC(getListElement(eres, "eb"))));
    sci::vector<double> ey(n, REAL(AS_NUMERIC(getListElement(eres, "ey"))));
    sci::vector<double> ez(n, REAL(AS_NUMERIC(getListElement(eres, "ez"))));
    sci::matrix<double>* enptr = sci::dnewcopy(*Qptr, REAL(AS_NUMERIC(getListElement(eres, "en"))));

    mapfit::phase_mstep(etotal, eb, ey, ez, *enptr, alpha, xi, *Qptr);

    SEXP ans = PROTECT(allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, ans_alpha);
    SET_VECTOR_ELT(ans, 1, ans_xi);
    SET_VECTOR_ELT(ans, 2, ans_Q);
    UNPROTECT(4);

    delete Qptr;
    delete enptr;

    return ans;
  }

  SEXP phfit_mstep_cf1(SEXP Rph, SEXP eres, SEXP Rdata) {
    int n = asInteger(getSlot(Rph, "size"));

    SEXP ans_alpha = PROTECT(allocVector(REALSXP, n));
    SEXP ans_xi = PROTECT(allocVector(REALSXP, n));

    sci::vector<double> alpha(n, REAL(AS_NUMERIC(ans_alpha)));
    alpha = REAL(AS_NUMERIC(getSlot(Rph, "alpha")));
    sci::vector<double> xi(n, REAL(AS_NUMERIC(ans_xi)));
    xi = REAL(AS_NUMERIC(getSlot(Rph, "xi")));
    sci::matrix<double>* Qptr = createMatrix(getSlot(Rph, "Q"));
    SEXP ans_Q = PROTECT(allocVector(REALSXP, Qptr->size));
    Qptr->ptr = REAL(AS_NUMERIC(ans_Q));

    double etotal = asReal(getListElement(eres, "etotal"));
    sci::vector<double> eb(n, REAL(AS_NUMERIC(getListElement(eres, "eb"))));
    sci::vector<double> ey(n, REAL(AS_NUMERIC(getListElement(eres, "ey"))));
    sci::vector<double> ez(n, REAL(AS_NUMERIC(getListElement(eres, "ez"))));
    sci::matrix<double>* enptr = sci::dnewcopy(*Qptr, REAL(AS_NUMERIC(getListElement(eres, "en"))));

    mapfit::phase_mstep(etotal, eb, ey, ez, *enptr, alpha, xi, *Qptr);
    mapfit::phase_bidiag_to_cf1(alpha, xi, *Qptr);

    SEXP ans = PROTECT(allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, ans_alpha);
    SET_VECTOR_ELT(ans, 1, ans_xi);
    SET_VECTOR_ELT(ans, 2, ans_Q);
    UNPROTECT(4);

    delete Qptr;
    delete enptr;

    return ans;
  }

  SEXP phfit_herlang_mstep(SEXP Rerl, SEXP eres, SEXP Rdata) {
    int n = asInteger(getSlot(Rerl, "size"));
    SEXP ans_mixrate = PROTECT(allocVector(REALSXP, n));
    SEXP ans_rate = PROTECT(allocVector(REALSXP, n));

    sci::vector<double> mixrate(n, REAL(AS_NUMERIC(ans_mixrate)));
    mixrate = REAL(AS_NUMERIC(getSlot(Rerl, "mixrate")));
    sci::vector<double> rate(n, REAL(AS_NUMERIC(ans_rate)));
    rate = REAL(AS_NUMERIC(getSlot(Rerl, "rate")));
    sci::vector<int> shape(n, INTEGER(AS_INTEGER(getSlot(Rerl, "shape"))));

    double etotal = asReal(getListElement(eres, "etotal"));
    sci::vector<double> eb(n, REAL(AS_NUMERIC(getListElement(eres, "eb"))));
    sci::vector<double> ew(n, REAL(AS_NUMERIC(getListElement(eres, "ew"))));

    mapfit::phase_erlang_mstep(etotal, eb, ew, mixrate, shape, rate);

    SEXP ans = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, ans_mixrate);
    SET_VECTOR_ELT(ans, 1, ans_rate);
    UNPROTECT(3);
    return ans;
  }

  SEXP mapfit_mstep_gen(SEXP Rmap, SEXP eres, SEXP Rdata) {
    int n = asInteger(getSlot(Rmap, "size"));

    SEXP ans_alpha = PROTECT(allocVector(REALSXP, n));
    sci::vector<double> alpha(n, REAL(AS_NUMERIC(ans_alpha)));
    alpha = REAL(AS_NUMERIC(getSlot(Rmap, "alpha")));

    sci::matrix<double>* D0ptr = createMatrix(getSlot(Rmap, "D0"));
    SEXP ans_D0 = PROTECT(allocVector(REALSXP, D0ptr->size));
    D0ptr->ptr = REAL(AS_NUMERIC(ans_D0));

    sci::matrix<double>* D1ptr = createMatrix(getSlot(Rmap, "D1"));
    SEXP ans_D1 = PROTECT(allocVector(REALSXP, D1ptr->size));
    D1ptr->ptr = REAL(AS_NUMERIC(ans_D1));

    sci::vector<double> eb(n, REAL(AS_NUMERIC(getListElement(eres, "eb"))));
    sci::vector<double> ez(n, REAL(AS_NUMERIC(getListElement(eres, "ez"))));
    sci::matrix<double>* en0ptr = sci::dnewcopy(*D0ptr, REAL(AS_NUMERIC(getListElement(eres, "en0"))));
    sci::matrix<double>* en1ptr = sci::dnewcopy(*D1ptr, REAL(AS_NUMERIC(getListElement(eres, "en1"))));

    mapfit::map_mstep(eb, ez, *en0ptr, *en1ptr, alpha, *D0ptr, *D1ptr);

    SEXP ans = PROTECT(allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, ans_alpha);
    SET_VECTOR_ELT(ans, 1, ans_D0);
    SET_VECTOR_ELT(ans, 2, ans_D1);
    UNPROTECT(4);

    delete D0ptr;
    delete D1ptr;
    delete en0ptr;
    delete en1ptr;

    return ans;
  }

  SEXP mapfit_hmm_erlang_mstep(SEXP Rerl, SEXP eres, SEXP Rdata) {
    int n = asInteger(getSlot(Rerl, "size"));
    SEXP ans_alpha = PROTECT(allocVector(REALSXP, n));
    SEXP ans_rate = PROTECT(allocVector(REALSXP, n));

    sci::vector<double> alpha(n, REAL(AS_NUMERIC(ans_alpha)));
    alpha = REAL(AS_NUMERIC(getSlot(Rerl, "alpha")));

    sci::vector<double> rate(n, REAL(AS_NUMERIC(ans_rate)));
    rate = REAL(AS_NUMERIC(getSlot(Rerl, "rate")));

    sci::matrix<double>* Pptr = createMatrix(getSlot(Rerl, "P"));
    SEXP ans_P = PROTECT(allocVector(REALSXP, Pptr->size));
    Pptr->ptr = REAL(AS_NUMERIC(ans_P));

    sci::vector<int> shape(n, INTEGER(AS_INTEGER(getSlot(Rerl, "shape"))));

    sci::vector<double> eb(n, REAL(AS_NUMERIC(getListElement(eres, "eb"))));
    sci::matrix<double>* enptr = sci::dnewcopy(*Pptr, REAL(AS_NUMERIC(getListElement(eres, "en"))));
    sci::vector<double> ew0(n, REAL(AS_NUMERIC(getListElement(eres, "ew0"))));
    sci::vector<double> ew1(n, REAL(AS_NUMERIC(getListElement(eres, "ew1"))));

    mapfit::hmm_erlang_mstep(eb, *enptr, ew0, ew1, alpha, *Pptr, shape, rate);

    SEXP ans = PROTECT(allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, ans_alpha);
    SET_VECTOR_ELT(ans, 1, ans_rate);
    SET_VECTOR_ELT(ans, 2, ans_P);
    UNPROTECT(4);

    delete Pptr;
    delete enptr;
    
    return ans;
  }
}

