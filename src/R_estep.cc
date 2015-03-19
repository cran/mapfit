
#include "R_common.h"

extern "C" {

  SEXP phfit_estep_gen_wtime(SEXP Rph, SEXP Rdata, SEXP eps, SEXP ufact) {
    int n = asInteger(getSlot(Rph, "size"));
    sci::vector<double> alpha(n, REAL(AS_NUMERIC(getSlot(Rph, "alpha"))));
    sci::vector<double> xi(n, REAL(AS_NUMERIC(getSlot(Rph, "xi"))));
    sci::matrix<double>* Qptr = createMatrix(getSlot(Rph, "Q"));
    sci::matrix<double>* Pptr = sci::dnewcopy(*Qptr);
    double qv = mexp::unif(*Qptr, *Pptr, asReal(ufact));
    int m = asInteger(getSlot(Rdata, "size"));
    sci::vector<double> tdat(m, REAL(AS_NUMERIC(getSlot(Rdata, "diff"))));
    sci::vector<double> wdat(m, REAL(AS_NUMERIC(getListElement(getSlot(Rdata, "data"), "weight"))));

    SEXP ans = PROTECT(allocVector(VECSXP, 6));
    SEXP ans_etotal = PROTECT(allocVector(REALSXP, 1));
    SEXP ans_eb = PROTECT(allocVector(REALSXP, n));
    SEXP ans_ey = PROTECT(allocVector(REALSXP, n));
    SEXP ans_ez = PROTECT(allocVector(REALSXP, n));
    SEXP ans_en = PROTECT(allocVector(REALSXP, Qptr->size));
    SEXP ans_llf = PROTECT(allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 0, ans_etotal);
    SET_VECTOR_ELT(ans, 1, ans_eb);
    SET_VECTOR_ELT(ans, 2, ans_ey);
    SET_VECTOR_ELT(ans, 3, ans_ez);
    SET_VECTOR_ELT(ans, 4, ans_en);
    SET_VECTOR_ELT(ans, 5, ans_llf);
    double *etotal = REAL(AS_NUMERIC(ans_etotal));
    sci::vector<double> eb(n, REAL(AS_NUMERIC(ans_eb)));
    sci::vector<double> ey(n, REAL(AS_NUMERIC(ans_ey)));
    sci::vector<double> ez(n, REAL(AS_NUMERIC(ans_ez)));
    sci::matrix<double>* enptr = sci::dnewcopy(*Qptr, REAL(AS_NUMERIC(ans_en)));
    double *llf = REAL(AS_NUMERIC(ans_llf));

    *llf = mapfit::phase_estep_wtime(alpha, xi, *Qptr, *Pptr, qv, tdat, wdat,
        *etotal, eb, ey, ez, *enptr, asReal(eps));

    UNPROTECT(7);

    delete Qptr;
    delete Pptr;
    delete enptr;

    return ans;
  }

  SEXP phfit_estep_gen_group(SEXP Rph, SEXP Rbaralpha, SEXP Rdata, SEXP gdatlast,
    SEXP eps, SEXP ufact) {
    int n = asInteger(getSlot(Rph, "size"));
    sci::vector<double> alpha(n, REAL(AS_NUMERIC(getSlot(Rph, "alpha"))));
    sci::vector<double> baralpha(n, REAL(AS_NUMERIC(Rbaralpha)));
    sci::vector<double> xi(n, REAL(AS_NUMERIC(getSlot(Rph, "xi"))));
    sci::vector<double> one(n, 1.0);
    sci::matrix<double>* Qptr = createMatrix(getSlot(Rph, "Q"));
    sci::matrix<double>* Pptr = sci::dnewcopy(*Qptr);
    double qv = mexp::unif(*Qptr, *Pptr, asReal(ufact));
    int m = asInteger(getSlot(Rdata, "size"));
    sci::vector<double> tdat(m, REAL(AS_NUMERIC(getListElement(getSlot(Rdata, "data"), "time"))));
    sci::vector<int> gdat(m, INTEGER(AS_INTEGER(getListElement(getSlot(Rdata, "data"), "counts"))));
    sci::vector<int> idat(m, INTEGER(AS_INTEGER(getListElement(getSlot(Rdata, "data"), "instant"))));

    SEXP ans = PROTECT(allocVector(VECSXP, 6));
    SEXP ans_etotal = PROTECT(allocVector(REALSXP, 1));
    SEXP ans_eb = PROTECT(allocVector(REALSXP, n));
    SEXP ans_ey = PROTECT(allocVector(REALSXP, n));
    SEXP ans_ez = PROTECT(allocVector(REALSXP, n));
    SEXP ans_en = PROTECT(allocVector(REALSXP, Qptr->size));
    SEXP ans_llf = PROTECT(allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 0, ans_etotal);
    SET_VECTOR_ELT(ans, 1, ans_eb);
    SET_VECTOR_ELT(ans, 2, ans_ey);
    SET_VECTOR_ELT(ans, 3, ans_ez);
    SET_VECTOR_ELT(ans, 4, ans_en);
    SET_VECTOR_ELT(ans, 5, ans_llf);
    double *etotal = REAL(AS_NUMERIC(ans_etotal));
    sci::vector<double> eb(n, REAL(AS_NUMERIC(ans_eb)));
    sci::vector<double> ey(n, REAL(AS_NUMERIC(ans_ey)));
    sci::vector<double> ez(n, REAL(AS_NUMERIC(ans_ez)));
    sci::matrix<double>* enptr = sci::dnewcopy(*Qptr, REAL(AS_NUMERIC(ans_en)));
    double *llf = REAL(AS_NUMERIC(ans_llf));

    *llf = mapfit::phase_estep_group_trunc(alpha, baralpha, xi, one, *Qptr, *Pptr,
        qv, tdat, gdat, asInteger(gdatlast), idat, *etotal, eb, ey, ez, *enptr, asReal(eps));

    UNPROTECT(7);

    delete Qptr;
    delete Pptr;
    delete enptr;

    return ans;
  }

  SEXP phfit_herlang_estep_wtime(SEXP Rerl, SEXP Rdata) {
    int n = asInteger(getSlot(Rerl, "size"));
    sci::vector<double> mixrate(n, REAL(AS_NUMERIC(getSlot(Rerl, "mixrate"))));
    sci::vector<int> shape(n, INTEGER(AS_INTEGER(getSlot(Rerl, "shape"))));
    sci::vector<double> rate(n, REAL(AS_NUMERIC(getSlot(Rerl, "rate"))));
    int m = asInteger(getSlot(Rdata, "size"));
    sci::vector<double> tdat(m, REAL(AS_NUMERIC(getSlot(Rdata, "diff"))));
    sci::vector<double> wdat(m, REAL(AS_NUMERIC(getListElement(getSlot(Rdata, "data"), "weight"))));

    SEXP ans = PROTECT(allocVector(VECSXP, 4));
    SEXP ans_etotal = PROTECT(allocVector(REALSXP, 1));
    SEXP ans_eb = PROTECT(allocVector(REALSXP, n));
    SEXP ans_ew = PROTECT(allocVector(REALSXP, n));
    SEXP ans_llf = PROTECT(allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 0, ans_etotal);
    SET_VECTOR_ELT(ans, 1, ans_eb);
    SET_VECTOR_ELT(ans, 2, ans_ew);
    SET_VECTOR_ELT(ans, 3, ans_llf);
    sci::vector<double> eb(n, REAL(AS_NUMERIC(ans_eb)));
    sci::vector<double> ew(n, REAL(AS_NUMERIC(ans_ew)));
    double *etotal = REAL(AS_NUMERIC(ans_etotal));
    double *llf = REAL(AS_NUMERIC(ans_llf));

    *llf = mapfit::phase_erlang_estep_wtime(mixrate, shape, rate, tdat, wdat, *etotal, eb, ew);

    UNPROTECT(5);
    return ans;
  }

  SEXP phfit_herlang_estep_group(SEXP Rerl, SEXP Rdata, SEXP gdatlast) {
    int n = asInteger(getSlot(Rerl, "size"));
    sci::vector<double> mixrate(n, REAL(AS_NUMERIC(getSlot(Rerl, "mixrate"))));
    sci::vector<int> shape(n, INTEGER(AS_INTEGER(getSlot(Rerl, "shape"))));
    sci::vector<double> rate(n, REAL(AS_NUMERIC(getSlot(Rerl, "rate"))));
    int m = asInteger(getSlot(Rdata, "size"));
    sci::vector<double> tdat(m, REAL(AS_NUMERIC(getListElement(getSlot(Rdata, "data"), "time"))));
    sci::vector<int> gdat(m, INTEGER(AS_INTEGER(getListElement(getSlot(Rdata, "data"), "counts"))));
    sci::vector<int> idat(m, INTEGER(AS_INTEGER(getListElement(getSlot(Rdata, "data"), "instant"))));

    SEXP ans = PROTECT(allocVector(VECSXP, 4));
    SEXP ans_etotal = PROTECT(allocVector(REALSXP, 1));
    SEXP ans_eb = PROTECT(allocVector(REALSXP, n));
    SEXP ans_ew = PROTECT(allocVector(REALSXP, n));
    SEXP ans_llf = PROTECT(allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 0, ans_etotal);
    SET_VECTOR_ELT(ans, 1, ans_eb);
    SET_VECTOR_ELT(ans, 2, ans_ew);
    SET_VECTOR_ELT(ans, 3, ans_llf);
    sci::vector<double> eb(n, REAL(AS_NUMERIC(ans_eb)));
    sci::vector<double> ew(n, REAL(AS_NUMERIC(ans_ew)));
    double *etotal = REAL(AS_NUMERIC(ans_etotal));
    double *llf = REAL(AS_NUMERIC(ans_llf));

    *llf = mapfit::phase_erlang_estep_group_trunc(mixrate, shape, rate,
        tdat, gdat, asInteger(gdatlast), idat, *etotal, eb, ew);

    UNPROTECT(5);
    return ans;
  }

  SEXP mapfit_estep_gen_group(SEXP Rmap, SEXP Rdata, SEXP eps, SEXP ufact) {
    int n = asInteger(getSlot(Rmap, "size"));
    sci::vector<double> alpha(n, REAL(AS_NUMERIC(getSlot(Rmap, "alpha"))));
    sci::vector<double> xi(n, 1.0);
    sci::matrix<double>* D0ptr = createMatrix(getSlot(Rmap, "D0"));
    sci::matrix<double>* D1ptr = createMatrix(getSlot(Rmap, "D1"));
    sci::matrix<double>* P0ptr = sci::dnewcopy(*D0ptr);
    sci::matrix<double>* P1ptr = sci::dnewcopy(*D1ptr);
    double qv = mapblas::unif(*D0ptr, *D1ptr, *P0ptr, *P1ptr, asReal(ufact));
    int m = asInteger(getSlot(Rdata, "size"));
    sci::vector<double> tdat(m, REAL(AS_NUMERIC(getListElement(getSlot(Rdata, "data"), "time"))));
    sci::vector<int> gdat(m, INTEGER(AS_INTEGER(getListElement(getSlot(Rdata, "data"), "counts"))));
    sci::vector<int> idat(m, INTEGER(AS_INTEGER(getListElement(getSlot(Rdata, "data"), "instant"))));

    SEXP ans = PROTECT(allocVector(VECSXP, 5));
    SEXP ans_eb = PROTECT(allocVector(REALSXP, n));
    SEXP ans_ez = PROTECT(allocVector(REALSXP, n));
    SEXP ans_en0 = PROTECT(allocVector(REALSXP, D0ptr->size));
    SEXP ans_en1 = PROTECT(allocVector(REALSXP, D1ptr->size));
    SEXP ans_llf = PROTECT(allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 0, ans_eb);
    SET_VECTOR_ELT(ans, 1, ans_ez);
    SET_VECTOR_ELT(ans, 2, ans_en0);
    SET_VECTOR_ELT(ans, 3, ans_en1);
    SET_VECTOR_ELT(ans, 4, ans_llf);
    double *llf = REAL(AS_NUMERIC(ans_llf));
    sci::vector<double> eb(n, REAL(AS_NUMERIC(ans_eb)));
    sci::vector<double> ez(n, REAL(AS_NUMERIC(ans_ez)));
    sci::matrix<double>* en0ptr = sci::dnewcopy(*D0ptr, REAL(AS_NUMERIC(ans_en0)));
    sci::matrix<double>* en1ptr = sci::dnewcopy(*D1ptr, REAL(AS_NUMERIC(ans_en1)));

    // *llf = mapfit::map_estep_group(alpha, xi, *D0ptr, *D1ptr, *P0ptr, *P1ptr,
    //     qv, tdat, gdat, idat, eb, ez, *en0ptr, *en1ptr, asReal(eps));
    *llf = mapfit::map_estep_groupNA(alpha, xi, *D0ptr, *D1ptr, *P0ptr, *P1ptr,
        qv, tdat, gdat, idat, eb, ez, *en0ptr, *en1ptr, asReal(eps));

    UNPROTECT(6);

    delete D0ptr;
    delete D1ptr;
    delete P0ptr;
    delete P1ptr;
    delete en0ptr;
    delete en1ptr;

    return ans;
  }

  SEXP mapfit_hmm_erlang_estep(SEXP Rerl, SEXP Rdata) {
    int n = asInteger(getSlot(Rerl, "size"));
    sci::vector<double> alpha(n, REAL(AS_NUMERIC(getSlot(Rerl, "alpha"))));
    sci::vector<double> xi(n, 1.0);
    sci::vector<int> shape(n, INTEGER(AS_INTEGER(getSlot(Rerl, "shape"))));
    sci::vector<double> rate(n, REAL(AS_NUMERIC(getSlot(Rerl, "rate"))));
    sci::matrix<double>* Pptr = createMatrix(getSlot(Rerl, "P"));

    int m = asInteger(getSlot(Rdata, "size"));
    sci::vector<double> tdat(m, REAL(AS_NUMERIC(getListElement(getSlot(Rdata, "data"), "time"))));

    SEXP ans = PROTECT(allocVector(VECSXP, 5));
    SEXP ans_eb = PROTECT(allocVector(REALSXP, n));
    SEXP ans_en = PROTECT(allocVector(REALSXP, Pptr->size));
    SEXP ans_ew0 = PROTECT(allocVector(REALSXP, n));
    SEXP ans_ew1 = PROTECT(allocVector(REALSXP, n));
    SEXP ans_llf = PROTECT(allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 0, ans_eb);
    SET_VECTOR_ELT(ans, 1, ans_en);
    SET_VECTOR_ELT(ans, 2, ans_ew0);
    SET_VECTOR_ELT(ans, 3, ans_ew1);
    SET_VECTOR_ELT(ans, 4, ans_llf);
    double *llf = REAL(AS_NUMERIC(ans_llf));
    sci::vector<double> eb(n, REAL(AS_NUMERIC(ans_eb)));
    sci::vector<double> ew0(n, REAL(AS_NUMERIC(ans_ew0)));
    sci::vector<double> ew1(n, REAL(AS_NUMERIC(ans_ew1)));
    sci::matrix<double>* enptr = sci::dnewcopy(*Pptr, REAL(AS_NUMERIC(ans_en)));

    *llf = mapfit::hmm_erlang_estep(alpha, xi, *Pptr, shape, rate, tdat, eb, *enptr, ew0, ew1);

    UNPROTECT(6);

    delete Pptr;
    delete enptr;

    return ans;
  }

  SEXP mapfit_estep_gmmpp(SEXP Rmap, SEXP Rdata, SEXP eps, SEXP divide) {
    int n = asInteger(getSlot(Rmap, "size"));
    sci::vector<double> alpha(n, REAL(AS_NUMERIC(getSlot(Rmap, "alpha"))));
    sci::vector<double> xi(n, 1.0);
    sci::matrix<double>* D0ptr = createMatrix(getSlot(Rmap, "D0"));
    sci::matrix<double>* D1ptr = createMatrix(getSlot(Rmap, "D1"));
    int m = asInteger(getSlot(Rdata, "size"));
    sci::vector<double> tdat(m, REAL(AS_NUMERIC(getListElement(getSlot(Rdata, "data"), "time"))));
    sci::vector<int> gdat(m, INTEGER(AS_INTEGER(getListElement(getSlot(Rdata, "data"), "counts"))));
    sci::vector<int> idat(m, INTEGER(AS_INTEGER(getListElement(getSlot(Rdata, "data"), "instant"))));

    SEXP ans = PROTECT(allocVector(VECSXP, 5));
    SEXP ans_eb = PROTECT(allocVector(REALSXP, n));
    SEXP ans_ez = PROTECT(allocVector(REALSXP, n));
    SEXP ans_en0 = PROTECT(allocVector(REALSXP, D0ptr->size));
    SEXP ans_en1 = PROTECT(allocVector(REALSXP, D1ptr->size));
    SEXP ans_llf = PROTECT(allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 0, ans_eb);
    SET_VECTOR_ELT(ans, 1, ans_ez);
    SET_VECTOR_ELT(ans, 2, ans_en0);
    SET_VECTOR_ELT(ans, 3, ans_en1);
    SET_VECTOR_ELT(ans, 4, ans_llf);
    double *llf = REAL(AS_NUMERIC(ans_llf));
    sci::vector<double> eb(n, REAL(AS_NUMERIC(ans_eb)));
    sci::vector<double> ez(n, REAL(AS_NUMERIC(ans_ez)));
    sci::matrix<double>* en0ptr = sci::dnewcopy(*D0ptr, REAL(AS_NUMERIC(ans_en0)));
    sci::matrix<double>* en1ptr = sci::dnewcopy(*D1ptr, REAL(AS_NUMERIC(ans_en1)));

    *llf = mapfit::gmmpp_estep(alpha, xi, *D0ptr, *D1ptr,
        tdat, gdat, idat, eb, ez, *en0ptr, *en1ptr, asInteger(divide), asReal(eps));

    UNPROTECT(6);

    delete D0ptr;
    delete D1ptr;
    delete en0ptr;
    delete en1ptr;

    return ans;
  }
}

