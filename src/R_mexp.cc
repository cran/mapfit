// R library

#include "R_common.h"

#include "poisson.h"

extern "C" {

  SEXP ctmc_st(SEXP RQ, SEXP Rite, SEXP Reps) {
    sci::matrix<double>* m = createMatrix(RQ);
    switch (m->type()) {
    case (DENSE): {
        const sci::dmatrix<double>& mm(dynamic_cast<const sci::dmatrix<double>&>(*m));
        SEXP v = PROTECT(allocVector(REALSXP, mm.nrow));
        sci::vector<double> x(mm.nrow, REAL(AS_NUMERIC(v)));
        mexp::ctmc_gth(mm, x);
        UNPROTECT(1);
        delete m;
        return v;
    }
    case (CSR): {
        const sci::csrmatrix<double>& mm(dynamic_cast<const sci::csrmatrix<double>&>(*m));
        SEXP v = PROTECT(allocVector(REALSXP, mm.nrow));
        sci::vector<double> x(mm.nrow, 1.0/mm.nrow);
        sci::vector<double> y(mm.nrow, REAL(AS_NUMERIC(v)));
        int info;
        info = mexp::ctmc_gs(mm, x, y, asInteger(Rite), asReal(Reps));
        if (info == -1) {
            Rf_warning("CTMCST: GS did not converge.");
        }
        UNPROTECT(1);
        delete m;
        return v;
    }
    case (CSC): {
        const sci::cscmatrix<double>& mm(dynamic_cast<const sci::cscmatrix<double>&>(*m));
        SEXP v = PROTECT(allocVector(REALSXP, mm.nrow));
        sci::vector<double> x(mm.nrow, 1.0/mm.nrow);
        sci::vector<double> y(mm.nrow, REAL(AS_NUMERIC(v)));
        int info;
        info = mexp::ctmc_gs(mm, x, y, asInteger(Rite), asReal(Reps));
        if (info == -1) {
            Rf_warning("CTMCST: GS did not converge.");
        }
        UNPROTECT(1);
        delete m;
        return v;
    }
    case (COO): {
        const sci::coomatrix<double>& mm(dynamic_cast<const sci::coomatrix<double>&>(*m));
        SEXP v = PROTECT(allocVector(REALSXP, mm.nrow));
        sci::vector<double> x(mm.nrow, 1.0/mm.nrow);
        sci::vector<double> y(mm.nrow, REAL(AS_NUMERIC(v)));
        int info;
        info = mexp::ctmc_gs(mm, x, y, asInteger(Rite), asReal(Reps));
        if (info == -1) {
            Rf_warning("CTMCST: GS did not converge.");
        }
        UNPROTECT(1);
        delete m;
        return v;
    }
    default:
        delete m;
        throw;
    }
  }

  SEXP marsolve(SEXP trans, SEXP Ralpha, SEXP RQ, SEXP Rx, SEXP Rite, SEXP Reps) {
    sci::matrix<double>* m = createMatrix(RQ);
    switch (m->type()) {
    case (DENSE): {
        const sci::dmatrix<double>& mm(dynamic_cast<const sci::dmatrix<double>&>(*m));
        sci::vector<double> x(mm.nrow, REAL(AS_NUMERIC(Rx)));
        SEXP v = PROTECT(allocVector(REALSXP, mm.nrow));
        sci::vector<double> y(mm.nrow, REAL(AS_NUMERIC(v)));
        if (asLogical(trans)) {
            sci::dgesv(sci::mat::T, asReal(Ralpha), mm, x, y);
        } else {
            sci::dgesv(sci::mat::N, asReal(Ralpha), mm, x, y);
        }
        UNPROTECT(1);
        delete m;
        return v;
    }
    case (CSR): {
        const sci::csrmatrix<double>& mm(dynamic_cast<const sci::csrmatrix<double>&>(*m));
        sci::vector<double> x(mm.nrow, REAL(AS_NUMERIC(Rx)));
        SEXP v = PROTECT(allocVector(REALSXP, mm.nrow));
        sci::vector<double> y(mm.nrow, REAL(AS_NUMERIC(v)));
        int info;
        if (asLogical(trans)) {
            info = sci::dgssolve(sci::mat::T, asReal(Ralpha), mm, x, y, asInteger(Rite), asReal(Reps));
        } else {
            info = sci::dgssolve(sci::mat::N, asReal(Ralpha), mm, x, y, asInteger(Rite), asReal(Reps));
        }
        if (info == -1) {
            Rf_warning("MSOLVE: GS did not converge.");
        }
        UNPROTECT(1);
        delete m;
        return v;
    }
    case (CSC): {
        const sci::cscmatrix<double>& mm(dynamic_cast<const sci::cscmatrix<double>&>(*m));
        sci::vector<double> x(mm.nrow, REAL(AS_NUMERIC(Rx)));
        SEXP v = PROTECT(allocVector(REALSXP, mm.nrow));
        sci::vector<double> y(mm.nrow, REAL(AS_NUMERIC(v)));
        int info;
        if (asLogical(trans)) {
            info = sci::dgssolve(sci::mat::T, asReal(Ralpha), mm, x, y, asInteger(Rite), asReal(Reps));
        } else {
            info = sci::dgssolve(sci::mat::N, asReal(Ralpha), mm, x, y, asInteger(Rite), asReal(Reps));
        }
        if (info == -1) {
            Rf_warning("MSOLVE: GS did not converge.");
        }
        UNPROTECT(1);
        delete m;
        return v;
    }
    case (COO): {
        const sci::coomatrix<double>& mm(dynamic_cast<const sci::coomatrix<double>&>(*m));
        sci::vector<double> x(mm.nrow, REAL(AS_NUMERIC(Rx)));
        SEXP v = PROTECT(allocVector(REALSXP, mm.nrow));
        sci::vector<double> y(mm.nrow, REAL(AS_NUMERIC(v)));
        int info;
        if (asLogical(trans)) {
            info = sci::dgssolve(sci::mat::T, asReal(Ralpha), mm, x, y, asInteger(Rite), asReal(Reps));
        } else {
            info = sci::dgssolve(sci::mat::N, asReal(Ralpha), mm, x, y, asInteger(Rite), asReal(Reps));
        }
        if (info == -1) {
            Rf_warning("MSOLVE: GS did not converge.");
        }
        UNPROTECT(1);
        delete m;
        return v;
    }
    default:
        delete m;
        throw;
    }
  }

  SEXP mpow_mat(SEXP RMA, SEXP Rm) {
    sci::matrix<double>* m = createMatrix(RMA);
    switch (m->type()) {
    case (DENSE): {
        const sci::dmatrix<double>& mm(dynamic_cast<const sci::dmatrix<double>&>(*m));
        int n = mm.nrow;
        SEXP ans = PROTECT(allocMatrix(REALSXP, n, n));
        sci::dmatrix<double> ME(n, n, REAL(ans));
        mexp::mpow(mm, ME, asInteger(Rm));
        UNPROTECT(1);
        delete m;
        return ans;
    }
    case (CSR): {
//        Rf_warning("mpow uses a dense matrix only.\n");
        const sci::dmatrix<double>& mm(dynamic_cast<const sci::csrmatrix<double>&>(*m));
        int n = mm.nrow;
        SEXP ans = PROTECT(allocMatrix(REALSXP, n, n));
        sci::dmatrix<double> ME(n, n, REAL(ans));
        mexp::mpow(mm, ME, asInteger(Rm));
        UNPROTECT(1);
        delete m;
        return ans;
    }
    case (CSC): {
//        Rf_warning("mpow uses a dense matrix only.\n");
        const sci::dmatrix<double>& mm(dynamic_cast<const sci::cscmatrix<double>&>(*m));
        int n = mm.nrow;
        SEXP ans = PROTECT(allocMatrix(REALSXP, n, n));
        sci::dmatrix<double> ME(n, n, REAL(ans));
        mexp::mpow(mm, ME, asInteger(Rm));
        UNPROTECT(1);
        delete m;
        return ans;
    }
    case (COO): {
//        Rf_warning("mpow uses a dense matrix only.\n");
        const sci::dmatrix<double>& mm(dynamic_cast<const sci::coomatrix<double>&>(*m));
        int n = mm.nrow;
        SEXP ans = PROTECT(allocMatrix(REALSXP, n, n));
        sci::dmatrix<double> ME(n, n, REAL(ans));
        mexp::mpow(mm, ME, asInteger(Rm));
        UNPROTECT(1);
        delete m;
        return ans;
    }
    default:
        delete m;
        throw;
    }
  }

  SEXP mexp_pade(SEXP Rn, SEXP RMA, SEXP eps) {
    int n = asInteger(Rn);
    sci::matrix<double>* MA = createMatrix(RMA);
    if (MA->type() != DENSE) {
        delete MA;
        Rf_error("Pade uses a dense matrix only.\n");
    }
    SEXP ans = PROTECT(allocMatrix(REALSXP, n, n));
    sci::dmatrix<double> ME(n, n, REAL(ans));
    mexp::mexp_pade(*MA, ME, asReal(eps));
    UNPROTECT(1);
    delete MA;
    return ans;
  }

  SEXP mexp_unif(SEXP Rn, SEXP RMA, SEXP Rt, SEXP eps, SEXP ufact, SEXP atol) {
    int n = asInteger(Rn);
    double t = asReal(Rt);
    sci::matrix<double>* MA = createMatrix(RMA);
    if (MA->type() != DENSE) {
        delete MA;
        Rf_error("mexp_unif uses a dense matrix only.\n");
    }
    sci::dmatrix<double> MP(n, n);
    SEXP ans = PROTECT(allocMatrix(REALSXP, n, n));
    sci::dmatrix<double> y(n, n, REAL(ans));
    sci::dmatrix<double> x(n, n);
    x.diag() = 1.0;

    double qv = mexp::unif(*MA, MP, asReal(ufact));
    int nmax = pois::rightbound(qv*t, asReal(eps));
    sci::range poi_range(0,nmax);
    sci::vector<double> poi(poi_range.size());
    double weight = pois::pmf(qv*t, poi_range, poi);
    mexp::mexp_unif(MP, qv, poi_range, poi, weight, x, y, asReal(atol));

    UNPROTECT(2);
    delete MA;
    return ans;
  }

  SEXP mexp_unifvec(SEXP trans, SEXP Rn, SEXP RMA, SEXP Rx, SEXP Rt,
    SEXP eps, SEXP ufact, SEXP atol) {
    int n = asInteger(Rn);
    double t = asReal(Rt);
    sci::matrix<double>* MA = createMatrix(RMA);
    sci::matrix<double>* MP = sci::dnewcopy(*MA);
    sci::vector<double> x(n, REAL(PROTECT(AS_NUMERIC(Rx))));
    SEXP ans = PROTECT(allocVector(REALSXP, n));
    sci::vector<double> y(n, REAL(ans));
    double qv = mexp::unif(*MA, *MP, asReal(ufact));
    int nmax = pois::rightbound(qv*t, asReal(eps));
    sci::range poi_range(0,nmax);
    sci::vector<double> poi(poi_range.size());
    double weight = pois::pmf(qv*t, poi_range, poi);
    if (asLogical(trans)) {
        mexp::mexp_unifvec(sci::mat::T, *MP, qv, poi_range, poi, weight, x, y, asReal(atol));
    } else {
        mexp::mexp_unifvec(sci::mat::N, *MP, qv, poi_range, poi, weight, x, y, asReal(atol));
    }
    UNPROTECT(2);
    delete MA;
    delete MP;
    return ans;
  }

  SEXP mexp_unifseq(SEXP trans, SEXP Rn, SEXP RMA, SEXP Rx, SEXP Rt,
    SEXP eps, SEXP ufact, SEXP atol) {
    int n = asInteger(Rn);
    int m = LENGTH(Rt);
    sci::matrix<double>* MA = createMatrix(RMA);
    sci::matrix<double>* MP = sci::dnewcopy(*MA);
    sci::vector<double> t(m, REAL(PROTECT(AS_NUMERIC(Rt))));
    sci::vector<double> x(n);
    x = REAL(AS_NUMERIC(Rx));

    SEXP ans = PROTECT(allocMatrix(REALSXP, n, m+1));
    sci::array< sci::vector<double> > y(m+1, sci::vector<double>(n, REAL(ans)));
    for (int k=0; k<=m; k++) {
        y[k].ptr = REAL(ans) + k*n;
    }
    double qv = mexp::unif(*MA, *MP, asReal(ufact));

    double tmax = sci::dmax(t);
    int nmax = pois::rightbound(qv*tmax, asReal(eps));
    sci::range poi_range(0,nmax);
    sci::range rows(1,n);
    sci::vector<double> poi(poi_range.size());
    double weight;
    y[0] = x;
    if (asLogical(trans)) {
        for (int k=1; k<=m; k++) {
            nmax = pois::rightbound(qv*t(k), asReal(eps));
            poi_range = sci::range(0, nmax);
            weight = pois::pmf(qv*t(k), poi_range, poi);
            x = mexp::mexp_unifvec(sci::mat::T, *MP, qv, poi_range, poi, weight,
                x, y[k], asReal(atol));
        }
    } else {
        for (int k=1; k<=m; k++) {
            nmax = pois::rightbound(qv*t(k), asReal(eps));
            poi_range = sci::range(0, nmax);
            weight = pois::pmf(qv*t(k), poi_range, poi);
            x = mexp::mexp_unifvec(sci::mat::N, *MP, qv, poi_range, poi, weight,
                x, y[k], asReal(atol));
        }
    }
    UNPROTECT(2);
    delete MA;
    delete MP;
    return ans;
  }

  SEXP mexp_kryvec(SEXP trans, SEXP Rn, SEXP RMA, SEXP Rx, SEXP Rt,
    SEXP Rksub, SEXP ite, SEXP tol, SEXP eps) {
    int n = asInteger(Rn);
    double t = asReal(Rt);
    int ksub = asInteger(Rksub);
    if (ksub > n) {
        ksub = n;
    }
    sci::spmatrix<double>* MA = createSpMatrix(RMA);
    sci::vector<double> x(n, REAL(PROTECT(AS_NUMERIC(Rx))));
    SEXP ans = PROTECT(allocVector(VECSXP, 2));
    SEXP ans1 = PROTECT(allocVector(REALSXP, n));
    sci::vector<double> y(n, REAL(ans1));
    SEXP ans2 = PROTECT(allocVector(REALSXP, 1));
    double *err = REAL(AS_NUMERIC(ans2));
    SET_VECTOR_ELT(ans, 0, ans1);
    SET_VECTOR_ELT(ans, 1, ans2);
    if (asLogical(trans)) {
        mexp::mexp_krylov_pade(sci::mat::T, *MA, t, x, y,
            ksub, *err, asInteger(ite), asReal(tol), asReal(eps));
    } else {
        mexp::mexp_krylov_pade(sci::mat::N, *MA, t, x, y,
            ksub, *err, asInteger(ite), asReal(tol), asReal(eps));
    }
    UNPROTECT(4);
    delete MA;
    return ans;
  }

}

