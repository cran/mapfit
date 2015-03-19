// R library

#include "R_common.h"

extern "C" {

  SEXP getListElement(SEXP list, const char *str) {
    SEXP elmt = R_NilValue, names = GET_NAMES(list);
    for (int i = 0; i < length(list); i++)
      if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
        elmt = VECTOR_ELT(list, i);
        break;
      }
    return elmt;
  }

  SEXP getSlot(SEXP obj, const char *str) {
    SEXP name;
    PROTECT(name = allocVector(STRSXP, 1));
    SET_STRING_ELT(name, 0, mkChar(str));
    SEXP result = getAttrib(obj, name);
    UNPROTECT(1);
    return result;
  }

  sci::matrix<double>* createMatrix(SEXP m) {
    SEXP cls = GET_CLASS(m);
    if (strcmp(CHAR(STRING_ELT(cls, 0)), "dgeMatrix") == 0) {
      // Dense
      int *dim = INTEGER(AS_INTEGER(getSlot(m, "Dim")));
      double *value = REAL(AS_NUMERIC(getSlot(m, "x")));
      sci::dmatrix<double> *m = new sci::dmatrix<double>(dim[0], dim[1]);
      *m = value;
      return m;
    } else {
      return createSpMatrix(m);
    }
  }

  sci::spmatrix<double>* createSpMatrix(SEXP m) {
    SEXP cls = GET_CLASS(m);
    if (strcmp(CHAR(STRING_ELT(cls, 0)), "dgRMatrix") == 0) {
      // CSR
      int *dim = INTEGER(AS_INTEGER(getSlot(m, "Dim")));
      int *rowptr = INTEGER(AS_INTEGER(getSlot(m, "p")));
      int *colind = INTEGER(AS_INTEGER(getSlot(m, "j")));
      SEXP Rvalue = getSlot(m, "x");
      sci::csrmatrix<double> *m = new sci::csrmatrix<double>(dim[0], dim[1],
        LENGTH(Rvalue), REAL(AS_NUMERIC(Rvalue)));
      m->rowptr = rowptr;
      m->rowptr += 1;
      m->colind = colind;
      m->colind += 1;
      return m;
    } else if (strcmp(CHAR(STRING_ELT(cls, 0)), "dgCMatrix") == 0) {
      // CSC
      int *dim = INTEGER(AS_INTEGER(getSlot(m, "Dim")));
      int *colptr = INTEGER(AS_INTEGER(getSlot(m, "p")));
      int *rowind = INTEGER(AS_INTEGER(getSlot(m, "i")));
      SEXP Rvalue = getSlot(m, "x");
      sci::cscmatrix<double> *m = new sci::cscmatrix<double>(dim[0], dim[1],
        LENGTH(Rvalue), REAL(AS_NUMERIC(Rvalue)));
      m->colptr = colptr;
      m->colptr += 1;
      m->rowind = rowind;
      m->rowind += 1;
      return m;
    } else if (strcmp(CHAR(STRING_ELT(cls, 0)), "dgTMatrix") == 0) {
      // COO
      int *dim = INTEGER(AS_INTEGER(getSlot(m, "Dim")));
      int *rowind = INTEGER(AS_INTEGER(getSlot(m, "i")));
      int *colind = INTEGER(AS_INTEGER(getSlot(m, "j")));
      SEXP Rvalue = getSlot(m, "x");
      sci::coomatrix<double> *m = new sci::coomatrix<double>(dim[0], dim[1],
        LENGTH(Rvalue), REAL(AS_NUMERIC(Rvalue)));
      m->rowind = rowind;
      m->rowind += 1;
      m->colind = colind;
      m->colind += 1;
      return m;
    } else {
      Rf_error("Cannot create C++matrix class from %s\n", CHAR(STRING_ELT(cls, 0)));
      return NULL;
    }
  }

}

sci::matrix<double>* ccMatrix(const sci::matrix<double>& m, SEXP& v) {
  switch (m.type()) {
  case (DENSE): {
    const sci::dmatrix<double>& mm(dynamic_cast<const sci::dmatrix<double>&>(m));
    v = PROTECT(allocVector(REALSXP, mm.nrow * mm.ncol));
    sci::dmatrix<double>* retm = new sci::dmatrix<double>(mm.nrow, mm.ncol, REAL(AS_NUMERIC(v)));
    return retm;
  }
  case (CSR): {
    const sci::csrmatrix<double>& mm(dynamic_cast<const sci::csrmatrix<double>&>(m));
    v = PROTECT(allocVector(REALSXP, mm.nnz));
    sci::csrmatrix<double>* retm = new sci::csrmatrix<double>(mm.nrow, mm.ncol, mm.nnz, REAL(AS_NUMERIC(v)));
    retm->rowptr = mm.rowptr;
    retm->rowptr += 1;
    retm->colind = mm.colind;
    retm->colind += 1;
    return retm;
  }
  case (CSC): {
    const sci::cscmatrix<double>& mm(dynamic_cast<const sci::cscmatrix<double>&>(m));
    v = PROTECT(allocVector(REALSXP, mm.nnz));
    sci::cscmatrix<double>* retm = new sci::cscmatrix<double>(mm.nrow, mm.ncol, mm.nnz, REAL(AS_NUMERIC(v)));
    retm->colptr = mm.colptr;
    retm->colptr += 1;
    retm->rowind = mm.rowind;
    retm->rowind += 1;
    return retm;
  }
  case (COO): {
    const sci::coomatrix<double>& mm(dynamic_cast<const sci::coomatrix<double>&>(m));
    v = PROTECT(allocVector(REALSXP, mm.nnz));
    sci::coomatrix<double>* retm = new sci::coomatrix<double>(mm.nrow, mm.ncol, mm.nnz, REAL(AS_NUMERIC(v)));
    retm->rowind = mm.rowind;
    retm->rowind += 1;
    retm->colind = mm.colind;
    retm->colind += 1;
    return retm;
  }
  default:
    throw;
  }      
}

