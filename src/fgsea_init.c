#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>

SEXP _signatureSearch_calcGseaStatCumulativeBatch(SEXP statsSEXP, SEXP gseaParamSEXP, SEXP pathwayScoresSEXP, SEXP pathwaysSizesSEXP, SEXP iterationsSEXP, SEXP seedSEXP);
SEXP _signatureSearch_calcGseaStatCumulative(SEXP statsSEXP, SEXP selectedStatsSEXP, SEXP gseaParamSEXP);
SEXP _signatureSearch_calcGseaStatBatchCpp(SEXP statsSEXP, SEXP selectedGenesSEXP, SEXP geneRanksSEXP);

R_CallMethodDef callMethods[]  = {
  {"_signatureSearch_calcGseaStatCumulativeBatch", (DL_FUNC) &_signatureSearch_calcGseaStatCumulativeBatch, 6},
  {"_signatureSearch_calcGseaStatCumulative", (DL_FUNC) &_signatureSearch_calcGseaStatCumulative, 3},
  {"_signatureSearch_calcGseaStatBatchCpp", (DL_FUNC) &_signatureSearch_calcGseaStatBatchCpp, 3},
  {NULL, NULL, 0}
};

void R_init_fgsea(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

