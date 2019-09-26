// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// expression
NumericVector expression(NumericVector p);
RcppExport SEXP _content_clust_expression(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(expression(p));
    return rcpp_result_gen;
END_RCPP
}
// small_in_big
bool small_in_big(IntegerVector i, int small_start, int small_end, int big_start, int big_end);
RcppExport SEXP _content_clust_small_in_big(SEXP iSEXP, SEXP small_startSEXP, SEXP small_endSEXP, SEXP big_startSEXP, SEXP big_endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type small_start(small_startSEXP);
    Rcpp::traits::input_parameter< int >::type small_end(small_endSEXP);
    Rcpp::traits::input_parameter< int >::type big_start(big_startSEXP);
    Rcpp::traits::input_parameter< int >::type big_end(big_endSEXP);
    rcpp_result_gen = Rcpp::wrap(small_in_big(i, small_start, small_end, big_start, big_end));
    return rcpp_result_gen;
END_RCPP
}
// have_overlap
bool have_overlap(IntegerVector i, int start1, int end1, int start2, int end2);
RcppExport SEXP _content_clust_have_overlap(SEXP iSEXP, SEXP start1SEXP, SEXP end1SEXP, SEXP start2SEXP, SEXP end2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type start1(start1SEXP);
    Rcpp::traits::input_parameter< int >::type end1(end1SEXP);
    Rcpp::traits::input_parameter< int >::type start2(start2SEXP);
    Rcpp::traits::input_parameter< int >::type end2(end2SEXP);
    rcpp_result_gen = Rcpp::wrap(have_overlap(i, start1, end1, start2, end2));
    return rcpp_result_gen;
END_RCPP
}
// substract_sorted
IntegerVector substract_sorted(IntegerVector a, IntegerVector i, int b_start, int b_end);
RcppExport SEXP _content_clust_substract_sorted(SEXP aSEXP, SEXP iSEXP, SEXP b_startSEXP, SEXP b_endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type b_start(b_startSEXP);
    Rcpp::traits::input_parameter< int >::type b_end(b_endSEXP);
    rcpp_result_gen = Rcpp::wrap(substract_sorted(a, i, b_start, b_end));
    return rcpp_result_gen;
END_RCPP
}
// find_children
IntegerVector find_children(IntegerVector content, int parent);
RcppExport SEXP _content_clust_find_children(SEXP contentSEXP, SEXP parentSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type content(contentSEXP);
    Rcpp::traits::input_parameter< int >::type parent(parentSEXP);
    rcpp_result_gen = Rcpp::wrap(find_children(content, parent));
    return rcpp_result_gen;
END_RCPP
}
// get_cell_set
IntegerVector get_cell_set(IntegerVector i, IntegerVector p, int cells, IntegerVector geneOrder, int start);
RcppExport SEXP _content_clust_get_cell_set(SEXP iSEXP, SEXP pSEXP, SEXP cellsSEXP, SEXP geneOrderSEXP, SEXP startSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type i(iSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type cells(cellsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type geneOrder(geneOrderSEXP);
    Rcpp::traits::input_parameter< int >::type start(startSEXP);
    rcpp_result_gen = Rcpp::wrap(get_cell_set(i, p, cells, geneOrder, start));
    return rcpp_result_gen;
END_RCPP
}
// get_sep_genes
IntegerVector get_sep_genes(IntegerVector i, IntegerVector p, IntegerVector geneOrder, IntegerVector child_genes);
RcppExport SEXP _content_clust_get_sep_genes(SEXP iSEXP, SEXP pSEXP, SEXP geneOrderSEXP, SEXP child_genesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type i(iSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type geneOrder(geneOrderSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type child_genes(child_genesSEXP);
    rcpp_result_gen = Rcpp::wrap(get_sep_genes(i, p, geneOrder, child_genes));
    return rcpp_result_gen;
END_RCPP
}
// process_oriented_graph
Rcpp::List process_oriented_graph(IntegerVector i, IntegerVector p, int cells, CharacterVector geneNames, IntegerVector geneOrder, IntegerVector content, int start);
RcppExport SEXP _content_clust_process_oriented_graph(SEXP iSEXP, SEXP pSEXP, SEXP cellsSEXP, SEXP geneNamesSEXP, SEXP geneOrderSEXP, SEXP contentSEXP, SEXP startSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type i(iSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type cells(cellsSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type geneNames(geneNamesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type geneOrder(geneOrderSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type content(contentSEXP);
    Rcpp::traits::input_parameter< int >::type start(startSEXP);
    rcpp_result_gen = Rcpp::wrap(process_oriented_graph(i, p, cells, geneNames, geneOrder, content, start));
    return rcpp_result_gen;
END_RCPP
}
// content_clust
Rcpp::List content_clust(IntegerVector i, IntegerVector p, int cells, CharacterVector geneNames, IntegerVector geneOrder);
RcppExport SEXP _content_clust_content_clust(SEXP iSEXP, SEXP pSEXP, SEXP cellsSEXP, SEXP geneNamesSEXP, SEXP geneOrderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type i(iSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type cells(cellsSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type geneNames(geneNamesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type geneOrder(geneOrderSEXP);
    rcpp_result_gen = Rcpp::wrap(content_clust(i, p, cells, geneNames, geneOrder));
    return rcpp_result_gen;
END_RCPP
}
// c_hello_world
void c_hello_world();
RcppExport SEXP _content_clust_c_hello_world() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    c_hello_world();
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_content_clust_expression", (DL_FUNC) &_content_clust_expression, 1},
    {"_content_clust_small_in_big", (DL_FUNC) &_content_clust_small_in_big, 5},
    {"_content_clust_have_overlap", (DL_FUNC) &_content_clust_have_overlap, 5},
    {"_content_clust_substract_sorted", (DL_FUNC) &_content_clust_substract_sorted, 4},
    {"_content_clust_find_children", (DL_FUNC) &_content_clust_find_children, 2},
    {"_content_clust_get_cell_set", (DL_FUNC) &_content_clust_get_cell_set, 5},
    {"_content_clust_get_sep_genes", (DL_FUNC) &_content_clust_get_sep_genes, 4},
    {"_content_clust_process_oriented_graph", (DL_FUNC) &_content_clust_process_oriented_graph, 7},
    {"_content_clust_content_clust", (DL_FUNC) &_content_clust_content_clust, 5},
    {"_content_clust_c_hello_world", (DL_FUNC) &_content_clust_c_hello_world, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_content_clust(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
