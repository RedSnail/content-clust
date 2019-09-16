#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]


//

// [[Rcpp::export]]
NumericVector SparseRowVar2(Eigen::SparseMatrix<double> mat,
                            NumericVector mu, 
                            int power){
  mat = mat.transpose();
  NumericVector allVars = no_init(mat.cols());
  for (int k=0; k<mat.outerSize(); ++k){
    double colSum = 0;
    int nZero = mat.rows();
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it) {
      nZero -= 1;
      colSum += pow(it.value() - mu[k], power);
    }
    colSum += pow(mu[k], power) * nZero;
    allVars[k] = colSum / (mat.rows() - 1);
  }
  return(allVars);
}


// [[Rcpp::export]]
NumericVector SparseRowVarStd(Eigen::SparseMatrix<double> mat,
                              NumericVector mu,
                              NumericVector sd,
                              double vmax){
  mat = mat.transpose();
  NumericVector allVars(mat.cols());
  for (int k=0; k<mat.outerSize(); ++k){
    if (sd[k] == 0) continue;
    double colSum = 0;
    int nZero = mat.rows();
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
    {
      nZero -= 1;
      colSum += pow(std::min(vmax, (it.value() - mu[k]) / sd[k]), 2);
    }
    colSum += pow((0 - mu[k]) / sd[k], 2) * nZero;
    allVars[k] = colSum / (mat.rows() - 1);
  }
  return(allVars);
}

// [[Rcpp::export]]
NumericVector SparseRowVarStd2(Eigen::SparseMatrix<double> mat,
                               NumericVector mu,
                               NumericVector var,
                               double varmax, 
                               int power){
  mat = mat.transpose();
  NumericVector allVars(mat.cols());
  for (int k=0; k<mat.outerSize(); ++k){
    if (var[k] == 0) continue;
    double colSum = 0;
    int nZero = mat.rows();
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
    {
      nZero -= 1;
      colSum += std::min(varmax*var[k], pow(it.value() - mu[k], power));
    }
    colSum += pow(mu[k], power) * nZero;
    allVars[k] = colSum/((mat.rows() - 1)*var[k]);
  }
  return(allVars);
}

// [[Rcpp::export]]
NumericVector Expression(NumericVector p) {
  NumericVector expression = no_init(p.length() - 1);
  for(int index = 0; index < (p.length() - 1); index++) {
    expression[index] = p[index + 1] - p[index];
  }
  
  return expression;
}

// [[Rcpp::export]]
bool OneInAnother(IntegerVector i, 
                  int start1, int end1, 
                  int start2, int end2) {
  bool ret = true;
  int big = start2;
  int small = start1;
  while (small < end1) {
    if(i[small] == i[big]) {
      small++;
      big++;
      continue;
    } else {
      big++;
    }
    
    if ((i[big] > i[small]) or (big >= end2)) {
      ret=false;
      break;
    }
  }
  
  return ret;
}


// [[Rcpp::export]]
bool HasOverlap(IntegerVector i, 
                int start1, int end1, 
                int start2, int end2) {
  bool ret = false;
  while(start1 < end1 and start2 < end2) {
    if(i[start1] == i[start2]) {
      ret = true;
      break;
    }
    
    if(i[start1] > i[start2]) {
      start2++;
    } else {
      start1++;
    }
  }
  
  return ret;
}


// [[Rcpp::export]]
IntegerVector SubstractSorted(IntegerVector a, IntegerVector b) {
  int* result = (int*) malloc(a.length()*sizeof(int));
  int result_index = 0;
  int j = 0;
  for (int i = 0; i < b.length(); i++) {
    while(a[j] <= b[i]) {
      if(a[j] != b[i]) {
        result[result_index] = a[j];
        result_index++;
      }
      j++;
    }
  }
  
  while(j < a.length()) {
    result[result_index] = a[j];
    j++;
    result_index++;
  }
  
  return IntegerVector(result, &(result[result_index]));
}

Rcpp::List ProcessOrientedGraph(IntegerVector i, IntegerVector p, int cells, 
                                CharacterVector geneNames, IntegerVector geneOrder,
                                IntegerVector content, int start) {
  printf("starting to process suborgraph\n");
  Rcpp::List outlist;
  Rcpp::StringVector list_names;
  IntegerVector child_vertices;
  for (int ind = start; ind < content.length(); ind++) {
    if (start == content[ind]) {
      child_vertices.push_back(ind);
    }
  }
  
  IntegerVector unclussified = NULL;
  
  if (start != -1) {
    unclussified = i[Range(p[geneOrder[start]], p[geneOrder[start] + 1] - 1)];
  } else {
    unclussified = Range(0, cells - 1);
  }
  
  printf("the num of cells of is %d\n", p[geneOrder[start] + 1] - p[geneOrder[start]] + 1);
  
  for(int index = 0; index < child_vertices.length(); index++) {
    bool separate = true;
    
    for(int j = 0; j < child_vertices.length(); j++) {
      if(index==j) {
        continue;
      }
      
      if (HasOverlap(i, 
                     p[geneOrder[(int) child_vertices[index]]],
                     p[geneOrder[(int) child_vertices[index]] + 1], 
                     p[geneOrder[(int) child_vertices[j]]],
                     p[geneOrder[(int) child_vertices[j]] + 1])) {
        separate = false;
        break;
      }
    }
    
    if (separate) {
      printf("captain, gene %d seems to be separate from brothers\n", geneOrder[(int) child_vertices[index]]);
      Rcpp::List child_sublist = ProcessOrientedGraph(i, p, cells, geneNames, geneOrder, content, (int) child_vertices[index]);
      Rcpp::String gene_name = geneNames[geneOrder[(int) child_vertices[index]]];
      list_names.push_back(gene_name);
      if (child_sublist.length() > 1) {
        outlist.push_back(child_sublist);  
      } else {
        outlist.push_back(child_sublist[0]);
      }
      if(start == -1) {
        printf("substracting\n");
        for(int p_index = p[geneOrder[(int) child_vertices[index]]]; p_index < p[geneOrder[(int) child_vertices[index]] + 1]; p_index++) {
          printf("    %d\n", i[p_index]);
        }
        
      }
      unclussified = SubstractSorted(unclussified, 
                                     i[Range(p[geneOrder[(int) child_vertices[index]]],
                                             p[geneOrder[(int) child_vertices[index]] + 1] - 1)]);
    }
  }
  
  if (unclussified.length() > 0) {
    outlist.push_back(unclussified);
    list_names.push_back("unclussified");
  }
  
  printf("suborgraph processing ended\n");
  outlist.names() = list_names;
  return outlist;
}

// [[Rcpp::export]]
Rcpp::List ContentClust(IntegerVector i, IntegerVector p, int cells,
                        CharacterVector geneNames, IntegerVector geneOrder) {
  IntegerVector content (geneOrder.length(), -1);
  for (int bigger_gene = 0; bigger_gene < geneOrder.length(); bigger_gene++) {
    for(int smaller_gene = bigger_gene+1; smaller_gene < geneOrder.length(); smaller_gene++) {
      if (OneInAnother(i, p[geneOrder[smaller_gene]], p[geneOrder[smaller_gene] + 1],
                       p[geneOrder[bigger_gene]], p[geneOrder[bigger_gene] + 1])) {
        printf("one in another\n");
        content[smaller_gene] = bigger_gene;
      }
    }
  }
  
  return ProcessOrientedGraph(i, p, cells, geneNames, geneOrder, content, -1);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
