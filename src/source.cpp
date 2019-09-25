/*
Copyright (C) 2019  Oleg Demiancheko

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>
*/

#include <Rcpp.h>

using namespace Rcpp;


// [[Rcpp::export]]
NumericVector expression(NumericVector p) {
  NumericVector expression = no_init(p.length() - 1);
  for(int index = 0; index < (p.length() - 1); index++) {
    expression[index] = p[index + 1] - p[index];
  }
  
  return expression;
}

// [[Rcpp::export]]
bool small_in_big(IntegerVector i, 
                  int small_start, int small_end, 
                  int big_start, int big_end) {
  bool ret = true;
  int big = big_start;
  int small = small_start;
  while (small < small_end) {
    if(i[small] == i[big]) {
      small++;
      big++;
      continue;
    } else {
      big++;
    }
    
    if ((i[big] > i[small]) or (big >= big_end)) {
      ret=false;
      break;
    }
  }
  
  return ret;
}


// [[Rcpp::export]]
bool have_overlap(IntegerVector i, 
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
IntegerVector substract_sorted(IntegerVector a, IntegerVector i, int b_start, int b_end) {
  int* result = (int*) malloc((a.length())*sizeof(int));
  int result_index = 0;
  int a_index = 0;
  for (int i_index = b_start; i_index < b_end; i_index++) {
    while(a[a_index] <= i[i_index]) {
      if(a[a_index] != i[i_index]) {
        result[result_index] = a[a_index];
        result_index++;
      }
      a_index++;
    }
  }
  
  while(a_index < a.length()) {
    result[result_index] = a[a_index];
    a_index++;
    result_index++;
  }
  
  return IntegerVector(result, &(result[result_index]));
}

// [[Rcpp::export]]
IntegerVector find_children(IntegerVector content, int parent) {
	IntegerVector ret;
	for (int index = parent+1; index < content.length(); index++) {
		if (parent == content[index]) {
			ret.push_back(index);
		}
	}

	return ret;
}

// [[Rcpp::export]]
IntegerVector get_cell_set(IntegerVector i, IntegerVector p, int cells, IntegerVector geneOrder, int start) {
	IntegerVector ret = NULL;

	if (start != -1) {
		ret = i[Range(p[geneOrder[start]], p[geneOrder[start] + 1] - 1)];
	} else {
		ret = Range(0, cells - 1);
	}

	return ret;
}

// [[Rcpp::export]]
IntegerVector get_sep_genes(IntegerVector i, IntegerVector p, IntegerVector geneOrder, IntegerVector child_genes) {
	IntegerVector ret;
	for(int i_ind = 0; i_ind < child_genes.length(); i_ind++) {
		bool sep = true;
		for(int j_ind = 0; j_ind < child_genes.length(); j_ind++) {
			if (i_ind == j_ind) {
				continue;
			}

			if (have_overlap(i, 
					 p[geneOrder[child_genes[i_ind]]], 
					 p[geneOrder[child_genes[i_ind]] + 1], 
					 p[geneOrder[child_genes[j_ind]]], 
					 p[geneOrder[child_genes[j_ind]] + 1])) {
				sep = false;
				break;
			}
		}
		if (sep) {
			ret.push_back(child_genes[i_ind]);
		}
	}

	return ret;
}

//[[Rcpp::export]]
Rcpp::List process_oriented_graph(IntegerVector i, IntegerVector p, int cells, 
                                  CharacterVector geneNames, IntegerVector geneOrder,
                                  IntegerVector content, int start) {
	printf("starting to process suborgraph\n");
	IntegerVector child_vertices = find_children(content, start);
  
	IntegerVector unclussified = get_cell_set(i, p, cells, geneOrder, start);
  
	printf("the num of cells of is %d\n", p[geneOrder[start] + 1] - p[geneOrder[start]]);

	IntegerVector sep_genes = get_sep_genes(i, p, geneOrder, child_vertices);
	
	Rcpp::List outlist;
	Rcpp::StringVector list_names;

	for (int sep_ind = 0; sep_ind < sep_genes.length(); sep_ind++) {
		printf("captain, gene %d seems to be separate from brothers\n", geneOrder[sep_genes[sep_ind]]);
		Rcpp::List child_sublist = process_oriented_graph(i, p, cells, geneNames, geneOrder, content, sep_genes[sep_ind]);
		Rcpp::String gene_name = geneNames[geneOrder[sep_genes[sep_ind]]];
		list_names.push_back(gene_name);
		if (child_sublist.length() > 1) {
			outlist.push_back(child_sublist);  
		} else {
			outlist.push_back(child_sublist[0]);
		}

		unclussified = substract_sorted(unclussified, i, 
                                      		p[geneOrder[sep_genes[sep_ind]]], 
                                      		p[geneOrder[sep_genes[sep_ind]] + 1]);
	}
  
	if (unclussified.length() > 0) {
		outlist.push_back(unclussified);
		list_names.push_back("unclussified");
	}
  
	printf("suborgraph processing ended\n");
	outlist.names() = list_names;
	return outlist;
}

//' @export
// [[Rcpp::export]]
Rcpp::List content_clust(IntegerVector i, IntegerVector p, int cells,
                         CharacterVector geneNames, IntegerVector geneOrder) {
  IntegerVector content (geneOrder.length(), -1);
  for (int bigger_gene = 0; bigger_gene < geneOrder.length(); bigger_gene++) {
    for(int smaller_gene = bigger_gene+1; smaller_gene < geneOrder.length(); smaller_gene++) {
      if (small_in_big(i, p[geneOrder[smaller_gene]], p[geneOrder[smaller_gene] + 1],
                       p[geneOrder[bigger_gene]], p[geneOrder[bigger_gene] + 1])) {
        printf("one in another\n");
        content[smaller_gene] = bigger_gene;
      }
    }
  }
  
  return process_oriented_graph(i, p, cells, geneNames, geneOrder, content, -1);
}

//' @export
// [[Rcpp::export]]
void c_hello_world() {
	printf("hello world, it's c, are you ok?");
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
