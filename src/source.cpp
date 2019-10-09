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
Rcpp::String concat(Rcpp::CharacterVector vec, Rcpp::String separator)  {
  std::vector<std::string> c_vec = Rcpp::as<std::vector<std::string>>(vec);
  std::string c_sep = separator.get_cstring();
  std::stringstream ret;
  for (int i = 0; i != (vec.length() - 1); i++) {
    ret << c_vec[i] << c_sep;
  }
  
  ret << c_vec[vec.length() - 1];
  
  return Rcpp::wrap(ret.str());
}


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
IntegerVector substract_sorted(IntegerVector a, IntegerVector i, IntegerVector b_starts, IntegerVector b_ends) {
	IntegerVector result;
  IntegerVector cur_vals = i[b_starts];
  int over = 0;
	int a_index = 0;
	while (over != b_starts.length()) {
	  int min_val = min(cur_vals);
	  for (int cur_val_i = 0; cur_val_i < cur_vals.length(); cur_val_i++) {
	    if (cur_vals[cur_val_i] == min_val) {
	      b_starts[cur_val_i] += 1;
	      if (b_starts[cur_val_i] == b_ends[cur_val_i]) {
	        cur_vals[cur_val_i] = 2147483647;
	        over++;
	      } else {
	        cur_vals[cur_val_i] = i[b_starts[cur_val_i]];
	      }
	    }
	  }
	  while ((a[a_index] <= min_val) && (a_index < a.length())) {
	    if (a[a_index] != min_val) {
	      result.push_back(a[a_index]);
	    }
	    a_index++;
	  }
	}
  
	while(a_index < a.length()) {
		result.push_back(a[a_index]);
		a_index++;
	}
  
	return result;
}

// [[Rcpp::export]]
IntegerVector find_children(IntegerVector content, IntegerVector parents) {
	IntegerVector ret;
	for (int index = min(parents)+1; index < content.length(); index++) {
	  for (int par_i = 0; par_i != parents.length(); par_i++) {
	    if (parents[par_i] == content[index]) {
	      ret.push_back(index);
	      break;
	    } 
	  }
	}

	return ret;
}

// [[Rcpp::export]]
IntegerVector get_cell_set(IntegerVector i, IntegerVector p, int cells, IntegerVector geneOrder, IntegerVector starts) {
  for (int ind = 0; ind != starts.length(); ind++) {
    if (starts[ind] == -1) {
      return Range(0, cells - 1);
    }
  }
  
  IntegerVector real_starts = geneOrder[starts];
  IntegerVector cur_ind = p[real_starts];
  IntegerVector ends = p[real_starts + 1];
  
  IntegerVector cur_vals = i[cur_ind];
  IntegerVector ret;
  int over = 0;
  
  while(over != starts.length()) {
    int min_val = min(cur_vals);
    
    for (int ind = 0; ind != starts.length(); ind++) {
      if (cur_vals[ind] == min_val) {
        cur_ind[ind] = cur_ind[ind] + 1;
        if (cur_ind[ind] == ends[ind]) {
          cur_vals[ind] = 2147483647;
          over++;
        } else {
          cur_vals[ind] = i[cur_ind[ind]];
        }
      }
    }
    
    ret.push_back(min_val);
  }
  
  // to optimise min search, it is possible to make a sorted vector of cur vals and modify it every time we modify curVals

  return ret;

}

// [[Rcpp::export]]
bool is_in(int x, IntegerVector a) {
  for(IntegerVector::iterator i = a.begin(); i != a.end(); i++) {
    if ((*i) == x) {
      return true;
    }
  }
  
  return false;
}

// [[Rcpp::export]]
Rcpp::List get_sep_genes(IntegerVector i, IntegerVector p, IntegerVector geneOrder, IntegerVector child_genes) {
  Rcpp::List outlist;
  
  while(child_genes.length() > 0) {
    IntegerVector conn_comp = {child_genes[0]};
    int last_unchecked_index = 0;
    child_genes.erase(child_genes.begin());
    while(last_unchecked_index != conn_comp.length()) {
      for(IntegerVector::iterator gene = child_genes.begin(); gene != child_genes.end();) {
        if(have_overlap(i, 
                        p[geneOrder[conn_comp[last_unchecked_index]]], 
                        p[geneOrder[conn_comp[last_unchecked_index]] + 1], 
                        p[geneOrder[*gene]], 
                        p[geneOrder[*gene] + 1])) {
          conn_comp.push_back(*gene);
          gene = child_genes.erase(gene);
        } else {
          gene++;
        }
      }
      last_unchecked_index++;
    }
    outlist.push_back(conn_comp);
    
  }
  
  return outlist;
}

//[[Rcpp::export]]
Rcpp::List process_oriented_graph(IntegerVector i, IntegerVector p, int cells, 
                                  CharacterVector geneNames, IntegerVector geneOrder,
                                  IntegerVector content, IntegerVector starts) {
	IntegerVector child_vertices = find_children(content, starts);
	IntegerVector unclassified = get_cell_set(i, p, cells, geneOrder, starts);

	Rcpp::List components = get_sep_genes(i, p, geneOrder, child_vertices);
	
	Rcpp::List outlist;
	Rcpp::StringVector list_names;

	for (int sep_ind = 0; sep_ind < components.length(); sep_ind++) {
		Rcpp::List child_sublist = process_oriented_graph(i, p, cells, geneNames, geneOrder, content, components[sep_ind]);
	  IntegerVector component = components[sep_ind];
	  IntegerVector real_inds = geneOrder[component];
	  CharacterVector names = geneNames[real_inds];
		Rcpp::String gene_name = concat(names.sort(), " ");
		list_names.push_back(gene_name);
		if (child_sublist.length() > 1) {
			outlist.push_back(child_sublist);  
		} else {
			outlist.push_back(child_sublist[0]);
		}

		IntegerVector start_ind = geneOrder[component];
		IntegerVector end_ind = start_ind + 1;
		unclassified = substract_sorted(unclassified, i, 
                                      		p[start_ind], 
                                      		p[end_ind]);

	}
  
	if (unclassified.length() > 0) {
		outlist.push_back(unclassified);
		list_names.push_back("unclassified");
	}
  
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
        content[smaller_gene] = bigger_gene;
      }
    }
  }
  IntegerVector start = {-1};
  return process_oriented_graph(i, p, cells, geneNames, geneOrder, content, start);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
