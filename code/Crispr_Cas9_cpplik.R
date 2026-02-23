
cpp_code <- '
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <math.h>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <cstdlib>
using namespace Rcpp;
using namespace std;
//[[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
arma::mat calculatePartials_children(arma::mat Partial1, arma::mat Partial2, arma::mat P1, arma::mat P2){
arma::mat PartialParent;
 PartialParent = (P1 * Partial1) % (P2 * Partial2);
return PartialParent;
}

//[[Rcpp::export]]
arma::mat calculatePartials_child(arma::mat Partial1, arma::mat P){
arma::mat PartialParent;
PartialParent = P * Partial1;
  return PartialParent;
}

//[[Rcpp::export]]
double At(double r, double t) {
  return std::exp(-r * t);
}

//[[Rcpp::export]]
double Bt(double l, double t) {
  return std::exp(-l * t);
}

//[[Rcpp::export]]
arma::vec Ct(arma::vec S, double at, double bt) {
  // Vectorized version (exactly equivalent)
  const double sumS = arma::accu(S);
  return S * (bt - at) / sumS;
}

//[[Rcpp::export]]
arma::mat P_mat(int k, List tree, arma::vec dis_node_root, arma::vec state, double t1, double t2, double r, double l, arma::vec S){
  arma::mat P = arma::mat(state.size() + 2, state.size() + 2, arma::fill::zeros);
  arma::vec edge_len = tree["edge.length"];
  arma::mat edge = tree["edge"];
  arma::uvec nodes_Cas9 = arma::find(dis_node_root >= (t1 + 1e-6) && dis_node_root <= (t2 + 1e-6))+1;
  //Rcout<<nodes_Cas9<<std::endl;

  //compute at and bt
  arma::vec t = edge_len(find(edge.col(1) == k));
  if(t(0) == 0){
    t(0) = 1e-5;
  }
  //Rcout<<"the branch length of the node:"<<std::endl;
  //Rcout<<t(0)<<std::endl;
  double at = At(r, t(0));
  double bt = Bt(l, t(0));
  //Rcout<<"at:"<<std::endl;
  //Rcout<<at<<std::endl;
  //Rcout<<"bt:"<<std::endl;
  //Rcout<<bt<<std::endl;

  // setting the value in specific sites in P
  P(0, 1) = 1 - bt;
  P(1, 1) = 1;
  //Rcout<<"step1:"<<std::endl;
  //Rcout<<P<<std::endl;
  //P.submat(2, 1, P.n_rows - 1, 1) = 1 - bt;
   int numRows = P.n_rows;
   int numCols = P.n_cols;
    // update the values of the second colum in P
    for (int i = 2; i < numRows; i++) {
      P(i, 1) = 1 - bt;
    }
    for (int i = 2; i < numRows; i++) {
      P(i, i) = bt;
    }
  //Rcout<<"step2:"<<std::endl;
  //Rcout<<P<<std::endl;

  //Rcout<<arma::any(nodes_Cas9 == k)<<std::endl;
  if (arma::any(nodes_Cas9 == k)) {
    arma::vec ct = Ct(S, at, bt);
    P(0,0) = at;
    for (int j = 2; j < numCols; j++) {
        P(0,j) = ct(j-2);
    }
  } else {
    P(0, 0) = bt;
  }
  //Rcout<<"step3:"<<std::endl;
  //Rcout<<P<<std::endl;
  return P;
}

//[[Rcpp::export]]
arma::mat calculatePartials(int k, List tree, arma::vec dis_node_root, arma::vec state, List tree_childlist, List Partials, double t1, double t2, double r, double l, arma::vec S)
{
  arma::mat ParentPartials;
  int pp=k;
  std::string pp_str= std::to_string(pp);
  //Rcout<<pp<< std::endl;
  arma::vec children_1 = tree_childlist[pp_str];
  //Rcout<<children_1<< std::endl;
  arma::ivec children(children_1.n_elem);
  for (arma::uword i = 0; i < children_1.n_elem; ++i) {
     children(i) = static_cast<int>(std::round(children_1(i))); 
  }
  //Rcout<<children<< std::endl;
  if (children.size() == 2) { // node with two children
    string c1= std::to_string(children[0]);
    //Rcout<< c1<< std::endl;
    string c2= std::to_string(children[1]);
    arma::mat PartialChild1 = Partials[c1];
    arma::mat PartialChild2 = Partials[c2];
    arma::mat P1 = P_mat(children[0], tree, dis_node_root, state, t1, t2, r, l, S);
    arma::mat P2 = P_mat(children[1], tree, dis_node_root, state, t1, t2, r, l, S);
    //cout<<c1<<c2<<PartialChild1<<PartialChild2<<P1<<P2;
    ParentPartials = calculatePartials_children(PartialChild1, PartialChild2, P1, P2);
  }
  if(children.size()>2){
    string c1 = std::to_string(children[0]);
    //Rcout<< c1<< std::endl;
    arma::mat PartialChild = Partials[c1];
    //Rcout<< PartialChild << std::endl;
    arma::mat P = P_mat(children[0], tree, dis_node_root, state, t1, t2, r, l, S);
    //Rcout<< P << std::endl;
    ParentPartials = calculatePartials_child(PartialChild, P);
    //Rcout<< ParentPartials << std::endl;

    for (int i = 1; i < children.size(); i++){
      //Rcout<< c1 << std::endl;
      string c1 = std::to_string(children[i]);
      arma::mat PartialChild = Partials[c1];
      P = P_mat(children[i], tree, dis_node_root, state, t1, t2, r, l, S);
      ParentPartials = ParentPartials % calculatePartials_child(PartialChild, P);
  }
}
  if(children.size()==1){// node with one children duo to scarring
    arma::mat PartialChild1 = Partials[std::to_string(children[0])];
    //Rcout<< PartialChild1 << std::endl;
    arma::mat P = P_mat(children[0], tree, dis_node_root, state, t1, t2, r, l, S);
    //Rcout<< P << std::endl;
    ParentPartials = calculatePartials_child(PartialChild1, P);
  }
return ParentPartials;

}
//[[Rcpp::export]]
double tree_likehood(List tree_Cas9, arma::vec tree_PostOrder, arma::vec dis_node_root, List tree_childlist, arma::mat state_matrix, arma::vec state, double t1, double t2, double r, double l, arma::vec S, bool loglik = true, bool linearS = true){
  const int N = state_matrix.n_rows;
  const int n_site = state_matrix.n_cols;
  const int Nnode = as<int>(tree_Cas9["Nnode"]);
  const int M = N + Nnode;
  const int nStates = state.size() + 2;

  // Extract edge/edge.length once (huge hot path in the old P_mat)
  arma::vec edge_len = tree_Cas9["edge.length"];
  arma::mat edge = tree_Cas9["edge"]; // columns: parent, child

  // Precompute branch length for each child node id (1..M)
  arma::vec brlen_by_child(M + 1, arma::fill::zeros);
  for (arma::uword e = 0; e < edge.n_rows; ++e) {
    const int child_id = (int)std::lround(edge(e, 1));
    if (child_id >= 1 && child_id <= M) {
      brlen_by_child(child_id) = edge_len(e);
    }
  }

  // Precompute Cas9 node membership once (P_mat checks this each time)
  std::vector<unsigned char> isCas9(M + 1, 0);
  const double lo = t1 + 1e-6;
  const double hi = t2 + 1e-6;
  for (int node = 1; node <= M && (node - 1) < (int)dis_node_root.n_elem; ++node) {
    const double d = dis_node_root(node - 1);
    if (d >= lo && d <= hi) isCas9[node] = 1;
  }

  // Preconvert children list (avoid string-key list access per node)
  std::vector< std::vector<int> > children_by_node(M + 1);
  Rcpp::CharacterVector cnames = tree_childlist.names();
  for (int i = 0; i < tree_childlist.size(); ++i) {
    if (cnames.size() == tree_childlist.size()) {
      const std::string name = Rcpp::as<std::string>(cnames[i]);
      const int node = std::atoi(name.c_str());
      if (node >= 1 && node <= M) {
        Rcpp::NumericVector ch = tree_childlist[i];
        std::vector<int> out;
        out.reserve(ch.size());
        for (int j = 0; j < ch.size(); ++j) {
          out.push_back((int)std::lround(ch[j]));
        }
        children_by_node[node] = std::move(out);
      }
    }
  }

  // Build fast lookup from observed state value -> column index in partial vector
  std::unordered_map<int, int> state_to_idx;
  state_to_idx.reserve(state.size() * 2);
  for (int j = 0; j < (int)state.n_elem; ++j) {
    state_to_idx[(int)std::lround(state(j))] = j + 2;
  }

  // Dense path precomputes full P matrices (O(S^2) memory/time per node).
  std::vector<arma::mat> P_by_child;

  // Linear-in-S path: exploit CRISPR model structure.
  // State indices: 0=unedited, 1=silenced, 2..=scar states.
  struct EdgeParams {
    double bt;
    double a0; // unedited->unedited prob: at if Cas9 active on this edge, else bt
    bool cas9;
    arma::vec ct; // length = (nStates-2) when cas9=true, empty otherwise
  };
  std::vector<EdgeParams> params_by_child;

  if (!linearS) {
    // Precompute transition matrices for each child node once per likelihood call.
    // This eliminates repeated P construction inside the (site x node) hot loop.
    P_by_child.resize(M + 1);
    for (int k = 1; k <= M; ++k) {
      arma::mat P = arma::mat(nStates, nStates, arma::fill::zeros);
      double t = brlen_by_child(k);
      if (t == 0.0) t = 1e-5;

      const double at = At(r, t);
      const double bt = Bt(l, t);

      P(0, 1) = 1 - bt;
      P(1, 1) = 1;
      for (int i = 2; i < nStates; i++) {
        P(i, 1) = 1 - bt;
        P(i, i) = bt;
      }

      if (isCas9[k]) {
        arma::vec ct = Ct(S, at, bt);
        P(0, 0) = at;
        for (int j = 2; j < nStates; j++) {
          P(0, j) = ct(j - 2);
        }
      } else {
        P(0, 0) = bt;
      }
      P_by_child[k] = std::move(P);
    }
  } else {
    // Precompute edge parameters (O(S) memory per node).
    params_by_child.resize(M + 1);
    for (int k = 1; k <= M; ++k) {
      double t = brlen_by_child(k);
      if (t == 0.0) t = 1e-5;
      const double at = At(r, t);
      const double bt = Bt(l, t);
      EdgeParams ep;
      ep.bt = bt;
      ep.cas9 = isCas9[k] != 0;
      ep.a0 = ep.cas9 ? at : bt;
      if (ep.cas9) {
        ep.ct = Ct(S, at, bt);
      }
      params_by_child[k] = std::move(ep);
    }
  }

  auto apply_P_fast = [&](const EdgeParams& ep, const arma::vec& x, arma::vec& y) {
    // y = P * x, using the model sparsity/structure.
    // P rows: parent states; columns: child states.
    y.zeros();
    const double x1 = x(1);

    // parent silenced -> child silenced w.p. 1
    y(1) = x1;

    // parent scar_i -> child silenced with (1-bt), stay scar_i with bt
    if (nStates > 2) {
      y.subvec(2, nStates - 1) = ep.bt * x.subvec(2, nStates - 1);
      y.subvec(2, nStates - 1) += (1.0 - ep.bt) * x1;
    }

    // parent unedited
    double dot_ct = 0.0;
    if (ep.cas9 && nStates > 2) {
      dot_ct = arma::dot(ep.ct, x.subvec(2, nStates - 1));
    }
    y(0) = (1.0 - ep.bt) * x1 + ep.a0 * x(0) + dot_ct;
  };

  double log_lik = 0.0;

  // Partials stored as columns indexed by node id (1..M). Column 0 unused.
  arma::mat partials(nStates, M + 1);
  arma::vec tmp1(nStates);
  arma::vec tmp2(nStates);

  for (int m = 0; m < n_site; m++) {
    partials.zeros();

    for (int i = 0; i < M; i++) {
      const int k = (int)tree_PostOrder(i);
      if (k <= N) {
        const double v = state_matrix(k - 1, m);
        int idx = 0;
        if (v == -1) {
          idx = 1;
        } else if (v > 0) {
          auto it = state_to_idx.find((int)std::lround(v));
          if (it != state_to_idx.end()) idx = it->second;
        }
        partials(idx, k) = 1.0;
      } else {
        const std::vector<int>& kids = children_by_node[k];
        if (kids.size() == 2) {
          const int c1 = kids[0];
          const int c2 = kids[1];
          if (!linearS) {
            const arma::mat& P1 = P_by_child[c1];
            const arma::mat& P2 = P_by_child[c2];
            partials.col(k) = (P1 * partials.col(c1)) % (P2 * partials.col(c2));
          } else {
            apply_P_fast(params_by_child[c1], partials.col(c1), tmp1);
            apply_P_fast(params_by_child[c2], partials.col(c2), tmp2);
            partials.col(k) = tmp1 % tmp2;
          }
        } else if (kids.size() == 1) {
          const int c1 = kids[0];
          if (!linearS) {
            partials.col(k) = P_by_child[c1] * partials.col(c1);
          } else {
            apply_P_fast(params_by_child[c1], partials.col(c1), tmp1);
            partials.col(k) = tmp1;
          }
        } else if (kids.size() > 2) {
          // Preserve behavior for multifurcations
          const int c1 = kids[0];
          if (!linearS) {
            arma::vec acc = P_by_child[c1] * partials.col(c1);
            for (size_t j = 1; j < kids.size(); ++j) {
              const int cj = kids[j];
              acc = acc % (P_by_child[cj] * partials.col(cj));
            }
            partials.col(k) = acc;
          } else {
            apply_P_fast(params_by_child[c1], partials.col(c1), tmp1);
            arma::vec acc = tmp1;
            for (size_t j = 1; j < kids.size(); ++j) {
              const int cj = kids[j];
              apply_P_fast(params_by_child[cj], partials.col(cj), tmp1);
              acc = acc % tmp1;
            }
            partials.col(k) = acc;
          }
        }
      }
    }

    const double root0 = partials(0, N + 1);
    log_lik += std::log(root0);
  }

  if (loglik == true) {
    return log_lik;
  }
  return std::exp(log_lik);
}
'

cpp_file <- file.path(tempdir(), "bilint_cas9_cpplik.cpp")
writeLines(cpp_code, cpp_file)
Rcpp::sourceCpp(cpp_file)
rm(cpp_code, cpp_file)
