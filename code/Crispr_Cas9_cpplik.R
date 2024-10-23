
sourceCpp(code='
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <math.h>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
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
    arma::vec ct(S.size());

    double sumS = 0.0;
    for (const auto& s : S) {
        sumS += s;
    }

    for (size_t i = 0; i < S.size(); ++i) {
        ct[i] = S[i] * (bt - at) / sumS;
    }

    return ct;
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

  arma::vec ct= Ct(S, at, bt);
  //Rcout<<arma::any(nodes_Cas9 == k)<<std::endl;
  if (arma::any(nodes_Cas9 == k)) {
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
  arma::ivec children(children_1.n_elem); // 创建一个与vec相同大小的arma::ivec对象
  for (arma::uword i = 0; i < children_1.n_elem; ++i) {
     children(i) = static_cast<int>(std::round(children_1(i))); // 将double值转换为int，并存储在intVec中
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
double tree_likehood(List tree_Cas9, arma::vec tree_PostOrder, arma::vec dis_node_root, List tree_childlist, arma::mat state_matrix, arma::vec state, double t1, double t2, double r, double l, arma::vec S, bool loglik = true){
  int N = state_matrix.n_rows;
  //Rcout << N << std::endl;
  int n_site = state_matrix.n_cols;
  int Nnode = tree_Cas9["Nnode"];
  int M = N + Nnode;
  double log_lik = 0.0;
  arma::uvec indices(1);
  for(int m=0; m < n_site; m++){
    //Rcout << m << std::endl;
    Rcpp::List Partials = Rcpp::List::create();
    for(int i=0; i < M; i++){
      int k = tree_PostOrder[i];
      if(k <= N){
        arma::vec Partial = arma::zeros(state.size() + 2);
        if(state_matrix(k-1, m)==0){
          indices(0) = 0;
        }
        if(state_matrix(k-1, m)==-1){
          indices(0) = 1;
        }
        if(state_matrix(k-1, m)>0){
          indices = find(state == state_matrix(k-1, m));
          indices += 2;
        }
        Partial.elem(indices).fill(1.0);
        Partials.push_back(Partial, std::to_string(k));
        //Rcout << "k=" << k << std::endl;
        //Rcout << Partial << std::endl;
      }else{
        arma::mat ParentPartials = calculatePartials(k, tree_Cas9, dis_node_root, state, tree_childlist, Partials, t1, t2, r, l, S);
        Partials.push_back(ParentPartials, std::to_string(k));
        //Rcout << "k=" << k << std::endl;
        //Rcout << ParentPartials << std::endl;
      }
    }
    arma::vec root_Partial = Partials[std::to_string(N+1)];
    //Rcout << root_Partial(0) << std::endl;
    log_lik = log_lik + log(root_Partial(0));
    //Rcout << "log_lik=" << log_lik << std::endl;
  }
  if(loglik==true){
    return log_lik;
  }else{
    return std::exp(log_lik);
  }
}
')
