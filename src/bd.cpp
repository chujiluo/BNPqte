#include "bd.h"

bool bd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, double sigma, 
        std::vector<size_t>& nv, std::vector<double>& pv, bool aug, rn& gen,
        std::vector<std::vector<double>>& mr_vec)
{
  tree::npv goodbots;  //nodes we could birth at (split on)
  double PBx = getpb(x,xi,pi,goodbots); //prob of a birth at x (CL add: just proposal)
  
  if(gen.uniform() < PBx) { //do birth or death
    
    //--------------------------------------------------
    //draw proposal
    tree::tree_p nx; //bottom node to split on
    size_t v,c; //splitting variable and cutpoint
    double pr; //part of metropolis ratio from proposal and prior
    bprop(x,xi,pi,goodbots,PBx,nx,v,c,pr,nv,pv,aug,gen);
    
    //--------------------------------------------------
    //compute sufficient statistics
    size_t nr,nl; //counts in proposed bots
    double syl, syr; //sum of y in proposed bots
    getsuff(x,nx,v,c,xi,di,nl,syl,nr,syr);
    
    //--------------------------------------------------
    //compute alpha
    double alpha=0.0, lalpha=0.0;
    double lhl, lhr, lht;
    if((nl>=5) && (nr>=5)) { //cludge?
      lhl = lh(nl,syl,sigma,pi.tau);
      lhr = lh(nr,syr,sigma,pi.tau);
      lht = lh(nl+nr,syl+syr,sigma,pi.tau);
      
      alpha=1.0;
      lalpha = log(pr) + (lhl+lhr-lht) + log(sigma);  //except log(pr), lalpha is liklihood ratio
      lalpha = std::min(0.0,lalpha);
    }
    
    //--------------------------------------------------
    //try metrop
    double mul,mur; //means for new bottom nodes, left and right
    double temp_mr; //temp Metropolis ratio for the potential newly added spliting variable
    double uu = gen.uniform();
    bool dostep = (alpha > 0) && (log(uu) < lalpha);
    if(dostep) {
      mul = drawnodemu(nl,syl,pi.tau,sigma,gen);
      mur = drawnodemu(nr,syr,pi.tau,sigma,gen);
      
      temp_mr = exp(lalpha);
      
      x.birthp(nx,v,c,mul,mur,temp_mr);  //CL add: x is the current tree, nx is the splitting node
      
      //update importance kernels
      nv[v]++;
      mr_vec[v].push_back(temp_mr);
      
      return true;
    } else {
      return false;
    }
  } else {
    //--------------------------------------------------
    //draw proposal
    double pr;  //part of metropolis ratio from proposal and prior
    tree::tree_p nx; //nog node to death at
    dprop(x,xi,pi,goodbots,PBx,nx,pr,gen);
    
    //--------------------------------------------------
    //compute sufficient statistics
    size_t nr,nl; //counts at bots of nx
    double syl, syr; //sum at bots of nx
    getsuff(x, nx->getl(), nx->getr(), xi, di, nl, syl, nr, syr);
    
    //--------------------------------------------------
    //compute alpha
    double lhl, lhr, lht;
    lhl = lh(nl,syl,sigma,pi.tau);
    lhr = lh(nr,syr,sigma,pi.tau);
    lht = lh(nl+nr,syl+syr,sigma,pi.tau);
    
    double lalpha = log(pr) + (lht - lhl - lhr) - log(sigma);
    lalpha = std::min(0.0,lalpha);
    
    //--------------------------------------------------
    //try metrop
    double mu;
    double temp_mr;
    
    if(log(gen.uniform()) < lalpha) {
      mu = drawnodemu(nl+nr,syl+syr,pi.tau,sigma,gen);
      
      size_t nx_getv = nx->getv();   //split var in proposed node
      
      nv[nx->getv()]--;
      
      x.deathp(nx,mu,temp_mr);
      
      //remove the Metropolis ratio
      auto mr_it = find(mr_vec[nx_getv].begin(), mr_vec[nx_getv].end(), temp_mr);
      int mr_index = distance(mr_vec[nx_getv].begin(), mr_it);
      mr_vec[nx_getv].erase(mr_vec[nx_getv].begin() + mr_index);
      
      return true;
    } else {
      return false;
    }
  }
}