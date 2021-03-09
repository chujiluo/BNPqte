#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"
#include "bart.h"
#include "rtnorm.h"
#include "heterbart.h"
#include "lambda.h"


#define TRDRAW(a, b) trdraw(a, b)
#define TEDRAW(a, b) tedraw(a, b)


RcppExport SEXP clbart(
    SEXP _in,            //number of observations in training data
    SEXP _ip,            //dimension of x
    SEXP _inp,           //number of observations in test data
    SEXP _ix,            //x, train,  pxn (transposed so rows are contiguous in memory)
    SEXP _iy,            //y, train,  nx1
    SEXP _ixp,           //x, test, pxnp (transposed so rows are contiguous in memory)
    SEXP _im,            //number of trees
    SEXP _inc,           //number of cut points
    SEXP _ind,           //number of kept draws (except for thinnning ..)
    SEXP _iburn,         //number of burn-in draws skipped
    SEXP _ipower,
    SEXP _ibase,
    SEXP _binaryOffset,
    SEXP _itau,
    SEXP _idart,         //dart prior: true(1)=yes, false(0)=no
    SEXP _ia,            //param a for sparsity prior
    SEXP _ib,            //param b for sparsity prior
    SEXP _irho,          //param rho for sparsity prior (default to p)
    SEXP _iaug,          //categorical strategy: true(1)=data augment false(0)=degenerate trees
    SEXP _inkeeptrain,
    SEXP _inkeeptest,
    SEXP _inkeeptreedraws,
    SEXP _inprintevery,
    SEXP _Xinfo
)
{
  //--------------------------------------------------
  //process args
  size_t n = Rcpp::as<int>(_in);
  size_t p = Rcpp::as<int>(_ip);
  size_t np = Rcpp::as<int>(_inp);
  Rcpp::NumericVector  xv(_ix);
  double *ix = &xv[0];
  Rcpp::IntegerVector  yv(_iy); // binary
  int *iy = &yv[0];
  Rcpp::NumericVector  xpv(_ixp);
  double *ixp = &xpv[0];
  size_t m = Rcpp::as<int>(_im);
  Rcpp::IntegerVector _nc(_inc);
  int *numcut = &_nc[0];
  size_t nd = Rcpp::as<int>(_ind);
  size_t burn = Rcpp::as<int>(_iburn);
  double mybeta = Rcpp::as<double>(_ipower);
  double alpha = Rcpp::as<double>(_ibase);
  // lbart does not currently employ the binaryOffset trick
  double binaryOffset = Rcpp::as<double>(_binaryOffset);
  double tau = Rcpp::as<double>(_itau);
  bool dart;
  if(Rcpp::as<int>(_idart)==1) dart=true;
  else dart=false;
  double a = Rcpp::as<double>(_ia);
  double b = Rcpp::as<double>(_ib);
  double rho = Rcpp::as<double>(_irho);
  bool aug;
  if(Rcpp::as<int>(_iaug)==1) aug=true;
  else aug=false;
  size_t nkeeptrain = Rcpp::as<int>(_inkeeptrain);
  size_t nkeeptest = Rcpp::as<int>(_inkeeptest);
  size_t nkeeptreedraws = Rcpp::as<int>(_inkeeptreedraws);
  size_t printevery = Rcpp::as<int>(_inprintevery);
  Rcpp::NumericMatrix Xinfo(_Xinfo);
  
  //return data structures (using Rcpp)
  Rcpp::NumericMatrix trdraw(nkeeptrain,n);
  Rcpp::NumericMatrix tedraw(nkeeptest,np);
  
  Rcpp::NumericMatrix varprb(nkeeptreedraws,p);
  Rcpp::IntegerMatrix varcnt(nkeeptreedraws,p);
  Rcpp::List mr_vecs(nkeeptreedraws);
  
  //random number generation
  arn gen;
  
  heterbart bm(m);
  
  if(Xinfo.size()>0) {
    xinfo _xi;
    _xi.resize(p);
    for(size_t i=0;i<p;i++) {
      _xi[i].resize(numcut[i]);
      //Rcpp::IntegerVector cutpts(Xinfo[i]);
      for(size_t j=0;j<numcut[i];j++) _xi[i][j]=Xinfo(i, j);
    }
    bm.setxinfo(_xi);
  }
  
  std::stringstream treess;  //string stream to write trees to
  treess.precision(10);
  treess << nkeeptreedraws << " " << m << " " << p << std::endl;
    
  printf("*****Into main of lbart\n");
  
  size_t skiptr,skipte,skiptreedraws;
  if(nkeeptrain) {skiptr=nd/nkeeptrain;}
  else skiptr = nd+1;
  if(nkeeptest) {skipte=nd/nkeeptest;}
  else skipte=nd+1;
  if(nkeeptreedraws) {skiptreedraws = nd/nkeeptreedraws;}
  else skiptreedraws=nd+1;
  
  //--------------------------------------------------
  //print args
  Rcpp::Rcout << "*****Data: " << "n, p, np: " << n << ", " << p << ", " << np << std::endl;
  printf("*****BinaryOffset: %lf\n",binaryOffset);
  printf("*****Number of Trees: %zu\n",m);
  printf("*****Prior: mybeta, alpha, tau: %lf,%lf,%lf\n", mybeta, alpha, tau);
  Rcpp::Rcout << "*****Dirichlet: sparse, a, b, rho, augment: " 
              << dart << ", " << a << ", " << b << ", " << rho << ", " << aug << std::endl;
  Rcpp::Rcout << "*****MCMC: (train) nskip, ndpost, keepevery: " << burn << ", " << nkeeptrain << ", " << skiptr << "\n"
              << "           (test) nskip, ndpost, keepevery: " << burn << ", " << nkeeptest << ", " << skipte << std::endl;
  
  //--------------------------------------------------
  //create logit latents
  //z = f(x) + eps, eps ~ N(0,lambda), f ~ BART
  double *z = new double[n]; //latent z's
  double *lambda = new double [n]; //latent lambda's
  //   double *yf = new double[n]; //??
  double *svec = new double[n]; //vector of standard dev for bart = sqrt(lambda)
  for(unsigned int i=0; i<n; i++) {
    if(iy[i]>0) z[i] = 1.0;
    else z[i]=-1.0;
    //iy[i]=z[i]; //iy is already +/- 1
    lambda[i] = 1.0;
    svec[i]=1.0; //square root of 1 is 1.
  }
  
  //--------------------------------------------------
  //set up BART model
  //heterbart bm(m);
  bm.setprior(alpha,mybeta,tau);
  bm.setdata(p,n,ix,z,numcut);
  bm.setdart(a,b,rho,aug,dart);
  
  // dart iterations
  std::vector<double> ivarprb (p,0.);
  std::vector<size_t> ivarcnt (p,0);
  
  //--------------------------------------------------
  //temporary storage
  //out of sample fit
  double* fhattest=0; 
  if(np) { fhattest = new double[np]; }
  
  //--------------------------------------------------
  //mcmc
  //size_t index;
  size_t trcnt=0; //count kept train draws
  size_t tecnt=0; //count kept test draws
  bool keeptest,keeptreedraw;
  
  time_t tp;
  int time1 = time(&tp);
  xinfo& xi = bm.getxinfo();
  
  for(size_t i=0;i<(nd+burn);i++) {
    if(i%printevery==0) 
      Rcpp::Rcout << "-------BART fit " << i << " out of " << (nd+burn) << std::endl;
    if(i==(burn/2)&&dart) bm.startdart();
    //draw bart
    bm.draw(svec,gen);
    
    for(size_t k=0; k<n; k++) {
      z[k]= iy[k]*rtnorm(iy[k]*bm.f(k), -iy[k]*binaryOffset, svec[k], gen);
      lambda[k]=draw_lambda_i(lambda[k], iy[k]*bm.f(k), 1000, 1, gen);
      svec[k] = sqrt(lambda[k]);
    }
    
    if(i>=burn) {
      if(nkeeptrain && (((i-burn+1) % skiptr) ==0)) {
        for(size_t k=0;k<n;k++) TRDRAW(trcnt,k)=bm.f(k);
        trcnt+=1;
      }
      keeptest = nkeeptest && (((i-burn+1) % skipte) ==0) && np;
      if(keeptest) {
        bm.predict(p,np,ixp,fhattest);
        for(size_t k=0;k<np;k++) TEDRAW(tecnt,k)=fhattest[k];
        tecnt+=1;
      }
      keeptreedraw = nkeeptreedraws && (((i-burn+1) % skiptreedraws) ==0);
      if(keeptreedraw) {
        for(size_t j=0;j<m;j++) {
          treess << bm.gettree(j);
        }
        
        ivarcnt=bm.getnv();
        ivarprb=bm.getpv();
        
        size_t k=(i-burn)/skiptreedraws;
        
        mr_vecs[k]=bm.getmrvec();
        
        for(size_t j=0;j<p;j++){
          varcnt(k,j)=ivarcnt[j];
          varprb(k,j)=ivarprb[j];
        }
      }
    }
  }
  int time2 = time(&tp);
  printf("Time elapsed: %ds\n",time2-time1);
  Rcpp::Rcout << "BART Finished!" << std::endl;
  
  //--------------------------------------------------
  if(fhattest) delete[] fhattest;
  delete[] z;
  delete[] lambda;
  delete[] svec;
  
  //--------------------------------------------------
  //return
  Rcpp::List ret;
  ret["yhat.train"]=trdraw;
  ret["yhat.test"]=tedraw;
  ret["varcount"]=varcnt;
  ret["varprob"]=varprb;
  ret["mr_vecs"]=mr_vecs;
  
  Rcpp::List xiret(xi.size());
  for(size_t i=0;i<xi.size();i++) {
    Rcpp::NumericVector vtemp(xi[i].size());
    std::copy(xi[i].begin(),xi[i].end(),vtemp.begin());
    xiret[i] = Rcpp::NumericVector(vtemp);
  }
  
  Rcpp::List treesL;
  treesL["cutpoints"] = xiret;
  treesL["trees"]=Rcpp::CharacterVector(treess.str());
  ret["treedraws"] = treesL;
  
  return ret;

    
}