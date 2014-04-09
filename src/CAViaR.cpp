#include "quantileVaR.h"

using namespace Rcpp;
//using namespace arma;

double calcAbsMean(const std::vector<double>& x, unsigned int first, unsigned int last);

class caviar{
  public:
  //---------------- Constructor --------------------------------------------------
	caviar(SEXP vX, SEXP prob){
		setY(vX);
		
		vBeta.resize(3, 0.5);
		RQ = 0.0;
		q = as<double>(prob);
		
	}
	
	//---------------- Setter ------------------------------------------------------
	void setY(SEXP vX){
	
		NumericVector vY(vX);
		y.resize(vY.size());
		iN = y.size();
		
		for(unsigned int i=0; i<iN; i++){
			y[i] = vY[i];
		}
		
		VAR.resize(iN);
		HIT.resize(iN);
	}

	//---------------- Methods for R ----------------------------------------------
	SEXP testInitials(NumericMatrix initMat, NumericVector Model){
		
		int N=initMat.nrow();
		int nBeta = initMat.ncol();
		vBeta.resize(nBeta);
		std::vector<double> RQscore;
		RQscore.resize(N);


		for(int i=0;i<N;i++){
		  for(short int j=0;j<nBeta;j++){
			vBeta[j]=initMat(i,j);
		  }
		  doCalc((int)Model[0]);
		  RQscore[i]=RQ;
		}

		return wrap(RQscore);
	}
	
	SEXP runModel(NumericVector BETA, NumericVector Model){
	//Limit the size of BETA
		vBeta = as< std::vector<double> >(BETA);
		doCalc((int)Model[0]); //Problem? 
		if(!::R_finite(RQ))RQ = 10000000; 
	return wrap(RQ);
	}
	
	//---------------- Run model and calculate HIT rate and RQ(score function) ------
	void doCalc(int Model){
  
	switch (Model){
		case 1: SAV();
			break;
		case 2: AS();
			break;
		case 3: IG();
			break;
		case 4: Adaptive();
			break;
		case 5: HAR();
			break;
		default:  	std::cout << "Unknown model, choose either:" << std::endl;
					std::cout << "1. Symetric Absolute Value" << std::endl;
					std::cout << "2. Asymetric Slope" << std::endl;
					std::cout << "3. Indirect GARCH" << std::endl;
					std::cout << "4. Adaptive" << std::endl;
					std::cout << "5. Quantile HAR" << std::endl;
					exit(0);
	  }
  
	RQ=0;
  
	for(unsigned int i=0;i<y.size();i++){//May start at i=1 since first VaR is constant
		y[i]<VAR[i]? HIT[i]=1 : HIT[i]=0;  //Do I realy need this?
		RQ+=(q-HIT[i])*(y[i]-VAR[i]);
	}
	
	}
  
	//---------------------Loop Symetric Abselute Value----------------------------------
	void SAV(){
  
    for(unsigned int i=1; i<iN; i++){
      VAR[i]=vBeta[0]+vBeta[1]*VAR[i-1]+vBeta[2]*::fabs(y[i-1]);
    }
	}
	//---------------------Asymetric Slope ----------------------------------
	void AS(){
  
    for(unsigned int i=1; i<iN; i++){
      VAR[i]=vBeta[0]+vBeta[1]*VAR[i-1]+vBeta[2]*y[i-1]*(y[i-1]>0)+
              vBeta[3]*y[i-1]*(y[i-1]<0);
    }
	}
	//--------------------- Indirect GARCH(1,1) ----------------------------------
	void IG(){
  
    for(unsigned int i=1; i<iN; i++){
      VAR[i]=sqrt(vBeta[0]+vBeta[1]*pow(VAR[i-1], 2)+vBeta[2]*pow(y[i-1], 2));
    }
	}
	//--------------------- Indirect GARCH(1,1) ----------------------------------
	void Adaptive(){
  
    for(unsigned int i=1; i<iN; i++){
      VAR[i]=VAR[i-1]+vBeta[0]*( 1/(1+exp(G*(y[i-1]-VAR[i-1]))) -q);
    }
	}
	//--------------------- Quantile HAR ----------------------------------
	void HAR(){
  
    for(unsigned int i=1; i<iN; i++){
      if(i>19) VAR[i]=vBeta[0]+vBeta[1]*calcAbsMean(y, i-20, i-1)+vBeta[2]*calcAbsMean(y, i-5, i-1)+vBeta[3]*::fabs(y[i-1]);
	  else if(i>4) VAR[i]=vBeta[0]+vBeta[1]*calcAbsMean(y, 19, 0)+vBeta[2]*calcAbsMean(y, i-5, i-1)+vBeta[3]*::fabs(y[i-1]);
	  else VAR[i]=vBeta[0]+vBeta[1]*calcAbsMean(y, 19, 0)+vBeta[2]*calcAbsMean(y, 4, 0)+vBeta[3]*::fabs(y[i-1]);
    }
	}
	//----------------------------------------------------------------------
  
	std::vector<double> y;
	std::vector<double> VAR;
	std::vector<double> HIT;
	std::vector<double> vBeta;
	double RQ;
	double q;
	static const int G=10;
	
	private:
		size_t iN;

};

//----------------------------- Rcpp Module ----------------------------------
RCPP_MODULE(cppMod){
    using namespace Rcpp ;
	                  
    class_<caviar>( "CAViaR" )
    // expose the default constructor
	.constructor<SEXP, SEXP>()
	.field( "y", &caviar::y )
	.field( "VAR", &caviar::VAR )
	.field( "HIT", &caviar::HIT )
	.field( "RQ", &caviar::RQ )
	.field( "Beta", &caviar::vBeta )
	.field( "q", &caviar::q )
    
    
    .method( "setY", &caviar::setY , "Assign data for analysis" )
	
	.method( "runModel", &caviar::runModel , "Runs a choosen mode, with choosen parameters" )
	.method( "testInitials", &caviar::testInitials , "get the message" )
	.method( "doCalc", &caviar::doCalc , "get the message" )
	
    //.method( "set", &World::set     , "set the message" )
    ;
}                     

SEXP fastDQ(SEXP ys_, SEXP Xs_, SEXP q_) {
	
	NumericVector ys(ys_); NumericMatrix Xs(Xs_); NumericVector q(q_);
    int n = Xs.nrow(), k = Xs.ncol();
    arma::mat X(Xs.begin(), n, k, false);
    arma::colvec y(ys.begin(), ys.size(), false);

    arma::colvec coef = arma::solve( X, y );

    double DQ = arma::as_scalar( coef.t()*X.t()*X*coef );
    DQ = DQ/(q[0]*(1-q[0]));


    return Rcpp::List::create(
        Rcpp::Named("coefficients") = coef,
        Rcpp::Named("DQ")       = DQ
    ) ;

}

double calcAbsMean(const std::vector<double>& x, unsigned int first, unsigned int last){

    double dOut=0;

    for(unsigned int i=first; i <= last; i++){
        dOut += (x[i]*x[i]);

    }

    return sqrt(dOut/(last-first+1));
}



