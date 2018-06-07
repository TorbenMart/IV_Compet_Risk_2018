#include <Rcpp.h>
using namespace Rcpp; 


// [[Rcpp::export]]

 double square_cpp(double x){
	return x*x;
}


// [[Rcpp::export]]

SEXP B_est_compet(NumericVector time, NumericVector status1, NumericVector status2, NumericVector stime,  NumericVector G, NumericVector X){
	int n=time.size();
	int k=stime.size();
	double muhat;
	NumericVector B1(k), dB1(k),status1_t(n), dN1(n), B2(k), dB2(k),status2_t(n), dN2(n), Gc(n);
	NumericMatrix dtmpmat1(k,n), tmpmat(k,n);
 

    muhat=0;
    for (int i=0;i<n;i++){
    	muhat+=G[i]/n;
    	status1_t[i]=(time[i]>=stime[0]) ? 1:0;
    	dN1[i]=(time[i]==stime[0]) ? 1:0;    	
    	status2_t[i]=(time[i]>=stime[0]) ? 1:0;
    	dN2[i]=(time[i]==stime[0]) ? 1:0;
    	}
    	double A1=0, A2=0;
    	for (int i=0;i<n;i++){
    		Gc[i]=G[i]-muhat;
    		dB1[0]+=Gc[i]*dN1[i];
    		A1+=Gc[i]*status1_t[i]*X[i];    		
    		dB2[0]+=Gc[i]*dN2[i];
    		A2+=Gc[i]*status2_t[i]*X[i];
    		dtmpmat1(0,i)=dN2[i];
    		tmpmat(0,i)=dN2[i];
    		}
    	
    		

    	dB1[0]=dB1[0]/A1;
    	B1[0]=dB1[0];
     	dB2[0]=dB2[0]/A2;
    	B2[0]=dB2[0];   
    	
    	for (int j=1;j<k;j++){
		    for (int i=0;i<n;i++){
    	status1_t[i]=(time[i]>=stime[j]) ? 1:0;
    	dN1[i]=(time[i]==stime[j]) ? status1[i]:0;
    	status2_t[i]=(time[i]>=stime[j]) ? 1:0;
    	dN2[i]=(time[i]==stime[j]) ? status2[i]:0;
    	}
    	double A1=0, A2=0;
    	for (int i=0;i<n;i++){
    		dB1[j]+=Gc[i]*exp(B1[j-1]*X[i]+B2[j-1]*X[i])*dN1[i];
    		A1+=Gc[i]*status1_t[i]*exp(B1[j-1]*X[i]+B2[j-1]*X[i])*X[i];
    		dB2[j]+=Gc[i]*exp(B1[j-1]*X[i]+B2[j-1]*X[i])*dN2[i];
    		A2+=Gc[i]*status2_t[i]*exp(B1[j-1]*X[i]+B2[j-1]*X[i])*X[i];
    		dtmpmat1(j,i)=exp(B1[j-1]*X[i])*dN2[i];
    		tmpmat(j,i)=tmpmat(j-1,i)+dtmpmat1(j,i);
    		}

    		
    	dB1[j]=dB1[j]/A1;
    	B1[j]=B1[j-1]+dB1[j];
    	dB2[j]=dB2[j]/A2;
    	B2[j]=B2[j-1]+dB2[j];


    }



  return List::create(
        Named("stime") = stime,
        Named("B1")= B1,
        Named("B2")= B2,
        Named("tmpmat")= tmpmat

    ) ;

}
