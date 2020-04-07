// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double kge(arma::vec x, arma::vec y) {
  mat sd_x, sd_y, mean_x, mean_y, r, a, b, kge;
  sd_x = arma::stddev(x);
  sd_y = arma::stddev(y);
  mean_x = arma::mean(x);
  mean_y = arma::mean(y);
  r = arma::cor(x, y);
  a = sd_x / sd_y;
  b = mean_x / mean_y;
  kge = 1 - sqrt(pow(r-1, 2)+ pow(a-1, 2) + pow(b-1, 2));
  return(kge(0,0));
}

// [[Rcpp::export]]
NumericVector convolve(arma::vec a, arma::vec b) {
  /* Define the numeric vectors and put in the data from R */
  vec xa(a);
  vec xb(b);

  /* Get the sizes of the vectors */
  int n_xa = xa.size(), n_xb = xb.size();

  /* The length of the overlap */
  int nab = n_xa + n_xb - 1;

  Rcpp::NumericVector xab(nab);
  for (int i = 0; i < n_xa; i++)
    for (int j = 0; j < n_xb; j++)
      xab[i + j] += xa[i] * xb[j];
  return xab;
}

// [[Rcpp::export]]
vec ttd_exponential(double x, arma::vec tw_ttd) {
  return((1 / x) * exp(-tw_ttd / x));
}

// [[Rcpp::export]]
List ttd_convolve(arma::vec tw, //Record of the input and output concentration (Weeks)
                  arma::vec Cin,  //Input (i.e., Precipitation) isotopic signature
                  arma::vec Cobs, //Output observed (i.e, streamflow) isotopic signature
                  const int k = 20,
                  const int u = 10,
                  const int n_runs = 10000,
                  const double TT_low = 0,
                  const double TT_up = 65,
                  const double thresh = 0.95){
  int b = Cin.size();       //number of records of Cin and Cout (e.g., 2 years biweekly, i.e., 52 Weeks)
  double Cobs_mn = mean(Cobs);   //mean of Cout (used for KGE at the end)
  double Cobs_var = var(Cobs);   //variance of Cout

  //extending cycles
  vec Cin_k = repmat(Cin,k,1);// ---Cin repeated k times
  int m = Cin_k.size();       //m -> total number of Weeks (i.e., m must eaual b*k)
  vec Cobs_k = repmat(Cobs,k,1);// ---Cobs repeated k times

  // 3. TTD extended u cycles
  vec tw_TTD = linspace(0, m*u, m*u);  //final length that the TTD will be extended

  // 4. MTT estimation
  // 4.1 Range of TTs to be avaluated
  mat M_z = zeros<mat>(n_runs,2);  // matrix of zeros for output (i.e, TT and KGE->model efficiency)
  // Random generation of n-runs TT between the limits of TTa and TTb to evaluated
  vec TT = (TT_up - TT_low) * randu<vec>(n_runs) + TT_low;
  M_z.col(0) = TT;        // save TT in firt column of predefined matrix of zeros

  // 4.2 Loop for conduct convolution for each TT
  vec KGEs = zeros<vec>(TT.size());
  for(int l = 0; l < n_runs; ++l) {
    vec TTD = (1/TT(l)) * exp(-tw_TTD / TT(l)); //exponential TTD
    //(The integarl of the trasfer function must be 1)
    TTD = TTD / sum(TTD);        //normalization to adjust for resoluton problems.
    vec Cmodfull = conv(Cin_k,TTD);  //returns full convolution (default) of the modeled output isotopic concentration
    int x = Cmodfull.size();
    vec Cmod = Cmodfull(span((k-1)*b, (k-1)*b + b - 1));
    KGEs(l) = kge(Cmod, Cobs);
  }
  M_z.col(1) = KGEs;

  //4.6 Selection of the TT that produces the maximun KGE

  double maxKGE = max(KGEs);                      //max KGE
  uvec maxKGEp = find(KGEs == maxKGE);      //position of the max KGE
  mat TT_maxNS = TT(maxKGEp);                        //TT of the max KGE

  // 5 Select parameters that yield KGE above threshold behavior

  double KGE_th =  thresh * maxKGE;   //define KGE threshols 0.95 * maxKGE
  uvec ind = find(KGEs >= KGE_th); //find indices of elements in M_z matrix whose KGEs are higher than the threshold
  mat Mat_R = zeros<mat>(ind.size(),2);
      Mat_R.col(0) = TT(ind);   //TTs which yield KGEs above threshold
      Mat_R.col(1) = KGEs(ind);   //KGE associated to each TT which yields KGEs above threshold

  // 6. MTT estimation for reanalysis

  // 6.1 Selection of new range of parameters within the selected trheshold for
  // reanalysis

  double TTa_re = min(Mat_R.col(0)); //define TT lower limit for reanalysis i.e., min of TT above threshold
  double TTb_re = max(Mat_R.col(0));  //define TT lower limit for reanalysis i.e., max of TT above threshold
  vec TT_re = (TTb_re - TTa_re) * randu<vec>(n_runs) + TTa_re; // Random generation of n-runs TT between the limits of TTa_re and TTb_re to evaluated

  mat  Mat_R_re = zeros<mat>(b, n_runs); // Mat of zeros for saving convolution n-runs times. b -> # of time steps which will be saved for each run

  vec KGEs_re = zeros<vec>(TT_re.size());
  for(int l = 0; l < n_runs; ++l) {
    vec TTD = (1/TT_re(l))*exp(-tw_TTD/TT_re(l)); //exponential TTD
    TTD = TTD / accu(TTD);        //normalization to adjust for resoluton problems (The integarl of the trasfer function must be 1).
    vec Cmodfull = conv(Cin_k, TTD);           //returns full convolution (default) of the modeled output isotopic concentration
    int x = Cmodfull.size();
    vec Cmod = Cmodfull(span((k-1)*b, (k-1)*b + b - 1));
    Mat_R_re.col(l) = Cmod;
    KGEs_re(l) = kge(Cmod, Cobs);
  }

  // 7. Results of the reanalysis (TT max, percs 5 & 95)

  double maxKGE_re = max(KGEs_re);                        //max KGE
  uvec maxKGEp_re = find(KGEs_re == maxKGE_re);      //position of the max KGE
  double TT_maxNS_re = as_scalar(TT_re(maxKGEp_re));  //TT of the max KGE
  vec P = {0.05, 0.95};
  vec TT_reP =  quantile(TT_re, P);
  mat Mz_re = zeros<mat>(TT_re.size(), 2);
  Mz_re.col(0) = TT_re;
  Mz_re.col(1) = KGEs_re;

  // 8. Estimation of 5 and 95 percentiles for CI band

  mat CI_P = quantile(Mat_R_re, P, 1);

  // 9. For plotting modeled Cout for the TT which brings max KGE

  vec TTD = (1 / TT_maxNS_re) * exp(-tw_TTD / TT_maxNS_re);   //exponential TTD
      TTD = TTD / accu(TTD);        //normalization to adjust for resoluton problems
  // (The integarl of the trasfer function must be 1).
  vec Cmodfull_ = conv(Cin_k, TTD);    //returns full convolution (default)
  vec Cmodfull_an_max = Cmodfull_(span((k-1)*b, (k-1)*b + b - 1));

  // 9.1 Statistics of the Cmodeled with max KGE

  double Cmod_mn = mean(Cmodfull_an_max);   //mean of Cmodeled
  double Cmod_var = var(Cmodfull_an_max);   //variance of Cmodeled
  double RMSE = sqrt(accu(pow(Cobs - Cmodfull_an_max,2)) / b); //Root Mean Square Error
  double MAE = accu(abs(Cobs - Cmodfull_an_max)) / b; // Mean Absolute Error
  double BIAS = accu(Cmodfull_an_max - Cobs) / b; //BIAS

  double Sr = accu(pow(Cmodfull_an_max - Cobs, 2));   //Sri -> (observed-modeled)^2 (vector)
  double St = accu(pow(Cobs - Cobs_mn, 2));   //Sti -> (modeled-mean)^2    (vector)

  double NSE = 1 - Sr / St;  //Nash Sutcliffe Efficiency

  // 9.3 For plotting the Cumulative distribution function

  vec CDF_exp = cumsum(TTD);

  //9.1.1 Save output Statistics (outpu_data)

  DataFrame output_data = DataFrame::create(Rcpp::Named("Obs_mean") = Cobs_mn,
                                  Rcpp::Named("Obs_var") = Cobs_var,
                                  Rcpp::Named("Mod_mean") = Cmod_mn,
                                  Rcpp::Named("Mod_var") = Cmod_var,
                                  Rcpp::Named("KGE_max") = maxKGE_re,
                                  Rcpp::Named("TT_max") = TT_maxNS_re,
                                  Rcpp::Named("TT_q05") = TT_reP(0),
                                  Rcpp::Named("TT_q95") = TT_reP(1),
                                  Rcpp::Named("NSE") = NSE,
                                  Rcpp::Named("RMSE") = RMSE,
                                  Rcpp::Named("MAE") = MAE,
                                  Rcpp::Named("BIAS") = BIAS);

  DataFrame Modeled_Full = DataFrame::create(Rcpp::Named("TTD") = M_z.col(0),
                                             Rcpp::Named("KGE") = M_z.col(1));

  DataFrame Modeled_CI = DataFrame::create(Rcpp::Named("TTD") = Mz_re.col(0),
                                             Rcpp::Named("KGE") = Mz_re.col(1));

  DataFrame CDF_graph = DataFrame::create(Rcpp::Named("Time_step_runs") = tw_TTD,
                                          Rcpp::Named("TTD") = TTD,
                                          Rcpp::Named("CDF") = CDF_exp);

  DataFrame Restuls_graph = DataFrame::create(Rcpp::Named("Time_step") = tw,
                                              Rcpp::Named("Input") = Cin,
                                              Rcpp::Named("Observed") = Cobs,
                                              Rcpp::Named("Modeled") = Cmodfull_an_max,
                                              Rcpp::Named("CI_q05") = CI_P.col(0),
                                              Rcpp::Named("CI_q95") = CI_P.col(1));

  return(
    Rcpp::List::create(
      Rcpp::Named("Model_stats") = output_data,
      Rcpp::Named("Modeled_full") = Modeled_Full,
      Rcpp::Named("TTD_fit") = CDF_graph,
      Rcpp::Named("Modeled") = Restuls_graph
    )
  );
}
