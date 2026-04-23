#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]


//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
arma::mat LSSA_arma(const arma::mat & dat, const arma::colvec & freq, int intercept) {

    arma::vec T = arma::conv_to< arma::vec >::from( dat.col(0) );
    arma::vec Y = arma::conv_to< arma::vec >::from( dat.col(1) );
    
    int n = T.size();
    int nf = freq.size();

    arma::vec theta(n);
    arma::vec chisq(nf);
    arma::mat coefs(nf, 4);

    arma::vec ones = Rcpp::rep(1.0, n);
    arma::vec model(n);

    if (intercept == 1) {
        for (int i = 0; i < nf; ++i) {
            theta = 2 * M_PI * freq[i] * T;
            arma::mat X = arma::join_rows(cos(theta), sin(theta), ones);
            arma::mat XtX = X.t() * X;
            arma::mat Xinv = arma::inv(XtX);
            arma::vec P = Xinv * (X.t() * Y);

            arma::vec Y0 = arma::conv_to< arma::vec >::from( P[0] * cos(theta) + P[1] * sin(theta) + P[2] );
            arma::vec sqdiff = Y - Y0;
            arma::vec sq = arma::square(sqdiff);
            double chisq_i = arma::sum(sq);

            chisq[i] = chisq_i;
            coefs(i , 0) = P[0];
            coefs(i , 1) = P[1];
            coefs(i , 2) = P[2];

        }
    } else {
        for (int i = 0; i < nf; ++i) {
            theta = 2 * M_PI * freq[i] * T;
            arma::mat X = arma::join_rows(cos(theta), sin(theta));
            arma::mat XtX = X.t() * X;
            arma::mat Xinv = arma::inv(XtX);
            arma::vec P = Xinv * (X.t() * Y);

            arma::vec Y0 = arma::conv_to< arma::vec >::from( P[0] * cos(theta) + P[1] * sin(theta) );
            arma::vec sqdiff = Y - Y0;
            arma::vec sq = arma::square(sqdiff);
            double chisq_i = arma::sum(sq);

            chisq[i] = chisq_i;
            coefs(i , 0) = P[0];
            coefs(i , 1) = P[1];
            coefs(i , 2) = 0;
        }
    }


    // double Ybar = ;
    arma::vec sqdiff0 = Y - arma::mean(Y);
    arma::vec sq0 = arma::square(sqdiff0);
    double chisq0 = arma::sum(sq0);
    arma::vec pow = (chisq0 - chisq) / chisq0;
    arma::mat pow_rss_coefs = arma::join_rows(pow, chisq, coefs);

    return pow_rss_coefs;
}

//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
arma::vec LSSA_resid_arma(const arma::mat & dat, double freq, int intercept) {

    arma::vec T = arma::conv_to< arma::vec >::from( dat.col(0) );
    arma::vec Y = arma::conv_to< arma::vec >::from( dat.col(1) );
    
    int n = T.size();
    arma::vec theta(n);
    arma::vec resid(n);
    arma::vec Y0 (n);

    if (intercept == 1) {
        arma::vec ones = Rcpp::rep(1.0, n);

        theta = 2 * M_PI * freq * T;
        arma::mat X = arma::join_rows(cos(theta), sin(theta), ones);
        arma::mat XtX = X.t() * X;
        arma::mat Xinv = arma::inv(XtX);
        arma::vec P = Xinv * (X.t() * Y);

        Y0 = arma::conv_to< arma::vec >::from( P[0] * cos(theta) + P[1] * sin(theta) + P[2] );
    } else {
        theta = 2 * M_PI * freq * T;
        arma::mat X = arma::join_rows(cos(theta), sin(theta));
        arma::mat XtX = X.t() * X;
        arma::mat Xinv = arma::inv(XtX);
        arma::vec P = Xinv * (X.t() * Y);

        Y0 = arma::conv_to< arma::vec >::from( P[0] * cos(theta) + P[1] * sin(theta) );

    }
    resid = Y - Y0;


    return resid;
}





// //' @useDynLib arkhaia
// //' @importFrom Rcpp sourceCpp
// // [[Rcpp::export]]
// arma::mat test(const arma::mat & dat, const arma::colvec & freq, int intercept) {

//     arma::vec T = arma::conv_to< arma::vec >::from( dat.col(0) );
//     arma::vec Y = arma::conv_to< arma::vec >::from( dat.col(1) );
    
//     int n = T.size();
//     int nf = freq.size();

//     arma::vec theta(n);
//     arma::vec chisq(nf);
//     arma::mat coefs(nf, 4);

//     arma::vec ones = Rcpp::rep(1.0, n);
//     arma::vec zeros = Rcpp::rep(0.0, n);

//     arma::vec model(n);

//     arma::mat X(n,3);

//     if (intercept == 1) {
//         for (int i = 0; i < nf; ++i) {
//             theta = 2 * M_PI * freq[i] * T;
//             X = arma::conv_to< arma::mat >::from( arma::join_rows(cos(theta), sin(theta), ones) );
//         }
//     } else {
//         for (int i = 0; i < nf; ++i) {
//             theta = 2 * M_PI * freq[i] * T;
//             X = arma::conv_to< arma::mat >::from( arma::join_rows(cos(theta), sin(theta), zeros) );
//         }
//     }


//     return X;
// }
