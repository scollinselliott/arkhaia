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



//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
arma::vec lambda_col_arma(const arma::mat & dat, arma::vec & lambda_grid, int omit_zero) {
    int nc = dat.n_cols;

    arma::vec lambdas = arma::vec(nc);

    if (omit_zero == 0) {
        for (int j = 0; j < nc; ++j) {

            arma::vec x = dat.col(j); 

            int n = x.size();
            arma::vec x_ = arma::unique(x) ;
            int nx = x_.size();
            arma::vec rat = arma::vec(nx);

            for (int i = 0; i < nx; ++i) {
                int counter = 0;
                for (int k = 0; k < n; k++) {
                    double xk = x[k];
                    if (xk <= x_[i]) {
                        counter++;
                    }
                }
                double ndouble = n;        
                rat[i] = counter/ndouble;
            }

            int nlam = lambda_grid.size();
            arma::vec y = arma::vec(nlam);

            for (int i = 0; i < nlam; ++i) {
                Rcpp::NumericVector x_numvec(nx);
                x_numvec = x_;
                arma::vec pois_ = Rcpp::dpois(x_numvec, lambda_grid[i]);
                arma::vec pois_ratio = arma::vec(nx);
                for (int k = 0; k < nx; k++) {
                    pois_ratio(k) = pois_[k];
                }
                pois_ratio = pois_ratio / pois_[nx - 1];

                arma::vec sqdiff = rat - pois_ratio;
                arma::vec sq = arma::square(sqdiff);
                double sumsquares = arma::sum(sq);

                y(i) = sumsquares;
            }

            lambdas(j) = lambda_grid[y.index_min() ];
        }
    } else {
        for (int j = 0; j < nc; ++j) {

            arma::vec x0 = dat.col(j); 
            int n0 = x0.size();

            int nz_counter = 0;
            for (int k = 0; k < n0; k++) {
                if (x0[k] > 0) {
                    nz_counter++;
                }
            }

            arma::vec x(nz_counter);

            nz_counter = 0;
            for (int k = 0; k < n0; k++) {
                if (x0[k] > 0) {
                    x[nz_counter] = x0[k];
                    nz_counter++;
                }
            }

            int n = x.size();
            arma::vec x_ = arma::unique(x) ;
            int nx = x_.size();
            arma::vec rat = arma::vec(nx);

            for (int i = 0; i < nx; ++i) {
                int counter = 0;
                for (int k = 0; k < n; k++) {
                    double xk = x[k];
                    if (xk <= x_[i]) {
                        counter++;
                    }
                }
                double ndouble = n;        
                rat[i] = counter/ndouble;
            }

            int nlam = lambda_grid.size();
            arma::vec y = arma::vec(nlam);

            for (int i = 0; i < nlam; ++i) {
                Rcpp::NumericVector x_numvec(nx);
                x_numvec = x_;
                arma::vec pois_ = Rcpp::dpois(x_numvec, lambda_grid[i]);
                arma::vec pois_ratio = arma::vec(nx);
                for (int k = 0; k < nx; k++) {
                    pois_ratio(k) = pois_[k];
                }
                pois_ratio = pois_ratio / pois_[nx - 1];

                arma::vec sqdiff = rat - pois_ratio;
                arma::vec sq = arma::square(sqdiff);
                double sumsquares = arma::sum(sq);

                y(i) = sumsquares;
            }

            lambdas(j) = lambda_grid[y.index_min() ];
        }
    }

    return lambdas;
}



//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
arma::mat trunc_pois_mat_arma(const arma::mat & dat, arma::vec & lambda_grid, int omit_zero) {
    int nc = dat.n_cols;
    int nr = dat.n_rows;

    arma::vec lambdas = lambda_col_arma(dat, lambda_grid, omit_zero);
    Rcpp::IntegerMatrix out(nr,nc);

    int u = 2 * dat.max();
    Rcpp::IntegerVector samp = Rcpp::Range(0,u);

    for (int i = 0; i < nc; ++i) {
    
        arma::vec x = dat.col(i);
        Rcpp::NumericVector probs = Rcpp::dpois(samp, lambdas[i]);

        arma::ivec resamp (nr);

        for (int j = 0; j < nr; ++j) {

            Rcpp::NumericVector probs2(u + 1);
            for (int k = (x[j]); k < u; ++k) {
                probs2(k) = probs(k);
            }
            double s = Rcpp::sum(probs2);
            probs2 = probs2 / s;
            int y = Rcpp::as<int>(Rcpp::sample(samp, 1, FALSE, probs2));
            out( j , i ) = y;
        }
    }
    return Rcpp::as<arma::mat>(out);
}



//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
double CR_arma(const arma::mat & X, double lambda) {

    arma::rowvec Xr = arma::sum(X,0);
    arma::colvec Xc = arma::sum(X,1);

    arma::mat outer = Xc * Xr;
    arma::mat Xe = outer / arma::accu(X);
    arma::mat X0 = arma::pow( (X / Xe), lambda) - 1;

    double out = arma::accu(X % X0) * (2/(lambda * (lambda + 1)));
    return out;
}



//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
double VB_arma(const arma::mat & X, double lambda) {
    int nc = X.n_cols;
    int nr = X.n_rows;

    double chi2 = CR_arma(X, lambda);
    double phi2 = chi2 / arma::accu(X);
    double Ephi2 = ((nr - 1) * (nc - 1) / (arma::accu(X) - 1));
    Rcpp::NumericVector a (1);
    a[0] = phi2 - Ephi2;
    if (a[0] > 0) {
        double denom = 0;
        if (nr < nc) {
            denom = nr - 1;
        } else {
            denom = nc - 1;
        }
        a[0] = a[0] / denom;
        a = Rcpp::sqrt(a);
    } else {
        a[0] = 0;
    }
    return a[0];
}



//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
Rcpp::NumericVector MC_pois_arma(const arma::mat & X, double lambda, arma::vec & lambda_grid, int omit_zero, int M) {
    Rcpp::NumericVector out (M);
    for (int i = 0; i < M; ++i) {
        arma::mat Y = trunc_pois_mat_arma(X, lambda_grid, omit_zero);
        out[i] = VB_arma(Y, lambda);
    }
    return out;
}


