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
arma::mat LSSA_LFI_arma(const arma::mat & x, int n_iter, int intercept) {

    int nr = x.n_rows;
    if (n_iter > nr) {
        n_iter = nr;
    }

    arma::mat coefs (n_iter, 3);
    arma::vec freqs (n_iter);
    arma::vec rss (n_iter);
    arma::vec aic (n_iter);

    arma::vec freq_input (999);
    for (int i = 0; i < 999; i++) {
        freq_input[i] = (1/1000) + 0.0005 * (i+1);
    }

    arma::mat first = LSSA_arma(x, freq_input, intercept);
    arma::vec power = first.col(0);
    int n_pow = power.n_elem;
    arma::vec diff (n_pow -1);
    for (int i = 0; i < (n_pow -1); i++) {
        diff[i] = power[i+1] - power[i];
    }

    int idx = 999;
    for (int i = 1; i < n_pow; i++) {
        if (diff[i-1] < 0) {
            if (i-1 < idx) {
                idx = i - 1;
            }
        }
    }

    double freq1 = freq_input[idx];
    double epsilon2 = first(idx,1); 

    freqs[0] = freq1;
    rss[0] = epsilon2;
    coefs(0,0) = first(idx,2);
    coefs(0,1) = first(idx,3);
    coefs(0,2) = first(idx,4);

    arma::vec resid = LSSA_resid_arma(x, freq1, intercept);

    int prev_freq_idx = idx;
    arma::mat x0 (nr , 2);
    x0.col(0) = x.col(0);
    x0.col(1) = resid;

    if (intercept == 1) {
        aic[0] = 2 * (3 + 2) + log(epsilon2 / nr) * nr + nr + nr * log(2 * M_PI);
    } else {
        aic[0] = 2 * (2 + 2) + log(epsilon2 / nr) * nr + nr + nr * log(2 * M_PI);
    }

    int end = 0;
    int limit = 0;

    if (n_iter > 1){
        for (int j = 1; j < n_iter; j++) {

            first = LSSA_arma(x0, freq_input, 0);
            power = first.col(0);
            n_pow = power.n_elem;
            for (int i = 1; i < n_pow; i++) {
                diff[i - 1] = power[i] - power[i - 1];
            }
            if (prev_freq_idx == (n_pow-1)) {
                end = 1;
            } else {
                end = 1;
                for (int i = prev_freq_idx; i < (n_pow-1); i++) {
                    if (diff[i] < 0) {
                        end = 0;
                    }
                }
            }

            if (end == 1) {
                limit = j;
                break;
            }

            idx = 999;
            for (int i = 1; i < n_pow; i++) {
                if (diff[i-1] < 0) {
                    if (i-1 < idx) {
                        if (i-1 > prev_freq_idx) {
                            idx = i - 1;
                        } 
                    }
                }
            }

            freq1 = freq_input[idx];
            epsilon2 = first(idx,1); 

            freqs[j] = freq1;
            rss[j] = epsilon2;
            coefs(j,0) = first(idx,2);
            coefs(j,1) = first(idx,3);
            coefs(j,2) = first(idx,4);

            resid = LSSA_resid_arma(x0, freq1, intercept);

            prev_freq_idx = idx;
            x0.col(1) = resid;

            if (intercept == 1) {
                aic[j] = 2 * (3 * (j+1) + 2) + log(epsilon2 / nr) * nr + nr + nr * log(2 * M_PI);
            } else {
                aic[j] = 2 * (3 * (j+1) + 1) + log(epsilon2 / nr) * nr + nr + nr * log(2 * M_PI);
            }
        }
    }

    // Rcpp::List out = Rcpp::List::create(Rcpp::Named("coefs") = coefs, Rcpp::Named("freqs") = freqs, Rcpp::Named("rss") = rss, Rcpp::Named("aic") = aic);

    int n_out;
    if (end == 1) {
        n_out = limit;
    } else {
        n_out = n_iter;
    }
    
    arma::mat out (n_out, 6);

    if (end == 1) {
        arma::mat coefs0 (limit, 3);
        arma::vec freqs0 (limit);
        arma::vec rss0 (limit);
        arma::vec aic0 (limit);
        for (int i = 0; i < limit; i++) {
            coefs0.row(i) = coefs.row(i);
            freqs0[i] = freqs[i];
            rss0[i] = rss[i];
            aic0[i] = aic[i];
        }
        out.col(0) = coefs0.col(0);
        out.col(1) = coefs0.col(1);
        out.col(2) = coefs0.col(2);
        out.col(3) = freqs0;
        out.col(4) = rss0;
        out.col(5) = aic0;
    } else  {
        out.col(0) = coefs.col(0);
        out.col(1) = coefs.col(1);
        out.col(2) = coefs.col(2);
        out.col(3) = freqs;
        out.col(4) = rss;
        out.col(5) = aic;
    }
        
    // out = Rcpp::List::create(Rcpp::Named("coefs") = coefs0, Rcpp::Named("freqs") = freqs0, Rcpp::Named("rss") = rss0, Rcpp::Named("aic") = aic0);

    return out;
}



//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
arma::mat LSSA_LFI_model_arma(const arma::mat & x, arma::vec t_, int n_iter, int intercept) {
    arma::mat L = LSSA_LFI_arma(x, n_iter, intercept);
    arma::mat coefs = L.cols(0,2);
    arma::vec freqs = L.col(3);
    int n = t_.n_elem;
    arma::mat out (n, 2);
    arma::vec y (n);

    arma::vec theta = t_ * 2 * M_PI * freqs[0];
    if (intercept == 1) {
        y = coefs(0,0) * cos(theta) + coefs(0,1) * sin(theta) + coefs(0,2);
    } else {
        y = coefs(0,0) * cos(theta) + coefs(0,1) * sin(theta);
    }

    if (n_iter > 1) {
        for (int i = 1; i < n_iter; i++) { 
            theta = t_ * 2 * M_PI * freqs[i];
            y = y + coefs(i,0) * cos(theta) + coefs(i,1) * sin(theta);
        }
    }
    
    out.col(0) = t_;
    out.col(1) = y;

    return(out);
}



//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
arma::mat partition_list_to_matrix_arma(Rcpp::List x) {
    int n = x.size();
    arma::vec ns (n);
    for (int i = 0; i < n; i++) {
        Rcpp::List xi = x[i];
        int nxi = xi.size();
        ns[i] = nxi;
    }
    int ncol = arma::max(ns);
    arma::mat out ( n , ncol );
    for (int i = 0; i < n; i++) {
        Rcpp::List r = x[i];
        int nj = ns[i];
        for (int j = 0; j < nj; j++) {
            arma::vec rj = r[j];
            int nk = rj.n_elem;
            for (int k = 0; k < nk; k++) {
                int m = rj[k];
                out(i , m - 1) = j;
            }
        }
    }
    return out;
}



//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
int LSSA_LFI_candidates_arma(arma::mat & x, arma::vec & start_idx, arma::vec & x_rows, arma::mat & sets, int n_iter, int intercept) {

    // // datasets
    // int n = x.size();
    // arma::vec x_nrows(n);
    // int n_obs = 0;
    // for (int i = 0; i < n; i++) {
    //     arma::mat M = x[i];
    //     x_nrows[i] = M.n_rows;
    //     n_obs = n_obs + M.n_rows;
    // }

    int n_obs = x.n_rows;

    // proceed by partitions
    int np = sets.n_rows;
    int nc = sets.n_cols;

    arma::vec aic_min (np);

    for (int P = 0; P < np; P++) {
        arma::rowvec A = sets.row(P);
        int nA = arma::max(A) + 1;

        arma::vec epsilon (n_iter); 
        int n_iter_max = n_iter;

        for (int Q = 0; Q < nA; Q++) {

            // nx nrow of combined data
            int nx = 0;
            for (int Q2 = 0; Q2 < nc; Q2++) {
                if (A[Q2] == Q) {
                    nx = nx + x_rows[Q2];
                }
            }
            arma::mat Y (nx,2);
            int k = 0;
            for (int Q2 = 0; Q2 < nc; Q2++) {
                if (A[Q2] == Q) {
                    int h = start_idx[Q2];
                    int nh = x_rows[Q2];
                    for (int hi = h; hi < (h + nh - 1); hi++) {
                        Y(k, 0) = x(hi, 0);
                        Y(k, 1) = x(hi, 1);
                        k = k + 1;
                    }
                }
            }

            arma::mat W = LSSA_LFI_arma(Y, n_iter, intercept);

            arma::vec rss = W.col(4);
            int n_rss = rss.n_elem;
            for (int i = 0; i < n_rss; i++) {
                epsilon[i] = epsilon[i] + rss[i];
            }

            if (n_rss < n_iter_max) {
                n_iter_max = n_rss;
            }
        }

        arma::vec epsilon2 (n_iter_max);
        for (int i = 0; i < n_iter_max; i++) {
            epsilon2[i] = epsilon[i];
        }

        arma::vec n_iter_out = arma::linspace(1, n_iter_max, n_iter_max);
        arma::vec aic = 2 * nA * (n_iter_out * 3 + 1 + intercept) + n_obs * log(epsilon2 / n_obs) + n_obs + n_obs * log(2 * M_PI);

        aic_min[P] =  arma::min(aic);

    }

    int out = arma::index_min(aic_min);

    return out;
}





//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
arma::mat LSSA_LFI_multi_arma(arma::mat & x, arma::vec & start_idx, arma::vec & x_rows, int n_iter, int intercept) {

    int n_obs = x.n_rows;

    int nA = start_idx.n_elem;

    arma::vec epsilon (n_iter); 
    int n_iter_max = n_iter;

    for (int Q = 0; Q < nA; Q++) {
        int h = start_idx[Q];
        int nh = x_rows[Q];

        arma::mat Y (nh,2);
        
        for (int hi = h; hi < (h + nh - 1); hi++) {
            Y(hi - h, 0) = x(hi, 0);
            Y(hi - h, 1) = x(hi, 1);
        }

        arma::mat W = LSSA_LFI_arma(Y, n_iter, intercept);

        arma::vec rss = W.col(4);
        int n_rss = rss.n_elem;
        for (int i = 0; i < n_rss; i++) {
            epsilon[i] = epsilon[i] + rss[i];
        }

        if (n_rss < n_iter_max) {
            n_iter_max = n_rss;
        }
    }

    arma::vec epsilon2 (n_iter_max);
    for (int i = 0; i < n_iter_max; i++) {
        epsilon2[i] = epsilon[i];
    }

    arma::vec n_iter_out = arma::linspace(1, n_iter_max, n_iter_max);
    arma::vec aic = 2 * nA * (n_iter_out * 3 + 1 + intercept) + n_obs * log(epsilon2 / n_obs) + n_obs + n_obs * log(2 * M_PI);

    arma::mat out (aic.n_elem , 2);
    out.col(0) = n_iter_out;
    out.col(1) = aic;

    return out;
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



