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
arma::mat LSSA_LFI_arma(const arma::mat & x, int n_iter, const arma::colvec & freq, int intercept) {

    int nr = x.n_rows;
    if (n_iter > nr) {
        n_iter = nr;
    }
    int nfreq = freq.n_elem;


    arma::mat coefs (n_iter, 3);
    arma::vec freqs (n_iter);
    arma::vec rss (n_iter);
    arma::vec aic (n_iter);

    // arma::vec freq_input (999);
    // for (int i = 0; i < 999; i++) {
    //     freq_input[i] = (1/1000) + 0.0005 * (i+1);
    // }

    arma::mat first = LSSA_arma(x, freq, intercept);
    arma::vec power = first.col(0);
    int n_pow = power.n_elem;
    arma::vec diff (n_pow -1);
    for (int i = 0; i < (n_pow -1); i++) {
        diff[i] = power[i+1] - power[i];
    }

    int idx = nfreq;
    for (int i = 1; i < n_pow; i++) {
        if (diff[i-1] < 0) {
            if (i-1 < idx) {
                idx = i - 1;
            }
        }
    }

    if (idx == nfreq ) {
        Rcpp::stop("no lowest peak in frequencies detected. Increase upper threshold of freqs");
    }

    // ERROR IS HERE

    double freq1 = freq[idx];
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
        aic[0] = 2 * (3 + 1) + log(epsilon2 / nr) * nr + nr + nr * log(2 * M_PI);
    }

    int end = 0;
    int limit = 0;

    if (n_iter > 1){
        for (int j = 1; j < n_iter; j++) {

            first = LSSA_arma(x0, freq, 0);
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

            freq1 = freq[idx];
            epsilon2 = first(idx,1); 

            freqs[j] = freq1;
            rss[j] = epsilon2;
            coefs(j,0) = first(idx,2);
            coefs(j,1) = first(idx,3);
            coefs(j,2) = first(idx,4);

            resid = LSSA_resid_arma(x0, freq1, 0);

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
arma::mat LSSA_LFI_model_arma(const arma::mat & x, arma::vec t_, int n_iter, const arma::colvec & freq, int intercept) {
    arma::mat L = LSSA_LFI_arma(x, n_iter, freq, intercept);
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
int LSSA_LFI_candidates_arma(arma::mat & x, arma::vec & start_idx, arma::vec & x_rows, arma::mat & sets, int n_iter, const arma::colvec & freq, int intercept) {

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

            arma::mat W = LSSA_LFI_arma(Y, n_iter, freq, intercept);

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
arma::mat LSSA_LFI_multi_arma(arma::mat & x, arma::vec & start_idx, arma::vec & x_rows, int n_iter, const arma::colvec & freq, int intercept) {

    int n_obs = x.n_rows;

    int nA = start_idx.n_elem;

    arma::vec epsilon (n_iter); 
    int n_iter_max = n_iter;

    for (int Q = 0; Q < nA; Q++) {
        int h = start_idx[Q];
        int nh = x_rows[Q];

        arma::mat Y (nh,2);
        
        for (int hi = h; hi < (h + nh); hi++) {
            Y(hi - h, 0) = x(hi, 0);
            Y(hi - h, 1) = x(hi, 1);
        }

        arma::mat W = LSSA_LFI_arma(Y, n_iter, freq, intercept);

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
    arma::vec aic = 2 * nA * (n_iter_out * 3 + 1 + intercept) + (n_obs * log(epsilon2 / n_obs)) + n_obs + (n_obs * log(2 * M_PI));

    arma::mat out (aic.n_elem , 3);
    out.col(0) = n_iter_out;
    out.col(1) = epsilon2;
    out.col(2) = aic;

    return out;
}



//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
arma::mat dtem_arma(const arma::mat & x, arma::vec & events, arma::vec & gammas,  arma::vec & lambdas, arma::vec & weights) {
    int nT = x.n_rows;
    int nE = events.n_elem;

    int ngam = gammas.n_elem;
    int nlam = lambdas.n_elem;

    arma::mat rss (ngam, nlam); 

    for (int g = 0; g < ngam; g++) {
        for (int l = 0; l < nlam; l++) {
            arma::vec xhat (x.n_rows);
            for (int i = 0; i < nT; i++) {
                for (int r = 0; r < nE; r++) {
                    if (events[r] < x(i, 0)) {
                        xhat[i] = xhat[i] + (x(i, 0) - events[r]) * gammas[g] * std::exp(-1 * lambdas[l] * (x(i, 0) - events[r])) ;
                    }
                }
            }
            arma::vec resid = xhat - x.col(1);
            arma::vec resid2 = arma::square(resid);
            resid2 = resid2 % weights;
            rss(g,l) = arma::accu(resid2);
        }
    }
    return rss;
}



//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
arma::cube dtem_intrcpt_arma(const arma::mat & x, arma::vec & events, arma::vec & gammas,  arma::vec & lambdas, arma::vec & intercept, arma::vec & weights) {
    int nT = x.n_rows;
    int nE = events.n_elem;

    int ngam = gammas.n_elem;
    int nlam = lambdas.n_elem;
    int nint = intercept.n_elem;

    arma::cube rss (ngam, nlam, nint); 

    for (int g = 0; g < ngam; g++) {
        for (int l = 0; l < nlam; l++) {
            for (int a = 0; a < nint; a++) {
                arma::vec xhat (x.n_rows);
                for (int i = 0; i < nT; i++) {
                    for (int r = 0; r < nE; r++) {
                        if (events[r] < x(i, 0)) {
                            xhat[i] = xhat[i] + (x(i, 0) - events[r]) * gammas[g] * std::exp(-1 * lambdas[l] * (x(i, 0) - events[r])) ;
                        }
                    }
                }
                xhat = xhat + intercept[a];
                arma::vec resid = xhat - x.col(1);
                arma::vec resid2 = arma::square(resid);
                resid2 = resid2 % weights;
                rss(g,l,a) = arma::accu(resid2);
            }
        }
    }
    return rss;
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



//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
arma::mat wass_R1_arma(arma::vec & X, arma::vec N, int scaled, double scaled_by, double incr) {
    int n_vec = N.n_elem; 

    arma::mat out (n_vec, n_vec);

    arma::vec idx_start (n_vec);
    arma::vec idx_end (n_vec);

    idx_start[0] = 0;
    idx_end[0] = N[0] - 1;
    for (int k = 1; k < n_vec; ++k) {
        idx_start[k] = idx_end[k-1] + 1;
        idx_end[k] = idx_start[k] + N[k] - 1;
    }
    
    for (int ki = 0; ki < (n_vec - 1); ++ki) {
       for (int kj = ki + 1; kj < n_vec; ++kj) {
            int idx11 = idx_start[ki];
            int idx12 = idx_end[ki];
            int idx21 = idx_start[kj];
            int idx22 = idx_end[kj];

            arma::colvec X1 = X.subvec(idx11,idx12);
            arma::colvec X2 = X.subvec(idx21,idx22);

            double minX1 = arma::min(X1);
            double minX2 = arma::min(X2);
            double maxX1 = arma::max(X1);
            double maxX2 = arma::max(X2);

            double minX = minX1;
            if (minX2 < minX1) {
                minX = minX2;
            }
            double maxX = maxX2;
            if (maxX2 > maxX1) {
                maxX = maxX2;
            }
            Rcpp::NumericVector len0v (1);
            len0v[0] = (maxX - minX) / incr;
            arma::vec len1v = Rcpp::round(len0v,0);
            int len1 = len1v[0];

            int n1 = X1.n_elem;
            int n2 = X2.n_elem;

            double nn1 = n1;
            double nn2 = n2;

            double t = minX;
            double summa = 0;
            for (int i = 0; i < len1; ++i) {
                t = minX + i * incr;
                int cnt1 = 0;
                int cnt2 = 0; 
                for (int j = 0; j < n1; ++j) {
                    if (X1[j] <= t) {
                        cnt1 = cnt1 + 1;
                    }
                }
                for (int j = 0; j < n2; ++j) {
                    if (X2[j] <= t) {
                        cnt2 = cnt2 + 1;
                    }
                }
                double f1 = cnt1 / nn1;
                double f2 = cnt2 / nn2;
                double f = f2 - f1;
                if (f < 0) {
                    f = f1 - f2;
                }
                summa = summa + f;
            }
            // Rcpp::NumericVector sqrtnn (1);
            // sqrtnn[0] = (nn1 * nn2) / (nn1 + nn2);
            // Rcpp::NumericVector sqrtn = Rcpp::sqrt(sqrtnn);
            double w = summa * incr;//* sqrtn[0];
            if (scaled == 1) {
                if (scaled_by < 0.000001) {
                    w = w / (maxX - minX);
                } else {
                    w = w / scaled_by;
                }
            } 
            out(kj,ki) = w;
        }
    }
    return out;
}





//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
double wass_R1_armaV2(Rcpp::NumericVector Y1, Rcpp::NumericVector Y2, int scaled, double scaled_by, double incr) {
    arma::vec X1 = Y1;
    arma::vec X2 = Y2;

    double minX1 = arma::min(X1);
    double minX2 = arma::min(X2);
    double maxX1 = arma::max(X1);
    double maxX2 = arma::max(X2);

    double minX = minX1;
    if (minX2 < minX1) {
        minX = minX2;
    }
    double maxX = maxX2;
    if (maxX2 > maxX1) {
        maxX = maxX2;
    }
    Rcpp::NumericVector len0v (1);
    len0v[0] = (maxX - minX) / incr;
    arma::vec len1v = Rcpp::round(len0v,0);
    int len1 = len1v[0];

    int n1 = X1.n_elem;
    int n2 = X2.n_elem;

    double nn1 = n1;
    double nn2 = n2;

    double t = minX;
    double summa = 0;
    for (int i = 0; i < len1; ++i) {
        t = minX + i * incr;
        int cnt1 = 0;
        int cnt2 = 0; 
        for (int j = 0; j < n1; ++j) {
            if (X1[j] <= t) {
                cnt1 = cnt1 + 1;
            }
        }
        for (int j = 0; j < n2; ++j) {
            if (X2[j] <= t) {
                cnt2 = cnt2 + 1;
            }
        }
        double f1 = cnt1 / nn1;
        double f2 = cnt2 / nn2;
        double f = f2 - f1;
        if (f < 0) {
            f = f1 - f2;
        }
        summa = summa + f;
    }
    // Rcpp::NumericVector sqrtnn (1);
    // sqrtnn[0] = (nn1 * nn2) / (nn1 + nn2);
    // Rcpp::NumericVector sqrtn = Rcpp::sqrt(sqrtnn);
    double w = summa * incr; // * sqrtn[0];
    if (scaled == 1) {
        if (scaled_by < 0.000001) {
            w = w / (maxX - minX);
        } else {
            w = w / scaled_by;
        }
    } 
    return w;
}







//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
arma::cube wass_R1_perm_arma(arma::vec & X, arma::vec N, int scaled, double scaled_by, double incr, int n_perm) {
    int n_vec = N.n_elem; 

    arma::cube out (n_vec, n_vec, n_perm);

    arma::vec idx_start (n_vec);
    arma::vec idx_end (n_vec);

    idx_start[0] = 0;
    idx_end[0] = N[0] - 1;
    for (int k = 1; k < n_vec; ++k) {
        idx_start[k] = idx_end[k-1] + 1;
        idx_end[k] = idx_start[k] + N[k] - 1;
    }

    for (int p = 0; p < n_perm ; ++p) {
        arma::vec Y = arma::shuffle(X);
        for (int ki = 0; ki < (n_vec - 1); ++ki) {
            for (int kj = ki + 1; kj < n_vec; ++kj) {
                int idx11 = idx_start[ki];
                int idx12 = idx_end[ki];
                int idx21 = idx_start[kj];
                int idx22 = idx_end[kj];

                arma::colvec X1 = Y.subvec(idx11,idx12);
                arma::colvec X2 = Y.subvec(idx21,idx22);

                double minX1 = arma::min(X1);
                double minX2 = arma::min(X2);
                double maxX1 = arma::max(X1);
                double maxX2 = arma::max(X2);

                double minX = minX1;
                if (minX2 < minX1) {
                    minX = minX2;
                }
                double maxX = maxX2;
                if (maxX2 > maxX1) {
                    maxX = maxX2;
                }
                Rcpp::NumericVector len0v (1);
                len0v[0] = (maxX - minX) / incr;
                arma::vec len1v = Rcpp::round(len0v,0);
                int len1 = len1v[0];

                int n1 = X1.n_elem;
                int n2 = X2.n_elem;

                double nn1 = n1;
                double nn2 = n2;

                double t = minX;
                double summa = 0;
                for (int i = 0; i < len1; ++i) {
                    t = minX + i * incr;
                    int cnt1 = 0;
                    int cnt2 = 0; 
                    for (int j = 0; j < n1; ++j) {
                        if (X1[j] <= t) {
                            cnt1 = cnt1 + 1;
                        }
                    }
                    for (int j = 0; j < n2; ++j) {
                        if (X2[j] <= t) {
                            cnt2 = cnt2 + 1;
                        }
                    }
                    double f1 = cnt1 / nn1;
                    double f2 = cnt2 / nn2;
                    double f = f2 - f1;
                    if (f < 0) {
                        f = f1 - f2;
                    }
                    summa = summa + f;
                }
                // Rcpp::NumericVector sqrtnn (1);
                // sqrtnn[0] = (nn1 * nn2) / (nn1 + nn2);
                // Rcpp::NumericVector sqrtn = Rcpp::sqrt(sqrtnn);
                double w = summa * incr; // * sqrtn[0];
                if (scaled == 1) {
                    if (scaled_by < 0.000001) {
                        w = w / (maxX - minX);
                    } else {
                        w = w / scaled_by;
                    }
                } 
                out(kj,ki,p) = w;
            }
        }
    }

    return out;
}






//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
double ks_Rcpp(Rcpp::NumericVector x1, Rcpp::NumericVector x2) {
    double ks = 0.0;
    int n1 = x1.length();
    int n2 = x2.length();
    double nn1 = n1;
    double nn2 = n2;
    int nz = n1 + n2;
    Rcpp::NumericVector z (nz);
    for (int i = 0; i < n1; ++i) {
        z[i] = x1[i];
    }
    for (int i = n1; i < nz; ++i) {
        z[i] = x2[i-n1];
    }
    Rcpp::NumericVector y = Rcpp::unique(z);
    int ny = y.length();
    Rcpp::NumericVector w (ny);
    for (int i = 0; i < ny; ++i) {
        int k = 0;
        double xk1 = 0.0;
        double xk2 = 0.0;
        for (int j = 0; j < n1; ++j) {
            if (x1[j] <= y[i]) {
              k = k + 1;
            }
        }
        xk1 = k / nn1;
        k = 0;
        for (int j = 0; j < n2; ++j) {
            if (x2[j] <= y[i]) {
              k = k + 1;
            }
        }
        xk2 = k / nn2;
        w[i] = xk1 - xk2;
    }
    w = Rcpp::abs(w);
    ks = Rcpp::max(w);
    // Rcpp::NumericVector sqrtnn (1);
    // sqrtnn[0] = (nn1 * nn2) / (nn1 + nn2);
    // Rcpp::NumericVector sqrtn = Rcpp::sqrt(sqrtnn);
    return ks ;
}




//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
Rcpp::List th_idt_ks_arma(const arma::vec & x1, const arma::vec & x2, double incr) {
    int n1 = x1.n_elem;
    int n2 = x2.n_elem;
    int nz = n1 + n2;
    arma::vec z (nz);
    for (int i = 0; i < n1; ++i) {
        z[i] = x1[i];
    }
    for (int i = n1; i < nz; ++i) {
        z[i] = x2[i-n1];
    }
    Rcpp::NumericVector len0v (1);
    double a0 = arma::min(z);
    double a1 = arma::max(z);
    len0v[0] = (a1 - a0) / incr;
    arma::vec len1v = Rcpp::round(len0v,0);
    int t_ncol = len1v[0];
    int h_nrow = t_ncol / 2;
    arma::vec t_ (t_ncol);
    arma::vec h_ (h_nrow);
    for (int i = 0; i < t_ncol; ++i) {
        t_[i] = a0 + i * incr;
    }
    for (int i = 0; i < h_nrow; ++i) {
        h_[i] = incr + i * incr;
    }
    arma::mat out (h_nrow , t_ncol, arma::fill::value(9999.0));
    for (int ti = 0; ti < t_ncol; ++ti) {
        for (int hi = 0; hi < h_nrow; ++hi) {
            if ((t_[ti] - h_[hi] >= a0) & (t_[ti] + h_[hi] <= a1)) {
                int n_x1a = 0;
                int n_x2a = 0;
                int n_x1p = 0;
                int n_x2p = 0;
                for (int i = 0; i < n1; ++i) {
                    if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti])) {
                        n_x1a = n_x1a + 1;
                    }
                    if ((t_[ti] <= x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                        n_x1p = n_x1p + 1;
                    }
                }
                for (int i = 0; i < n2; ++i) {
                    if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti])) {
                        n_x2a = n_x2a + 1;
                    }
                    if ((t_[ti] <= x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                        n_x2p = n_x2p + 1;
                    }
                }
                if (((n_x1a > 0) & (n_x2a > 0)) & ((n_x1p > 0) & (n_x2p > 0))) {
                    Rcpp::NumericVector x1a (n_x1a);
                    Rcpp::NumericVector x2a (n_x2a);
                    Rcpp::NumericVector x1p (n_x1p);
                    Rcpp::NumericVector x2p (n_x2p);
                    int k = 0;
                    for (int i = 0; i < n1; ++i) {
                        if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti])) {
                            x1a[k] = x1[i];
                            k = k + 1;
                        }
                    }            
                    k = 0;
                    for (int i = 0; i < n1; ++i) {
                        if ((t_[ti] <= x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                            x1p[k] = x1[i];
                            k = k + 1;
                        }
                    }    
                    k = 0;
                    for (int i = 0; i < n2; ++i) {
                        if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti])) {
                            x2a[k] = x2[i];
                            k = k + 1;
                        }
                    }    
                    k = 0;
                    for (int i = 0; i < n2; ++i) {
                        if ((t_[ti] <= x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                            x2p[k] = x2[i];
                            k = k + 1;
                        }
                    }
                    double bef = ks_Rcpp(x1a,x2a);
                    double aft = ks_Rcpp(x1p,x2p);
                    // Rcpp::Rcout << "t : " << t_[ti] << "\n";
                    // Rcpp::Rcout << "h : " << h_[hi] << "\n";
                    // Rcpp::Rcout << "bef : " << bef << "\n";
                    // Rcpp::Rcout << "aft : " << aft << "\n\n";

                    out(hi,ti) = bef - aft;
                }
            }
        }
    }
    return Rcpp::List::create(out, t_, h_);
}





//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
arma::cube th_idt_ks_perm_arma(const arma::vec & x1_input, const arma::vec & x2_input, double incr, int n_perm) {
    int n1 = x1_input.n_elem;
    int n2 = x2_input.n_elem;
    int nz = n1 + n2;
    arma::vec z (nz);
    for (int i = 0; i < n1; ++i) {
        z[i] = x1_input[i];
    }
    for (int i = n1; i < nz; ++i) {
        z[i] = x2_input[i-n1];
    }
    Rcpp::NumericVector len0v (1);
    double a0 = arma::min(z);
    double a1 = arma::max(z);
    len0v[0] = (a1 - a0) / incr;
    arma::vec len1v = Rcpp::round(len0v,0);
    int t_ncol = len1v[0];
    int h_nrow = t_ncol / 2;
    arma::vec t_ (t_ncol);
    arma::vec h_ (h_nrow);
    for (int i = 0; i < t_ncol; ++i) {
        t_[i] = a0 + i * incr;
    }
    for (int i = 0; i < h_nrow; ++i) {
        h_[i] = incr + i * incr;
    }
    arma::cube out (h_nrow , t_ncol, n_perm, arma::fill::value(9999.0));

    for (int p = 0; p < n_perm ; ++p) {
        arma::vec y = arma::shuffle(z);
        arma::vec x1 (n1);
        arma::vec x2 (n2);
        for (int k = 0; k < n1; ++k) {
            x1[k] = y[k];
        }
        for (int k = 0; k < n2; ++k) {
            x2[k] = y[k+n1];
        }
        for (int ti = 0; ti < t_ncol; ++ti) {
            for (int hi = 0; hi < h_nrow; ++hi) {
                if ((t_[ti] - h_[hi] >= a0) & (t_[ti] + h_[hi] <= a1)) {
                    int n_x1a = 0;
                    int n_x2a = 0;
                    int n_x1p = 0;
                    int n_x2p = 0;
                    for (int i = 0; i < n1; ++i) {
                        if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti])) {
                            n_x1a = n_x1a + 1;
                        }
                        if ((t_[ti] <= x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                            n_x1p = n_x1p + 1;
                        }
                    }
                    for (int i = 0; i < n2; ++i) {
                        if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti])) {
                            n_x2a = n_x2a + 1;
                        }
                        if ((t_[ti] <= x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                            n_x2p = n_x2p + 1;
                        }
                    }
                    if (((n_x1a > 0) & (n_x2a > 0)) & ((n_x1p > 0) & (n_x2p > 0))) {
                        Rcpp::NumericVector x1a (n_x1a);
                        Rcpp::NumericVector x2a (n_x2a);
                        Rcpp::NumericVector x1p (n_x1p);
                        Rcpp::NumericVector x2p (n_x2p);
                        int k = 0;
                        for (int i = 0; i < n1; ++i) {
                            if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti])) {
                                x1a[k] = x1[i];
                                k = k + 1;
                            }
                        }            
                        k = 0;
                        for (int i = 0; i < n1; ++i) {
                            if ((t_[ti] <= x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                                x1p[k] = x1[i];
                                k = k + 1;
                            }
                        }    
                        k = 0;
                        for (int i = 0; i < n2; ++i) {
                            if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti])) {
                                x2a[k] = x2[i];
                                k = k + 1;
                            }
                        }    
                        k = 0;
                        for (int i = 0; i < n2; ++i) {
                            if ((t_[ti] <= x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                                x2p[k] = x2[i];
                                k = k + 1;
                            }
                        }
                        double bef = ks_Rcpp(x1a,x2a);
                        double aft = ks_Rcpp(x1p,x2p);
                        out(hi,ti,p) = bef - aft;
                    }
                }
            }
        }
    }
    return out;
}




//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
Rcpp::List th_idt_wass_arma(const arma::vec & x1, const arma::vec & x2, int scaled, double scaled_by, double quad_by, double incr) {
    int n1 = x1.n_elem;
    int n2 = x2.n_elem;
    int nz = n1 + n2;
    arma::vec z (nz);
    for (int i = 0; i < n1; ++i) {
        z[i] = x1[i];
    }
    for (int i = n1; i < nz; ++i) {
        z[i] = x2[i-n1];
    }
    Rcpp::NumericVector len0v (1);
    double a0 = arma::min(z);
    double a1 = arma::max(z);
    len0v[0] = (a1 - a0) / incr;
    arma::vec len1v = Rcpp::round(len0v,0);
    int t_ncol = len1v[0];
    int h_nrow = t_ncol / 2;
    arma::vec t_ (t_ncol);
    arma::vec h_ (h_nrow);
    for (int i = 0; i < t_ncol; ++i) {
        t_[i] = a0 + i * incr;
    }
    for (int i = 0; i < h_nrow; ++i) {
        h_[i] = incr + i * incr;
    }
    arma::mat out (h_nrow , t_ncol, arma::fill::value(9999.0));
    for (int ti = 0; ti < t_ncol; ++ti) {
        for (int hi = 0; hi < h_nrow; ++hi) {
            if ((t_[ti] - h_[hi] >= a0) & (t_[ti] + h_[hi] <= a1)) {
                int n_x1a = 0;
                int n_x2a = 0;
                int n_x1p = 0;
                int n_x2p = 0;
                for (int i = 0; i < n1; ++i) {
                    if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti])) {
                        n_x1a = n_x1a + 1;
                    }
                    if ((t_[ti] <= x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                        n_x1p = n_x1p + 1;
                    }
                }
                for (int i = 0; i < n2; ++i) {
                    if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti])) {
                        n_x2a = n_x2a + 1;
                    }
                    if ((t_[ti] <= x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                        n_x2p = n_x2p + 1;
                    }
                }
                if (((n_x1a > 0) & (n_x2a > 0)) & ((n_x1p > 0) & (n_x2p > 0))) {
                    Rcpp::NumericVector x1a (n_x1a);
                    Rcpp::NumericVector x2a (n_x2a);
                    Rcpp::NumericVector x1p (n_x1p);
                    Rcpp::NumericVector x2p (n_x2p);
                    int k = 0;
                    for (int i = 0; i < n1; ++i) {
                        if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti])) {
                            x1a[k] = x1[i];
                            k = k + 1;
                        }
                    }            
                    k = 0;
                    for (int i = 0; i < n1; ++i) {
                        if ((t_[ti] <= x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                            x1p[k] = x1[i];
                            k = k + 1;
                        }
                    }    
                    k = 0;
                    for (int i = 0; i < n2; ++i) {
                        if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti])) {
                            x2a[k] = x2[i];
                            k = k + 1;
                        }
                    }    
                    k = 0;
                    for (int i = 0; i < n2; ++i) {
                        if ((t_[ti] <= x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                            x2p[k] = x2[i];
                            k = k + 1;
                        }
                    }
                    double sc_by = 0;
                    if (scaled_by < 0) {
                        sc_by = h_[hi];
                    } else {
                        sc_by = scaled_by;
                    }
                    double bef = wass_R1_armaV2(x1a,x2a,scaled,sc_by,quad_by);
                    double aft = wass_R1_armaV2(x1p,x2p,scaled,sc_by,quad_by);
                    // Rcpp::Rcout << "t : " << t_[ti] << "\n";
                    // Rcpp::Rcout << "h : " << h_[hi] << "\n";
                    // Rcpp::Rcout << "bef : " << bef << "\n";
                    // Rcpp::Rcout << "aft : " << aft << "\n\n";

                    out(hi,ti) = bef - aft;
                }
            }
        }
    }
    return Rcpp::List::create(out, t_, h_);
}




//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
arma::cube th_idt_wass_perm_arma(const arma::vec & x1_input, const arma::vec & x2_input, int scaled, double scaled_by, double quad_by, double incr, int n_perm) {
    int n1 = x1_input.n_elem;
    int n2 = x2_input.n_elem;
    int nz = n1 + n2;
    arma::vec z (nz);
    for (int i = 0; i < n1; ++i) {
        z[i] = x1_input[i];
    }
    for (int i = n1; i < nz; ++i) {
        z[i] = x2_input[i-n1];
    }
    Rcpp::NumericVector len0v (1);
    double a0 = arma::min(z);
    double a1 = arma::max(z);
    len0v[0] = (a1 - a0) / incr;
    arma::vec len1v = Rcpp::round(len0v,0);
    int t_ncol = len1v[0];
    int h_nrow = t_ncol / 2;
    arma::vec t_ (t_ncol);
    arma::vec h_ (h_nrow);
    for (int i = 0; i < t_ncol; ++i) {
        t_[i] = a0 + i * incr;
    }
    for (int i = 0; i < h_nrow; ++i) {
        h_[i] = incr + i * incr;
    }
    arma::cube out (h_nrow , t_ncol, n_perm, arma::fill::value(9999.0));

    for (int p = 0; p < n_perm ; ++p) {
        arma::vec y = arma::shuffle(z);
        arma::vec x1 (n1);
        arma::vec x2 (n2);
        for (int k = 0; k < n1; ++k) {
            x1[k] = y[k];
        }
        for (int k = 0; k < n2; ++k) {
            x2[k] = y[k+n1];
        }
        for (int ti = 0; ti < t_ncol; ++ti) {
            for (int hi = 0; hi < h_nrow; ++hi) {
                if ((t_[ti] - h_[hi] >= a0) & (t_[ti] + h_[hi] <= a1)) {
                    int n_x1a = 0;
                    int n_x2a = 0;
                    int n_x1p = 0;
                    int n_x2p = 0;
                    for (int i = 0; i < n1; ++i) {
                        if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti])) {
                            n_x1a = n_x1a + 1;
                        }
                        if ((t_[ti] <= x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                            n_x1p = n_x1p + 1;
                        }
                    }
                    for (int i = 0; i < n2; ++i) {
                        if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti])) {
                            n_x2a = n_x2a + 1;
                        }
                        if ((t_[ti] <= x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                            n_x2p = n_x2p + 1;
                        }
                    }
                    if (((n_x1a > 0) & (n_x2a > 0)) & ((n_x1p > 0) & (n_x2p > 0))) {
                        Rcpp::NumericVector x1a (n_x1a);
                        Rcpp::NumericVector x2a (n_x2a);
                        Rcpp::NumericVector x1p (n_x1p);
                        Rcpp::NumericVector x2p (n_x2p);
                        int k = 0;
                        for (int i = 0; i < n1; ++i) {
                            if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti])) {
                                x1a[k] = x1[i];
                                k = k + 1;
                            }
                        }            
                        k = 0;
                        for (int i = 0; i < n1; ++i) {
                            if ((t_[ti] <= x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                                x1p[k] = x1[i];
                                k = k + 1;
                            }
                        }    
                        k = 0;
                        for (int i = 0; i < n2; ++i) {
                            if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti])) {
                                x2a[k] = x2[i];
                                k = k + 1;
                            }
                        }    
                        k = 0;
                        for (int i = 0; i < n2; ++i) {
                            if ((t_[ti] <= x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                                x2p[k] = x2[i];
                                k = k + 1;
                            }
                        }
                        double sc_by = 0;
                        if (scaled_by < 0) {
                           sc_by = h_[hi];
                        } else {
                            sc_by = scaled_by;
                        }
                        double bef = wass_R1_armaV2(x1a,x2a,scaled,sc_by,quad_by);
                        double aft = wass_R1_armaV2(x1p,x2p,scaled,sc_by,quad_by);
                        out(hi,ti,p) = bef - aft;
                    }
                }
            }
        }
    }
    return out;
}





// With sample size added



//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
Rcpp::List th_idt_ks_stat_arma(const arma::vec & x1, const arma::vec & x2, double incr) {
    int n1 = x1.n_elem;
    int n2 = x2.n_elem;
    int nz = n1 + n2;
    arma::vec z (nz);
    for (int i = 0; i < n1; ++i) {
        z[i] = x1[i];
    }
    for (int i = n1; i < nz; ++i) {
        z[i] = x2[i-n1];
    }
    Rcpp::NumericVector len0v (1);
    double a0 = arma::min(z);
    double a1 = arma::max(z);
    len0v[0] = (a1 - a0) / incr;
    arma::vec len1v = Rcpp::round(len0v,0);
    int t_ncol = len1v[0];
    int h_nrow = t_ncol / 2;
    arma::vec t_ (t_ncol);
    arma::vec h_ (h_nrow);
    for (int i = 0; i < t_ncol; ++i) {
        t_[i] = a0 + i * incr;
    }
    for (int i = 0; i < h_nrow; ++i) {
        h_[i] = incr + i * incr;
    }
    arma::mat out (h_nrow , t_ncol, arma::fill::value(9999.0));
    for (int ti = 0; ti < t_ncol; ++ti) {
        for (int hi = 0; hi < h_nrow; ++hi) {
            if ((t_[ti] - h_[hi] >= a0) & (t_[ti] + h_[hi] <= a1)) {
                int n_x1a = 0;
                int n_x2a = 0;
                int n_x1p = 0;
                int n_x2p = 0;
                for (int i = 0; i < n1; ++i) {
                    if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti])) {
                        n_x1a = n_x1a + 1;
                    }
                    if ((t_[ti] <= x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                        n_x1p = n_x1p + 1;
                    }
                }
                for (int i = 0; i < n2; ++i) {
                    if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti])) {
                        n_x2a = n_x2a + 1;
                    }
                    if ((t_[ti] <= x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                        n_x2p = n_x2p + 1;
                    }
                }
                if (((n_x1a > 0) & (n_x2a > 0)) & ((n_x1p > 0) & (n_x2p > 0))) {
                    Rcpp::NumericVector x1a (n_x1a);
                    Rcpp::NumericVector x2a (n_x2a);
                    Rcpp::NumericVector x1p (n_x1p);
                    Rcpp::NumericVector x2p (n_x2p);
                    int k = 0;
                    for (int i = 0; i < n1; ++i) {
                        if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti])) {
                            x1a[k] = x1[i];
                            k = k + 1;
                        }
                    }            
                    k = 0;
                    for (int i = 0; i < n1; ++i) {
                        if ((t_[ti] <= x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                            x1p[k] = x1[i];
                            k = k + 1;
                        }
                    }    
                    k = 0;
                    for (int i = 0; i < n2; ++i) {
                        if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti])) {
                            x2a[k] = x2[i];
                            k = k + 1;
                        }
                    }    
                    k = 0;
                    for (int i = 0; i < n2; ++i) {
                        if ((t_[ti] <= x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                            x2p[k] = x2[i];
                            k = k + 1;
                        }
                    }
                    double nn1a = x1a.length();
                    double nn2a = x2a.length();
                    double nn1p = x1p.length();
                    double nn2p = x2p.length();
                    Rcpp::NumericVector sqrtnna (1);
                    sqrtnna[0] = (nn1a * nn2a) / (nn1a + nn2a);
                    Rcpp::NumericVector sqrtna = Rcpp::sqrt(sqrtnna);
                    Rcpp::NumericVector sqrtnnp (1);
                    sqrtnnp[0] = (nn1p * nn2p) / (nn1p + nn2p);
                    Rcpp::NumericVector sqrtnp = Rcpp::sqrt(sqrtnnp);
                    
                    double bef = ks_Rcpp(x1a,x2a) * sqrtna[0];
                    double aft = ks_Rcpp(x1p,x2p) * sqrtnp[0];

                    out(hi,ti) = bef - aft;
                }
            }
        }
    }
    return Rcpp::List::create(out, t_, h_);
}





//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
arma::cube th_idt_ks_stat_perm_arma(const arma::vec & x1_input, const arma::vec & x2_input, double incr, int n_perm) {
    int n1 = x1_input.n_elem;
    int n2 = x2_input.n_elem;
    int nz = n1 + n2;
    arma::vec z (nz);
    for (int i = 0; i < n1; ++i) {
        z[i] = x1_input[i];
    }
    for (int i = n1; i < nz; ++i) {
        z[i] = x2_input[i-n1];
    }
    Rcpp::NumericVector len0v (1);
    double a0 = arma::min(z);
    double a1 = arma::max(z);
    len0v[0] = (a1 - a0) / incr;
    arma::vec len1v = Rcpp::round(len0v,0);
    int t_ncol = len1v[0];
    int h_nrow = t_ncol / 2;
    arma::vec t_ (t_ncol);
    arma::vec h_ (h_nrow);
    for (int i = 0; i < t_ncol; ++i) {
        t_[i] = a0 + i * incr;
    }
    for (int i = 0; i < h_nrow; ++i) {
        h_[i] = incr + i * incr;
    }
    arma::cube out (h_nrow , t_ncol, n_perm, arma::fill::value(9999.0));

    for (int p = 0; p < n_perm ; ++p) {
        arma::vec y = arma::shuffle(z);
        arma::vec x1 (n1);
        arma::vec x2 (n2);
        for (int k = 0; k < n1; ++k) {
            x1[k] = y[k];
        }
        for (int k = 0; k < n2; ++k) {
            x2[k] = y[k+n1];
        }
        for (int ti = 0; ti < t_ncol; ++ti) {
            for (int hi = 0; hi < h_nrow; ++hi) {
                if ((t_[ti] - h_[hi] >= a0) & (t_[ti] + h_[hi] <= a1)) {
                    int n_x1a = 0;
                    int n_x2a = 0;
                    int n_x1p = 0;
                    int n_x2p = 0;
                    for (int i = 0; i < n1; ++i) {
                        if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti])) {
                            n_x1a = n_x1a + 1;
                        }
                        if ((t_[ti] <= x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                            n_x1p = n_x1p + 1;
                        }
                    }
                    for (int i = 0; i < n2; ++i) {
                        if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti])) {
                            n_x2a = n_x2a + 1;
                        }
                        if ((t_[ti] <= x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                            n_x2p = n_x2p + 1;
                        }
                    }
                    if (((n_x1a > 0) & (n_x2a > 0)) & ((n_x1p > 0) & (n_x2p > 0))) {
                        Rcpp::NumericVector x1a (n_x1a);
                        Rcpp::NumericVector x2a (n_x2a);
                        Rcpp::NumericVector x1p (n_x1p);
                        Rcpp::NumericVector x2p (n_x2p);
                        int k = 0;
                        for (int i = 0; i < n1; ++i) {
                            if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti])) {
                                x1a[k] = x1[i];
                                k = k + 1;
                            }
                        }            
                        k = 0;
                        for (int i = 0; i < n1; ++i) {
                            if ((t_[ti] <= x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                                x1p[k] = x1[i];
                                k = k + 1;
                            }
                        }    
                        k = 0;
                        for (int i = 0; i < n2; ++i) {
                            if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti])) {
                                x2a[k] = x2[i];
                                k = k + 1;
                            }
                        }    
                        k = 0;
                        for (int i = 0; i < n2; ++i) {
                            if ((t_[ti] <= x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                                x2p[k] = x2[i];
                                k = k + 1;
                            }
                        }
                        double nn1a = x1a.length();
                        double nn2a = x2a.length();
                        double nn1p = x1p.length();
                        double nn2p = x2p.length();
                        Rcpp::NumericVector sqrtnna (1);
                        sqrtnna[0] = (nn1a * nn2a) / (nn1a + nn2a);
                        Rcpp::NumericVector sqrtna = Rcpp::sqrt(sqrtnna);
                        Rcpp::NumericVector sqrtnnp (1);
                        sqrtnnp[0] = (nn1p * nn2p) / (nn1p + nn2p);
                        Rcpp::NumericVector sqrtnp = Rcpp::sqrt(sqrtnnp);
                        double bef = ks_Rcpp(x1a,x2a) * sqrtna[0];
                        double aft = ks_Rcpp(x1p,x2p) * sqrtnp[0];
                        out(hi,ti,p) = bef - aft;
                    }
                }
            }
        }
    }
    return out;
}




//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
Rcpp::List th_idt_wass_stat_arma(const arma::vec & x1, const arma::vec & x2, int scaled, double scaled_by, double quad_by, double incr) {
    int n1 = x1.n_elem;
    int n2 = x2.n_elem;
    int nz = n1 + n2;
    arma::vec z (nz);
    for (int i = 0; i < n1; ++i) {
        z[i] = x1[i];
    }
    for (int i = n1; i < nz; ++i) {
        z[i] = x2[i-n1];
    }
    Rcpp::NumericVector len0v (1);
    double a0 = arma::min(z);
    double a1 = arma::max(z);
    len0v[0] = (a1 - a0) / incr;
    arma::vec len1v = Rcpp::round(len0v,0);
    int t_ncol = len1v[0];
    int h_nrow = t_ncol / 2;
    arma::vec t_ (t_ncol);
    arma::vec h_ (h_nrow);
    for (int i = 0; i < t_ncol; ++i) {
        t_[i] = a0 + i * incr;
    }
    for (int i = 0; i < h_nrow; ++i) {
        h_[i] = incr + i * incr;
    }
    arma::mat out (h_nrow , t_ncol, arma::fill::value(9999.0));
    for (int ti = 0; ti < t_ncol; ++ti) {
        for (int hi = 0; hi < h_nrow; ++hi) {
            if ((t_[ti] - h_[hi] >= a0) & (t_[ti] + h_[hi] <= a1)) {
                int n_x1a = 0;
                int n_x2a = 0;
                int n_x1p = 0;
                int n_x2p = 0;
                for (int i = 0; i < n1; ++i) {
                    if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti])) {
                        n_x1a = n_x1a + 1;
                    }
                    if ((t_[ti] <= x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                        n_x1p = n_x1p + 1;
                    }
                }
                for (int i = 0; i < n2; ++i) {
                    if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti])) {
                        n_x2a = n_x2a + 1;
                    }
                    if ((t_[ti] <= x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                        n_x2p = n_x2p + 1;
                    }
                }
                if (((n_x1a > 0) & (n_x2a > 0)) & ((n_x1p > 0) & (n_x2p > 0))) {
                    Rcpp::NumericVector x1a (n_x1a);
                    Rcpp::NumericVector x2a (n_x2a);
                    Rcpp::NumericVector x1p (n_x1p);
                    Rcpp::NumericVector x2p (n_x2p);
                    int k = 0;
                    for (int i = 0; i < n1; ++i) {
                        if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti])) {
                            x1a[k] = x1[i];
                            k = k + 1;
                        }
                    }            
                    k = 0;
                    for (int i = 0; i < n1; ++i) {
                        if ((t_[ti] <= x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                            x1p[k] = x1[i];
                            k = k + 1;
                        }
                    }    
                    k = 0;
                    for (int i = 0; i < n2; ++i) {
                        if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti])) {
                            x2a[k] = x2[i];
                            k = k + 1;
                        }
                    }    
                    k = 0;
                    for (int i = 0; i < n2; ++i) {
                        if ((t_[ti] <= x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                            x2p[k] = x2[i];
                            k = k + 1;
                        }
                    }
                    double sc_by = 0;
                    if (scaled_by < 0) {
                        sc_by = h_[hi];
                    } else {
                        sc_by = scaled_by;
                    }
                    double nn1a = x1a.length();
                    double nn2a = x2a.length();
                    double nn1p = x1p.length();
                    double nn2p = x2p.length();
                    Rcpp::NumericVector sqrtnna (1);
                    sqrtnna[0] = (nn1a * nn2a) / (nn1a + nn2a);
                    Rcpp::NumericVector sqrtna = Rcpp::sqrt(sqrtnna);
                    Rcpp::NumericVector sqrtnnp (1);
                    sqrtnnp[0] = (nn1p * nn2p) / (nn1p + nn2p);
                    Rcpp::NumericVector sqrtnp = Rcpp::sqrt(sqrtnnp);
                    double bef = wass_R1_armaV2(x1a,x2a,scaled,sc_by,quad_by) * sqrtna[0];
                    double aft = wass_R1_armaV2(x1p,x2p,scaled,sc_by,quad_by) * sqrtna[0];

                    out(hi,ti) = bef - aft;
                }
            }
        }
    }
    return Rcpp::List::create(out, t_, h_);
}




//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
arma::cube th_idt_wass_stat_perm_arma(const arma::vec & x1_input, const arma::vec & x2_input, int scaled, double scaled_by, double quad_by, double incr, int n_perm) {
    int n1 = x1_input.n_elem;
    int n2 = x2_input.n_elem;
    int nz = n1 + n2;
    arma::vec z (nz);
    for (int i = 0; i < n1; ++i) {
        z[i] = x1_input[i];
    }
    for (int i = n1; i < nz; ++i) {
        z[i] = x2_input[i-n1];
    }
    Rcpp::NumericVector len0v (1);
    double a0 = arma::min(z);
    double a1 = arma::max(z);
    len0v[0] = (a1 - a0) / incr;
    arma::vec len1v = Rcpp::round(len0v,0);
    int t_ncol = len1v[0];
    int h_nrow = t_ncol / 2;
    arma::vec t_ (t_ncol);
    arma::vec h_ (h_nrow);
    for (int i = 0; i < t_ncol; ++i) {
        t_[i] = a0 + i * incr;
    }
    for (int i = 0; i < h_nrow; ++i) {
        h_[i] = incr + i * incr;
    }
    arma::cube out (h_nrow , t_ncol, n_perm, arma::fill::value(9999.0));

    for (int p = 0; p < n_perm ; ++p) {
        arma::vec y = arma::shuffle(z);
        arma::vec x1 (n1);
        arma::vec x2 (n2);
        for (int k = 0; k < n1; ++k) {
            x1[k] = y[k];
        }
        for (int k = 0; k < n2; ++k) {
            x2[k] = y[k+n1];
        }
        for (int ti = 0; ti < t_ncol; ++ti) {
            for (int hi = 0; hi < h_nrow; ++hi) {
                if ((t_[ti] - h_[hi] >= a0) & (t_[ti] + h_[hi] <= a1)) {
                    int n_x1a = 0;
                    int n_x2a = 0;
                    int n_x1p = 0;
                    int n_x2p = 0;
                    for (int i = 0; i < n1; ++i) {
                        if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti])) {
                            n_x1a = n_x1a + 1;
                        }
                        if ((t_[ti] <= x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                            n_x1p = n_x1p + 1;
                        }
                    }
                    for (int i = 0; i < n2; ++i) {
                        if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti])) {
                            n_x2a = n_x2a + 1;
                        }
                        if ((t_[ti] <= x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                            n_x2p = n_x2p + 1;
                        }
                    }
                    if (((n_x1a > 0) & (n_x2a > 0)) & ((n_x1p > 0) & (n_x2p > 0))) {
                        Rcpp::NumericVector x1a (n_x1a);
                        Rcpp::NumericVector x2a (n_x2a);
                        Rcpp::NumericVector x1p (n_x1p);
                        Rcpp::NumericVector x2p (n_x2p);
                        int k = 0;
                        for (int i = 0; i < n1; ++i) {
                            if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti])) {
                                x1a[k] = x1[i];
                                k = k + 1;
                            }
                        }            
                        k = 0;
                        for (int i = 0; i < n1; ++i) {
                            if ((t_[ti] <= x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                                x1p[k] = x1[i];
                                k = k + 1;
                            }
                        }    
                        k = 0;
                        for (int i = 0; i < n2; ++i) {
                            if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti])) {
                                x2a[k] = x2[i];
                                k = k + 1;
                            }
                        }    
                        k = 0;
                        for (int i = 0; i < n2; ++i) {
                            if ((t_[ti] <= x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                                x2p[k] = x2[i];
                                k = k + 1;
                            }
                        }
                        double sc_by = 0;
                        if (scaled_by < 0) {
                           sc_by = h_[hi];
                        } else {
                            sc_by = scaled_by;
                        }
                        double nn1a = x1a.length();
                        double nn2a = x2a.length();
                        double nn1p = x1p.length();
                        double nn2p = x2p.length();
                        Rcpp::NumericVector sqrtnna (1);
                        sqrtnna[0] = (nn1a * nn2a) / (nn1a + nn2a);
                        Rcpp::NumericVector sqrtna = Rcpp::sqrt(sqrtnna);
                        Rcpp::NumericVector sqrtnnp (1);
                        sqrtnnp[0] = (nn1p * nn2p) / (nn1p + nn2p);
                        Rcpp::NumericVector sqrtnp = Rcpp::sqrt(sqrtnnp);
                        double bef = wass_R1_armaV2(x1a,x2a,scaled,sc_by,quad_by) * sqrtna[0];
                        double aft = wass_R1_armaV2(x1p,x2p,scaled,sc_by,quad_by) * sqrtnp[0];
                        out(hi,ti,p) = bef - aft;
                    }
                }
            }
        }
    }
    return out;
}






















// idf functions









//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
Rcpp::List idt_ks_arma(const arma::vec & x1, const arma::vec & x2, double incr) {
    int n1 = x1.n_elem;
    int n2 = x2.n_elem;
    int nz = n1 + n2;
    arma::vec z (nz);
    for (int i = 0; i < n1; ++i) {
        z[i] = x1[i];
    }
    for (int i = n1; i < nz; ++i) {
        z[i] = x2[i-n1];
    }
    Rcpp::NumericVector len0v (1);
    double a0 = arma::min(z);
    double a1 = arma::max(z);
    len0v[0] = (a1 - a0) / incr;
    arma::vec len1v = Rcpp::round(len0v,0);
    int t_ncol = len1v[0];
    int h_nrow = t_ncol / 2;
    arma::vec t_ (t_ncol);
    arma::vec h_ (h_nrow);
    for (int i = 0; i < t_ncol; ++i) {
        t_[i] = a0 + i * incr;
    }
    for (int i = 0; i < h_nrow; ++i) {
        h_[i] = incr + i * incr;
    }
    arma::mat out (h_nrow , t_ncol, arma::fill::value(9999.0));
    for (int ti = 0; ti < t_ncol; ++ti) {
        for (int hi = 0; hi < h_nrow; ++hi) {
            if ((t_[ti] - h_[hi] >= a0) & (t_[ti] + h_[hi] <= a1)) {
                int n_x1a = 0;
                int n_x2a = 0;
                for (int i = 0; i < n1; ++i) {
                    if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                        n_x1a = n_x1a + 1;
                    }
                }
                for (int i = 0; i < n2; ++i) {
                    if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                        n_x2a = n_x2a + 1;
                    }
                }
                if ((n_x1a > 0) & (n_x2a > 0)) {
                    Rcpp::NumericVector x1a (n_x1a);
                    Rcpp::NumericVector x2a (n_x2a);
                    int k = 0;
                    for (int i = 0; i < n1; ++i) {
                        if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                            x1a[k] = x1[i];
                            k = k + 1;
                        }
                    }            
                    k = 0;
                    for (int i = 0; i < n2; ++i) {
                        if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                            x2a[k] = x2[i];
                            k = k + 1;
                        }
                    }    
                    // Rcpp::Rcout << "t : " << t_[ti] << "\n";
                    // Rcpp::Rcout << "h : " << h_[hi] << "\n";
                    // Rcpp::Rcout << "bef : " << bef << "\n";
                    // Rcpp::Rcout << "aft : " << aft << "\n\n";

                    out(hi,ti) = ks_Rcpp(x1a,x2a);
                }
            }
        }
    }
    return Rcpp::List::create(out, t_, h_);
}





//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
arma::cube idt_ks_perm_arma(const arma::vec & x1_input, const arma::vec & x2_input, double incr, int n_perm) {
    int n1 = x1_input.n_elem;
    int n2 = x2_input.n_elem;
    int nz = n1 + n2;
    arma::vec z (nz);
    for (int i = 0; i < n1; ++i) {
        z[i] = x1_input[i];
    }
    for (int i = n1; i < nz; ++i) {
        z[i] = x2_input[i-n1];
    }
    Rcpp::NumericVector len0v (1);
    double a0 = arma::min(z);
    double a1 = arma::max(z);
    len0v[0] = (a1 - a0) / incr;
    arma::vec len1v = Rcpp::round(len0v,0);
    int t_ncol = len1v[0];
    int h_nrow = t_ncol / 2;
    arma::vec t_ (t_ncol);
    arma::vec h_ (h_nrow);
    for (int i = 0; i < t_ncol; ++i) {
        t_[i] = a0 + i * incr;
    }
    for (int i = 0; i < h_nrow; ++i) {
        h_[i] = incr + i * incr;
    }
    arma::cube out (h_nrow , t_ncol, n_perm, arma::fill::value(9999.0));

    for (int p = 0; p < n_perm ; ++p) {
        arma::vec y = arma::shuffle(z);
        arma::vec x1 (n1);
        arma::vec x2 (n2);
        for (int k = 0; k < n1; ++k) {
            x1[k] = y[k];
        }
        for (int k = 0; k < n2; ++k) {
            x2[k] = y[k+n1];
        }
        for (int ti = 0; ti < t_ncol; ++ti) {
            for (int hi = 0; hi < h_nrow; ++hi) {
                if ((t_[ti] - h_[hi] >= a0) & (t_[ti] + h_[hi] <= a1)) {
                    int n_x1a = 0;
                    int n_x2a = 0;
                    for (int i = 0; i < n1; ++i) {
                        if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                            n_x1a = n_x1a + 1;
                        }
                    }
                    for (int i = 0; i < n2; ++i) {
                        if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                            n_x2a = n_x2a + 1;
                        }
                    }
                    if ((n_x1a > 0) & (n_x2a > 0)) {
                        Rcpp::NumericVector x1a (n_x1a);
                        Rcpp::NumericVector x2a (n_x2a);
                        int k = 0;
                        for (int i = 0; i < n1; ++i) {
                            if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                                x1a[k] = x1[i];
                                k = k + 1;
                            }
                        }
                        k = 0;
                        for (int i = 0; i < n2; ++i) {
                            if ((t_[ti] - h_[hi] < x2[i] ) & (x2[i] < t_[ti]+ h_[hi])) {
                                x2a[k] = x2[i];
                                k = k + 1;
                            }
                        }   
                        // Rcpp::Rcout << "t : " << t_[ti] << "\n";
                        // Rcpp::Rcout << "h : " << h_[hi] << "\n";     
                        out(hi,ti,p) = ks_Rcpp(x1a,x2a);
                    }
                }
            }
        }
    }
    // Rcpp::Rcout << "h : " << out << "\n";   
    return out;
}




//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
Rcpp::List idt_wass_arma(const arma::vec & x1, const arma::vec & x2, int scaled, double scaled_by, double quad_by, double incr) {
    int n1 = x1.n_elem;
    int n2 = x2.n_elem;
    int nz = n1 + n2;
    arma::vec z (nz);
    for (int i = 0; i < n1; ++i) {
        z[i] = x1[i];
    }
    for (int i = n1; i < nz; ++i) {
        z[i] = x2[i-n1];
    }
    Rcpp::NumericVector len0v (1);
    double a0 = arma::min(z);
    double a1 = arma::max(z);
    len0v[0] = (a1 - a0) / incr;
    arma::vec len1v = Rcpp::round(len0v,0);
    int t_ncol = len1v[0];
    int h_nrow = t_ncol / 2;
    arma::vec t_ (t_ncol);
    arma::vec h_ (h_nrow);
    for (int i = 0; i < t_ncol; ++i) {
        t_[i] = a0 + i * incr;
    }
    for (int i = 0; i < h_nrow; ++i) {
        h_[i] = incr + i * incr;
    }
    arma::mat out (h_nrow , t_ncol, arma::fill::value(9999.0));
    for (int ti = 0; ti < t_ncol; ++ti) {
        for (int hi = 0; hi < h_nrow; ++hi) {
            if ((t_[ti] - h_[hi] >= a0) & (t_[ti] + h_[hi] <= a1)) {
                int n_x1a = 0;
                int n_x2a = 0;
                for (int i = 0; i < n1; ++i) {
                    if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                        n_x1a = n_x1a + 1;
                    }
                }
                for (int i = 0; i < n2; ++i) {
                    if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                        n_x2a = n_x2a + 1;
                    }
                }
                if (((n_x1a > 0) & (n_x2a > 0))) {
                    Rcpp::NumericVector x1a (n_x1a);
                    Rcpp::NumericVector x2a (n_x2a);
                    int k = 0;
                    for (int i = 0; i < n1; ++i) {
                        if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                            x1a[k] = x1[i];
                            k = k + 1;
                        }
                    }            
                    k = 0;
                    for (int i = 0; i < n2; ++i) {
                        if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                            x2a[k] = x2[i];
                            k = k + 1;
                        }
                    }    
                    k = 0;
                    double sc_by = 0;
                    if (scaled_by < 0) {
                        sc_by = h_[hi];
                    } else {
                        sc_by = scaled_by;
                    }
                    // Rcpp::Rcout << "t : " << t_[ti] << "\n";
                    // Rcpp::Rcout << "h : " << h_[hi] << "\n";
                    // Rcpp::Rcout << "bef : " << bef << "\n";
                    // Rcpp::Rcout << "aft : " << aft << "\n\n";

                    out(hi,ti) = wass_R1_armaV2(x1a,x2a,scaled,sc_by,quad_by);
                }
            }
        }
    }
    return Rcpp::List::create(out, t_, h_);
}




//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
arma::cube idt_wass_perm_arma(const arma::vec & x1_input, const arma::vec & x2_input, int scaled, double scaled_by, double quad_by, double incr, int n_perm) {
    int n1 = x1_input.n_elem;
    int n2 = x2_input.n_elem;
    int nz = n1 + n2;
    arma::vec z (nz);
    for (int i = 0; i < n1; ++i) {
        z[i] = x1_input[i];
    }
    for (int i = n1; i < nz; ++i) {
        z[i] = x2_input[i-n1];
    }
    Rcpp::NumericVector len0v (1);
    double a0 = arma::min(z);
    double a1 = arma::max(z);
    len0v[0] = (a1 - a0) / incr;
    arma::vec len1v = Rcpp::round(len0v,0);
    int t_ncol = len1v[0];
    int h_nrow = t_ncol / 2;
    arma::vec t_ (t_ncol);
    arma::vec h_ (h_nrow);
    for (int i = 0; i < t_ncol; ++i) {
        t_[i] = a0 + i * incr;
    }
    for (int i = 0; i < h_nrow; ++i) {
        h_[i] = incr + i * incr;
    }
    arma::cube out (h_nrow , t_ncol, n_perm, arma::fill::value(9999.0));

    for (int p = 0; p < n_perm ; ++p) {
        arma::vec y = arma::shuffle(z);
        arma::vec x1 (n1);
        arma::vec x2 (n2);
        for (int k = 0; k < n1; ++k) {
            x1[k] = y[k];
        }
        for (int k = 0; k < n2; ++k) {
            x2[k] = y[k+n1];
        }
        for (int ti = 0; ti < t_ncol; ++ti) {
            for (int hi = 0; hi < h_nrow; ++hi) {
                if ((t_[ti] - h_[hi] >= a0) & (t_[ti] + h_[hi] <= a1)) {
                    int n_x1a = 0;
                    int n_x2a = 0;
                    for (int i = 0; i < n1; ++i) {
                        if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                            n_x1a = n_x1a + 1;
                        }
                    }
                    for (int i = 0; i < n2; ++i) {
                        if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                            n_x2a = n_x2a + 1;
                        }
                    }
                    if (((n_x1a > 0) & (n_x2a > 0))) {
                        Rcpp::NumericVector x1a (n_x1a);
                        Rcpp::NumericVector x2a (n_x2a);
                        int k = 0;
                        for (int i = 0; i < n1; ++i) {
                            if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                                x1a[k] = x1[i];
                                k = k + 1;
                            }
                        }             
                        k = 0;
                        for (int i = 0; i < n2; ++i) {
                            if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                                x2a[k] = x2[i];
                                k = k + 1;
                            }
                        }    
                        double sc_by = 0;
                        if (scaled_by < 0) {
                           sc_by = h_[hi];
                        } else {
                            sc_by = scaled_by;
                        }
                        double bef = wass_R1_armaV2(x1a,x2a,scaled,sc_by,quad_by);
                        out(hi,ti,p) = bef;
                    }
                }
            }
        }
    }
    return out;
}





// With sample size added



//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
Rcpp::List idt_ks_stat_arma(const arma::vec & x1, const arma::vec & x2, double incr) {
    int n1 = x1.n_elem;
    int n2 = x2.n_elem;
    int nz = n1 + n2;
    arma::vec z (nz);
    for (int i = 0; i < n1; ++i) {
        z[i] = x1[i];
    }
    for (int i = n1; i < nz; ++i) {
        z[i] = x2[i-n1];
    }
    Rcpp::NumericVector len0v (1);
    double a0 = arma::min(z);
    double a1 = arma::max(z);
    len0v[0] = (a1 - a0) / incr;
    arma::vec len1v = Rcpp::round(len0v,0);
    int t_ncol = len1v[0];
    int h_nrow = t_ncol / 2;
    arma::vec t_ (t_ncol);
    arma::vec h_ (h_nrow);
    for (int i = 0; i < t_ncol; ++i) {
        t_[i] = a0 + i * incr;
    }
    for (int i = 0; i < h_nrow; ++i) {
        h_[i] = incr + i * incr;
    }
    arma::mat out (h_nrow , t_ncol, arma::fill::value(9999.0));
    for (int ti = 0; ti < t_ncol; ++ti) {
        for (int hi = 0; hi < h_nrow; ++hi) {
            if ((t_[ti] - h_[hi] >= a0) & (t_[ti] + h_[hi] <= a1)) {
                int n_x1a = 0;
                int n_x2a = 0;
                for (int i = 0; i < n1; ++i) {
                    if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                        n_x1a = n_x1a + 1;
                    }
                }
                for (int i = 0; i < n2; ++i) {
                    if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                        n_x2a = n_x2a + 1;
                    }
                }
                if (((n_x1a > 0) & (n_x2a > 0))) {
                    Rcpp::NumericVector x1a (n_x1a);
                    Rcpp::NumericVector x2a (n_x2a);
                    int k = 0;
                    for (int i = 0; i < n1; ++i) {
                        if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                            x1a[k] = x1[i];
                            k = k + 1;
                        }
                    }             
                    k = 0;
                    for (int i = 0; i < n2; ++i) {
                        if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                            x2a[k] = x2[i];
                            k = k + 1;
                        }
                    }    
                    double nn1a = x1a.length();
                    double nn2a = x2a.length();
                    Rcpp::NumericVector sqrtnna (1);
                    sqrtnna[0] = (nn1a * nn2a) / (nn1a + nn2a);
                    Rcpp::NumericVector sqrtna = Rcpp::sqrt(sqrtnna);

                    double bef = ks_Rcpp(x1a,x2a) * sqrtna[0];

                    out(hi,ti) = bef;
                }
            }
        }
    }
    return Rcpp::List::create(out, t_, h_);
}





//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
arma::cube idt_ks_stat_perm_arma(const arma::vec & x1_input, const arma::vec & x2_input, double incr, int n_perm) {
    int n1 = x1_input.n_elem;
    int n2 = x2_input.n_elem;
    int nz = n1 + n2;
    arma::vec z (nz);
    for (int i = 0; i < n1; ++i) {
        z[i] = x1_input[i];
    }
    for (int i = n1; i < nz; ++i) {
        z[i] = x2_input[i-n1];
    }
    Rcpp::NumericVector len0v (1);
    double a0 = arma::min(z);
    double a1 = arma::max(z);
    len0v[0] = (a1 - a0) / incr;
    arma::vec len1v = Rcpp::round(len0v,0);
    int t_ncol = len1v[0];
    int h_nrow = t_ncol / 2;
    arma::vec t_ (t_ncol);
    arma::vec h_ (h_nrow);
    for (int i = 0; i < t_ncol; ++i) {
        t_[i] = a0 + i * incr;
    }
    for (int i = 0; i < h_nrow; ++i) {
        h_[i] = incr + i * incr;
    }
    arma::cube out (h_nrow , t_ncol, n_perm, arma::fill::value(9999.0));

    for (int p = 0; p < n_perm ; ++p) {
        arma::vec y = arma::shuffle(z);
        arma::vec x1 (n1);
        arma::vec x2 (n2);
        for (int k = 0; k < n1; ++k) {
            x1[k] = y[k];
        }
        for (int k = 0; k < n2; ++k) {
            x2[k] = y[k+n1];
        }
        for (int ti = 0; ti < t_ncol; ++ti) {
            for (int hi = 0; hi < h_nrow; ++hi) {
                if ((t_[ti] - h_[hi] >= a0) & (t_[ti] + h_[hi] <= a1)) {
                    int n_x1a = 0;
                    int n_x2a = 0;
                    for (int i = 0; i < n1; ++i) {
                        if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                            n_x1a = n_x1a + 1;
                        }
                    }
                    for (int i = 0; i < n2; ++i) {
                        if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                            n_x2a = n_x2a + 1;
                        }
                    }
                    if ((n_x1a > 0) & (n_x2a > 0)) {
                        Rcpp::NumericVector x1a (n_x1a);
                        Rcpp::NumericVector x2a (n_x2a);
                        int k = 0;
                        for (int i = 0; i < n1; ++i) {
                            if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                                x1a[k] = x1[i];
                                k = k + 1;
                            }
                        }
                        k = 0;
                        for (int i = 0; i < n2; ++i) {
                            if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                                x2a[k] = x2[i];
                                k = k + 1;
                            }
                        }    
                        double nn1a = x1a.length();
                        double nn2a = x2a.length();
                        Rcpp::NumericVector sqrtnna (1);
                        sqrtnna[0] = (nn1a * nn2a) / (nn1a + nn2a);
                        Rcpp::NumericVector sqrtna = Rcpp::sqrt(sqrtnna);
                        double bef = ks_Rcpp(x1a,x2a) * sqrtna[0];
                        out(hi,ti,p) = bef;
                    }
                }
            }
        }
    }
    return out;
}




//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
Rcpp::List idt_wass_stat_arma(const arma::vec & x1, const arma::vec & x2, int scaled, double scaled_by, double quad_by, double incr) {
    int n1 = x1.n_elem;
    int n2 = x2.n_elem;
    int nz = n1 + n2;
    arma::vec z (nz);
    for (int i = 0; i < n1; ++i) {
        z[i] = x1[i];
    }
    for (int i = n1; i < nz; ++i) {
        z[i] = x2[i-n1];
    }
    Rcpp::NumericVector len0v (1);
    double a0 = arma::min(z);
    double a1 = arma::max(z);
    len0v[0] = (a1 - a0) / incr;
    arma::vec len1v = Rcpp::round(len0v,0);
    int t_ncol = len1v[0];
    int h_nrow = t_ncol / 2;
    arma::vec t_ (t_ncol);
    arma::vec h_ (h_nrow);
    for (int i = 0; i < t_ncol; ++i) {
        t_[i] = a0 + i * incr;
    }
    for (int i = 0; i < h_nrow; ++i) {
        h_[i] = incr + i * incr;
    }
    arma::mat out (h_nrow , t_ncol, arma::fill::value(9999.0));
    for (int ti = 0; ti < t_ncol; ++ti) {
        for (int hi = 0; hi < h_nrow; ++hi) {
            if ((t_[ti] - h_[hi] >= a0) & (t_[ti] + h_[hi] <= a1)) {
                int n_x1a = 0;
                int n_x2a = 0;
                for (int i = 0; i < n1; ++i) {
                    if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                        n_x1a = n_x1a + 1;
                    }
                }
                for (int i = 0; i < n2; ++i) {
                    if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                        n_x2a = n_x2a + 1;
                    }
                }
                if (((n_x1a > 0) & (n_x2a > 0))) {
                    Rcpp::NumericVector x1a (n_x1a);
                    Rcpp::NumericVector x2a (n_x2a);
                    int k = 0;
                    for (int i = 0; i < n1; ++i) {
                        if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                            x1a[k] = x1[i];
                            k = k + 1;
                        }
                    } 
                    k = 0;
                    for (int i = 0; i < n2; ++i) {
                        if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                            x2a[k] = x2[i];
                            k = k + 1;
                        }
                    }    
                    double sc_by = 0;
                    if (scaled_by < 0) {
                        sc_by = h_[hi];
                    } else {
                        sc_by = scaled_by;
                    }
                    double nn1a = x1a.length();
                    double nn2a = x2a.length();
                    Rcpp::NumericVector sqrtnna (1);
                    sqrtnna[0] = (nn1a * nn2a) / (nn1a + nn2a);
                    Rcpp::NumericVector sqrtna = Rcpp::sqrt(sqrtnna);
                    double bef = wass_R1_armaV2(x1a,x2a,scaled,sc_by,quad_by) * sqrtna[0];

                    out(hi,ti) = bef;
                }
            }
        }
    }
    return Rcpp::List::create(out, t_, h_);
}




//' @useDynLib arkhaia
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
arma::cube idt_wass_stat_perm_arma(const arma::vec & x1_input, const arma::vec & x2_input, int scaled, double scaled_by, double quad_by, double incr, int n_perm) {
    int n1 = x1_input.n_elem;
    int n2 = x2_input.n_elem;
    int nz = n1 + n2;
    arma::vec z (nz);
    for (int i = 0; i < n1; ++i) {
        z[i] = x1_input[i];
    }
    for (int i = n1; i < nz; ++i) {
        z[i] = x2_input[i-n1];
    }
    Rcpp::NumericVector len0v (1);
    double a0 = arma::min(z);
    double a1 = arma::max(z);
    len0v[0] = (a1 - a0) / incr;
    arma::vec len1v = Rcpp::round(len0v,0);
    int t_ncol = len1v[0];
    int h_nrow = t_ncol / 2;
    arma::vec t_ (t_ncol);
    arma::vec h_ (h_nrow);
    for (int i = 0; i < t_ncol; ++i) {
        t_[i] = a0 + i * incr;
    }
    for (int i = 0; i < h_nrow; ++i) {
        h_[i] = incr + i * incr;
    }
    arma::cube out (h_nrow , t_ncol, n_perm, arma::fill::value(9999.0));

    for (int p = 0; p < n_perm ; ++p) {
        arma::vec y = arma::shuffle(z);
        arma::vec x1 (n1);
        arma::vec x2 (n2);
        for (int k = 0; k < n1; ++k) {
            x1[k] = y[k];
        }
        for (int k = 0; k < n2; ++k) {
            x2[k] = y[k+n1];
        }
        for (int ti = 0; ti < t_ncol; ++ti) {
            for (int hi = 0; hi < h_nrow; ++hi) {
                if ((t_[ti] - h_[hi] >= a0) & (t_[ti] + h_[hi] <= a1)) {
                    int n_x1a = 0;
                    int n_x2a = 0;
                    for (int i = 0; i < n1; ++i) {
                        if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                            n_x1a = n_x1a + 1;
                        }
                    }
                    for (int i = 0; i < n2; ++i) {
                        if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                            n_x2a = n_x2a + 1;
                        }
                    }
                    if (((n_x1a > 0) & (n_x2a > 0)) ) {
                        Rcpp::NumericVector x1a (n_x1a);
                        Rcpp::NumericVector x2a (n_x2a);
                        int k = 0;
                        for (int i = 0; i < n1; ++i) {
                            if ((t_[ti] - h_[hi] < x1[i]) & (x1[i] < t_[ti] + h_[hi])) {
                                x1a[k] = x1[i];
                                k = k + 1;
                            }
                        } 
                        k = 0;
                        for (int i = 0; i < n2; ++i) {
                            if ((t_[ti] - h_[hi] < x2[i]) & (x2[i] < t_[ti] + h_[hi])) {
                                x2a[k] = x2[i];
                                k = k + 1;
                            }
                        }    
                        double sc_by = 0;
                        if (scaled_by < 0) {
                           sc_by = h_[hi];
                        } else {
                            sc_by = scaled_by;
                        }
                        double nn1a = x1a.length();
                        double nn2a = x2a.length();
                        Rcpp::NumericVector sqrtnna (1);
                        sqrtnna[0] = (nn1a * nn2a) / (nn1a + nn2a);
                        Rcpp::NumericVector sqrtna = Rcpp::sqrt(sqrtnna);
                        double bef = wass_R1_armaV2(x1a,x2a,scaled,sc_by,quad_by) * sqrtna[0];
                        out(hi,ti,p) = bef ;
                    }
                }
            }
        }
    }
    return out;
}

