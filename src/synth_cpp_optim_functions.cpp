//synth_cpp_optim_functions.cpp

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//C++ version of the chunkmult function from kernlab
arma::mat chunkmultCpp(arma::mat Z, double csize, arma::vec colscale) {
  int n = Z.n_rows;
  int m = Z.n_cols;
  arma::vec d = sqrt(colscale);
  double nchunks = ceil(m / csize);
  //the results matix
  arma::mat res = arma::zeros<arma::mat>(n, n);

  for (int i = 1; i <= nchunks; i++) {
    double lowerb = (i - 1) * csize + 1;
    double upperb = i * csize;
    if (upperb > m) upperb = m;
    arma::mat buffer = trans(Z.cols(lowerb - 1, upperb - 1));
    arma::mat bufferd = d(arma::span(lowerb - 1, upperb - 1));
    buffer = buffer / bufferd;
    res = res + buffer.t() * buffer;
  }
  return res;
}

//Helper Function to convert to NumericVector, apply pmax, with bound,
//and return and arma vector
arma::vec pmax_arma(arma::vec x, double bound) {
  int length_x = x.n_elem;
  for (int i = 0; i < length_x; i++) {
    if (x(i) < bound) x(i) = bound;
  }
  return x;
}

//ipopCpp --For further documentation see the ipopCpp package on
//github
List ipopCpp(arma::vec c, arma::mat H, arma::mat A, arma::vec b,
	     arma::vec l, arma::vec u, arma::vec r,
	     int sigf = 7, int maxiter = 40, double margin = 0.05,
	     double bound = 10) {

  //The number of rows in H
  int n = H.n_rows;
  int cols = H.n_cols;
  double inv_tol = 1.0e-025;

  //Check for a decomposed matrix
  int smw = 0;
  if (n == cols) {
    smw = 0;
  }
  if (n > cols) {
    smw = 1;
  }
  if (n < cols) {
    smw = 1;
    n = cols;
    H = H.t();
  }

  //The number of rows in A
  int m = A.n_rows;

  //primal vector
  arma::vec primal(n);
  primal.zeros();

  int nrows_c = c.n_rows;
  //Make sure vectors have the right dimensions
  if (n != nrows_c) {
    stop("H and c ar incompatible");
  }
  int ncols_A = A.n_cols;
  if (n != ncols_A) {
    stop("A and c are incompatible");
  }
  int nrows_b = b.n_rows;
  if (m != nrows_b) {
    stop("A and b are incompatible");
  }
  int nrows_u = u.n_rows;
  if (n != nrows_u) {
    stop("u is incompatible with H");
  }
  int nrows_l = l.n_rows;
  if (n != nrows_l) {
    stop("l is incompatible with H");
  }

  //Update n
  n = A.n_cols;
  //Get the diagonals of H
  arma::vec H_diag = H.diag();

  arma::mat H_x;
  if (smw == 0) {
    H_x = H;
  } else if (smw == 1) {
    H_x = H.t();
  }

  //b+1 and c+1
  double b_plus_1 = max(svd(b)) + 1;
  double c_plus_1 = max(svd(c)) + 1;

  //starting point of the optimization
  int smwn = 0;
  if (smw == 0) {
    //replace the diagonals of H_x with h_diag + 1
    H_x.diag() = H_diag + 1;
  } else {
    smwn = H.n_cols;
  }

  //Generate and identity matrix with m columns
  arma::mat H_y = arma::eye<arma::mat>(m, m);
  arma::vec c_x = c;
  arma::vec c_y = b;
  // solve the system [-H.x A' A H.y] [x, y] = [c.x c.y]
  //initialize vec x and vec y and other variables
  arma::vec x;
  arma::vec y;
  arma::mat AP;
  arma::mat smwinner;
  arma::mat smwa1;   arma::mat smwc1; arma::mat smwa2; arma::mat smwc2;
  if (smw == 0 ) {
    //Note: matix indices START AT ZER0.
    AP = arma::zeros<arma::mat>(m+n,m+n);
    //Add in -H_x for the first n x n matrix of AP 0 to n for
    //both rows and columns
    AP(arma::span(0, n - 1),  arma::span(0, n - 1)) = -1 * H_x;
    //The last m rows and the first n columns will be set to A
    AP(arma::span(n, n + m - 1), arma::span(0, n - 1)) = A;
    //The first n columsn and the last m colums will set to A'
    AP(arma::span(0, n - 1), arma::span(n, n + m - 1)) = A.t();
    //The last m rows and the last m columns set to H_y
    AP(arma::span(n, n + m - 1), arma::span(n, n + m - 1)) = H_y;
    //Solve s_tmp = AP^{-1} %*% c(c.x,c.y)
    arma::vec cvec = join_vert(c_x, c_y);  // stack the c_x and c_y vectors
    // Solve
    arma::mat AP_pinv = pinv(AP, inv_tol);
    arma::vec s_tmp = AP_pinv *  cvec;
    //arma::vec s_tmp  = arma::solve(AP, cvec);

    // Get x and y
    x = s_tmp.rows(0, n - 1);
    y = s_tmp.rows(n, n + m - 1);

  } else {
    arma::mat V = arma::eye<arma::mat>(smwn, smwn);
    arma::mat crossH = cross(H, H);
    arma::mat smwinner = V + chol(H.t() * H);
    arma::mat smwa1 = A.t();
    arma::mat smwc1 = c_x;
    // arma::mat smwa2 = smwa1 -
    //   (H * arma::solve(smwinner, arma::solve(smwinner.t(), (H.t() * smwa1))));
    arma::mat smwinner_pinv = arma::pinv(smwinner, inv_tol);
    arma::mat smwa2 = smwa1 -
      (smwinner_pinv * (smwinner_pinv.t() * (H.t() * smwa1)));
    // arma::mat smwc2 = smwc1 -
    //   (H * arma::solve(smwinner, arma::solve(smwinner.t(), (H.t() * smwc1))));
    arma::mat smwc2 = smwc1 -
      (H * (smwinner_pinv * (smwinner_pinv.t() * (H.t() * smwc1))));
    // y = arma::solve(A * smwa2 + H_y, c_y + A * smwc2);
    y = arma::pinv(A * smwa2 + H_y, inv_tol) * c_y + A * smwc2;
    x = smwa2 * y - smwc2;
  }

  //Create the slack conditions. Note to create the slack conditions, we
  //need to use Rcpp syntactic sugar pmax() function and thus must convert
  //x and y to Rcpp NumericVectors. This will be done with the above pmax_arma()
  //function
  arma::vec g = pmax_arma(abs(x - l), bound);
  arma::vec z = pmax_arma(abs(x), bound);
  arma::vec t = pmax_arma(abs(u - x), bound);
  arma::vec s = pmax_arma(abs(x), bound);
  arma::vec v = pmax_arma(abs(y), bound);
  arma::vec w = pmax_arma(abs(y), bound);
  arma::vec p = pmax_arma(abs(r - w), bound);
  arma::vec q = pmax_arma(abs(y), bound);
  //Define mu -- a double
  double mu = as_scalar((z.t() * g + v.t() * w + s.t() * t + p.t() * q) / (2 * (m + n)));

  //flag variables for the optimization
  double sigfig = 0;
  // double old_sigfig = 0;
  int counter = 0;
  double alfa = 1;

  //initialize other variables
  arma::mat rho; arma::mat nu; arma::mat tau; arma::mat alpha; arma::mat sigma;
  arma::mat beta;
  // The while look for optimization
  while (counter < maxiter) {
    // update the counter
    counter++;
    //central path (predictor)
    arma::mat H_dot_x;
    if (smw == 0) {
      H_dot_x = H * x;
    } else if (smw == 1) {
      H_dot_x = H * (H.t() * x);
    }
    rho =  b - A * x + w;
    nu  = l - x + g;
    tau = u - x - t;
    alpha = r - w - p;
    sigma = c - A.t() * y - z + s + H_dot_x;
    beta = y + q - v;
    arma::vec gamma_z = -1 * z;
    arma::vec gamma_w = -1 * w;
    arma::vec gamma_s = -1 * s;
    arma::vec gamma_q = -1 * q;

    //instrumentation
    arma::mat x_dot_H_dot_x = x.t() * H_dot_x;
    //In the kernlab::ipop() function, they calculate
    // the primal and dual infeasibility in each iterate for
    //reporting. This program will only calculate them one
    //and report below.
    // double primal_infeasibility = max(svd(join_vert(join_vert(join_vert(rho, tau), alpha), nu))) / b_plus_1;
    // double dual_infeasibility = max(svd(join_vert(sigma, beta))) / c_plus_1;
    arma::mat primal_obj = c.t() * x + 0.5 * x_dot_H_dot_x;
    arma::mat dual_obj = b.t() * y - 0.5 * x_dot_H_dot_x + l.t() * z - u.t() * s - r.t() * q;
    //old_sigfig = sigfig;  // not used again in ipop() from kernlab
    //Note that if the object is a matrix, the armadillo as_scalar() funcion will return a
    //double.
    sigfig = as_scalar(-1*log10(abs(primal_obj - dual_obj) / (abs(primal_obj) + 1)));



    if (sigfig <= 0) sigfig = 0;
    if (sigfig >= sigf) break;

    //Some more intermediate variables (the hat section)
    arma::mat hat_beta = beta - v % gamma_w / w;
    arma::mat hat_alpha = alpha - p * gamma_q / q;
    arma::mat hat_nu = nu + g % gamma_z / z;
    arma::mat hat_tau = tau - t % gamma_s / s;

    // The diagonal terms
    arma::mat d = z / g + s / t;
    arma::mat e = 1 / (v / w + q / p);

    //The initialization before the big cholesky
    if (smw == 0) {
      H_x.diag() = H_diag + d;
    }
    H_y.diag() = e;
    c_x = sigma - z % hat_nu / g - s % hat_tau / t;
    c_y = rho - e % (hat_beta - q % hat_alpha / p);



    // and solve the system [-H.x A' A H.y] [delta.x, delta.y] <- [c.x c.y]
    //Add in -H_x for the first n x n matrix of AP 0 to n for
    //both rows and columns
    arma::vec delta_x;
    arma::vec delta_y;
    if (smw == 0) {
      AP(arma::span(0, n - 1),  arma::span(0, n - 1)) = -1 * H_x;
      //The last m rows and the last m columns set to H_y
      AP(arma::span(n, n + m - 1), arma::span(n, n + m - 1)) = H_y;
      //Solve s_tmp = AP^{-1} %*% c(c.x,c.y)
      arma::vec cvec = join_vert(c_x, c_y);  // stack the c_x and c_y vectors
      // Solve
      //arma::vec s1_tmp  = arma::solve(AP, cvec);
      arma::vec s1_tmp = arma::pinv(AP, inv_tol) * cvec;
      delta_x = s1_tmp.rows(0, n - 1);
      delta_y = s1_tmp.rows(n, n + m - 1);
    } else {
      arma::mat V = arma::eye(smwn, smwn);
      smwinner = chol(V + chunkmultCpp(H.t(), 2000, d));
      arma::mat smwa1 = A.t();
      smwa1  = smwa1 / d;
      arma::mat smwc1 = c_x / d;
      //Use pinv to allow for tolerance in matrix inversion
      arma::mat smwinner_pinv = arma::pinv(smwinner, inv_tol);
      arma::mat smwa2 = A.t() - (smwinner_pinv * (smwinner_pinv.t() * (H.t() * smwa1)));
      // arma::mat smwa2 = A.t() -
      // 	(H * arma::solve(smwinner,
      // 			 arma::solve(smwinner.t(),
      // 				     H.t() * smwa1,
      // 				     )
      // 			 ));
      smwa2 = smwa2 / d;
      arma::mat smwc2 = c_x - (smwinner_pinv * (smwinner_pinv.t() * (H.t() * smwc1)));
      smwc2 = smwc2 / d;

      // arma::mat smwc2 = (c_x -
      // 			 (H * arma::solve(smwinner,
      // 					  arma::solve(smwinner.t(),
      // 						      H.t() * smwc1
      // 						      )
      // 					  )
      // 			  )) / d;
      delta_y = arma::pinv(A * smwa2 + H_y, inv_tol) * (c_y + A * smwc2);
      //delta_y = arma::solve(A * smwa2 + H_y, c_y + A * smwc2);
      delta_x = smwa2 * delta_y - smwc2;
    }

    //backsubstitution
    arma::vec delta_w  = -1*e % (hat_beta - q % hat_alpha / p + delta_y);
    arma::vec delta_s = s % (delta_x - hat_tau) / t;
    arma::vec delta_z = z % (hat_nu - delta_x) / g;
    arma::vec delta_q = q % (delta_w - hat_alpha) / p;
    arma::vec delta_v = v % (gamma_w - delta_w) / w;
    arma::vec delta_p = p % (gamma_q - delta_q) / q;
    arma::vec delta_g = g % (gamma_z - delta_z) / z;
    arma::vec delta_t = t % (gamma_s - delta_s) / s;
    // compute update step now (sebastian's trick)
    double alfa_num = -1 * (1 - margin);
    arma::vec alfa_denom = join_vert(delta_g / g, delta_w / w);
    alfa_denom = join_vert(alfa_denom, delta_t / t);
    alfa_denom = join_vert(alfa_denom, delta_p / p);
    alfa_denom = join_vert(alfa_denom, delta_z / z);
    alfa_denom = join_vert(alfa_denom, delta_s / s);
    alfa_denom = join_vert(alfa_denom, delta_q / q);
    double alfa_denom_double = min(alfa_denom);
    if (-1 < alfa_denom_double) alfa_denom_double = -1;
    alfa = alfa_num / alfa_denom_double;
    double newmu = as_scalar((z.t() * g + v.t() * w + s.t() * t + p.t() * q)/ (2 * (m + n)));
    newmu = mu * pow(((alfa - 1) / (alfa + 10)),2);
    gamma_z = mu / g - z - delta_z % delta_g / g;
    gamma_w = mu / v - w - delta_w % delta_v / v;
    gamma_s = mu / t - s - delta_s % delta_t / t;
    gamma_q = mu / p - q - delta_q % delta_p / p;
    // Some more intermediate variables (the hat section)
    hat_beta = beta - v % gamma_w / w;
    hat_alpha = alpha - p % gamma_q / q;
    hat_nu = nu + g % gamma_z / z;
    hat_tau = tau - t % gamma_s / s;
    // //initialization before the big cholesky
    c_x = sigma - z % hat_nu / g - s % hat_tau / t;
    c_y = rho - e % (hat_beta - q % hat_alpha / p);

    // and solve the system [-H.x A' A H.y] [delta.x, delta.y] <- [c.x c.y]
    if (smw == 0) {
      //Add in -H_x for the first n x n matrix of AP 0 to n for
      //both rows and columns
      AP(arma::span(0, n - 1),  arma::span(0, n - 1)) = -1 * H_x;
      //The last m rows and the last m columns set to H_y
      AP(arma::span(n, n + m - 1), arma::span(n, n + m - 1)) = H_y;
      // Solve
      // arma::vec s1_tmp  = arma::solve(AP, join_vert(c_x, c_y));
      arma::vec s1_tmp  = arma::pinv(AP, inv_tol) * join_vert(c_x, c_y);
      delta_x = s1_tmp.rows(0, n - 1);
      delta_y = s1_tmp.rows(n, n + m - 1);
    } else if (smw == 1) {
      smwc1 = c_x / d;
      arma::mat smwinner_pinv = arma::pinv(smwinner, inv_tol);
      // smwc2 = (c_x -
      // 	       (H * arma::solve(smwinner, arma::solve(smwinner.t(), H.t() * smwc1)))) / d;
      smwc2 = (c_x -
	       (H * (smwinner_pinv * (smwinner_pinv.t() * (H.t() * smwc1))))) / d;
      // delta_y = arma::solve(A * smwa2 + H_y, c_y + A * smwc2);
      delta_y = arma::pinv(A * smwa2 + H_y, inv_tol) * (c_y + A * smwc2);
      delta_x = smwa2 * delta_y - smwc2;
    }

    //backsubtitution
    delta_w  = -1*e % (hat_beta - q % hat_alpha / p + delta_y);
    delta_s = s % (delta_x - hat_tau) / t;
    delta_z = z % (hat_nu - delta_x) / g;
    delta_q = q % (delta_w - hat_alpha) / p;
    delta_v = v % (gamma_w - delta_w) / w;
    delta_p = p % (gamma_q - delta_q) / q;
    delta_g = g % (gamma_z - delta_z) / z;
    delta_t = t % (gamma_s - delta_s) / s;
    // compute the updates
    alfa_num = -1 * (1 - margin);
    alfa_denom = join_vert(delta_g / g, delta_w / w);
    alfa_denom = join_vert(alfa_denom, delta_t / t);
    alfa_denom = join_vert(alfa_denom, delta_p / p);
    alfa_denom = join_vert(alfa_denom, delta_z / z);
    alfa_denom = join_vert(alfa_denom, delta_s / s);
    alfa_denom = join_vert(alfa_denom, delta_q / q);
    alfa_denom_double = min(alfa_denom);
    if (-1 < alfa_denom_double) alfa_denom_double = -1;
    alfa = alfa_num / alfa_denom_double;
    x += delta_x * alfa;
    g += delta_g * alfa;
    w += delta_w * alfa;
    t += delta_t * alfa;
    p += delta_p * alfa;
    y += delta_y * alfa;
    z += delta_z * alfa;
    v += delta_v * alfa;
    s += delta_s * alfa;
    q += delta_q * alfa;

    mu = newmu;

  }




  //primal and dual infeasiblity
  double primal_infeasibility = max(svd(join_vert(join_vert(join_vert(rho, tau), alpha), nu))) / b_plus_1;
    double dual_infeasibility = max(svd(join_vert(sigma, beta))) / c_plus_1;

  //the convergence
  std::string convergence;
  if ((sigfig > sigf) & (counter < maxiter)) {
    convergence  = "converged";
  } else {
    if ((primal_infeasibility > 10e5) & (dual_infeasibility > 10e5)) {
      convergence = "primal and dual infeasible";
    }
    if (primal_infeasibility > 10e5) {
      convergence = "primal infeasible";
    }
    if (dual_infeasibility > 10e5) {
      convergence = "dual infeasible";
    } else {
      convergence = "slow convergence, change bound?";
    }
  }

  //return a list with the results
  List ret;
  ret["primal"] = x;
  ret["dual"] = as_scalar(y);
  ret["convergence"] = convergence;
  return ret;
}

//' To get the solution for the weights from the weights on each
//' predictor variable and the X scaled matrices for the treated unit
//' and control units
//'
//' @title solution for region weights
//' @param solution_v the weighte on each predictor
//' @param X0_scaled scaled matrix of control units
//' @param X1_scaled scaled matrix for the treated unit
//' @return A vector with the solution with the regional weights
//' @author Chandler Lutz
//'
// [[Rcpp::export]]
arma::vec solution_w_cpp(arma::vec solution_v, arma::mat X0_scaled, arma::mat X1_scaled) {

  //Identity matrix with length V
  arma::mat V_mat = arma::diagmat(solution_v);

  //H matrix
  arma::mat H_mat = X0_scaled.t() * V_mat * X0_scaled;

  //c vector
  //arma::mat c_mat = -1 * (X1_scaled.t() * V_mat * X0_scaled);

  arma::vec c_vec = -1 * vectorise(X1_scaled.t() * V_mat * X0_scaled);

  //A matrix -- row unit vector with length c
  arma::mat A_mat = arma::ones<arma::mat>(1, c_vec.n_elem);

  //b vector
  arma::vec b_vec = arma::eye<arma::vec>(1, 1);

  //l vector
  arma::vec l_vec = arma::zeros<arma::vec>(c_vec.n_elem, 1);

  //u vector
  arma::vec u_vec = arma::ones<arma::vec>(c_vec.n_elem, 1);

  //r vector
  arma::vec r_vec = arma::zeros<arma::vec>(1, 1);

  Rcpp::List ipop_results = ipopCpp(c_vec, H_mat, A_mat,
				    b_vec, l_vec, u_vec,
				    r_vec, 5, 1000,
				    0.0005, 10);

  //The the solution to the optimization problem which are the
  //weights in the synthetic control process
  arma::vec solution_out = ipop_results["primal"];

  // std::string conv_out = ipop_results["convergence"];
  // double dual_out = ipop_results["dual"];
  // Rcout << conv_out << std::endl;
  // Rcout << dual_out << std::endl;

  return solution_out;
}

//' The Optimization function for the synthetic control implemented
//' using Rcpp
//'
//' @title Synth Optimization func for predictor variable weights
//' @param variables_v The weights on the Synth predictor variables
//' @param X0_scaled scaled matrix of control units
//' @param X1_scaled scaled matrix for the treated unit
//' @param Z0 matrix of outcome observations for the control units
//' @param Z1 matrix of outcome observations for the treated unit
//' @return the loss from the parameter optimizaiton
//' @author Chandler Lutz
//'
// [[Rcpp::export]]
double fn_v_cpp(arma::vec variables_v,
		arma::mat X0_scaled, arma::mat X1_scaled,
		arma::mat Z0, arma::mat Z1) {

  //Absolute value of starting values
  variables_v = abs(variables_v);

  //The weights for each region
  arma::vec solution_w_temp = solution_w_cpp(variables_v, X0_scaled, X1_scaled);

  // The loss for v (the weights on each predictor variable)
  double loss_v = arma::as_scalar(arma::trans(Z1 - (Z0 * solution_w_temp)) *
				  (Z1 - (Z0 * solution_w_temp)));

  loss_v /= Z0.n_rows;

  return loss_v;

}
