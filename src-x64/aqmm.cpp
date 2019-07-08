#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List C_hfun_aqmm(NumericVector y, NumericMatrix X, NumericMatrix B, NumericMatrix Z, NumericVector weights, NumericVector theta_x, NumericVector v, NumericMatrix u, NumericVector Phi_inv, NumericMatrix Sigma_inv, int M, int N, IntegerVector ni, int p, int Q, int H, int s, double tau, double omega){

	double fidelity = 0;
	double penalty = 0;
	NumericVector res(N);
	NumericVector Bv(N);
	NumericVector Zu(N);
	NumericVector w(N);
	double val = 0;
	double aa = (tau - 1)*omega;
	double bb = tau*omega;
	double Ajj = 0;
	double bvec = 0;
	double cvec = 0;

	// Penalty random effects
	double peni = 0;
	for(int i = 0; i < M; ++i){
		peni = 0;
		for(int j = 0; j < Q; ++j){
			for(int k = 0; k < Q; ++k){
				peni += u(i,k)*Sigma_inv(k,j)*u(i,j);
			}
		}
		penalty += peni*weights[i];
	}

	// Penalty smooth terms
	
	for(int i = 0; i < H; ++i){
		penalty += pow(v[i], 2)*Phi_inv[i];
	}

	// Bv
	
	for(int i = 0; i < N; ++i){
		for(int j = 0; j < H; ++j){
			Bv[i] += B(i,j) * v(j);
		}
	}
	
	// Zu
	
	int start = 0;
	int stop = 0;
	for(int i = 0; i < M; ++i){
		stop += ni[i];
		for(int j = start; j < stop; ++j){
			for(int k = 0; k < Q; ++k){
				Zu[j] += Z(j,k) * u(i,k);
			}
			w[j] = weights[i];
		}
		start += ni[i];
	}

	// y - X*theta_x - Zu - Bv
	
	for (int i = 0; i < N; ++i){
		double xb = 0;
		for (int j = 0; j < p; ++j){
			xb += X(i,j)*theta_x[j];
			}
		res[i] = y[i] - xb - Bv[i] - Zu[i];
	}
	
	// fidelity
	
	for(int i = 0; i < N; ++i){
		int s = 0;
		if(res[i] <= aa){
			s = -1;
		} else if(res[i] >= bb){
			s = 1;
		}
		Ajj = (1 - pow(s, 2))/omega;
		bvec = s*((2*tau - 1)*s + 1);
		cvec = 0.5*(1-2*tau)*omega*s - 0.5*(1-2*tau+2*pow(tau,2))*omega*pow(s,2);
		
		fidelity += (pow(res[i], 2)*Ajj + bvec*(res[i]) + cvec)*w[i];
	}

	val = fidelity + penalty;
	//Rcout << "val = " << val << std::endl;

	List ans;
	ans["val"] = val;
	ans["resid"] = res;
	return ans;
}

// [[Rcpp::export]]
List C_hfunD_aqmm(NumericVector y, NumericMatrix X, NumericMatrix B, NumericMatrix Z, NumericVector weights, NumericVector theta_x, NumericVector v, NumericMatrix u, NumericVector Phi_inv, NumericMatrix Sigma_inv, int M, int N, IntegerVector ni, int p, int Q, int H, int s, double tau, double omega){

	double fidelity = 0;
	double penalty = 0;
	NumericVector res(N);
	NumericVector Bv(N);
	NumericVector Zu(N);
	NumericVector w(N);
	NumericVector gradient(H+M*Q);
	NumericVector phiv(H);
	NumericVector sigmau(M*Q);
	NumericMatrix hessian(H+M*Q,H+M*Q);
	double val = 0;
	double aa = (tau - 1)*omega;
	double bb = tau*omega;
	NumericVector A(N);
	double bvec = 0;
	double cvec = 0;

	// Penalty random effects
	double peni = 0;
	for(int i = 0; i < M; ++i){
		peni = 0;
		for(int j = 0; j < Q; ++j){
			for(int k = 0; k < Q; ++k){
				peni += u(i,k)*Sigma_inv(k,j)*u(i,j);
				sigmau[j + (Q*i)] += 2*Sigma_inv(j,k)*u(i,k)*weights[i];
			}
		}
		penalty += peni*weights[i];
	}
	
	// Penalty smooth terms
	
	for(int i = 0; i < H; ++i){
		penalty += pow(v[i], 2)*Phi_inv[i];
		phiv[i] = 2*v[i]*Phi_inv[i];
	}

	// Bv
	
	for(int i = 0; i < N; ++i){
		for(int j = 0; j < H; ++j){
			Bv[i] += B(i,j) * v(j);
		}
	}
	
	// Zu
	
	int start = 0;
	int stop = 0;
	for(int i = 0; i < M; ++i){
		stop += ni[i];
		for(int j = start; j < stop; ++j){
			for(int k = 0; k < Q; ++k){
				Zu[j] += Z(j,k) * u(i,k);
			}
			w[j] = weights[i];
		}
		start += ni[i];
	}

	// y - X*theta_x - Zu - Bv
	
	for (int i = 0; i < N; ++i){
		double xb = 0;
		for (int j = 0; j < p; ++j){
			xb += X(i,j)*theta_x[j];
			}
		res[i] = y[i] - xb - Bv[i] - Zu[i];
	}
	
	// fidelity
	
	for(int i = 0; i < N; ++i){
		int s = 0;
		if(res[i] <= aa){
			s = -1;
		} else if(res[i] >= bb){
			s = 1;
		}
		A[i] = (1 - pow(s, 2))/omega;
		bvec = s*((2*tau - 1)*s + 1);
		cvec = 0.5*(1-2*tau)*omega*s - 0.5*(1-2*tau+2*pow(tau,2))*omega*pow(s,2);
		
		fidelity += (pow(res[i], 2)*A[i] + bvec*(res[i]) + cvec)*w[i];
		res[i] = (2*A[i]*res[i] + bvec)*w[i];
	}
	
	// gradient
	
	for(int j = 0; j < H; ++j){
		for(int i = 0; i < N; ++i){
			gradient[j] += -B(i,j)*res[i];
		}
		gradient[j] += phiv[j];
	}

	start = 0;
	stop = 0;
	for(int i = 0; i < M; ++i){
		stop += ni[i];
		for(int k = 0; k < Q; ++k){
			for(int j = start; j < stop; ++j){
				gradient[H + k + (Q*i)] += -Z(j,k)*res[j];
			}
		gradient[H + k + (Q*i)] += sigmau[k + (Q*i)];
		}
		start += ni[i];
	}
	

	// return value, gradient, Hessian
	
	val = fidelity + penalty;
	List ans;
	ans["val"] = val;
	ans["gradient"] = gradient;
	ans["weights"] = A;
	return ans;
}