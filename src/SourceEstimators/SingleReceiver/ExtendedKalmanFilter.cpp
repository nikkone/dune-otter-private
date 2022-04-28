// This code is written by Sverre Søbstad Løvskar as a part of his masters thesis named: 
// "Positioning of periodic acoustic emitters using an omnidirectional hydrophone on an AUV platform"
// The code has been modified by Nikolai Lauvås

#include "ExtendedKalmanFilter.h"

using namespace std;
using namespace arma;
namespace SourceEstimators
{
  //! Insert explanation on task behaviour here.
  //! @author Nikolai Lauvås
  namespace SingleReceiver
  {
ExtendedKalmanFilter::ExtendedKalmanFilter() {

	//vec P_diag = {250000.0, 160000.0, 1.0};
	vec Q_diag = {0.0000, 0.0000, 0.00001};
	vec R_diag = {1.5*3.8706, 2.1638e6};

	F = eye<mat>(3,3);
	Q = diagmat(Q_diag);
	R = diagmat(R_diag);
	R_rows = R.n_rows; // Matlab: size(R,1)
	x = {0.0,0.0,0.0};
	//P = diagmat(P_diag);
	P = { 	{66458, -26820, 0},
			{-26820, 12116, 0},
			{0, 	0, 		0}};
	k = 0;
	initialized = false;
	pre_init_data = 0;

}
        
ExtendedKalmanFilter::ExtendedKalmanFilter(vec X_0, mat P_0, mat Qm, mat Rm) {
	F = eye<mat>(3,3);
	Q = Qm;
	R = Rm;
	x = X_0;
	P = P_0;
	k = 0;
	R_rows = R.n_rows; // Matlab: size(R,1)
	initialized = false;
	pre_init_data = 0;
}
int ExtendedKalmanFilter::set_process_noise(mat Qm) {
	Q = Qm;
	return 1;
}
int ExtendedKalmanFilter::set_measurement_noise(mat Rm) {
	R = Rm;
	return 1;
}
int ExtendedKalmanFilter::get_steps() {
	return k;
}
int ExtendedKalmanFilter::initialize(vec pos) {
	x = pos;
	initialized = true;
	return 0;
}


void ExtendedKalmanFilter::measurementStep(vec ranging, vec position_previous, vec position_current, double rtoa) {
	if (!initialized) {
		if (rtoa == 0) {
			rtoa = 0.00000001;
		}
		//position_current(2) += ((double) std::rand() / RAND_MAX)/100;
		if (pre_init_data == 0) {
			posi = zeros<mat>(3,5);
			ToA = zeros<vec>(5);

			ToA(pre_init_data) = 0;
			posi.col(pre_init_data) = position_previous;
			ToA(pre_init_data) = 0;
			pre_init_data++;
			posi.col(pre_init_data) = position_current;
			ToA(pre_init_data) = ToA(pre_init_data-1) + rtoa;
			pre_init_data++;
		} else {
			posi.col(pre_init_data) = position_current;
			ToA(pre_init_data) = ToA(pre_init_data-1) + rtoa;
			pre_init_data++;
		}
		if (pre_init_data > 5-1) {
			try {
				batchSolve();
				initialized = true;
			} catch (...) {
				cout << "Ah, nope.. exception.." << endl;
			}
		}
		return;
	}
// Start calculate Jacobian of measurement model at X_k
	vec distance1 = x-position_previous; // X_e-X_rx0
	vec distance2 = x-position_current;  // X_e-X_rx1

	double r1 = norm(distance1,2); // Find euclidean norm (p-norm, p=2). ||X_e-X_rx0||
	double r2 = norm(distance2,2); // Find euclidean norm (p-norm, p=2). ||X_e-X_rx1||

	vec h(R_rows); // Measurement equation
	h(0) = r2 - r1; // h is eq (2.16) in masters
	mat H = (distance2.t()/r2) - (distance1.t())/r1; // The Jacobian of h found in eq (2.17)
	/* If two rows in R matrix, then the SNR is also used, so the new measurement function
	   is added to h and its Jacobian H */
	if (R_rows == 2) {
		h(1) = r2; // Eq (2.19) in masters
		H = join_vert(H,distance2.t()/r2); // Exends the Jacobian with eq (2.20)
	}

//Compute Kalman gain 
	// K = P_predict*H'/(H*P_predict*H' + R);
	// K = B*A^(-1) => KA=BI=B => A^T*K^T = B^T is in for that can be solved by solve() function (Ax=B)
	mat A = H*P*H.t() + R;
	mat B = P*H.t();
	K = solve(A.t(),B.t()); // K' = A' \ B'
	K = K.t(); // (K^T)^T=K
	
//Measurement update
	x = x + K*(ranging-h); // x_est = x_predict + K*(z-h);

// Compute error covariance for updated estimate
	P = (eye<mat>(3,3) - K*H)*P; // P_est = (eye(3) - K*H)*P_predict;


}
int ExtendedKalmanFilter::predictionStep(void) {
	if(initialized) {
		//Predict Step / Project ahead
		x = F*x; // x_predict = F*x_est;
		P = F*P*F.t() + Q; // P_predict = F*P_est*F' + Q;
		return ++k;
	}
	return k;
}

// TODO: Is this strictly a part of the Kalman Filter?
bool ExtendedKalmanFilter::batchSolve() {
	cout << "Batch solve time!" << endl;

	cout << "Posi = " << endl << posi << endl;
	cout << "ToA = " << endl << ToA << endl;
	mat M = zeros<mat>(3,pre_init_data); // Equivalent to [A B C] in masters eq (2.5)
	vec D = zeros<vec>(pre_init_data);

	cout << "Starting loop!" << endl;
	for (unsigned int m = 2; m < posi.n_cols; m++) {
    
	    double ddm = ToA(m);
	    double dd2 = ToA(1);
		
		M.col(m) = (2*(posi.col(m) - posi.col(0)) / ddm) - (2*(posi.col(1)-posi.col(0)) / dd2);
	    cout << "M.Col" << endl << M.col(m) << endl;
	    double off1 = arma::sum(posi.col(0)%posi.col(0));
	    double off2 = arma::sum(posi.col(1)%posi.col(1));
	    double offm = arma::sum(posi.col(m)%posi.col(m));
	    
	    D(m) = ddm - dd2 + (off1-offm)/ddm - (off1-off2)/dd2;
	   	cout << "Iter " << m << " ok" << endl;
	}
    cout << "merging matrices" << endl;
    cout << "Extracting from " << 2 << " to " << posi.n_cols-1 << endl;
	M = M.cols(2,posi.n_cols-1);
	D = -D.rows(2,posi.n_cols-1);
	cout << "attempting solution" << endl;
	cout << "M = " << endl << M << endl;
	cout << "D = " << endl << D << endl;
	if (solve(x,M.t(),D)) {
		cout 	<< "Solution: " << endl 
				<< "\tx(0) = " << (int) x(0) << endl
				<< "\tx(1) = " << (int) x(1) << endl
				<< "\tx(2) = " << (int) x(2) << endl;
		x(2) = -2;
		return true;
	} 
	cout << "Nope. Couldn't solve it." << endl;
	return false;
}

bool ExtendedKalmanFilter::isInitialized(void) {
	return initialized;
}
  }
}