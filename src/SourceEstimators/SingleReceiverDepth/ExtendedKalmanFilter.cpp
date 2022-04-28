// This code is based on the work performed by Sverre Søbstad Løvskar as a part of his masters thesis named: 
// "Positioning of periodic acoustic emitters using an omnidirectional hydrophone on an AUV platform"
// The code has been modified by Nikolai Lauvås

#include "ExtendedKalmanFilter.h"

using namespace std;
using namespace arma;
namespace SourceEstimators
{
  namespace SingleReceiverDepth
  {
ExtendedKalmanFilter::ExtendedKalmanFilter() {

	/*//vec P_diag = {250000.0, 160000.0, 1.0};
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
	k = 0;*/
	initialized = false;
	//pre_init_data = 0;

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
	//pre_init_data = 0;
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


void ExtendedKalmanFilter::measurementStep(vec z, vec position_previous, vec position_current) {
// Start calculate Jacobian of measurement model at X_k
	if(initialized) {
		vec distance1 = x-position_previous; // X_e-X_rx0
		vec distance2 = x-position_current;  // X_e-X_rx1

		double r1 = norm(distance1,2); // Find euclidean norm (p-norm, p=2). ||X_e-X_rx0||
		double r2 = norm(distance2,2); // Find euclidean norm (p-norm, p=2). ||X_e-X_rx1||

		vec h(R_rows); // Measurement equation
		//std::cout << "R" << R << std::endl;
		h(0) = r2 - r1; // h is eq (2.16) in masters
		mat H = (distance2.t()/r2) - (distance1.t())/r1; // The Jacobian found in eq (2.17)
		/* If two rows in R matrix, then the SNR is also used, so the new measurement function
		is added to h and its Jacobian H */
		if (R_rows > 1) {
			h(1) = r2; // Eq (2.19) in masters
			H = join_vert(H,distance2.t()/r2); // Exends the Jacobian with eq (2.20)
			std::cout << "using SNR"<< std::endl;
		}
		if (R_rows > 2) {
			h(2) = x(2);
			arma::rowvec Hdepth= {0.0,0.0,1.0};
			H = join_vert(H,Hdepth);
			std::cout << "using Depth"<< std::endl;
		}

	//Compute Kalman gain 
		// K = P_predict*H'/(H*P_predict*H' + R);
		// K = B*A^(-1) => KA=BI=B => A^T*K^T = B^T is in for that can be solved by solve() function (Ax=B)
		mat A = H*P*H.t() + R;
		mat B = P*H.t();
		K = solve(A.t(),B.t()); // K' = A' \ B'
		K = K.t(); // (K^T)^T=K
		
	//Measurement update
		x = x + K*(z-h); // x_est = x_predict + K*(z-h);

	// Compute error covariance for updated estimate
		P = (eye<mat>(3,3) - K*H)*P; // P_est = (eye(3) - K*H)*P_predict;
	}
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

bool ExtendedKalmanFilter::isInitialized(void) {
	return initialized;
}
  }
}