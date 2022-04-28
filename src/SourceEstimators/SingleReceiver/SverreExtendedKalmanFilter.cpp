#include "SverreExtendedKalmanFilter.h"

using namespace std;
using namespace arma;
namespace SourceEstimators
{
  //! Insert explanation on task behaviour here.
  namespace SingleReceiver
  {
SverreExtendedKalmanFilter::SverreExtendedKalmanFilter() {

	//vec P_diag = {250000.0, 160000.0, 1.0};
	vec Q_diag = {0.0000, 0.0000, 0.00001};
	vec R_diag = {1.5*3.8706, 2.1638e6};

	F = eye<mat>(3,3);
	Q = diagmat(Q_diag);
	R = diagmat(R_diag);
	x = {0.0,0.0,0.0};
	//P = diagmat(P_diag);
	P = { 	{66458, -26820, 0},
			{-26820, 12116, 0},
			{0, 	0, 		0}};
	k = 0;
	initialized = false;
	pre_init_data = 0;

}
        
SverreExtendedKalmanFilter::SverreExtendedKalmanFilter(vec X_0, mat P_0, mat Qm, mat Rm) {
	F = eye<mat>(3,3);
	Q = Qm;
	R = Rm;
	x = X_0;
	P = P_0;
	k = 0;
	initialized = false;
	pre_init_data = 0;
}
int SverreExtendedKalmanFilter::set_process_noise(mat Qm) {
	Q = Qm;
	return 1;
}
int SverreExtendedKalmanFilter::set_measurement_noise(mat Rm) {
	R = Rm;
	return 1;
}
int SverreExtendedKalmanFilter::steps() {
	return k;
}
int SverreExtendedKalmanFilter::initialize(vec pos) {
	x = pos;
	initialized = true;
	return 0;


}

int SverreExtendedKalmanFilter::update(vec ranging, vec position_previous, vec position_current, double rtoa) {

	if (!initialized) {
		if (rtoa == 0) {
			rtoa = 0.00000001;
		}
		//position_current(2) += ((double) std::rand() / RAND_MAX)/100;
		if (pre_init_data == 0) {
			posi = zeros<mat>(3,5);
			ToA = zeros<vec>(5);

			ToA(pre_init_data) = 0; // WhY?
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
		return k;
	}

	int n = R.n_rows; // Matlab: size(R,1)
	//int n = R.n_cols; // Matlab: size(R,2)
	//cout << "Kalman::Update(): Predicting step" << endl << flush;
	//Predict Step
	x = F*x; // x_predict = F*x_est;
	P = F*P*F.t() + Q; // P_predict = F*P_est*F' + Q;

	//cout << "Kalman::Update(): Calc vectors" << endl << flush;
	vec distance1 = x-position_previous; // dist_meas1 = x_predict - pos1;
	vec distance2 = x-position_current; // dist_meas2 = x_predict - pos2;

	double r1 = norm(distance1,2);
	double r2 = norm(distance2,2);

	vec h(n);
	h(0) = r2 - r1; // h(1,1) = norm(dist_meas2) - norm(dist_meas1);

	//cout << "r1 =" << r1 << "  r2 =" << r2 << endl;
	mat H = (distance2.t()/r2) - (distance1.t())/r1; // H = ((dist_meas2')/norm(dist_meas2)) - ((dist_meas1')/norm(dist_meas1));
	//cout << "H(1) =" << endl << H << endl;
	if (n == 2) {				// if n == 2
		h(1) = r2; 				// h(2,1) = norm(dist_meas2);
		H = join_vert(H,distance2.t()/r2); // H = [H; (dist_meas2')/norm(dist_meas2)];
	}
	//cout << "H =" << endl << H << endl;
	//cout << "Kalman::Update(): Measurement Update" << endl << flush;
	//Measurement update
	//Compute Kalman gain 
	// K = P_predict*H'/(H*P_predict*H' + R);
	mat A = H*P*H.t() + R;
	//cout << "P =" << endl << P << endl;
	mat B = P*H.t();
	//cout << "A' =" <<endl<< A.t() << endl;
	//cout << "B' =" <<endl<< B.t() << endl;
	mat K = solve(A.t(),B.t()); // K = A' \ B'
	K = K.t();
	

	//cout << "New K: " << endl << K << endl << flush;
	//cout << "RangeDiff: " << endl <<  ranging << endl << flush;
	//cout << "h vect: " << endl << h << endl << flush;
	x = x + K*(ranging-h); // x_est = x_predict + K*(z-h);
	P = (eye<mat>(3,3) - K*H)*P; // P_est = (eye(3) - K*H)*P_predict;

	return ++k;
}


bool SverreExtendedKalmanFilter::batchSolve() {
	cout << "Sverre Batch solve time!" << endl;

	cout << "Posi = " << endl << posi << endl;
	cout << "ToA = " << endl << ToA << endl;
	mat M = zeros<mat>(3,pre_init_data);
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
  }
}