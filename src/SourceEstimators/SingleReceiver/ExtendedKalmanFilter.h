// This code is written by Sverre Søbstad Løvskar as a part of his masters thesis named: 
// "Positioning of periodic acoustic emitters using an omnidirectional hydrophone on an AUV platform"
// The code has been modified by Nikolai Lauvås


#ifndef EXTENDEDKALMANFILTER_PLAIN_H
#define EXTENDEDKALMANFILTER_PLAIN_H

#include <iostream>
#include <armadillo>
namespace SourceEstimators
{
  //! Insert explanation on task behaviour here.
  //! @author Nikolai Lauvås
  namespace SingleReceiver
  {
class ExtendedKalmanFilter {
    arma::mat F; // System / Plant
    arma::mat Q; // Process Noise
    arma::mat R; // Measurement Noise
	int R_rows; // Rows in R
    unsigned long k; // Steps
    bool initialized;

    arma::mat posi;
    arma::vec ToA;
    int pre_init_data;
public:
    ExtendedKalmanFilter();
    ExtendedKalmanFilter(arma::vec X, arma::mat P, arma::mat Q, arma::mat R);
    arma::mat x; // State Estimation
    arma::mat P; // Error Coef Matrix
    arma::mat K; // Kalman Gain

    int initialize(arma::vec pos);
    bool isInitialized(void);
    int predictionStep(void);
    void measurementStep(arma::vec ranging, arma::vec position_previous, arma::vec position_current, double rtoa);
    int get_steps();
    int set_process_noise(arma::mat Qm);
    int set_measurement_noise(arma::mat Rm);
    bool batchSolve();
};
  }
}

#endif //EXTENDEDKALMANFILTER_H