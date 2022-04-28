// This code is written by Sverre Søbstad Løvskar as a part of his masters thesis named: 
// "Positioning of periodic acoustic emitters using an omnidirectional hydrophone on an AUV platform"
// The code has been modified by Nikolai Lauvås


#ifndef EXTENDEDKALMANFILTER_DEPTH_H
#define EXTENDEDKALMANFILTER_DEPTH_H

#include <iostream>
#include <armadillo>
namespace SourceEstimators
{
  namespace SingleReceiverDepth
  {
    class ExtendedKalmanFilter {
        arma::mat F; // System / Plant
        arma::mat Q; // Process Noise
        arma::mat R; // Measurement Noise
        int R_rows; // Rows in R
        unsigned long k; // Steps
        bool initialized;
    public:
        ExtendedKalmanFilter();
        ExtendedKalmanFilter(arma::vec X, arma::mat P, arma::mat Q, arma::mat R);
        arma::vec x; // State Estimation
        arma::mat P; // Error Coef Matrix
        arma::mat K; // Kalman Gain

        int initialize(arma::vec pos);
        bool isInitialized(void);
        int predictionStep(void);
        void measurementStep(arma::vec z, arma::vec position_previous, arma::vec position_current);
        int get_steps();
        int set_process_noise(arma::mat Qm);
        int set_measurement_noise(arma::mat Rm);
        bool batchSolve();
    };
  }
}

#endif //EXTENDEDKALMANFILTER_H