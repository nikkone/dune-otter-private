#ifndef SVEXTENDEDKALMANFILTER_PLAIN_H
#define SVEXTENDEDKALMANFILTER_PLAIN_H

#include <iostream>
#include <armadillo>
namespace SourceEstimators
{
  //! Insert explanation on task behaviour here.
  namespace SingleReceiver
  {
class SverreExtendedKalmanFilter {
    arma::mat F; // System / Plant
    arma::mat Q; // Process Noise
    arma::mat R; // Measurement Noise
    unsigned long k; // Steps
    bool initialized;

    arma::mat posi;
    arma::vec ToA;
    int pre_init_data;
public:
    SverreExtendedKalmanFilter();
    SverreExtendedKalmanFilter(arma::vec X, arma::mat P, arma::mat Q, arma::mat R);
    arma::mat x; // State Estimation
    arma::mat P; // Error Coef Matrix

    int initialize(arma::vec pos);
    int update(arma::vec RangeDiff, arma::vec Pos1, arma::vec Pos2, double td);
    int steps();
    int set_process_noise(arma::mat Qm);
    int set_measurement_noise(arma::mat Rm);
    bool batchSolve();
};
  }
}

#endif //EXTENDEDKALMANFILTER_H