#ifndef ALGEBRAICSOLVER_DEPTH_H
#define ALGEBRAICSOLVER_DEPTH_H

#include <armadillo>
namespace SourceEstimators
{
  namespace SingleReceiverDepth
  {
    class AlgebraicSolver {
        arma::mat posi;
        arma::vec ToA;
        double depth;
        uint8_t measurements;
        uint8_t neededMeasurements;
    public:
        AlgebraicSolver();
        arma::vec x; // State Estimation

        bool addMeasurement(arma::vec z, arma::vec position_previous, arma::vec position_current);
        bool solve();
    };
  }
}

#endif //ALGEBRAICSOLVER_H