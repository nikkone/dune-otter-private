#ifndef DUNE_SOURCEESTIMATORS_XKF_KF
#define DUNE_SOURCEESTIMATORS_XKF_KF

#include <DUNE/DUNE.hpp>
using DUNE_NAMESPACES;
  namespace SourceEstimators
  {
      namespace FishTagXKF
      {
        class KalmanFilter
        {
          public:
          bool active;

          double dt; // Timestep

          uint8_t nx;     // Num states
          uint8_t ny;     // Num outputs
          uint8_t nyk;    // Num time varying outputs
          uint8_t nu;     // Num inputs
          uint8_t nd;     // Num unknown inputs - noise
          double qq_cov;  // Process noise Covariance - position
          double rr_cov;  // Measurement noise Covariance - range
          double rz_cov;  // Measurement noise Covariance - depth
          
          DUNE::Math::Matrix xHat;    // State estimates
          DUNE::Math::Matrix PHat;    // State Covariance estimates
          DUNE::Math::Matrix A;       // State transition matrix
          DUNE::Math::Matrix B;       // Input Matrix
          DUNE::Math::Matrix C;       // Observation Matrix
          DUNE::Math::Matrix Ck;      // Time varying observation matrix
          DUNE::Math::Matrix D;       // unknown input matrix
          DUNE::Math::Matrix Q;       // Process noise covariance matrix
          DUNE::Math::Matrix R;       // Measurement noise covariance matrix
          DUNE::Math::Matrix Rk;      // Time varying Measurement noise covariance matrix
          DUNE::Math::Matrix innov;   // Time varying innovation vector
          DUNE::Math::Matrix yk;      // Measurement vector
          DUNE::Math::Matrix ykest;   // Estimate of the measurements

          void initialize(uint8_t dimension_measurement, double (&initPosition)[3]);
          void predict();
          void update();
        };
      }
  }
#endif //DUNE_SOURCEESTIMATORS_XKF_KF_DATASTRUCTURES