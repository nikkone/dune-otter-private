#include "KalmanFilter.hpp"
namespace SourceEstimators
{
  namespace FishTagXKF
  {
    void KalmanFilter::initialize(uint8_t dimension_measurement, double (&initPosition)[3]) {
      active = 0;
      nx = 3;
      nu = 0;
      ny = dimension_measurement;
      nyk = ny; //TODO: Do something with this
      nd = 3; //TODO: remove this hardcoding

      A.resizeAndFill(nx, nx,0.0);
      C.resizeAndFill(ny, nx,0.0);
      Ck.resizeAndFill(nyk, nx,0.0);
      PHat.resizeAndFill(nx, nx,0.0);
      Q.resizeAndFill(nx, nx,0.0);
      //R.resizeAndFill(ny, ny,0.0);
      Rk.resizeAndFill(nyk, nyk,0.0);
      xHat.resizeAndFill(nx,1,0.0);
      innov.resizeAndFill(nyk,1,0.0);
      yk.resizeAndFill(nyk,1,0.0);
      ykest.resizeAndFill(nyk,1,0.0);

      if(nu != 0)
      {
          B.resizeAndFill(nx, nu, 0.0);
      }
      if(nd != 0)
      {
          D.resizeAndFill(nx, nd, 0.0);
      }
      // Set the state transition matrix (does not change)
      double A_v[] = {1, 0, 0,
                      0, 1, 0,
                      0, 0, 1};

      double Q_v[] = {qq_cov, 0, 0,
                      0, qq_cov, 0,
                      0, 0, qq_cov/10};

      // // Initial position of Fish tag in [m]
      double x0_v[] = {initPosition[0], initPosition[1], initPosition[2]};
      double P0_v[] = {100, 0, 0,
                          0, 100, 0,
                          0, 0, 10};

      double D_v[] = {dt*1, 0, 0,
                      0, dt*1, 0,
                      0, 0, dt*1};

      A.fill(nx, nx, A_v);
      PHat.fill(nx, nx, P0_v);
      Q.fill(nx, nx, Q_v);
      xHat.fill(nx, 1, x0_v);
      D.fill(nx, nd, D_v);
      // // B.fill(nx, nu, B_v);

      }


      void
      KalmanFilter::update()
      {
        Matrix K = PHat*transpose(Ck)*inverse(Ck*PHat*transpose(Ck) + Rk);
        Matrix I(nx); // create identity matrix

        innov = yk - ykest;
        xHat = xHat + K*innov;
        PHat = (I - K*Ck)*PHat;
      } // End of xkf_update_stage3 function

      void KalmanFilter::predict() {
        xHat = A*xHat;
        PHat = A*PHat*transpose(A) + D*Q*transpose(D);
      }
    }          
  }