#include "EKF.hpp"
namespace SourceEstimators
{
  namespace FishTagXKF
  {
    EKF::EKF() {
    }

void EKF::constructCandR(DUNE::Math::Matrix receiverPositions, DUNE::Math::Matrix RDOA, double tagDepth, bool validCombination[3], Matrix xHatInn)
    {
      //>> Get the state estimate from previous stage
      Matrix x_est = xHatInn;

      //>> Compute C and R matrices
      double qN = x_est.element(0,0);
      double qE = x_est.element(1,0);
      double qD = x_est.element(2,0);

      double d1 = std::sqrt(std::pow((receiverPositions(0,0) - qN),2) + std::pow((receiverPositions(1,0) - qE),2) + std::pow((receiverPositions(2,0) - qD),2));
      double d2 = std::sqrt(std::pow((receiverPositions(0,1) - qN),2) + std::pow((receiverPositions(1,1) - qE),2) + std::pow((receiverPositions(2,1) - qD),2));
      double d3 = std::sqrt(std::pow((receiverPositions(0,2) - qN),2) + std::pow((receiverPositions(1,2) - qE),2) + std::pow((receiverPositions(2,2) - qD),2));

      double c11 = (receiverPositions(0,1) - qN)/(d2) - (receiverPositions(0,0) - qN)/(d1);
      double c12 = (receiverPositions(1,1) - qE)/(d2) - (receiverPositions(1,0) - qE)/(d1);
      double c13 = (receiverPositions(2,1) - qD)/(d2) - (receiverPositions(2,0) - qD)/(d1);

      double c21 = (receiverPositions(0,2) - qN)/(d3) - (receiverPositions(0,1) - qN)/(d2);
      double c22 = (receiverPositions(1,2) - qE)/(d3) - (receiverPositions(1,1) - qE)/(d2);
      double c23 = (receiverPositions(2,2) - qD)/(d3) - (receiverPositions(2,1) - qD)/(d2);

      double c31 = (receiverPositions(0,0) - qN)/(d1) - (receiverPositions(0,2) - qN)/(d3);
      double c32 = (receiverPositions(1,0) - qE)/(d1) - (receiverPositions(1,2) - qE)/(d3);
      double c33 = (receiverPositions(2,0) - qD)/(d1) - (receiverPositions(2,2) - qD)/(d3);

      double c41 = 0;
      double c42 = 0;
      double c43 = 1;

      double C_v[] = {c11, c12, c13,
                      c21, c22, c23,
                      -c31, -c32, -c33,
                      c41, c42, c43};

      // Stage 3 Covariance matrix - new with cross Covariance terms
      /*double R_v[] = {2*rr_cov, rr_cov, rr_cov, 0,
                      rr_cov, 2*rr_cov, rr_cov, 0,
                      rr_cov, rr_cov, 2*rr_cov, 0,
                      0, 0, 0, rz_cov};*/

      C.fill(ny, nx, C_v);
      //R.fill(ny, ny, R_v);

      // Compute number of measurements. Plus 1 is for depth measurement.
      // This can be done this function is called only when atleast one measurement is available.

      uint8_t m12 = validCombination[0];
      uint8_t m23 = validCombination[1];
      uint8_t m31 = validCombination[2];

      nyk = m12 + m23 + m31 + 1;
      //debug(DTR("Number of measurements: %d"), nyk);

      // Resize matrices
      Ck.resizeAndFill(nyk, nx, 0.0);
      Rk.resizeAndFill(nyk, nyk, rr_cov);
      yk.resizeAndFill(nyk,1,0.0);
      ykest.resizeAndFill(nyk,1,0.0);
      innov.resizeAndFill(nyk, 1, 0.0);

      // Compute error from stage 3 and stage 2 estimate
      Matrix error_stage3 = C*(xHat - xHatInn);
      //debug(DTR("resizeAndFill successful"));
      uint8_t temp_ind = 0;
      // Construct Ck, Rk and yk
      if(m12 == 1)
      {
        double rangeDiff1 = RDOA(0,0);//speed_of_sound_in_water*(r1->unix_milisecond_time - r2->unix_milisecond_time)/1000;
        double rangeEst1 = d1 - d2;
        //debug(DTR("R12: Range difference computed %f"), rangeDiff1);
        //debug(DTR("R12: Extracted first row from C"));
        Ck.put(temp_ind, 0, C.row(0));
        //debug(DTR("R12: Put first row into Ck"));
        Rk(temp_ind, temp_ind) = 2*rr_cov;
        //debug(DTR("R12: Update Rk"));
        yk(temp_ind, 0) = rangeDiff1;
        ykest(temp_ind, 0) = rangeEst1 + error_stage3(0,0);
        temp_ind++;
        //debug(DTR("Measurement update from R12"));
      }

      if(m23 == 1)
      {
        double rangeDiff2 = RDOA(1,0);//speed_of_sound_in_water*(r2->unix_milisecond_time - r3->unix_milisecond_time)/1000;
        double rangeEst2 = d2 - d3;
        //debug(DTR("R23: Range difference computed %f"), rangeDiff2);
        Ck.put(temp_ind, 0, C.row(1));
        Rk(temp_ind, temp_ind) = 2*rr_cov;
        yk(temp_ind, 0) = rangeDiff2;
        ykest(temp_ind, 0) = rangeEst2 + error_stage3(1,0);
        temp_ind++;
        //debug(DTR("Measurement update from R23"));
      }

      if(m31 == 1)
      {
        double rangeDiff3 = RDOA(2,0);//speed_of_sound_in_water*(r3->unix_milisecond_time - r1->unix_milisecond_time)/1000;
        double rangeEst3 = d3 - d1;
        //debug(DTR("R31: Range difference computed %f"), rangeDiff3);
        Ck.put(temp_ind, 0, C.row(2));
        Rk(temp_ind, temp_ind) = 2*rr_cov;
        yk(temp_ind, 0) = -rangeDiff3;
        ykest(temp_ind, 0) = -rangeEst3 + error_stage3(2,0);
        temp_ind++;
        //debug(DTR("Measurement update from R31"));
      }

      Matrix R_temp = 0.0*Rk;
      Rk.put(temp_ind, 0, R_temp.row(temp_ind));
      Rk.put(0, temp_ind, R_temp.column(temp_ind));

      Ck.put(temp_ind, 0, C.row(3));
      Rk(temp_ind, temp_ind) = rz_cov;
      yk(temp_ind, 0) = tagDepth;
      ykest(temp_ind, 0) = qD + error_stage3(3,0);
      //debug(DTR("Update of C and R matrices successful"));
    } // End xkf_constructCandR_stage 3()
  }
}