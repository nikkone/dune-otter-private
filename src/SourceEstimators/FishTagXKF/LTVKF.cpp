#include "LTVKF.hpp"
namespace SourceEstimators
{
  namespace FishTagXKF
  {
    LTVKF::LTVKF() {
    }

void LTVKF::constructCandR(DUNE::Math::Matrix receiverPositions, DUNE::Math::Matrix RDOA, double tagDepth, bool validCombination[3], double m_dr)
    {
        //>> Compute C and R matrices
        double Czq_v[] = {-(receiverPositions(0,0) - receiverPositions(0,2)), -(receiverPositions(1,0) - receiverPositions(1,2)), 0.0,
                            -(receiverPositions(0,1) - receiverPositions(0,2)), -(receiverPositions(1,1) - receiverPositions(1,2)), 0.0,
                                            0,                0, 0.5};

        // TODO: Change the way covariance matrix is initialized
        double R_v[] = {0.1, 0, 0,
                        0, 0.1, 0,
                        0, 0, 0.2};

        C.fill(ny, nx, Czq_v);
        R.fill(ny, ny, R_v);

        // Compute number of measurements. Plus 1 is for depth measurement.
        // This can be done this function is called only when atleast one measurement is available.

        uint8_t m23 = validCombination[1];
        uint8_t m31 = validCombination[2];

        nyk = m23 + m31 + 1;
        //debug(DTR("Number of measurements: %d"), nyk);

        // Resize matrices
        Ck.resizeAndFill(nyk, nx, 0.0);
        Rk.resizeAndFill(nyk, nyk, 0.0);
        yk.resizeAndFill(nyk,1,0.0);
        ykest.resizeAndFill(nyk,1,0.0);
        innov.resizeAndFill(nyk, 1, 0.0);

        //debug(DTR("resizeAndFill successful"));
        uint8_t temp_ind = 0;
        Matrix Yest =  C*xHat;
        double d1_temp = 0.0, d2_temp = 0.0;
        // Construct Ck, Rk and yk
        if(m31 == 1)
        {
            double d1 = -RDOA(2,0);
            double z1 = 0.5*(d1*d1 - receiverPositions(0,0)*receiverPositions(0,0) - receiverPositions(1,0)*receiverPositions(1,0) + (receiverPositions(0,2))*(receiverPositions(0,2)) + (receiverPositions(1,2))*(receiverPositions(1,2)));
            //debug(DTR("R13: Measurement computed %f"), z1);
            Ck.put(temp_ind, 0, C.row(0));
            double z1_cov = (1/2)*pow(2*rr_cov,2) + (d1 + m_dr)*(d1 + m_dr)*2*rr_cov;
            Rk(temp_ind, temp_ind) = z1_cov;
            yk(temp_ind, 0) = z1 + m_dr*d1;
            ykest(temp_ind, 0) = Yest(temp_ind, 0);
            temp_ind++;
            //debug(DTR("Measurement update from R13"));
            d1_temp = d1;
        }

        if(m23 == 1)
        {
            double d2 = RDOA(1, 0);
            double z2 = 0.5*(d2*d2 - receiverPositions(0,1)*receiverPositions(0,1) - receiverPositions(1,1)*receiverPositions(1,1) + (receiverPositions(0,2))*(receiverPositions(0,2)) + (receiverPositions(1,2))*(receiverPositions(1,2)));
            //debug(DTR("R23: Measurement computed %f"), z2);
            Ck.put(temp_ind, 0, C.row(1));
            double z2_cov = (1/2)*pow(2*rr_cov,2) + (d2 + m_dr)*(d2 + m_dr)*2*rr_cov;
            Rk(temp_ind, temp_ind) = z2_cov;
            yk(temp_ind, 0) = z2 + m_dr*d2;
            ykest(temp_ind, 0) = Yest(temp_ind, 0);
            temp_ind++;
            //debug(DTR("Measurement update from R23"));

            d2_temp = d2;
        }

        if(nyk == 3)
        {
            Rk(0,1) = 0.5*pow(rr_cov,2) + (d1_temp*d2_temp + m_dr*(d1_temp + d2_temp) + m_dr*m_dr)*rr_cov;
            Rk(1,0) = 0.5*pow(rr_cov,2) + (d1_temp*d2_temp + m_dr*(d1_temp + d2_temp) + m_dr*m_dr)*rr_cov;
        }


        Ck.put(temp_ind, 0, C.row(2));
        Rk(temp_ind, temp_ind) = rz_cov;
        yk(temp_ind, 0) = tagDepth/2; // Why divided by two?
        ykest(temp_ind, 0) = Yest(temp_ind, 0);

    } // End xkf_constructCandR_stage2()

  }
}