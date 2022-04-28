#include "LeastSquares.hpp"
namespace SourceEstimators
{
  namespace FishTagXKF
  {
    LeastSquares::LeastSquares() {
        xHat = Matrix(3,1);
    }

    double
    LeastSquares::resolveRAmbiguity(double R1, double R2)
    {
        double R_temp;
        double max_range = 700; // TODO: Move this elsewhere

        if((R1 > 0.0) && (R1 < max_range))
        {
            if((R2 > 0.0) && (R2 < max_range))
            {
            //printf("FishTagXKF2: resolveRAmbiguity: Cannot resolve ambiguity [Valid] %f %f", R1, R2);
            R_temp = R1; // choose one of them TODO: some trick here will help
            }
            else
            {
            //debug(DTR("FishTagXKF2: resolveRAmbiguity: R1 valid! %f %f"), R1, R2);
            R_temp = R1;
            }
        }
        else
        {
            if((R2 > 0.0) && (R2 < max_range))
            {
            //debug(DTR("FishTagXKF2: resolveRAmbiguity: R2 valid! %f %f"), R1, R2);
            R_temp = R2;
            }
            else
            {
            //debug(DTR("FishTagXKF2: resolveRAmbiguity: Cannot resolve ambiguity [Invalid] %f %f"), R1, R2);
            R_temp = 0;
            }
        }
        return R_temp;
    } // End resolveRAmbiguity function
    
    uint8_t
    LeastSquares::update(DUNE::Math::Matrix receiverPositions, DUNE::Math::Matrix RDOA, double tagDepth)
    {
        uint8_t isSuccess = 0;

        // >> Construct Least squares and measurement matrices
        double d1 = -RDOA(2,0);
        double d2 = RDOA(1,0);

        double Czq_v[] = {-(receiverPositions(0,0) - receiverPositions(0,2)), -(receiverPositions(1,0) - receiverPositions(1,2)), 0.0,
                            -(receiverPositions(0,1) - receiverPositions(0,2)), -(receiverPositions(1,1) - receiverPositions(1,2)), 0.0,
                                            0,                0, 0.5};


        Matrix Czq(Czq_v, 3, 3);
        double l_v[] = {d1,d2,0};
        Matrix l(l_v, 3, 1);

        double z_v[]  = {(d1*d1 - receiverPositions(0,0)*receiverPositions(0,0) - receiverPositions(1,0)*receiverPositions(1,0) + receiverPositions(0,2)*receiverPositions(0,2) + (receiverPositions(1,2))*(receiverPositions(1,2))),
                         (d2*d2 - receiverPositions(0,1)*receiverPositions(0,1) - receiverPositions(1,1)*receiverPositions(1,1) + receiverPositions(0,2)*receiverPositions(0,2) + (receiverPositions(1,2))*(receiverPositions(1,2))),
                         tagDepth};

        Matrix z(z_v, 3, 1);

        Matrix invCzq = inverse(transpose(Czq)*Czq)*transpose(Czq);
        Matrix c = invCzq*l;
        Matrix w = 0.5*invCzq*z;

        // temp variables
        double ctc = c(0,0)*c(0,0) + c(1,0)*c(1,0) + c(2,0)*c(2,0);
        double ptc = receiverPositions(0,2)*c(0,0) +  receiverPositions(1,2)*c(1,0);
        double wtc = w(0,0)*c(0,0) + w(1,0)*c(1,0) + w(2,0)*c(2,0);
        double ptw = receiverPositions(0,2)*w(0,0) +  receiverPositions(1,2)*w(1,0);
        double wtw = w(0,0)*w(0,0) + w(1,0)*w(1,0) + w(2,0)*w(2,0);
        double ptp = (receiverPositions(0,2))*(receiverPositions(0,2)) + (receiverPositions(1,2))*(receiverPositions(1,2));

        // quadratic equation coefficients for computation of d3 = m_dr, p = p_r
        double aa = 1 - ctc; // 1 - c'c
        double bb = 2*(ptc - wtc); // 2(p'c - w'c)
        double cc = 2*ptw - wtw - ptp; //2p'w - w'w - p'p
        // Compute solution for quadratic term
        double R1, R2;
        double max_range = 700;
        Matrix fp_ls(3,1);
        if ((ctc == 1) || ((bb*bb - 4*aa*cc) <= 0.0))
        {
            if(ctc == 1) {
                R1 = -cc/bb;
            } else {
                R1 = -bb/(2*aa);
            }

            if((R1 > 0.0) && (R1 <= max_range)) {
                // Unique solution
                fp_ls = (R1*c + w);
                m_dr = R1;
                isSuccess = 1;
            }
            // Invalid solution
        } else {
            double s = sqrt(bb*bb - 4*aa*cc);
            R1 = (-bb + s)/(2*aa);
            R2 = (-bb - s)/(2*aa);
            //m_rlog.precision(15);
            //m_rlog << Clock::getSinceEpochMsec() << "," << R1 << "," << R2 << std::endl;
            R1 = resolveRAmbiguity(R1, R2);
            // Compute two candidate solutions - required to resolve ambiguity and select a particular R
            if(R1 == 0) {
                isSuccess = 0;
            } else {
                fp_ls = (R1*c + w);
                m_dr = R1;
                isSuccess = 1;
            }
        }

        // Compute position of the source for logging
        if(isSuccess)
        {
            xHat = fp_ls;
        }
        return isSuccess;
    } // End of LeastSquaresEstimate function

  }
}