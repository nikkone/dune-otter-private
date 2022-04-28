

#include "XKF.hpp"
namespace SourceEstimators
{
  namespace FishTagXKF
  {
    XKF::XKF() {
      //initializeXKF();
    }

    //! Initialize Extended Kalman filter
    //! @param [in] initPosition[0] North coordinate in reference frame given in meters.
    //! @param [in] initPosition[1] East coordinate in reference frame given in meters.
    //! @param [in] initPosition[2] Height in reference frame given in meters.
    void
    XKF::initialize(double (&initPosition)[3])
    {
      stage2.initialize(3, initPosition);
      stage3.initialize(4, initPosition);
      initialized = true;
    } // End of initializeXKF function


    void
    XKF::predict()
    {
      stage2.predict();
      stage3.predict();
    }

    void XKF::update(DUNE::Math::Matrix receiverPositions, DUNE::Math::Matrix RDOA, double tagDepth, bool validCombination[3]) {

        if(validCombination[1] && validCombination[2])
        {
          bool isSuccess = stage1.update(receiverPositions, RDOA, tagDepth);

          if(isSuccess && !stage2.active)
          {
              stage3.active = 1;
              stage2.active = 1;
          }
        }
        //>> Step 2.2: Measurement Update - stage - 2
        if(stage2.active && (validCombination[1] && validCombination[2])) // conditions for stage 2 update
        {
          stage2.constructCandR(receiverPositions, RDOA, tagDepth, validCombination, stage1.m_dr);
          stage2.update();
        }
        if(stage3.active)
        {
          stage3.constructCandR(receiverPositions, RDOA, tagDepth, validCombination, stage2.xHat);
          stage3.update();
        }
    }

    bool XKF::isInitialized(void) {
      return initialized;
    }

    void XKF::setTimestep(double timestep) {
      stage2.dt=timestep;
      stage3.dt=timestep;
    }

    void XKF::setDiagonalCovarianceR(double rr_cov) {
      stage3.rr_cov = rr_cov;
      stage2.rr_cov = rr_cov;
    }

    void XKF::setDiagonalCovarianceQ(double qq_cov) {
      stage3.qq_cov = qq_cov;
      stage2.qq_cov = qq_cov;
    }

    void XKF::setVarianceRZ(double rz_var) {
      stage3.rz_cov = rz_var;
      stage2.rz_cov = rz_var;
    }

    void XKF::setActive(bool activate) {
      active=activate;
    }

    bool XKF::isActive(void) {
      return active;
    }
  }
}