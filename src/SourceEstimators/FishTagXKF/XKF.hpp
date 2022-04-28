
#ifndef DUNE_SOURCEESTIMATORS_XKF
#define DUNE_SOURCEESTIMATORS_XKF

#include <DUNE/DUNE.hpp>
#include "LTVKF.hpp"
#include "EKF.hpp"
#include "LeastSquares.hpp"
namespace SourceEstimators
{
  namespace FishTagXKF
  {
    class XKF {
      public:
        LeastSquares stage1;
        LTVKF stage2;
        EKF stage3;

        XKF();

        void initialize(double (&initPosition)[3]);
        void predict();
        void update(DUNE::Math::Matrix receiverPositions, DUNE::Math::Matrix RDOA, double tagDepth, bool validCombination[3]);


        bool isInitialized(void);
        void setReferenceCoordinate(double lat, double lon, double hae);
        void setInitialPosition(double north, double east, double down);
        void setTimestep(double timestep);
        void setDiagonalCovarianceR(double rr_cov);
        void setDiagonalCovarianceQ(double qq_cov);
        void setVarianceRZ(double rz_var);
        void setActive(bool activate);
        bool isActive(void);
      private:
        bool initialized;
        bool active;
    };
  }
}
#endif //DUNE_SOURCEESTIMATORS_XKF