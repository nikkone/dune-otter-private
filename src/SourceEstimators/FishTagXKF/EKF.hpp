#ifndef DUNE_SOURCEESTIMATORS_XKF_EKF
#define DUNE_SOURCEESTIMATORS_XKF_EKF

#include <DUNE/DUNE.hpp>
#include "KalmanFilter.hpp"
namespace SourceEstimators
{
  namespace FishTagXKF
  {
    class EKF: public FishTagXKF::KalmanFilter{
        public:
        EKF();
        void constructCandR(DUNE::Math::Matrix receiverPositions, DUNE::Math::Matrix RDOA, double tagDepth, bool validCombination[3], Matrix xHatInn);

    };
  }
}
#endif //DUNE_SOURCEESTIMATORS_XKF_EKF
