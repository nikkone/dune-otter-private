#ifndef DUNE_SOURCEESTIMATORS_XKF_LTVKF
#define DUNE_SOURCEESTIMATORS_XKF_LTVKF
#include <DUNE/DUNE.hpp>
#include "KalmanFilter.hpp"
namespace SourceEstimators
{
  namespace FishTagXKF
  {
    class LTVKF: public FishTagXKF::KalmanFilter{
        public:
        LTVKF();
        void constructCandR(DUNE::Math::Matrix receiverPositions, DUNE::Math::Matrix RDOA, double tagDepth, bool validCombination[3], double m_dr);
        
    };
  }
}
#endif //DUNE_SOURCEESTIMATORS_XKF_LTVKF