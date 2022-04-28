
#ifndef DUNE_SOURCEESTIMATORS_XKF_LEASTSQUARES
#define DUNE_SOURCEESTIMATORS_XKF_LEASTSQUARES
#include <DUNE/DUNE.hpp>
namespace SourceEstimators
{
  namespace FishTagXKF
  {
    class LeastSquares{
      public:
        //! Reference coordinate system
        DUNE::Math::Matrix xHat;
        double m_dr;

        LeastSquares();
        void initialize();
        uint8_t update(DUNE::Math::Matrix receiverPositions, DUNE::Math::Matrix RDOA, double tagDepth);

      private:
        double resolveRAmbiguity(double R1, double R2);
    };
  }
}
#endif //DUNE_SOURCEESTIMATORS_XKF_LEASTSQUARES