#include "AlgebraicSolution.hpp"

using namespace std;
using namespace arma;
namespace SourceEstimators
{
  namespace SingleReceiverDepth
  {
    AlgebraicSolver::AlgebraicSolver() {
      posi = arma::zeros<mat>(3,5);
      ToA = arma::zeros<vec>(5);
      measurements = 0;
      neededMeasurements = 5;
    }
    bool AlgebraicSolver::addMeasurement(vec z, vec position_previous, vec position_current) {
      depth = z(2);
      double rtoa = z(0);
      if (rtoa == 0) {
          rtoa = 0.00000001;
      }
      //position_current(2) += ((double) std::rand() / RAND_MAX)/100;
      if (measurements == 0) {
          ToA(measurements) = 0;
          posi.col(measurements) = position_previous;
          ToA(measurements) = 0;
          measurements++;
          posi.col(measurements) = position_current;
          ToA(measurements) = ToA(measurements-1) + rtoa;
          measurements++;
      } else {
          //std::cout << "Measurements" << ((measurements % neededMeasurements) -1) % neededMeasurements << std::endl;
          posi.col(measurements % neededMeasurements) = position_current;
          int prevMeasurement = ((measurements % neededMeasurements) -1) % neededMeasurements;
          if(prevMeasurement < 0 ) prevMeasurement=neededMeasurements-1;
          ToA(measurements % neededMeasurements) = rtoa + ToA(prevMeasurement);
          measurements++;
      }
      if (measurements > 5-1) {
          try {
              solve();
              return true;
          } catch (...) {
              //cout << "Ah, nope.. exception.." << endl;
              return false;
          }
      }
      return false;
        
    }
    bool AlgebraicSolver::solve() {
        //cout << " AlgebraicSolver Batch solve time!" << endl;
        //cout << "Posi = " << endl << posi << endl;
        //cout << "ToA = " << endl << ToA << endl;
        //double depth = 3.0;
        arma::mat M = arma::zeros<arma::mat>(3,neededMeasurements+1);
        arma::vec temp={0.0,0.0,1.0};
        arma::vec D = arma::zeros<arma::vec>(neededMeasurements+1);
        M.col(neededMeasurements) = temp;
        D(neededMeasurements) = depth;
        //cout << "Starting loop!" << endl;
        for (unsigned int m = 2; m < posi.n_cols; m++) {
        
            double ddm = ToA(m);
            double dd2 = ToA(1);
            
            M.col(m) = (2*(posi.col(m) - posi.col(0)) / ddm) - (2*(posi.col(1)-posi.col(0)) / dd2);
            //cout << "M.Col" << endl << M.col(m) << endl;
            double off1 = arma::sum(posi.col(0)%posi.col(0));
            double off2 = arma::sum(posi.col(1)%posi.col(1));
            double offm = arma::sum(posi.col(m)%posi.col(m));
            
            D(m) = ddm - dd2 + (off1-offm)/ddm - (off1-off2)/dd2;
            //cout << "Iter " << m << " ok" << endl;
        }
        //cout << "Extracting from " << 2 << " to " << posi.n_cols-1 << endl;
        M = M.cols(2,neededMeasurements);
        D = -D.rows(2,neededMeasurements);
        //cout << "M = " << endl << M << endl;
        //cout << "D = " << endl << D << endl;
        if (arma::solve(x,M.t(),D)) {
            return true;
        } 
        return false;
    }
  }
}