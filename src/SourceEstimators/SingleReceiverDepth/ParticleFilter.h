// This code is written by Sverre Søbstad Løvskar as a part of his masters thesis named: 
// "Positioning of periodic acoustic emitters using an omnidirectional hydrophone on an AUV platform"


#ifndef PARTICLEFILTER_DEPTH_H
#define PARTICLEFILTER_DEPTH_H
namespace SourceEstimators
{
  namespace SingleReceiverDepth
  {

        #include <iostream>
        #include <armadillo>
        #include <cmath>

        class ParticleFilter {
        public:
            arma::vec x; // State Estimatate

            ParticleFilter(void);
            //! Constructor
            //! @param [in] n particles one direction
            ParticleFilter(int particle_width_inn, double max_range_inn, double sigma_rd_inn, double sigma_r_inn);

            int initParticles(arma::vec center, int n);

            int update(arma::vec RangeDiff, arma::vec Pos1, arma::vec Pos2);
            int get_steps();
            arma::vec get_sigma();
            void set_sigma(double sigma_rd_inn, double sigma_r_inn);
        private:
            unsigned long k; // Steps
            bool initialized;
            int n_particles;
            int width;
            double max_range; //in meters 
            arma::mat particles;
            arma::vec _center;
            double dist_scale_r;
            double dist_scale_rd;
            double var_r;
            double var_rd;
            double sigma_r;
            double sigma_rd;
            double exp_denominator_r;
            double exp_denominator_rd;
        };
  }
}
#endif //PARTICLEFILTER_H   
