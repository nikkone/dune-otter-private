// This code is written by Sverre Søbstad Løvskar as a part of his masters thesis named: 
// "Positioning of periodic acoustic emitters using an omnidirectional hydrophone on an AUV platform"


#include <iostream>
#include <armadillo>
#include <cmath>
#include "ParticleFilter.h"
namespace SourceEstimators
{
  namespace SingleReceiverDepth
  {
    using namespace std;
    using namespace arma;


    ParticleFilter::ParticleFilter(void) {
        k = 0;
    }
    ParticleFilter::ParticleFilter(int particle_width_inn, double max_range_inn, double sigma_rd_inn, double sigma_r_inn) : width(particle_width_inn), max_range(max_range_inn) {
        set_sigma(sigma_rd_inn, sigma_r_inn);
        n_particles = particle_width_inn*particle_width_inn;
        k = 0;    
    }

    int ParticleFilter::initParticles(vec center, int n)     {
        width = n;
        n_particles = n*n;   
        particles.zeros(3,n_particles);
        
        _center = center;

        int i = 0;
        double x_inc = (2*max_range)/(n-1);
        double y_inc = (2*max_range)/(n-1);
        double x_coord = center(0) - max_range;
        
        // Initialize particles in grid
        for(int x_iter=0; x_iter < n; x_iter++) {
            double y_coord = center(1) - max_range;
            for(int y_iter=0; y_iter < n; y_iter++) {
                vec pos = {x_coord,y_coord,-4};
                particles.col(i) = pos;
                
                i++;
                y_coord += y_inc;
            }
            x_coord += x_inc;
        }
        initialized = true;
        // std::cout << "init = [";
        // for(int i = 0; i < n_particles;i++) {
        //     vec p = particles.col(i) - _center;
        //     std::cout << (int)p(0) << "," << (int)p(1) << "," << (int)p(2) << "; ";
        // }
        // std::cout << "];" << std::endl;
        return n;
    }

    int ParticleFilter::update(vec ranging, vec pos_prev, vec pos_curr) {
        // Step 1: Initiate a set of N particles spread across some space.
        if (!initialized) {
            initParticles(pos_prev, width);
        } 
        
        std::vector<double> P(n_particles);
        std::vector<double> cumsum(n_particles);
        
        {
            double rangeDiff = ranging(0);
            double range = ranging(1);
            
            double y_r;  // h_k_1
            double y_rd; // h_k_2
            double sum_P_r = 0;
            double sum_P_rd = 0; 
            double max_P_r = 0;
            double max_P_rd = 0;
            std::vector<double> P_r(n_particles); //r=range Taken from SNR
            std::vector<double> P_rd(n_particles); //rd=range Taken from TDOA

            // Step 2: Calculate particle weights from the probability of each particle given some measurement.
            for(int i = 0; i<n_particles; i++) {
                y_r = norm(particles.col(i) - pos_curr,2); // ||X_p_i-X_rx1|| eq (2.26) in master
                P_r[i] = dist_scale_r * exp(-(pow(range - y_r,2.0)/exp_denominator_r)); // Eq (2.28) in master
                sum_P_r += P_r[i];  
                if (P_r[i] > max_P_r) max_P_r = P_r[i];

                y_rd = norm(particles.col(i) - pos_curr,2) - norm(particles.col(i) - pos_prev,2);  // ||X_p_i-X_rx1|| - ||X_p_i-X_rx0|| eq (2.25) in masters
                P_rd[i] = dist_scale_rd * exp(-pow(rangeDiff - y_rd,2.0)/exp_denominator_rd); // Eq (2.27) in master
                sum_P_rd += P_rd[i];
                if (P_rd[i] > max_P_rd) max_P_rd = P_rd[i];
                
            }
            // Step 3: Normalize weights
            if(max_P_r > 0) {
                for (int i = 0; i < n_particles; i++) {
                    P_r[i] = P_r[i]/sum_P_r;
                    //cout << "P_r[" << 100+i << "]: " << P_r[i] << endl;
                }
            }
            if(max_P_rd > 0) {
                for (int i = 0; i < n_particles; i++) {
                    P_rd[i] = P_rd[i]/sum_P_rd;
                    //cout << "P_rd[" << 100+i << "]: " << P_rd[i] << endl;
                }
            }
            
            for (int i = 0; i < n_particles; i++) {
                P[i] = (P_r[i] + P_rd[i])/2;
                cumsum[i] = P[i];
            }
        }
        

        // Step 4: Re-sample particles based on their weights
        for (int i = 1; i < n_particles; i++) {
            cumsum[i] += cumsum[i-1];
            //cout << "Cumsum[" << i << "]: " << cumsum[i] << endl;
        }
        
        double ranval;
        for (int i = 0; i < n_particles; i++) {
            ranval = ((double) std::rand() / RAND_MAX); // u_i. Use of rand() not really random
            int index = 0;

            // Reassign particles (last eq on page 28)      
            for (int j = 0; j < n_particles; j++) {
                if (cumsum[j] > ranval) {
                    index = j;
                    break;
                }
            }
            
            particles.col(i) = particles.col(index) + randn<vec>(3);
            
        }
    /*
        std::ofstream logOutStream;
        logOutStream.open(("log/particles.log"), std::fstream::app);
        if (logOutStream.good()) {
            logOutStream << "s"<<k<<" = [";
            for(int i = 0; i < n_particles;i++) {
                vec p = particles.col(i) - _center;
                logOutStream << (int)p(0) << "," << (int)p(1) << "," << (int)p(2) << "; ";
            }
            logOutStream << "];" << std::endl;
            logOutStream.close();
        }   
    */  

        // Step 5: Extract estimate from the posteriori particles (e.g. mean of particle).
        x = mean(particles,1); //mean of rows dim=1
        return ++k;
    }
    int ParticleFilter::get_steps() {
        return k;
    }
    vec ParticleFilter::get_sigma() {
        vec sigma = {sigma_rd,sigma_r};
        return sigma;
    }
    void ParticleFilter::set_sigma(double sigma_rd_inn, double sigma_r_inn) {
        // Sigma values determine how sharp distribution is. Low values for precise (low variance) and high values for imprecise (high variance) measurements
        sigma_r =  sigma_r_inn; //Corresponds to sigma_rd^2 in master
        sigma_rd = sigma_rd_inn; //Corresponds to sigma_rd^2 in master
        dist_scale_r = 1/sqrt(2.0*M_PI*sigma_r);
        dist_scale_rd = 1/sqrt(2.0*M_PI*sigma_rd);
        exp_denominator_r = 2*sigma_r;
        exp_denominator_rd = 2*sigma_rd;
    }
  }
}