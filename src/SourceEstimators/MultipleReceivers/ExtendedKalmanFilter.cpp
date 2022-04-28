#include "ExtendedKalmanFilter.hpp"
namespace SourceEstimators
{
  namespace MultipleReceivers
  {
    EKFilter::EKFilter(uint8_t states, uint8_t outputs,uint8_t outputs_k, uint8_t inputs, uint8_t noise_inputs, double timestep) :nx(states), ny(outputs),nyk(outputs_k),nu(inputs),nd(noise_inputs), dt(timestep) {
      qq_cov = 0.1;//FishPos Cov
      rr_cov = 0.1;//ToA Cov
      rz_cov = 0.2;//Depth Cov
      max_time_shift_ms = 500;
      speed_of_sound_in_water = 1485;
      position[0] = 100;
      position[1] = 100;
      position[2] = 7;
      
      
//      position = {100, 100, 7};
      //constructRmatrix();
      initializeEKF();
    }

//! Initialize Extended Kalman filter
void EKFilter::initializeEKF(void) {
    // Set the size of Kalman filter matrices
       A=Matrix(nx, nx,0.0);
       C=Matrix(ny, nx,0.0);
      Ck=Matrix(nyk, nx,0.0);
    PHat=Matrix(nx, nx,0.0);
       Q=Matrix(nx, nx,0.0);
       //R=Matrix(ny, ny,0.0); // Not used, only Rk used
      Rk=Matrix(nyk, nyk,0.0);
    xHat=Matrix(nx,1,0.0);
   innov=Matrix(nyk,1,0.0);
      yk=Matrix(nyk,1,0.0);
   ykest=Matrix(nyk,1,0.0);

    if(nu != 0)
    {
        B.resizeAndFill(nx, nu, 0.0);
    }
    if(nd != 0)
    {
        D.resizeAndFill(nx, nd, 0.0);
    }
    // Set the state transition matrix (does not change)
    double A_v[] = {1, 0, 0,
                    0, 1, 0,
                    0, 0, 1};

    double Q_v[] = {qq_cov, 0, 0,
                    0, qq_cov, 0,
                    0, 0, qq_cov/10.0};

    double x0_v[] = {position[0], position[1], position[2]};
    double P0_v[] = {1, 0, 0,
                        0, 1, 0,
                        0, 0, 1};

    //D_v: dt is timestep of random walk modell. Set equal to timestep of the filter
    double D_v[] = {dt*1, 0, 0,
                    0, dt*1, 0,
                    0, 0, dt*1};

    A.fill(nx, nx, A_v);
    
    PHat.fill(nx, nx, P0_v);
    Q.fill(nx, nx, Q_v);
    xHat.fill(nx, 1, x0_v);
    D.fill(nx, nd, D_v);
    //std::cout << A << std::endl;

    initialized = 1;
}

/*void EKFilter::constructRmatrix(void) {
    double R_v[] = {2*rr_cov, 0, 0, 0,
                    0, 2*rr_cov, 0, 0,
                    0, 0, 2*rr_cov, 0,
                    0, 0, 0, rz_cov};
    R.fill(ny, ny, R_v);
}*/

//! Builds Jacobian matrix for available measurements, and prepares Rk, Yk and Ykest
void EKFilter::constructCmatrix(ReceiverData* r1, ReceiverData* r2, ReceiverData* r3) {
    //>> Get the state estimate
    Matrix x_est = xHat;
    //std::cout << x_est << std::endl;

    //Set estimated position as tag position Q
    double qN = x_est.element(0,0);
    double qE = x_est.element(1,0);
    double qD = x_est.element(2,0);

    // Compute distance between estimated position and receiver positions providing tag detections.
    double d1 = std::sqrt(std::pow((r1->N - qN),2) + std::pow((r1->E - qE),2) + std::pow((r1->D - qD),2));
    double d2 = std::sqrt(std::pow((r2->N - qN),2) + std::pow((r2->E - qE),2) + std::pow((r2->D - qD),2));
    double d3 = std::sqrt(std::pow((r3->N - qN),2) + std::pow((r3->E - qE),2) + std::pow((r3->D - qD),2));

    /*m_d1 = d1;
    m_d2 = d2;
    m_d3 = d3;*/

    //printf(DTR("ekf_constructCandR: range d1: %f, d2: %f, d3: %f"), d1, d2, d3);
    // Calculates the Jacobian of the measurement function?

/*    
    double c11 = (r2->N - qN)/(d2) - (r1->N - qN)/(d1);
    double c12 = (r2->E - qE)/(d2) - (r1->E - qE)/(d1);
    double c13 = (r2->D - qD)/(d2) - (r1->D - qD)/(d1);

    double c21 = (r3->N - qN)/(d3) - (r2->N - qN)/(d2);
    double c22 = (r3->E - qE)/(d3) - (r2->E - qE)/(d2);
    double c23 = (r3->D - qD)/(d3) - (r2->D - qD)/(d2);

    double c31 = (r1->N - qN)/(d1) - (r3->N - qN)/(d3);
    double c32 = (r1->E - qE)/(d1) - (r3->E - qE)/(d3);
    double c33 = (r1->D - qD)/(d1) - (r3->D - qD)/(d3);
*/
    
    double c11 = (qN - r2->N)/(d2) - (qN - r1->N)/(d1);
    double c12 = (qE - r2->E)/(d2) - (qE - r1->E)/(d1);
    double c13 = (qD - r2->D)/(d2) - (qD - r1->D)/(d1);

    double c21 = (qN - r3->N)/(d3) - (qN - r2->N)/(d2);
    double c22 = (qE - r3->E)/(d3) - (qE - r2->E)/(d2);
    double c23 = (qD - r3->D)/(d3) - (qD - r2->D)/(d2);

    double c31 = (qN - r1->N)/(d1) - (qN - r3->N)/(d3);
    double c32 = (qE - r1->E)/(d1) - (qE - r3->E)/(d3);
    double c33 = (qD - r1->D)/(d1) - (qD - r3->D)/(d3);

    double c41 = 0;
    double c42 = 0;
    double c43 = 1;

    double C_v[] = {c11, c12, c13,
                    c21, c22, c23,
                    c31, c32, c33,
                    c41, c42, c43};

    C.fill(ny, nx, C_v);

    // Compute number of measurements. Plus 1 is for depth measurement.
    // This can be done this function is called only when atleast one measurement is available.

    uint8_t m12 = ((r1->isMeasure || r2->isMeasure) && (timeShiftCorrect(r1, r2)));
    uint8_t m23 = ((r2->isMeasure || r3->isMeasure) && (timeShiftCorrect(r2, r3)));
    uint8_t m31 = ((r3->isMeasure || r1->isMeasure) && (timeShiftCorrect(r3, r1)));

    nyk = m12 + m23 + m31 + 1;
    /*std::cout << "dist3 " << nyk << " " << (timeShiftCorrect(r1, r2)) << std::endl;
    if(m12) {
      std::cout << "true " << nyk << " " << (timeShiftCorrect(r1, r2)) << std::endl;
    }*/
    //printf(DTR("Number of measurements: %d"), nyk);

    // Resize matrices (Because varying amounts of measurements available.)
    Ck.resizeAndFill(nyk, nx, 0.0);
    Rk.resizeAndFill(nyk, nyk, 0.0);
    yk.resizeAndFill(nyk,1,0.0);
    ykest.resizeAndFill(nyk,1,0.0);
    innov.resizeAndFill(nyk, 1, 0.0);

    //printf(DTR("resizeAndFill successful"));
    uint8_t temp_ind = 0;
    // Construct Ck, Rk and yk
    if(m12 == 1)
    {
        double rangeDiff1 = speed_of_sound_in_water*(r2->unix_milisecond_time - r1->unix_milisecond_time)/1000;
        double rangeEst1 = d2 - d1;
        //std::cout << "dist1 " << rangeDiff1 << " " << rangeEst1 << std::endl;
        //printf(DTR("R12: Range difference computed %f"), rangeDiff1);
        Matrix temp_C1 = C.row(0);
        //printf(DTR("R12: Extracted first row from C"));
        Ck.put(temp_ind, 0, temp_C1);
        //printf(DTR("R12: Put first row into Ck"));
        Rk(temp_ind, temp_ind) = 2*rr_cov;
        //printf(DTR("R12: Update Rk"));
        yk(temp_ind, 0) = rangeDiff1;
        ykest(temp_ind, 0) = rangeEst1;
        temp_ind++;
        //printf(DTR("Measurement update from R12"));
        r1->isMeasure = 0; r2->isMeasure = 0;
    }

    if(m23 == 1)
    {
        double rangeDiff2 = speed_of_sound_in_water*(r3->unix_milisecond_time - r2->unix_milisecond_time)/1000;
        double rangeEst2 = d3 - d2;
        //printf(DTR("R23: Range difference computed %f"), rangeDiff2);
        //std::cout << "dist2 " << rangeDiff2 << " " << rangeEst2 << std::endl;
        Ck.put(temp_ind, 0, C.row(1));
        Rk(temp_ind, temp_ind) = 2*rr_cov;
        yk(temp_ind, 0) = rangeDiff2;
        ykest(temp_ind, 0) = rangeEst2;
        temp_ind++;
        //printf(DTR("Measurement update from R23"));
        r2->isMeasure = 0; r3->isMeasure = 0;
    }

    if(m31 == 1)
    {
        double rangeDiff3 = speed_of_sound_in_water*(r1->unix_milisecond_time - r3->unix_milisecond_time)/1000;
        double rangeEst3 = d1 - d3;
        //printf(DTR("R31: Range difference computed %f"), rangeDiff3);
        //std::cout << "dist3 " << rangeDiff3 << " " << rangeEst3 << std::endl;
        Ck.put(temp_ind, 0, C.row(2));
        Rk(temp_ind, temp_ind) = 2*rr_cov;
        yk(temp_ind, 0) = rangeDiff3;
        ykest(temp_ind, 0) = rangeEst3;
        temp_ind++;
        //printf(DTR("Measurement update from R31"));
        r3->isMeasure = 0; r1->isMeasure = 0;
    }

    // Takes average of depth from used tags
    double depth = m12*(r1->sensor_data + r2->sensor_data)/2 +
                    m23*(r2->sensor_data + r3->sensor_data)/2 +
                    m31*(r3->sensor_data + r1->sensor_data)/2;
    //printf(DTR("ekf_constructCandR:: Depth before %f"), depth);
    depth = depth/(m12 + m23 + m31);
    //printf(DTR("ekf_constructCandR:: Depth after %f by %d"), depth, (m12 + m23 + m31));

    Ck.put(temp_ind, 0, C.row(3));
    Rk(temp_ind, temp_ind) = rz_cov;
    yk(temp_ind, 0) = depth;
    ykest(temp_ind, 0) = qD;
    //printf(DTR("Update of C and R matrices successful"));
    } // End ekf_constructCandR()

    void
    EKFilter::predict()
    {

    //TODO: Manage situations when inputs are present
    xHat = A*xHat;
    PHat = A*PHat*transpose(A) + D*Q*transpose(D);

    //printf(DTR("Predicted state: %f %f %f"), xHat.element(0,0), xHat.element(1,0), xHat.element(2,0));
    //printf(DTR("Predict successful"));
}

void EKFilter::update() {
//std::cout << "Nikkofilter" << std::endl;
    //std::cout << "A"<< std::endl<< A << std::endl;
    //std::cout << "B"<< std::endl<< B << std::endl; 
    //std::cout << "C"<< std::endl<< C << std::endl; 
    //std::cout << "D"<< std::endl<< D << std::endl; 
    //std::cout << "Q"<< std::endl<< Q << std::endl; 
    //std::cout << "R"<< std::endl<< R << std::endl; 
    //std::cout << "Ck"<< std::endl<< Ck << std::endl; 
    //std::cout << "Rk"<< std::endl<< Rk << std::endl; 
    //std::cout << "yk"<< std::endl<< yk << std::endl; 
    Matrix K = PHat*transpose(Ck)*inverse(Ck*PHat*transpose(Ck) + Rk);
    Matrix I(nx); // create identity matrix

    innov = yk - ykest;
    xHat = xHat + K*innov;
    PHat = (I - K*Ck)*PHat;
    //m_norm_innov = innov.norm_p(2);

    // Make covariance matrix symmetric
    // PHat = 0.5*(PHat + transpose(PHat));
    //printf(DTR("Estimated state: %f %f %f"), xHat.element(0,0), xHat.element(1,0), xHat.element(2,0));
    // //printf(DTR("Innovation: %f %f %f %f"), innov(0,0), innov(1,0), innov(2,0), innov(3,0));
    //printf(DTR("Update successful"));
}

//! Check if time shift between messages from 3 different sensors is within acceptable limits defined by user
bool EKFilter::timeShiftCorrect(ReceiverData* r1, ReceiverData* r2) {
    double t12 = std::abs(r1->unix_milisecond_time - r2->unix_milisecond_time);

    //printf(DTR("Timeshifts: %f"), t12);

    if((t12 <= max_time_shift_ms))
        return true;
    else
        return false;
    } // End function timeShiftCorrect
  }
}