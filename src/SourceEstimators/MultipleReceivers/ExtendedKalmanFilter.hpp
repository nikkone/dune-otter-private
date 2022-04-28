#include <DUNE/DUNE.hpp>

namespace SourceEstimators
{
  //! Insert short task description here.
  //!
  //! Insert explanation on task behaviour here.
  //! @author Nikolai Lauvås
  namespace MultipleReceivers
  {
    using DUNE_NAMESPACES;

    struct ReceiverData
    {
      //bool valid;
      bool isMeasure; // Signifies if there if the data is new

      //int serial_no;
      //long unix_timestamp;
      //long milisecond_timestamp;
      double unix_milisecond_time;
      //std::string code_type;
      //int tag_id;
      int sensor_data; // assumed depth
      //int snr;
      //int memory_number;

      // Used for position of receiver
      double lat;
      double lon;
      // Used by kalman filter to calculate displacement from receiver
      double N;
      double E;
      double D;

      //int fix;
      //float hdop;
      //int satellites;
      //uint8_t checksum;
      //std::string msg;
    };
    class EKFilter
    {
      public:

        // Originally arguments in task
        //double m_d1, m_d2, m_d3;
        //ReceiverData m_r1, m_r2, m_r3;
        double speed_of_sound_in_water;
        //double m_norm_innov;
        double position[3];
        int max_time_shift_ms;

        EKFilter(uint8_t states, uint8_t outputs,uint8_t outputs_k, uint8_t inputs, uint8_t noise_inputs, double timestep);

        void constructCmatrix(ReceiverData* r1, ReceiverData* r2, ReceiverData* r3);
        void predict();
        void update();
        bool timeShiftCorrect(ReceiverData* r1, ReceiverData* r2);
        bool isInitialized(void) {
          return initialized;
        }
        bool isActive(void) {
          return active;
        }
        bool setActive(bool activate) {
          active=activate;
          return active;
        }
        Matrix getEstimate(void) {
          return xHat;
        }
      //private:
        void initializeEKF(void);
        //void constructRmatrix(void);
        Matrix xHat;    // State estimates
        Matrix PHat;    // State Covariance estimates
        Matrix A;       // State transition matrix
        Matrix B;       // Input Matrix
        Matrix C;       // Observation Matrix
        Matrix Ck;      // Time varying observation matrix
        Matrix D;       // unknown input matrix
        Matrix Q;       // Process noise covariance matrix
        //Matrix R;       // Measurement noise covariance matrix
        Matrix Rk;      // Time varying Measurement noise covariance matrix
        Matrix innov;   // Time varying innovation vector
        Matrix yk;      // Measurement vector
        Matrix ykest;   // Estimate of the measurements        
        uint8_t nx;     // Num states
        uint8_t ny;     // Num outputs
        uint8_t nyk;    // Num time varying outputs
        uint8_t nu;     // Num inputs
        uint8_t nd;     // Num unknown inputs - noise
        double dt;      // Timestep

        double qq_cov;  // Process noise Covariance - position
        double rr_cov;  // Measurement noise Covariance - range
        double rz_cov;  // Measurement noise Covariance - depth

        bool active;
        bool initialized;
    };
  }
}