//***************************************************************************
// Copyright 2013-2022 Norwegian University of Science and Technology (NTNU)*
// Department of Engineering Cybernetics (ITK)                              *
//***************************************************************************
// This file was developed for use in DUNE: Unified Navigation Environment. *
//                                                                          *
// Commercial Licence Usage                                                 *
// Licencees holding valid commercial DUNE licences may use this file in    *
// accordance with the commercial licence agreement provided with the       *
// Software or, alternatively, in accordance with the terms contained in a  *
// written agreement between you and Universidade do Porto. For licensing   *
// terms, conditions, and further information contact lsts@fe.up.pt.        *
//                                                                          *
// European Union Public Licence - EUPL v.1.1 Usage                         *
// Alternatively, this file may be used under the terms of the EUPL,        *
// Version 1.1 only (the "Licence"), appearing in the file LICENCE.md       *
// included in the packaging of this file. You may not use this work        *
// except in compliance with the Licence. Unless required by applicable     *
// law or agreed to in writing, software distributed under the Licence is   *
// distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF     *
// ANY KIND, either express or implied. See the Licence for the specific    *
// language governing permissions and limitations at                        *
// http://ec.europa.eu/idabc/eupl.html.                                     *
//***************************************************************************
// Author: Nikolai Lauvås                                                   *
//***************************************************************************

// DUNE headers.
#include <DUNE/DUNE.hpp>
#include <Eigen/Core>

// CPP STD headers
#include <algorithm>  // std::count
#include <iterator>

// External Libraries
//#include <OpenFilterPack/AlgebraicSolution.hpp>
//#include <OpenFilterPack/ExtendedKalmanFilter.hpp>

  namespace SourceEstimators
  {
      namespace MultipleReceiverXKF2
      {
        using DUNE_NAMESPACES;

        const uint8_t c_receivers = 4; // TODO: Make everything scalable according to this number so x receivers can be used

        //! %Task arguments.
        struct Arguments
        {
          //! Reference coordinate position (degrees)
          std::vector<double> reference;
          //! Initial position of the fish tag (NED) for EKF
          std::vector<double> position;
          //! Sensor serial numbers to use
          std::vector<uint32_t> receiver_serial;
          //! Fish Position Covariance
          double qq_cov;
          //! Time of Arrival Covariance
          double rr_cov;
          //! Depth measurement Covariance
          double rz_cov;
          //! Initial value for speed of sound in water [m/s]
          double speed_of_sound_in_water;
          //! Maximum allowed time [ms] shift between receivers' messages
          double max_time_shift_ms;

// Kalman Filter
          //! Extended Kalman filter - Qm
          std::vector<double> ekf_Qm;      
          //! Extended Kalman filter - Rm
          std::vector<double> ekf_Rm; 
          //! Extended Kalman filter - P0
          std::vector<double> ekf_P0;
          //! Extended Kalman filter - x0
          std::vector<double> ekf_x0;
        };

        struct Task: public DUNE::Tasks::Periodic
        {

          //! Task arguments
          Arguments m_args;
          //! Output files for logging
          std::ofstream m_o[c_receivers];
          //! Storage for most recent fish tag detection per receiver.
          IMC::TBRFishTag m_tagDetection[c_receivers];
          //! To keep track of registrations that are used.
          bool m_newDetection[c_receivers];
          //! Value to compare against when desciding if a pair of measurements are valid
          double m_max_rdoa;
          //! Reference coordinate used to calculate NED frame
          double m_refCoord[3];
          //! Position in NED used to initialize the estimator
          double m_initPosition[3];
          //! Position estimate filter
          //OFP::ExtendedKalmanFilter<double, 3, 3, 0> m_ekf;
          //! Constructor.
          //! @param[in] name task name.
          //! @param[in] ctx context.
          Task(const std::string& name, Tasks::Context& ctx):
            DUNE::Tasks::Periodic(name, ctx)
          {
            param("Receiver Serial Numbers", m_args.receiver_serial)
            .description("Receiver Serial Numbers")
            .size(c_receivers)
            .defaultValue("632, 634, 631");

            param("Max time shift [ms]", m_args.max_time_shift_ms)
            .description("Maximum allowed time [ms] shift between receivers' messages")
            .defaultValue("500");

            param("Speed of sound in water [m/s]", m_args.speed_of_sound_in_water)
            .description("Speed of sound in water [m/s]")
            .defaultValue("1485");

            param("FishPos Cov", m_args.qq_cov)
            .description("Fish Position Covariance")
            .defaultValue("0");

            param("ToA Cov", m_args.rr_cov)
            .description("Time of Arrival Covariance")
            .defaultValue("0");

            param("Depth Cov", m_args.rz_cov)
            .description("Depth measurement Covariance")
            .defaultValue("0");

            param("Reference Coordinate", m_args.reference)
            .units(Units::Degree)
            .size(2)
            .description("Origin of the reference coordinate system");

            param("Initial Position", m_args.position)
            .size(3)
            .description("Initial Fish tag position in NED from reference coordinate frame");
// Kalman Filter Parameters
            param("x0", m_args.ekf_x0)
            .size(3)
            .description("Initial X value for the extended Kalman filter")
            .defaultValue("0.0, 0.0, 0.0");

            param("P0", m_args.ekf_P0)
            .size(9)
            .description("Initial P matrix value for the extended Kalman filter, first row")
            .defaultValue("66458, -26820, 0, -26820, 12116, 0, 0, 0, 0");

            param("Rm", m_args.ekf_Rm)
            .description("Initial X value for the extended Kalman filter")
            .defaultValue("1.5*3.8706, 0, 0, 2.1638e6");

            param("Qm", m_args.ekf_Qm)
            .size(9)
            .description("Initial X value for the extended Kalman filter")
            .defaultValue("0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00001");
            // Setup processing of IMC messages
            // Initialize messages.
            bind<IMC::TBRFishTag>(this);
          }

          //! Update internal state with new parameter values.
          void
          onUpdateParameters(void)
          {
            m_refCoord[0] = Math::Angles::radians(m_args.reference[0]);
            m_refCoord[1] = Math::Angles::radians(m_args.reference[1]);
            m_refCoord[2] = 0.0;

            m_max_rdoa = m_args.max_time_shift_ms/1000*m_args.speed_of_sound_in_water;
          }

          //! Reserve entity identifiers.
          void
          onEntityReservation(void)
          {
          }

          //! Resolve entity names.
          void
          onEntityResolution(void)
          {
          }

          //! Acquire resources.
          void
          onResourceAcquisition(void)
          {
          }

          //! Initialize resources.
          void
          onResourceInitialization(void)
          {
        // m_ekf
          /*m_ekf.A << 1.0, 0.0, 0.0,
                    0.0, 1.0, 0.0,
                    0.0, 0.0, 1.0;
          m_ekf.Q = Eigen::Map<Eigen::Matrix<double, 3, 3> >(m_args.ekf_Qm.data());
          m_ekf.R = Eigen::Map<Eigen::Matrix<double, 3, 3> >(m_args.ekf_Rm.data());
          m_ekf.PHat = Eigen::Map<Eigen::Matrix<double, 3, 3> >(m_args.ekf_P0.data());
          m_ekf.xHat = Eigen::Map<Eigen::Matrix<double, 3, 1> >(m_args.ekf_x0.data());
          //! Set timer for periodic part of filter
          m_ekf.dt = getFrequency();*/

          /*m_ekf.h = [this](Eigen::Matrix<double, 3, 1> x, std::vector<std::tuple<double, double, double>> NED) {
            // Find euclidean norm (p-norm, p=2) between measurements and estimated tag position
            Eigen::Matrix<double, 3, 1> distance1 = x;//-this->pos_previous;//z.block(0,0,3,1); // X_e-X_rx0
            Eigen::Matrix<double, 3, 1> distance2 = x;//-this->pos_current;//z.block(3,0,3,1);  // X_e-X_rx1
            double r1 = distance1.norm();//  ||X_e-X_rx0||
            double r2 = distance2.norm();// ||X_e-X_rx1||
            
            // Calculate estimated measurements
            Eigen::Matrix<double, 3, 1> ykest;
            ykest(0) = r2 - r1; // h is eq (2.16) in masters
            ykest(1) = r2; // Eq (2.19) in masters
            ykest(2) = x(2); // Depth estimate
            return ykest;
          };*/
/*
          m_ekf.calculateJacobian = [this](Eigen::Matrix<double, 3, 1> x) {
            // Find euclidean norm (p-norm, p=2) between measurements and estimated tag position
            Eigen::Matrix<double, 3, 1> distance1 = x-this->pos_previous;//z.block(0,0,3,1); // X_e-X_rx0
            Eigen::Matrix<double, 3, 1> distance2 = x-this->pos_current;//z.block(3,0,3,1);  // X_e-X_rx1
            double r1 = distance1.norm();//  ||X_e-X_rx0||
            double r2 = distance2.norm();// ||X_e-X_rx1||

            // Calculate Jacobian with RDOA, SNR and Depth
            Eigen::Matrix<double, 3, 3> C;       // Observation matrix
            C.row(0) = (distance2/r2) - (distance1/r1); // Eq (2.18)
            C.row(1) = (distance2/r2); // Exends the Jacobian with eq (2.20)
            Eigen::Matrix<double, 1, 3> Hdepth= {0.0,0.0,1.0};
            C.row(2) = Hdepth;  

            return C;
          };
        */
            setEntityState(IMC::EntityState::ESTA_NORMAL, Status::CODE_ACTIVE);
          }

          //! Release resources.
          void onResourceRelease(void)
          {
          }

          void consume(const IMC::TBRFishTag* msg) {

            // Iterate through receivers
            for(unsigned i=0;i<c_receivers;i++) {
              if (m_args.receiver_serial[i] == msg->serial_no)
              {
                debug(DTR("Message from R%d arrived"), i);
                m_tagDetection[i] = *msg;
                m_newDetection[i] = true;
                return;
              }
            }
          }

          //! Check if the RDOA indicates a time shift larger than accepted
          //! @param [in] RDOA Range Difference of Arrival 
          //! @return Boolean representing accepted/not accepted
          bool timeShiftCorrect(const double RDOA)
          {
            if((std::abs(RDOA) <= m_max_rdoa))
              return true;
            else
              return false;
          }

          //! Turns the latitude and longtitude of the input to a NED representation with refCoord as origin.
          //! @param [in] input Tag detection to take lat/lon [rad] from 
          //! @param [in] refCoord Reference coordinate in {lat[rad], lon [rad], elevation [m]} 
          //! @param [out] output NED frame representation of input in {North, East, Down} [meters] relative to the reference coordinate
          void toNEDframe(const IMC::TBRFishTag &input, const double refCoord[3], std::tuple<double, double, double> &output)
          {
            WGS84::displacement(refCoord[0], refCoord[1], refCoord[2], input.lat, input.lon, 0.0, &(std::get<0>(output)), &(std::get<1>(output)), &(std::get<2>(output)));
          }

          void findFishPosition(const IMC::TBRFishTag tagData[c_receivers]) {


            static Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, c_receivers,1> RDOA2;
            static Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, c_receivers,3> NED2;
            RDOA2.resize(0,1);
            NED2.resize(0,3);
            std::tuple<double, double, double> tempNED;

            std::vector<double> RDOA; // In meters
            unsigned used_counter = 0;
            //For loop to create unique order invariant permutations (Each receiver is combined with another only once.)
            std::vector<bool> used(c_receivers, false); // Initializes all to false
            for(uint8_t outerreceiver = 0;outerreceiver<c_receivers;++outerreceiver) {
              for(uint8_t receiver = outerreceiver+1;receiver<c_receivers;++receiver) {
                if( !(used[outerreceiver] && used[receiver]) && (m_newDetection[outerreceiver] || m_newDetection[receiver]) ) {
                  double tempRDOA = m_args.speed_of_sound_in_water*((static_cast<double>(tagData[receiver].unix_timestamp) - tagData[outerreceiver].unix_timestamp) * 1000.0 + (static_cast<double>(tagData[receiver].millis) - tagData[outerreceiver].millis))/1000.0;
                  if(timeShiftCorrect(tempRDOA)) {
                    inf("%d, %d", outerreceiver, receiver);
                    RDOA.push_back(tempRDOA);
                    RDOA2.resize(used_counter+1,1); // +1 because used_counter starts from 0.
                    RDOA2.row(used_counter) << tempRDOA;
                    used[outerreceiver] = true;
                    used[receiver] = true;
                    used_counter++;
                  }
                }
              }
            }
            if(RDOA2.rows() < 1) {
              inf("No valid RDOA2 combinations");
              return;
            }
            /*if(RDOA.empty()) {
              inf("No valid combinations");
              return;
            }*/

NED2.resize(used_counter+1,3);

            // Create NED representation of used detections, and remove them from the new detections
            std::vector<std::tuple<double, double, double>> NED(c_receivers);
            std::vector<double> depth; 

            for(std::vector<bool>::iterator it = used.begin(); it != used.end(); it++) {
              if(*it) { // Enter if receiver is used
                inf("%lu", it - used.begin());
                toNEDframe(tagData[it - used.begin()], m_refCoord, NED[it - used.begin()]);

toNEDframe(tagData[it - used.begin()], m_refCoord, tempNED);
NED2.row(it - used.begin()) << std::get<0>(tempNED), std::get<1>(tempNED),std::get<2>(tempNED);

                m_newDetection[it - used.begin()] = false;
                depth.push_back(tagData[it - used.begin()].trans_data*0.392); //0.392 from S256 data spec?
              }
            }
            // FOR DEBUGGING
            for(std::vector<bool>::iterator it = used.begin(); it != used.end(); it++) {
              if(*it) { // Enter if receiver is used
                debug("%f, %f, %f", std::get<0>(NED[it - used.begin()]), std::get<1>(NED[it - used.begin()]),std::get<2>(NED[it - used.begin()]));
              }
            } 
            for(std::vector<double>::iterator it = RDOA.begin(); it != RDOA.end(); it++) {
              inf("RDOA %f", *it);
            }
inf("NED2 POST");
std::cout << NED2 << std::endl;
inf("RDOA2");
std::cout << RDOA2 << std::endl;

            //RDOA2.resize(0,0);
            //NED2.resize(0,0);
            inf("End findFishPos");
          }
/*
          void findFishPosition(const IMC::TBRFishTag tagData[c_receivers]) {
              uint8_t new_detections=0;
              uint8_t last_valid_detection = 0;
              //Find last valid receiver and count valid receivers
              for(uint8_t receiver = 0;receiver<c_receivers;receiver++) {
                  if(m_newDetection[receiver]) {
                      new_detections++;
                      last_valid_detection = receiver;
                  }
              }
              if(new_detections > 1) {
                spew("New detections: %u, Last valid detection: %u", new_detections, last_valid_detection);
                // Build RDOA vector
                double *RDOA = new double(new_detections-1);
                uint8_t combinations = 0;
                for(uint8_t receiver = 0;receiver<last_valid_receiver;receiver++) {
                  if(m_newDetection[receiver]) {
                      double tempRDOA = m_args.speed_of_sound_in_water*((static_cast<double>(tagData[receiver].unix_timestamp) - tagData[last_valid_receiver].unix_timestamp) * 1000.0 + (static_cast<double>(tagData[receiver].millis) - tagData[last_valid_receiver].millis))/1000.0;
                      if(timeShiftCorrect(tempRDOA)) {
                        RDOA[combinations] = tempRDOA;
                        combinations++;
                        m_newDetection[receiver] = false;
                      }
                  }
                }

                delete [] RDOA;
              }
          }*/
          //! Main loop.
          void
          task(void)
          {
              findFishPosition(m_tagDetection);
              inf("other side");
          } // End of task function
        };
      }
    }

DUNE_TASK
