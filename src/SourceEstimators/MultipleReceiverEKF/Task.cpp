//***************************************************************************
// Copyright 2013-2021 Norwegian University of Science and Technology (NTNU)*
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
// Author: Nikolai Lauvås (Based on the works of:                           *
//         Praveen Jain and Artur Zolich                                    *
//***************************************************************************

// DUNE headers.
#include <DUNE/DUNE.hpp>
#include <Eigen/Core>
#include <algorithm>  // std::count
#include <iterator>

  namespace SourceEstimators
  {
      namespace MultipleReceiverEKF
      {
        using DUNE_NAMESPACES;

        const uint8_t c_receivers = 3; // TODO: Make everything scalable according to this number so x receivers can be used

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

            // Setup processing of IMC messages
            // Initialize messages.
            bind<IMC::TBRFishTag>(this);
          }

          //! Update internal state with new parameter values.
          void
          onUpdateParameters(void)
          {
            m_max_rdoa = m_args.max_time_shift_ms/1000*m_args.speed_of_sound_in_water;

            m_refCoord[0] = Math::Angles::radians(m_args.reference[0]);
            m_refCoord[1] = Math::Angles::radians(m_args.reference[1]);
            m_refCoord[2] = 0.0;

            m_initPosition[0] = m_args.position[0];
            m_initPosition[1] = m_args.position[1];
            m_initPosition[2] = m_args.position[2];

            debug(DTR("Reference Cooridnate: %f %f"), m_args.reference[0], m_args.reference[1]);
            debug(DTR("Fish Position: %f %f %f"), m_args.position[0], m_args.position[1], m_args.position[2]);
/*
            for (unsigned i = 0; i<c_receivers;i++) {
              m_o[i].open("/home/nikolai/ExperimentDataSet/Set4/S" + std::to_string(i) + "3output.txt", std::ios::out | std::ios::trunc);
              if(m_o[i].is_open())
              {
                debug(DTR("Stage %d open!"), i);
              } 
            }
*/
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
                //m_xkf.setActive(true);
                return;
              }
            }

            // This only gets executed if no if clause is fullfilled
            debug(DTR("consume: Message from non registered carrier with ID %d"), msg->getSource());
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
          void toNEDframe(const IMC::TBRFishTag &input, const double refCoord[3], double (&output)[3])
          {
            double input_d[3] = {input.lat, input.lon, 0.0};
            toNEDframe(input_d,refCoord, output);
          }

          //! Turns the latitude and longtitude of the input to a NED representation with refCoord as origin.
          //! @param [in] input Location in WGS84 {lat[rad], lon [rad], elevation [m]} 
          //! @param [in] refCoord Reference coordinate in WGS84 {lat[rad], lon [rad], elevation [m]} 
          //! @param [out] output NED frame representation of input in {North, East, Down} [meters] relative to the reference coordinate
          void toNEDframe(const double input[3], const double refCoord[3], double (&output)[3])
          {
            WGS84::displacement(refCoord[0], refCoord[1], refCoord[2], input[0], input[1], input[2], &(output[0]), &(output[1]), &(output[2]));
          }

          //! Takes a NED frame position and transforms it to a WGS84 lat/lon/elevation position
          //! @param [in] input NED frame position to transform {North, East, Down} [meters] relative to the reference coordinate
          //! @param [in] refCoord Reference coordinate in WGS84 {lat[rad], lon [rad], elevation [m]} 
          //! @param [out] output Input position converted to WGS84 coordinates {lat[rad], lon [rad], elevation [m]} 
          void fromNEDframe(const double input[3], const double refCoord[3], double (&output)[3]) {
            output[0] = refCoord[0];
            output[1] = refCoord[1];
            output[2] = refCoord[2];
            WGS84::displace(input[0], input[1], input[2], &(output[0]), &(output[1]), &(output[2]));
          }
/*
          //! Function responsible for loggin results from the three stages of the estimator
          void logFishPosition()
          {
            double result[3][3] = {{m_xkf.stage1.xHat(0,0), m_xkf.stage1.xHat(1,0), m_xkf.stage1.xHat(2,0)},
                                   {m_xkf.stage2.xHat(0,0), m_xkf.stage2.xHat(1,0), m_xkf.stage2.xHat(2,0)},
                                   {m_xkf.stage3.xHat(0,0), m_xkf.stage3.xHat(1,0), m_xkf.stage3.xHat(2,0)}
                                   };
            double latlon[3];
            //Log timestamps and estimates for stage 1-3
            for (unsigned i = 0; i<3;i++) {
              m_o[i].precision(15);
              
              m_o[i] << m_tagDetection[0].unix_timestamp << m_tagDetection[0].millis << "," << m_tagDetection[1].unix_timestamp << m_tagDetection[1].millis << "," << m_tagDetection[2].unix_timestamp << m_tagDetection[2].millis << ","; 
              fromNEDframe(result[i], m_refCoord, latlon);
              m_o[i] << DUNE::Math::Angles::degrees(latlon[0]) << "," << DUNE::Math::Angles::degrees(latlon[1]) << "," << latlon[2];
            }

            // Log Stage 1: LS Estimate
            m_o[0] << "," << m_xkf.stage1.m_dr << "," << Clock::getSinceEpochMsec() << std::endl;

            // Log Stage 2: LTV KF Estimate
            m_o[1] << "," << m_xkf.stage2.innov.norm_p(2) << "," << m_xkf.stage2.PHat.norm_p(2) << "," << m_xkf.stage2.PHat.trace() << "," << Clock::getSinceEpochMsec() << std::endl;

            // Log Stage 3: Linearized KF Estimate
            m_o[2] << "," << m_xkf.stage3.innov.norm_p(2) << "," << m_xkf.stage3.PHat.norm_p(2) << "," << m_xkf.stage3.PHat.trace() << "," << Clock::getSinceEpochMsec() << std::endl;
            
            // Send output to Neptus/DUNE log
            IMC::RemoteSensorInfo tagPosition;
            tagPosition.lat = latlon[0];
            tagPosition.lon = latlon[1];
            tagPosition.alt = -latlon[2];
            tagPosition.data = std::to_string(m_tagDetection[0].unix_timestamp) + std::to_string(m_tagDetection[0].millis) + "," + std::to_string(m_tagDetection[1].unix_timestamp) + std::to_string(m_tagDetection[1].millis) + "," + std::to_string(m_tagDetection[2].unix_timestamp) + std::to_string(m_tagDetection[2].millis);
            // + "," + m_xkf.stage3.innov.norm_p(2) + "," << m_xkf.stage3.PHat.norm_p(2) + "," << m_xkf.stage3.PHat.trace();
            tagPosition.id = "FishTagXKF";
            dispatch(tagPosition);
          }
*/

          //! Function that prepares the measurements for use in the estimator
          //! @param [in] receiver1 The most recent detection made by a receiver
          //! @param [in] receiver2 The most recent detection made by a receiver
          //! @param [in] receiver3 The most recent detection made by a receiver
          void findFishPosition(const IMC::TBRFishTag receiver[c_receivers]) {
            // Step 0: Ensure Initialized Kalman Filter
            /*if(!m_xkf.isInitialized())
            {
              m_xkf.initialize(m_initPosition);
            }*/

            // Step 1: Predict
            //m_xkf.predict();

            // Step 2: Update (checks for validity and available measurements)
            double RDOA_i[] = {m_args.speed_of_sound_in_water*((static_cast<double>(receiver[0].unix_timestamp) - receiver[1].unix_timestamp) * 1000.0 + (static_cast<double>(receiver[0].millis) - receiver[1].millis))/1000.0,
                               m_args.speed_of_sound_in_water*((static_cast<double>(receiver[1].unix_timestamp) - receiver[2].unix_timestamp) * 1000.0 + (static_cast<double>(receiver[1].millis) - receiver[2].millis))/1000.0,
                               m_args.speed_of_sound_in_water*((static_cast<double>(receiver[2].unix_timestamp) - receiver[0].unix_timestamp) * 1000.0 + (static_cast<double>(receiver[2].millis) - receiver[0].millis))/1000.0};

            DUNE::Math::Matrix RDOA(RDOA_i,3,1);
            bool newMeasurement[c_receivers] = {m_newDetection[0], m_newDetection[1], m_newDetection[2]};
            bool validCombination[c_receivers] = {((m_newDetection[0] || m_newDetection[1]) && (timeShiftCorrect(RDOA(0,0)))), 
                                                  ((m_newDetection[1] || m_newDetection[2]) && (timeShiftCorrect(RDOA(1,0)))),
                                                  ((m_newDetection[2] || m_newDetection[0]) && (timeShiftCorrect(RDOA(2,0))))
                                                  };
            std::vector<bool> validCombination_vec(validCombination, validCombination+c_receivers);
            if(validCombination[0] || validCombination[1] || validCombination[2])
            {

              //>> Step 2.1: Update System Matrices and Covariance
              double NED1[3];
              double NED2[3];
              double NED3[3];
              toNEDframe(receiver[0], m_refCoord, NED1);
              toNEDframe(receiver[1], m_refCoord, NED2);
              toNEDframe(receiver[2], m_refCoord, NED3);

              Eigen::Matrix<double, 3,3> receiverPositions_e;
              receiverPositions_e << NED1[0], NED2[0],NED3[0],
                                     NED1[1], NED2[1],NED3[1],
                                     NED1[2], NED2[2],NED3[2];
              Eigen::Matrix<double, 3,1> tagDepth_e;
              tagDepth_e << receiver[0].trans_data*0.392, receiver[1].trans_data*0.392, receiver[2].trans_data*0.392;

              double receiverPositions_i[] = {NED1[0], NED2[0],NED3[0],
                                              NED1[1], NED2[1],NED3[1],
                                              NED1[2], NED2[2],NED3[2]};

              double tagDepth_i[] = {receiver[0].trans_data*0.392, receiver[1].trans_data*0.392, receiver[2].trans_data*0.392}; 
              DUNE::Math::Matrix tagDepth(tagDepth_i, 3, 1);

              DUNE::Math::Matrix receiverPositions(receiverPositions_i, 3,3);

std::cout << "ValidCombinations" << std::endl << validCombination[0]*1 << std::endl;
std::cout << "ValidCombinations" << std::endl << validCombination[1]*1 << std::endl;
std::cout << "ValidCombinations" << std::endl << validCombination[2]*1 << std::endl;
std::cout << "ValidCombinations_t" << std::endl << std::count(validCombination_vec.begin(), validCombination_vec.end(), true) << std::endl;

std::cout << "receiverPositions" << std::endl << receiverPositions << std::endl;
std::cout << "receiverPositions_e" << std::endl << receiverPositions_e << std::endl;
            

              //m_xkf.update(receiverPositions, RDOA, tagDepth, newMeasurement, validCombination);
///////////
{
      double qN = 10;
      double qE = 10;
      double qD = 5;
Eigen::Matrix<double, 3, 1> q;
q << qN, qE, qD;

Eigen::Matrix<double, 3, 3> delta = receiverPositions_e.colwise()-q;
Eigen::Matrix<double, 3, 1> d = delta.colwise().norm();;
      double d1 = std::sqrt(std::pow((receiverPositions(0,0) - qN),2) + std::pow((receiverPositions(1,0) - qE),2) + std::pow((receiverPositions(2,0) - qD),2));
      double d2 = std::sqrt(std::pow((receiverPositions(0,1) - qN),2) + std::pow((receiverPositions(1,1) - qE),2) + std::pow((receiverPositions(2,1) - qD),2));
      double d3 = std::sqrt(std::pow((receiverPositions(0,2) - qN),2) + std::pow((receiverPositions(1,2) - qE),2) + std::pow((receiverPositions(2,2) - qD),2));

//std::cout << "d" << std::endl << d << std::endl;
//std::cout << "d1" << std::endl << d1 << std::endl;
//std::cout << "d2" << std::endl << d2 << std::endl;
//std::cout << "d3" << std::endl << d3 << std::endl;
//std::cout << "delta" << std::endl << delta << std::endl;
//std::cout << "receiverPositions(0,0) - qN" << std::endl << receiverPositions(0,0) - qN << std::endl;
//std::cout << "receiverPositions(1,0) - qE" << std::endl << receiverPositions(1,0) - qE << std::endl;
//std::cout << "receiverPositions(2,0) - qD" << std::endl << receiverPositions(2,0) - qD << std::endl;

//delta.col(1)/d.row1-delta.col(0)/d.row(0);

      double c11 = (receiverPositions(0,1) - qN)/(d2) - (receiverPositions(0,0) - qN)/(d1);
      double c12 = (receiverPositions(1,1) - qE)/(d2) - (receiverPositions(1,0) - qE)/(d1);
      double c13 = (receiverPositions(2,1) - qD)/(d2) - (receiverPositions(2,0) - qD)/(d1);

//std::cout << "delta.col(1)/d.row1-delta.col(0)/d.row(0)" << std::endl << delta(0,1)/d(1,0)-delta(0,0)/d(0,0) << std::endl;
//std::cout << "c11" << std::endl << c11 << std::endl;

      double c21 = (receiverPositions(0,2) - qN)/(d3) - (receiverPositions(0,1) - qN)/(d2);
      double c22 = (receiverPositions(1,2) - qE)/(d3) - (receiverPositions(1,1) - qE)/(d2);
      double c23 = (receiverPositions(2,2) - qD)/(d3) - (receiverPositions(2,1) - qD)/(d2);

      double c31 = (receiverPositions(0,0) - qN)/(d1) - (receiverPositions(0,2) - qN)/(d3);
      double c32 = (receiverPositions(1,0) - qE)/(d1) - (receiverPositions(1,2) - qE)/(d3);
      double c33 = (receiverPositions(2,0) - qD)/(d1) - (receiverPositions(2,2) - qD)/(d3);

      double c41 = 0;
      double c42 = 0;
      double c43 = 1;

      double C_v[] = {c11, c12, c13,
                      c21, c22, c23,
                      -c31, -c32, -c33, // Hvorfor - her?
                      c41, c42, c43};
DUNE::Math::Matrix C(C_v, 4,3);
Eigen::Matrix<double, 4, 3> c_e;
c_e <<   delta(0,1)/d(1,0)-delta(0,0)/d(0,0),    delta(1,1)/d(1,0)-delta(1,0)/d(0,0),    delta(2,1)/d(1,0)-delta(2,0)/d(0,0),
         delta(0,2)/d(2,0)-delta(0,1)/d(1,0),    delta(1,2)/d(2,0)-delta(1,1)/d(1,0),    delta(2,2)/d(2,0)-delta(2,1)/d(1,0),
       -(delta(0,0)/d(0,0)-delta(0,2)/d(2,0)), -(delta(1,0)/d(0,0)-delta(1,2)/d(2,0)), -(delta(2,0)/d(0,0)-delta(2,2)/d(2,0)),
         0,0,1;

std::cout << "C" << std::endl << C << std::endl;
std::cout << "c_e" << std::endl << c_e << std::endl;
}
            }
///////////

      /*        if(validCombination[0]) {
                m_newDetection[0]= false;
                m_newDetection[1]= false;
              }
              if(validCombination[1]) {
                m_newDetection[1]= false;
                m_newDetection[2]= false;
              }
              if(validCombination[2]) {
                m_newDetection[2]= false;
                m_newDetection[0]= false;
              }                            
            }*/
            //logFishPosition();
          }

          //! Main loop.
          void
          task(void)
          {
            findFishPosition(m_tagDetection);
            /*if(m_xkf.isActive())
            {
              findFishPosition(m_tagDetection[0], m_tagDetection[1], m_tagDetection[2]);

            } // End of if*/
          } // End of task function
        };
      }
    }

DUNE_TASK
