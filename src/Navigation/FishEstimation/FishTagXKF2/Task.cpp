//***************************************************************************
// Copyright 2007-2016 Universidade do Porto - Faculdade de Engenharia      *
// Laboratório de Sistemas e Tecnologia Subaquática (LSTS)                  *
//***************************************************************************
// This file is part of DUNE: Unified Navigation Environment.               *
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
// Author: Praveen Jain and Artur Zolich                                                  *
//***************************************************************************

// DUNE headers.
#include <DUNE/DUNE.hpp>
#include <iterator>

namespace
{
  namespace Navigation
  {
    namespace FishEstimation
    {
      namespace FishTagXKF2
      {
        using DUNE_NAMESPACES;

        struct ReceiverData
        {
          bool valid;
          bool isMeasure;

          int serial_no;
          long unix_timestamp;
          long milisecond_timestamp;
          double unix_milisecond_time;
          std::string code_type;
          int tag_id;
          double sensor_data; // assumed depth
          int snr;
          int memory_number;

          double lat;
          double lon;
          double N;
          double E;
          double D;

          int fix;
          float hdop;
          int satellites;
          uint8_t checksum;
          std::string msg;
        };

        struct RefCoordinate
        {
          double lat;
          double lon;
          double hae;
        };

        struct FishPosition
        {
          double lat;
          double lon;
          double hae;
          double n;
          double e;
          double d;
        };

        struct EKFilter
        {
          bool isEKFActive;
          bool isEKFInitialized;

          uint8_t nx;     // Num states
          uint8_t ny;     // Num outputs
          uint8_t nyk;    // Num time varying outputs
          uint8_t nu;     // Num inputs
          uint8_t nd;     // Num unknown inputs - noise
          double qq_cov;  // Process noise Covariance - position
          double rr_cov;  // Measurement noise Covariance - range
          double rz_cov;  // Measurement noise Covariance - depth

          Matrix xHat;    // State estimates
          Matrix PHat;    // State Covariance estimates
          Matrix A;       // State transition matrix
          Matrix B;       // Input Matrix
          Matrix C;       // Observation Matrix
          Matrix Ck;      // Time varying observation matrix
          Matrix D;       // unknown input matrix
          Matrix Q;       // Process noise covariance matrix
          Matrix R;       // Measurement noise covariance matrix
          Matrix Rk;      // Time varying Measurement noise covariance matrix
          Matrix innov;   // Time varying innovation vector
          Matrix yk;      // Measurement vector
          Matrix ykest;   // Estimate of the measurements
        };

        //! %Task arguments.
        struct Arguments
        {
          //! Receiver 1 carrier
          std::string rec_1_vehicle_name;
          //! Receiver 2 carrier
          std::string rec_2_vehicle_name;
          //! Receiver 3 carrier
          std::string rec_3_vehicle_name;

          //! Maximum allowed time [ms] shift between receivers' messages
          int max_time_shift_ms;
          int speed_of_sound_in_water; // speed of sound in water [m/s]
          double dt;

          //! Reference coordinate position (degrees)
          std::vector<double> reference;

          //! Initial position of the fish tag (NED) for EKF
          std::vector<double> position;

          uint8_t num_states;

          bool isXKFInitialized;
          bool isXKFActive;
        };

        struct Task: public DUNE::Tasks::Periodic
        {

          //! Task arguments
          Arguments m_args;

          //! Receiver 1 carrier system ID.
          unsigned int m_r1_name;
          //! Receiver 2 carrier system ID.
          unsigned int m_r2_name;
          //! Receiver 3 carrier system ID.
          unsigned int m_r3_name;

          //! Structure that hold the parsed fish tag data
          ReceiverData m_r1, m_r2, m_r3;

          //! Structure to hold the EKF (LTV-KF) data
          EKFilter m_kf_stage2, m_kf_stage3;
          double m_dr; // set by the LeastSquaresEstimate function
          FishPosition m_fp_stage1, m_fp_stage2, m_fp_stage3;

          //! Output file only for use in simulations
          std::ofstream m_o1, m_o2, m_o3, m_rlog;

          //! Reference coordinate system
          RefCoordinate m_refCoord;

          //! Constructor.
          //! @param[in] name task name.
          //! @param[in] ctx context.
          Task(const std::string& name, Tasks::Context& ctx):
            DUNE::Tasks::Periodic(name, ctx)
          {

            param("Receiver 1 Carrier Name", m_args.rec_1_vehicle_name)
            .description("Receiver 1 Carrier Name")
            .defaultValue("alfa-07");

            param("Receiver 2 Carrier Name", m_args.rec_2_vehicle_name)
            .description("Receiver 2 Carrier Name")
            .defaultValue("alfa-07");

            param("Receiver 3 Carrier Name", m_args.rec_3_vehicle_name)
            .description("Receiver 3 Carrier Name")
            .defaultValue("alfa-07");

            param("Max time shift [ms]", m_args.max_time_shift_ms)
            .description("Maximum allowed time [ms] shift between receivers' messages")
            .defaultValue("500");

            param("Speed of sound in water [m/s]", m_args.speed_of_sound_in_water)
            .description("Speed of sound in water [m/s]")
            .defaultValue("1484");

            param("FishPos Cov", m_kf_stage2.qq_cov)
            .description("Fish Position Covariance")
            .defaultValue("0");

            param("ToA Cov", m_kf_stage2.rr_cov)
            .description("Time of Arrival Covariance")
            .defaultValue("0");

            param("Depth Cov", m_kf_stage2.rz_cov)
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
            bind<IMC::DevDataText>(this);
            //bind<IMC::EstimatedState>(this);
          }

          //! Update internal state with new parameter values.
          void
          onUpdateParameters(void)
          {
            m_args.dt = 1/getFrequency();

            m_refCoord.lat = Math::Angles::radians(m_args.reference[0]);
            m_refCoord.lon = Math::Angles::radians(m_args.reference[1]);
            m_refCoord.hae = 0.0;

            // m_refCoord.lat = m_args.reference[0];
            // m_refCoord.lon = m_args.reference[1];
            // m_refCoord.hae = 0.0;

            m_args.num_states = 3;

            debug(DTR("Reference Cooridnate: %f %f"), m_args.reference[0], m_args.reference[1]);
            debug(DTR("Fish Position: %f %f %f"), m_args.position[0], m_args.position[1], m_args.position[2]);

            // m_r1_name = resolveSystemName(m_args.rec_1_vehicle_name);
            // m_r2_name = resolveSystemName(m_args.rec_2_vehicle_name);
            // m_r3_name = resolveSystemName(m_args.rec_3_vehicle_name);
            // inf(DTR("Carrier of receiver 1 name is %s, with ID %d"), m_args.rec_1_vehicle_name.c_str(), m_r1_name);
            // inf(DTR("Carrier of receiver 2 name is %s, with ID %d"), m_args.rec_2_vehicle_name.c_str(), m_r2_name);
            // inf(DTR("Carrier of receiver 3 name is %s, with ID %d"), m_args.rec_3_vehicle_name.c_str(), m_r3_name);

            //>> Used for simulation purposes only
            m_r1_name = 632; // Duckling 1
            m_r2_name = 634; // Duckling 2
            m_r3_name = 631; // Duckling 3

            m_o1.open("/home/nikolai/ExperimentDataSet/Set4/S1output.txt", std::ios::out | std::ios::trunc);
            if(m_o1.is_open())
            {
              debug(DTR("Stage 1 open!"));
            }

            m_o2.open("/home/nikolai/ExperimentDataSet/Set4/S2output.txt", std::ios::out | std::ios::trunc);
            if(m_o2.is_open())
            {
              debug(DTR("Stage 2 open!"));
            }

            m_o3.open("/home/nikolai/ExperimentDataSet/Set4/S3output.txt", std::ios::out | std::ios::trunc);
            if(m_o3.is_open())
            {
              debug(DTR("Stage 3 open!"));
            }

            m_rlog.open("/home/nikolai/ExperimentDataSet/Set4/Routput.txt", std::ios::out | std::ios::trunc);
            if(m_rlog.is_open())
            {
              debug(DTR("Stage 3 open!"));
            }
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
          void
          onResourceRelease(void)
          {
          }

          uint8_t nmeaCheckSum(std::string sentence)
          {
            uint8_t check = 0;

            std::size_t star = sentence.find('*');
            if(star != std::string::npos)
            {
              // iterate over the string, XOR each byte with the total sum:
              for (unsigned int c = 1; c < star; c++)
              {
                //  debug(DTR("CHK: %c"), sentence.at(c));
                check = uint8_t(check ^ (uint8_t)sentence.at(c));
              }
            }
            debug(DTR("Checksum %x"), check);
            // return the result
            return check;
          }

          std::vector<std::string>
          splitString( std::string s )
          {
            std::replace( s.begin(), s.end(), '$', ' ' );
            std::replace( s.begin(), s.end(), '*', ' ' );
            std::replace( s.begin(), s.end(), ',', ' ' );
            std::size_t found = s.find("  ");
            if(found != std::string::npos)
            {
              s.replace( found, 2, " " );
              debug(DTR("Found double space"));
            }

            debug(DTR("MSG %s"), s.c_str());

            std::istringstream stream( s );
            std::vector<std::string> result;
            for( ;; )
            {
              std::string word;
              if( !( stream >> word ) )
              {
                break;
              }
              result.push_back(word);
            }
            return result;
          }

          //>> Function: parseFakeFishTagData used for simulation purposes only
          ReceiverData
          parseFakeFishTagData(std::string msg)
          {
            ReceiverData data;

            debug(DTR("Message: %s"), msg.c_str());
            std::vector<std::string> data_fields = splitString(msg);

            data.valid = false;

            if(data_fields.size() == 10){
              data.msg = msg.c_str();
              data.serial_no = std::atoi(data_fields[0].c_str());
              data.unix_milisecond_time = std::atof(data_fields[1].c_str());
              debug(DTR("Unix ms time %f"), data.unix_milisecond_time);
              data.sensor_data = std::atof(data_fields[5].c_str()); // assumed depth
              debug(DTR("Depth %d"), data.sensor_data);
              data.lat = std::atof(data_fields[2].c_str());
              debug(DTR("Decimal Lat %f"), data.lat);
              data.lon = std::atof(data_fields[3].c_str());
              debug(DTR("Decimal Long %f"), data.lon);
              data.valid = true;
            }
            else
            {
              inf(DTR("Wrong number of fields in the message, received %d"), (int)data_fields.size());
            }

            return data;
          } // End parseFakeFishTagData

          //Function: Parse the Fish tag data
          ReceiverData
          parseFishTagData(std::string msg)
          {
            ReceiverData data;

            debug(DTR("Message: %s"), msg.c_str());
            std::vector<std::string> data_fields = splitString(msg);

            data.valid = false;

            //for(int i=0; i<17; i++){
            //  inf(DTR("%d: %s"), i, data_fields[i].c_str());
            //}
            if(data_fields.size() == 18){  // was 17
                  data.msg = msg.c_str();
                  data.serial_no = std::atoi(data_fields[1].c_str());

                  data.unix_timestamp = std::atol(data_fields[2].c_str());
                  data.milisecond_timestamp = std::atol(data_fields[3].c_str());
                  data.unix_milisecond_time = data.unix_timestamp*1000 + data.milisecond_timestamp;
                  debug(DTR("Unix ms time %f"), data.unix_milisecond_time);

                  data.code_type = data_fields[4];
                  data.tag_id = std::atoi(data_fields[5].c_str());
                  data.sensor_data = std::atof(data_fields[6].c_str()) * 0.392; // assumed depth
                  debug(DTR("Depth %d"), data.sensor_data);

                  data.snr = std::atoi(data_fields[7].c_str());
                  data.memory_number = std::atoi(data_fields[9].c_str()); // was 8 and so on

                  double tmpLat = std::atof(data_fields[10].c_str());
                  int tmpDeg = (int)(tmpLat / 100.0);
                  double tmpMin = tmpLat - tmpDeg*100.0;
                  data.lat = Angles::convertDMSToDecimal(tmpDeg, tmpMin);
                  if(data_fields[11]=="S")
                    data.lat = -1.0 * data.lat;

                  debug(DTR("Decimal Lat %f"), data.lat);
                  data.lat = Angles::radians(data.lat);

                  double tmpLon = std::atof(data_fields[12].c_str());
                  tmpDeg = (int)(tmpLon / 100.0);
                  tmpMin = tmpLon - tmpDeg*100.0;
                  data.lon = Angles::convertDMSToDecimal(tmpDeg, tmpMin);
                  if(data_fields[13]=="W")
                    data.lon = -1.0 * data.lon;

                  debug(DTR("Decimal Long %f"), data.lon);
                  data.lon = Angles::radians(data.lon);
                  data.fix = std::atoi(data_fields[14].c_str());
                  data.satellites = std::atoi(data_fields[15].c_str());
                  data.hdop = std::atoi(data_fields[16].c_str());
                  data.checksum = strtol(data_fields[17].c_str(), NULL, 16);

                  debug(DTR("Data CHK %s, %d"), data_fields[17].c_str(), data.checksum);
                  uint8_t checksum = nmeaCheckSum(msg);
                  if(checksum == data.checksum)
                  {
                    data.valid = true;
                    debug(DTR("Checksum OK"));
                  }
                  else
                  {
                    data.valid = false;
                    inf(DTR("Checksum Failed"));
                  }
                }
                else
                {
                  inf(DTR("Wrong number of fields in the message, received %d"), (int)data_fields.size());
                }

            return data;

          } // End Function: parseFishTagData

          //! Check if time shift between messages from 3 different sensors is within acceptable limits defined by user
          bool timeShiftCorrect(ReceiverData* r1, ReceiverData* r2)
          {
            double t12 = std::abs(r1->unix_milisecond_time - r2->unix_milisecond_time);

            debug(DTR("Timeshifts: %f"), t12);

            if((t12 <= m_args.max_time_shift_ms))
              return true;
            else
              return false;
          } // End function timeShiftCorrect

          std::string
          atoi(int a)
          {
            std::stringstream ss;
            ss << a;
            std::string str = ss.str();
            return str;
          }

          void
          createNEDframe(ReceiverData *r1, ReceiverData *r2, ReceiverData *r3)
          {
            WGS84::displacement(m_refCoord.lat, m_refCoord.lon, 0.0, r1->lat, r1->lon, 0.0, &(r1->N), &(r1->E), &(r1->D));
            WGS84::displacement(m_refCoord.lat, m_refCoord.lon, 0.0, r2->lat, r2->lon, 0.0, &(r2->N), &(r2->E), &(r2->D));
            WGS84::displacement(m_refCoord.lat, m_refCoord.lon, 0.0, r3->lat, r3->lon, 0.0, &(r3->N), &(r3->E), &(r3->D));
          }

/*          void
          consume(const IMC::EstimatedState* msg)
          {
            inf("Estimated arrived from %d", msg->getSource());
          }*/

          //! Message from FishTag Receivers
          void
          consume(const IMC::DevDataText* msg)
          {
            if((m_r1_name == msg->getSource()) || (m_r2_name == msg->getSource()) || (m_r3_name == msg->getSource()))
            {
              ReceiverData tmp;
              tmp = parseFishTagData(msg->value);

              IMC::RemoteSensorInfo tagPosition;
              tagPosition.lat = tmp.lat;
              tagPosition.lon = tmp.lon;
              tagPosition.alt = -tmp.sensor_data;
              tagPosition.data = tmp.msg;

              if(tmp.valid)
              {
                if (m_r1_name == msg->getSource())
                {
                  debug(DTR("consume: Message from R1 arrived"));
                  m_r1 = tmp;
                  m_r1.isMeasure = 1;
                  tagPosition.id = "Receiver_1";
                  dispatch(tagPosition);
                  m_args.isXKFActive = 1;
                }
                else if (m_r2_name == msg->getSource())
                {
                  debug(DTR("consume: Message from R2 arrived"));
                  m_r2 = tmp;
                  m_r2.isMeasure = 1;
                  tagPosition.id = "Receiver_2";
                  dispatch(tagPosition);
                  m_args.isXKFActive = 1;
                }
                else if (m_r3_name == msg->getSource())
                {
                  debug(DTR("consume: Message from R3 arrived"));
                  m_r3 = tmp;
                  m_r3.isMeasure = 1;
                  tagPosition.id = "Receiver_3";
                  dispatch(tagPosition);
                  m_args.isXKFActive = 1;
                }
                else
                {
                  inf(DTR("consume: Message from non registered carrier with ID %d"), msg->getSource());
                }
              }
              else
              {
                inf(DTR("Invalid"));
              }
            }
          } // End of consume DevDataText


          //! Initialize Extended Kalman filter
          void
          initializeXKF(void)
          {
            m_kf_stage2.isEKFActive = 0;
            m_kf_stage3.isEKFActive = 0;
            //>> Initialization of Stage 2: LTV Kalman filter
            m_kf_stage2.nx = 3;
            m_kf_stage2.nu = 0;
            m_kf_stage2.ny = 3;
            m_kf_stage2.nyk = m_kf_stage2.ny; //TODO: Do something with this
            m_kf_stage2.nd = 3; //TODO: remove this hardcoding

            m_kf_stage2.A.resizeAndFill(m_kf_stage2.nx, m_kf_stage2.nx,0.0);
            m_kf_stage2.C.resizeAndFill(m_kf_stage2.ny, m_kf_stage2.nx,0.0);
            m_kf_stage2.Ck.resizeAndFill(m_kf_stage2.nyk, m_kf_stage2.nx,0.0);
            m_kf_stage2.PHat.resizeAndFill(m_kf_stage2.nx, m_kf_stage2.nx,0.0);
            m_kf_stage2.Q.resizeAndFill(m_kf_stage2.nx, m_kf_stage2.nx,0.0);
            m_kf_stage2.R.resizeAndFill(m_kf_stage2.ny, m_kf_stage2.ny,0.0);
            m_kf_stage2.Rk.resizeAndFill(m_kf_stage2.nyk, m_kf_stage2.nyk,0.0);
            m_kf_stage2.xHat.resizeAndFill(m_kf_stage2.nx,1,0.0);
            m_kf_stage2.innov.resizeAndFill(m_kf_stage2.nyk,1,0.0);
            m_kf_stage2.yk.resizeAndFill(m_kf_stage2.nyk,1,0.0);
            m_kf_stage2.ykest.resizeAndFill(m_kf_stage2.nyk,1,0.0);

            if(m_kf_stage2.nu != 0)
            {
              m_kf_stage2.B.resizeAndFill(m_kf_stage2.nx, m_kf_stage2.nu, 0.0);
            }
            if(m_kf_stage2.nd != 0)
            {
              m_kf_stage2.D.resizeAndFill(m_kf_stage2.nx, m_kf_stage2.nd, 0.0);
            }

            // Set the state transition matrix (does not change)
            double A_v[] = {1, 0, 0,
                            0, 1, 0,
                            0, 0, 1};

            double Q_v[] = {m_kf_stage2.qq_cov, 0, 0,
                            0, m_kf_stage2.qq_cov, 0,
                            0, 0, m_kf_stage2.qq_cov/10};

            double x0_v[] = {m_args.position[0], m_args.position[1], m_args.position[2]};
            double P0_v[] = {100, 0, 0,
                             0, 100, 0,
                             0, 0, 10};

            double D_v[] = {m_args.dt*1, 0, 0,
                            0, m_args.dt*1, 0,
                            0, 0, m_args.dt*1};

            m_kf_stage2.A.fill(m_kf_stage2.nx, m_kf_stage2.nx, A_v);
            m_kf_stage2.PHat.fill(m_kf_stage2.nx, m_kf_stage2.nx, P0_v);
            m_kf_stage2.Q.fill(m_kf_stage2.nx, m_kf_stage2.nx, Q_v);
            m_kf_stage2.xHat.fill(m_kf_stage2.nx, 1, x0_v);
            m_kf_stage2.D.fill(m_kf_stage2.nx, m_kf_stage2.nd, D_v);

            //>> Initialization of stage 3 EKF
            m_kf_stage3.nx = 3;
            m_kf_stage3.nu = 0;
            m_kf_stage3.ny = 4;
            m_kf_stage3.nyk = m_kf_stage3.ny; //TODO: Do something with this
            m_kf_stage3.nd = 3; //TODO: remove this hardcoding

            m_kf_stage3.A.resizeAndFill(m_kf_stage3.nx, m_kf_stage3.nx,0.0);
            m_kf_stage3.C.resizeAndFill(m_kf_stage3.ny, m_kf_stage3.nx,0.0);
            m_kf_stage3.Ck.resizeAndFill(m_kf_stage3.nyk, m_kf_stage3.nx,0.0);
            m_kf_stage3.PHat.resizeAndFill(m_kf_stage3.nx, m_kf_stage3.nx,0.0);
            m_kf_stage3.Q.resizeAndFill(m_kf_stage3.nx, m_kf_stage3.nx,0.0);
            m_kf_stage3.R.resizeAndFill(m_kf_stage3.ny, m_kf_stage3.ny,0.0);
            m_kf_stage3.Rk.resizeAndFill(m_kf_stage3.nyk, m_kf_stage3.nyk,0.0);
            m_kf_stage3.xHat.resizeAndFill(m_kf_stage3.nx,1,0.0);
            m_kf_stage3.innov.resizeAndFill(m_kf_stage3.nyk,1,0.0);
            m_kf_stage3.yk.resizeAndFill(m_kf_stage3.nyk,1,0.0);
            m_kf_stage3.ykest.resizeAndFill(m_kf_stage3.nyk,1,0.0);

            if(m_kf_stage3.nu != 0)
            {
              m_kf_stage3.B.resizeAndFill(m_kf_stage3.nx, m_kf_stage3.nu, 0.0);
            }
            if(m_kf_stage3.nd != 0)
            {
              m_kf_stage3.D.resizeAndFill(m_kf_stage3.nx, m_kf_stage3.nd, 0.0);
            }
            // Set the state transition matrix (does not change)
            double A_v1[] = {1, 0, 0,
                            0, 1, 0,
                            0, 0, 1};

            double Q_v1[] = {m_kf_stage2.qq_cov, 0, 0,
                            0, m_kf_stage2.qq_cov, 0,
                            0, 0, m_kf_stage2.qq_cov/10};

            // // Initial position of Fish tag in [m]
            double x0_v1[] = {m_args.position[0], m_args.position[1], m_args.position[2]};
            double P0_v1[] = {100, 0, 0,
                             0, 100, 0,
                             0, 0, 10};

            double D_v1[] = {m_args.dt*1, 0, 0,
                            0, m_args.dt*1, 0,
                            0, 0, m_args.dt*1};

            m_kf_stage3.A.fill(m_kf_stage3.nx, m_kf_stage3.nx, A_v1);
            m_kf_stage3.PHat.fill(m_kf_stage3.nx, m_kf_stage3.nx, P0_v1);
            m_kf_stage3.Q.fill(m_kf_stage3.nx, m_kf_stage3.nx, Q_v1);
            m_kf_stage3.xHat.fill(m_kf_stage3.nx, 1, x0_v1);
            m_kf_stage3.D.fill(m_kf_stage3.nx, m_kf_stage3.nd, D_v1);
            // // m_kf_stage3.B.fill(m_kf_stage3.nx, m_kf_stage3.nu, B_v);

            // Compute fish position of stage 2 for logging
            m_fp_stage2.n = m_kf_stage2.xHat(0,0);
            m_fp_stage2.e = m_kf_stage2.xHat(1,0);
            m_fp_stage2.d = m_kf_stage2.xHat(2,0);

            m_fp_stage2.lat = m_refCoord.lat; m_fp_stage2.lon = m_refCoord.lon; m_fp_stage2.hae = m_refCoord.hae;
            WGS84::displace(m_fp_stage2.n, m_fp_stage2.e, m_fp_stage2.d, &(m_fp_stage2.lat), &(m_fp_stage2.lon), &(m_fp_stage2.hae));

            // Compute fish position of stage 3 for logging
            m_fp_stage3.n = m_kf_stage3.xHat(0,0);
            m_fp_stage3.e = m_kf_stage3.xHat(1,0);
            m_fp_stage3.d = m_kf_stage3.xHat(2,0);

            m_fp_stage3.lat = m_refCoord.lat; m_fp_stage3.lon = m_refCoord.lon; m_fp_stage3.hae = m_refCoord.hae;
            WGS84::displace(m_fp_stage3.n, m_fp_stage3.e, m_fp_stage3.d, &(m_fp_stage3.lat), &(m_fp_stage3.lon), &(m_fp_stage3.hae));
          } // End of initializeXKF function

          double
          resolveRAmbiguity(double R1, double R2)
          {
            double R;
            double max_range = 700; // TODO: Move this elsewhere

            if((R1 > 0.0) && (R1 < max_range))
            {
              if((R2 > 0.0) && (R2 < max_range))
              {
                debug(DTR("FishTagXKF2: resolveRAmbiguity: Cannot resolve ambiguity [Valid] %f %f"), R1, R2);
                R = R1; // choose one of them TODO: some trick here will help
              }
              else
              {
                debug(DTR("FishTagXKF2: resolveRAmbiguity: R1 valid! %f %f"), R1, R2);
                R = R1;
              }
            }
            else
            {
              if((R2 > 0.0) && (R2 < max_range))
              {
                debug(DTR("FishTagXKF2: resolveRAmbiguity: R2 valid! %f %f"), R1, R2);
                R = R2;
              }
              else
              {
                debug(DTR("FishTagXKF2: resolveRAmbiguity: Cannot resolve ambiguity [Invalid] %f %f"), R1, R2);
                R = 0;
              }
            }
            return R;
          } // End resolveRAmbiguity function

          uint8_t
          LeastSquaresEstimate(ReceiverData* r1, ReceiverData* r2, ReceiverData* r3)
          {
            uint8_t isSuccess = 0;

            // >> Construct Least squares and measurement matrices
            double d1 = ((r1->unix_milisecond_time - r3->unix_milisecond_time)*m_args.speed_of_sound_in_water)/1000.0;
            double d2 = ((r2->unix_milisecond_time - r3->unix_milisecond_time)*m_args.speed_of_sound_in_water)/1000.0;
            double depth = (m_r1.isMeasure*r1->sensor_data + m_r2.isMeasure*r2->sensor_data + m_r3.isMeasure*r3->sensor_data)/(m_r1.isMeasure + m_r2.isMeasure + m_r3.isMeasure);

            debug(DTR("Depth %d %f"), r1->sensor_data, depth);

            double Czq_v[] = {-(r1->N - r3->N), -(r1->E - r3->E), 0.0,
                              -(r2->N - r3->N), -(r2->E - r3->E), 0.0,
                                             0,                0, 0.5};

            Matrix Czq(Czq_v, 3, 3);
            double l_v[] = {d1,d2,0};
            Matrix l(l_v, 3, 1);

            double z_v[]  = {(d1*d1 - r1->N*r1->N - r1->E*r1->E + (r3->N)*(r3->N) + (r3->E)*(r3->E)),
                             (d2*d2 - r2->N*r2->N - r2->E*r2->E + (r3->N)*(r3->N) + (r3->E)*(r3->E)),
                             depth};

            Matrix z(z_v, 3, 1);

            Matrix invCzq = inverse(transpose(Czq)*Czq)*transpose(Czq);
            Matrix c = invCzq*l;
            Matrix w = 0.5*invCzq*z;

            // temp variables
            double ctc = c(0,0)*c(0,0) + c(1,0)*c(1,0) + c(2,0)*c(2,0);
            double ptc = r3->N*c(0,0) +  r3->E*c(1,0);
            double wtc = w(0,0)*c(0,0) + w(1,0)*c(1,0) + w(2,0)*c(2,0);
            double ptw = r3->N*w(0,0) +  r3->E*w(1,0);
            double wtw = w(0,0)*w(0,0) + w(1,0)*w(1,0) + w(2,0)*w(2,0);
            double ptp = (r3->N)*(r3->N) + (r3->E)*(r3->E);

            // quadratic equation coefficients for computation of d3 = m_dr, p = p_r
            double aa = 1 - ctc; // 1 - c'c
            double bb = 2*(ptc - wtc); // 2(p'c - w'c)
            double cc = 2*ptw - wtw - ptp; //2p'w - w'w - p'p
            debug(DTR("FishTagXKF2::LeastSquaresEstimate:: aa: %f, bb: %f, cc: %f, Czq det: %f invCzq det: %f"), aa,bb,cc, Czq.det(), invCzq.det());
            // Compute solution for quadratic term
            double R1, R2;
            double max_range = 700;
            Matrix fp_ls(3,1);
            if ((ctc == 1) || ((bb*bb - 4*aa*cc) <= 0.0))
            {
              if(ctc == 1)
              {
                R1 = -cc/bb;
              }
              else
              {
                R1 = -bb/(2*aa);
              }

              if((R1 > 0.0) && (R1 <= max_range))
              {
                fp_ls = (R1*c + w);
                m_dr = R1;
                debug(DTR("FishTagXKF2: LeastSquaresEstimate: Unique Solution %f"), R1);
                isSuccess = 1;
              }
              else
              {
                debug(DTR("FishTagXKF2: LeastSquaresEstimate: Unique Solution, invalid m_dr %f"), R1);
              }
            }
            else
            {
              debug(DTR("FishTagXKF2: LeastSquaresEstimate: Two possible solutions"));
              double s = sqrt(bb*bb - 4*aa*cc);
              R1 = (-bb + s)/(2*aa);
              R2 = (-bb - s)/(2*aa);
              m_rlog.precision(15);
              m_rlog << Clock::getSinceEpochMsec() << "," << R1 << "," << R2 << std::endl;
              R1 = resolveRAmbiguity(R1, R2);
              // Compute two candidate solutions - required to resolve ambiguity and select a particular R
              if(R1 == 0)
              {
                isSuccess = 0;
              }
              else
              {
                  fp_ls = (R1*c + w);
                  m_dr = R1;
                  isSuccess = 1;
              }
            }

            // Compute position of the source for logging
            if(isSuccess)
            {
              m_fp_stage1.n = fp_ls(0,0);
              m_fp_stage1.e = fp_ls(1,0);
              m_fp_stage1.d = fp_ls(2,0);

              m_fp_stage1.lat = m_refCoord.lat; m_fp_stage1.lon = m_refCoord.lon; m_fp_stage1.hae = m_refCoord.hae;
              WGS84::displace(fp_ls(0,0), fp_ls(1,0), fp_ls(2,0), &(m_fp_stage1.lat), &(m_fp_stage1.lon), &(m_fp_stage1.hae));
            }
            return isSuccess;
          } // End of LeastSquaresEstimate function

          void
          xkf_constructCandR_stage3(ReceiverData* r1, ReceiverData* r2, ReceiverData* r3)
          {
            //>> Get the state estimate
            Matrix x_est = m_kf_stage2.xHat;

            //>> Compute C and R matrices
            double qN = x_est.element(0,0);
            double qE = x_est.element(1,0);
            double qD = x_est.element(2,0);

            double d1 = std::sqrt(std::pow((r1->N - qN),2) + std::pow((r1->E - qE),2) + std::pow((r1->D - qD),2));
            double d2 = std::sqrt(std::pow((r2->N - qN),2) + std::pow((r2->E - qE),2) + std::pow((r2->D - qD),2));
            double d3 = std::sqrt(std::pow((r3->N - qN),2) + std::pow((r3->E - qE),2) + std::pow((r3->D - qD),2));

            // m_d1 = d1;
            // m_d2 = d2;
            // m_d3 = d3;

            debug(DTR("xkf_constructCandR_stage3: range d1: %f, d2: %f, d3: %f"), d1, d2, d3);
            double c11 = (r2->N - qN)/(d2) - (r1->N - qN)/(d1);
            double c12 = (r2->E - qE)/(d2) - (r1->E - qE)/(d1);
            double c13 = (r2->D - qD)/(d2) - (r1->D - qD)/(d1);

            double c21 = (r3->N - qN)/(d3) - (r2->N - qN)/(d2);
            double c22 = (r3->E - qE)/(d3) - (r2->E - qE)/(d2);
            double c23 = (r3->D - qD)/(d3) - (r2->D - qD)/(d2);

            double c31 = (r1->N - qN)/(d1) - (r3->N - qN)/(d3);
            double c32 = (r1->E - qE)/(d1) - (r3->E - qE)/(d3);
            double c33 = (r1->D - qD)/(d1) - (r3->D - qD)/(d3);

            double c41 = 0;
            double c42 = 0;
            double c43 = 1;

            double C_v[] = {c11, c12, c13,
                            c21, c22, c23,
                            -c31, -c32, -c33,
                            c41, c42, c43};

            // Stage 3 Covariance matrix - new with cross Covariance terms
            double R_v[] = {2*m_kf_stage2.rr_cov, m_kf_stage2.rr_cov, m_kf_stage2.rr_cov, 0,
                            m_kf_stage2.rr_cov, 2*m_kf_stage2.rr_cov, m_kf_stage2.rr_cov, 0,
                            m_kf_stage2.rr_cov, m_kf_stage2.rr_cov, 2*m_kf_stage2.rr_cov, 0,
                            0, 0, 0, m_kf_stage2.rz_cov};

            m_kf_stage3.C.fill(m_kf_stage3.ny, m_kf_stage3.nx, C_v);
            m_kf_stage3.R.fill(m_kf_stage3.ny, m_kf_stage3.ny, R_v);

            // Compute number of measurements. Plus 1 is for depth measurement.
            // This can be done this function is called only when atleast one measurement is available.

            uint8_t m12 = ((m_r1.isMeasure || m_r2.isMeasure) && (timeShiftCorrect(&m_r1, &m_r2)));
            uint8_t m23 = ((m_r2.isMeasure || m_r3.isMeasure) && (timeShiftCorrect(&m_r2, &m_r3)));
            uint8_t m31 = ((m_r3.isMeasure || m_r1.isMeasure) && (timeShiftCorrect(&m_r3, &m_r1)));

            m_kf_stage3.nyk = m12 + m23 + m31 + 1;
            debug(DTR("Number of measurements: %d"), m_kf_stage3.nyk);

            // Resize matrices
            m_kf_stage3.Ck.resizeAndFill(m_kf_stage3.nyk, m_kf_stage3.nx, 0.0);
            m_kf_stage3.Rk.resizeAndFill(m_kf_stage3.nyk, m_kf_stage3.nyk, m_kf_stage2.rr_cov);
            m_kf_stage3.yk.resizeAndFill(m_kf_stage3.nyk,1,0.0);
            m_kf_stage3.ykest.resizeAndFill(m_kf_stage3.nyk,1,0.0);
            m_kf_stage3.innov.resizeAndFill(m_kf_stage3.nyk, 1, 0.0);

            // Compute error from stage 3 and stage 3 estimate
            Matrix error_stage3 = m_kf_stage3.C*(m_kf_stage3.xHat - m_kf_stage2.xHat);
            debug(DTR("resizeAndFill successful"));
            uint8_t temp_ind = 0;
            // Construct Ck, Rk and yk
            if(m12 == 1)
            {
              double rangeDiff1 = m_args.speed_of_sound_in_water*(r1->unix_milisecond_time - r2->unix_milisecond_time)/1000;
              double rangeEst1 = d1 - d2;
              debug(DTR("R12: Range difference computed %f"), rangeDiff1);
              debug(DTR("R12: Extracted first row from C"));
              m_kf_stage3.Ck.put(temp_ind, 0, m_kf_stage3.C.row(0));
              debug(DTR("R12: Put first row into Ck"));
              m_kf_stage3.Rk(temp_ind, temp_ind) = 2*m_kf_stage2.rr_cov;
              debug(DTR("R12: Update Rk"));
              m_kf_stage3.yk(temp_ind, 0) = rangeDiff1;
              m_kf_stage3.ykest(temp_ind, 0) = rangeEst1 + error_stage3(0,0);
              temp_ind++;
              debug(DTR("Measurement update from R12"));
              m_r1.isMeasure = 0; m_r2.isMeasure = 0;
            }

            if(m23 == 1)
            {
              double rangeDiff2 = m_args.speed_of_sound_in_water*(r2->unix_milisecond_time - r3->unix_milisecond_time)/1000;
              double rangeEst2 = d2 - d3;
              debug(DTR("R23: Range difference computed %f"), rangeDiff2);
              m_kf_stage3.Ck.put(temp_ind, 0, m_kf_stage3.C.row(1));
              m_kf_stage3.Rk(temp_ind, temp_ind) = 2*m_kf_stage2.rr_cov;
              m_kf_stage3.yk(temp_ind, 0) = rangeDiff2;
              m_kf_stage3.ykest(temp_ind, 0) = rangeEst2 + error_stage3(1,0);
              temp_ind++;
              debug(DTR("Measurement update from R23"));
              m_r2.isMeasure = 0; m_r3.isMeasure = 0;
            }

            if(m31 == 1)
            {
              double rangeDiff3 = m_args.speed_of_sound_in_water*(r3->unix_milisecond_time - r1->unix_milisecond_time)/1000;
              double rangeEst3 = d3 - d1;
              debug(DTR("R31: Range difference computed %f"), rangeDiff3);
              m_kf_stage3.Ck.put(temp_ind, 0, m_kf_stage3.C.row(2));
              m_kf_stage3.Rk(temp_ind, temp_ind) = 2*m_kf_stage2.rr_cov;
              m_kf_stage3.yk(temp_ind, 0) = -rangeDiff3;
              m_kf_stage3.ykest(temp_ind, 0) = -rangeEst3 + error_stage3(2,0);
              temp_ind++;
              debug(DTR("Measurement update from R31"));
              m_r3.isMeasure = 0; m_r1.isMeasure = 0;
            }

            double depth = m12*(r1->sensor_data + r2->sensor_data)/2 +
                           m23*(r2->sensor_data + r3->sensor_data)/2 +
                           m31*(r3->sensor_data + r1->sensor_data)/2;
            debug(DTR("xkf_constructCandR_stage3:: Depth before %f"), depth);
            depth = depth/(m12 + m23 + m31);
            debug(DTR("xkf_constructCandR_stage3:: Depth after %f by %d"), depth, (m12 + m23 + m31));

            Matrix R_temp = 0.0*m_kf_stage3.Rk;
            m_kf_stage3.Rk.put(temp_ind, 0, R_temp.row(temp_ind));
            m_kf_stage3.Rk.put(0, temp_ind, R_temp.column(temp_ind));

            m_kf_stage3.Ck.put(temp_ind, 0, m_kf_stage3.C.row(3));
            m_kf_stage3.Rk(temp_ind, temp_ind) = m_kf_stage2.rz_cov;
            m_kf_stage3.yk(temp_ind, 0) = depth;
            m_kf_stage3.ykest(temp_ind, 0) = qD + error_stage3(3,0);
            debug(DTR("Update of C and R matrices successful"));
          } // End xkf_constructCandR_stage 3()

          void
          xkf_constructCandR_stage2(ReceiverData* r1, ReceiverData* r2, ReceiverData* r3)
          {
            //>> Compute C and R matrices
            double Czq_v[] = {-(r1->N - r3->N), -(r1->E - r3->E), 0.0,
                              -(r2->N - r3->N), -(r2->E - r3->E), 0.0,
                                             0,                0, 0.5};

            // TODO: Change the way covariance matrix is initialized
            double R_v[] = {0.1, 0, 0,
                            0, 0.1, 0,
                            0, 0, 0.2};

            m_kf_stage2.C.fill(m_kf_stage2.ny, m_kf_stage2.nx, Czq_v);
            m_kf_stage2.R.fill(m_kf_stage2.ny, m_kf_stage2.ny, R_v);

            // Compute number of measurements. Plus 1 is for depth measurement.
            // This can be done this function is called only when atleast one measurement is available.

            uint8_t m23 = ((m_r2.isMeasure || m_r3.isMeasure) && (timeShiftCorrect(&m_r2, &m_r3)));
            uint8_t m31 = ((m_r3.isMeasure || m_r1.isMeasure) && (timeShiftCorrect(&m_r3, &m_r1)));

            m_kf_stage2.nyk = m23 + m31 + 1;
            debug(DTR("Number of measurements: %d"), m_kf_stage2.nyk);

            // Resize matrices
            m_kf_stage2.Ck.resizeAndFill(m_kf_stage2.nyk, m_kf_stage2.nx, 0.0);
            m_kf_stage2.Rk.resizeAndFill(m_kf_stage2.nyk, m_kf_stage2.nyk, 0.0);
            m_kf_stage2.yk.resizeAndFill(m_kf_stage2.nyk,1,0.0);
            m_kf_stage2.ykest.resizeAndFill(m_kf_stage2.nyk,1,0.0);
            m_kf_stage2.innov.resizeAndFill(m_kf_stage2.nyk, 1, 0.0);

            debug(DTR("resizeAndFill successful"));
            uint8_t temp_ind = 0;
            Matrix Yest =  m_kf_stage2.C*m_kf_stage2.xHat;
						double d1_temp, d2_temp;
            // Construct Ck, Rk and yk
            if(m31 == 1)
            {
              double r13_temp = (r1->unix_milisecond_time - r3->unix_milisecond_time)/1000;
              double d1 = m_args.speed_of_sound_in_water*r13_temp;
              double z1 = 0.5*(d1*d1 - r1->N*r1->N - r1->E*r1->E + (r3->N)*(r3->N) + (r3->E)*(r3->E));
              debug(DTR("R13: Measurement computed %f"), z1);
              m_kf_stage2.Ck.put(temp_ind, 0, m_kf_stage2.C.row(0));
              double z1_cov = (1/2)*pow(2*m_kf_stage2.rr_cov,2) + (d1 + m_dr)*(d1 + m_dr)*2*m_kf_stage2.rr_cov;
              m_kf_stage2.Rk(temp_ind, temp_ind) = z1_cov;
              m_kf_stage2.yk(temp_ind, 0) = z1 + m_dr*d1;
              m_kf_stage2.ykest(temp_ind, 0) = Yest(temp_ind, 0);
              temp_ind++;
              debug(DTR("Measurement update from R13"));
              // m_r1.isMeasure = 0; m_r3.isMeasure = 0;
              d1_temp = d1;
            }

            if(m23 == 1)
            {
              double d2 = m_args.speed_of_sound_in_water*(r2->unix_milisecond_time - r3->unix_milisecond_time)/1000;
              double z2 = 0.5*(d2*d2 - r2->N*r2->N - r2->E*r2->E + (r3->N)*(r3->N) + (r3->E)*(r3->E));
              debug(DTR("R23: Measurement computed %f"), z2);
              m_kf_stage2.Ck.put(temp_ind, 0, m_kf_stage2.C.row(1));
              double z2_cov = (1/2)*pow(2*m_kf_stage2.rr_cov,2) + (d2 + m_dr)*(d2 + m_dr)*2*m_kf_stage2.rr_cov;
              m_kf_stage2.Rk(temp_ind, temp_ind) = z2_cov;
              m_kf_stage2.yk(temp_ind, 0) = z2 + m_dr*d2;
              m_kf_stage2.ykest(temp_ind, 0) = Yest(temp_ind, 0);
              temp_ind++;
              debug(DTR("Measurement update from R23"));
              // m_r2.isMeasure = 0; m_r3.isMeasure = 0;
              d2_temp = d2;
            }

            if(m_kf_stage2.nyk == 3)
            {
              m_kf_stage2.Rk(0,1) = 0.5*pow(m_kf_stage2.rr_cov,2) + (d1_temp*d2_temp + m_dr*(d1_temp + d2_temp) + m_dr*m_dr)*m_kf_stage2.rr_cov;
              m_kf_stage2.Rk(1,0) = 0.5*pow(m_kf_stage2.rr_cov,2) + (d1_temp*d2_temp + m_dr*(d1_temp + d2_temp) + m_dr*m_dr)*m_kf_stage2.rr_cov;
            }

            double depth = m23*(r2->sensor_data + r3->sensor_data)/2 + m31*(r3->sensor_data + r1->sensor_data)/2;
            debug(DTR("xkf_constructCandR_stage2:: Depth before %f"), depth);
            depth = depth/(m23 + m31);
            debug(DTR("xkf_constructCandR_stage2:: Depth after %f by %d"), depth, (m23 + m31));

            m_kf_stage2.Ck.put(temp_ind, 0, m_kf_stage2.C.row(2));
            m_kf_stage2.Rk(temp_ind, temp_ind) = m_kf_stage2.rz_cov;
            m_kf_stage2.yk(temp_ind, 0) = depth/2;
            m_kf_stage2.ykest(temp_ind, 0) = Yest(temp_ind, 0);
            debug(DTR("xkf_constructCandR_stage2:: Update of C and R matrices successful"));

          } // End xkf_constructCandR_stage2()

          void
          xkf_predict()
          {
            //>> Stage 2 LTV KF motion update
            //TODO: Manage situations when inputs are present
            m_kf_stage2.xHat = m_kf_stage2.A*m_kf_stage2.xHat;
            m_kf_stage2.PHat = m_kf_stage2.A*m_kf_stage2.PHat*transpose(m_kf_stage2.A) + m_kf_stage2.D*m_kf_stage2.Q*transpose(m_kf_stage2.D);

            m_kf_stage3.xHat = m_kf_stage3.A*m_kf_stage3.xHat;
            m_kf_stage3.PHat = m_kf_stage3.A*m_kf_stage3.PHat*transpose(m_kf_stage3.A) + m_kf_stage3.D*m_kf_stage3.Q*transpose(m_kf_stage3.D);

            // debug(DTR("Predicted state: %f %f %f"), m_kf_stage2.xHat.element(0,0), m_kf_stage2.xHat.element(1,0), m_kf_stage2.xHat.element(2,0));
            // debug(DTR("Predict successful"));
          }

          void
          xkf_update_stage3()
          {
            Matrix K = m_kf_stage3.PHat*transpose(m_kf_stage3.Ck)*inverse(m_kf_stage3.Ck*m_kf_stage3.PHat*transpose(m_kf_stage3.Ck) + m_kf_stage3.Rk);
            Matrix I(m_kf_stage3.nx); // create identity matrix

            m_kf_stage3.innov = m_kf_stage3.yk - m_kf_stage3.ykest;
            m_kf_stage3.xHat = m_kf_stage3.xHat + K*m_kf_stage3.innov;
            m_kf_stage3.PHat = (I - K*m_kf_stage3.Ck)*m_kf_stage3.PHat;
            // m_norm_innov = m_kf_stage3.innov.norm_p(2);

            debug(DTR("Update successful"));

            // Compute fish position of stage 3 for logging
            m_fp_stage3.n = m_kf_stage3.xHat(0,0);
            m_fp_stage3.e = m_kf_stage3.xHat(1,0);
            m_fp_stage3.d = m_kf_stage3.xHat(2,0);

            m_fp_stage3.lat = m_refCoord.lat; m_fp_stage3.lon = m_refCoord.lon; m_fp_stage3.hae = m_refCoord.hae;
            WGS84::displace(m_fp_stage3.n, m_fp_stage3.e, m_fp_stage3.d, &(m_fp_stage3.lat), &(m_fp_stage3.lon), &(m_fp_stage3.hae));
          } // End of xkf_update_stage3 function

          void
          xkf_update_stage2()
          {
            Matrix K = m_kf_stage2.PHat*transpose(m_kf_stage2.Ck)*inverse(m_kf_stage2.Ck*m_kf_stage2.PHat*transpose(m_kf_stage2.Ck) + m_kf_stage2.Rk);
            Matrix I(m_kf_stage2.nx); // create identity matrix

            m_kf_stage2.innov = m_kf_stage2.yk - m_kf_stage2.ykest;
            m_kf_stage2.xHat = m_kf_stage2.xHat + K*m_kf_stage2.innov;
            m_kf_stage2.PHat = (I - K*m_kf_stage2.Ck)*m_kf_stage2.PHat;
            // m_norm_innov = m_kf_stage2.innov.norm_p(2);
            //
            // // Make covariance matrix symmetric
            // // m_kf_stage2.PHat = 0.5*(m_kf_stage2.PHat + transpose(m_kf_stage2.PHat));
            // debug(DTR("Estimated state: %f %f %f"), m_kf_stage2.xHat.element(0,0), m_kf_stage2.xHat.element(1,0), m_kf_stage2.xHat.element(2,0));
            // // debug(DTR("Innovation: %f %f %f %f"), m_kf_stage2.innov(0,0), m_kf_stage2.innov(1,0), m_kf_stage2.innov(2,0), m_kf_stage2.innov(3,0));
            debug(DTR("Update successful"));

            // Compute fish position of stage 2 for logging
            m_fp_stage2.n = m_kf_stage2.xHat(0,0);
            m_fp_stage2.e = m_kf_stage2.xHat(1,0);
            m_fp_stage2.d = m_kf_stage2.xHat(2,0);

            m_fp_stage2.lat = m_refCoord.lat; m_fp_stage2.lon = m_refCoord.lon; m_fp_stage2.hae = m_refCoord.hae;
            WGS84::displace(m_fp_stage2.n, m_fp_stage2.e, m_fp_stage2.d, &(m_fp_stage2.lat), &(m_fp_stage2.lon), &(m_fp_stage2.hae));
          } // End of xkf_update_stage2 function

          void logFishPosition()
          {
            double time_log = Clock::getSinceEpochMsec();
            //>> Log Stage 1: LS Estimate
            m_o1.precision(15);
            m_o1 << m_r1.unix_milisecond_time << "," << m_r2.unix_milisecond_time << "," << m_r3.unix_milisecond_time << "," << DUNE::Math::Angles::degrees(m_fp_stage1.lat) << "," << DUNE::Math::Angles::degrees(m_fp_stage1.lon) << "," << m_fp_stage1.d << "," << m_dr << "," << time_log << std::endl;

            //>> Log Stage 2: LTV KF Estimate
            m_o2.precision(15);
            m_o2 << m_r1.unix_milisecond_time << "," << m_r2.unix_milisecond_time << "," << m_r3.unix_milisecond_time << "," << DUNE::Math::Angles::degrees(m_fp_stage2.lat) << "," << DUNE::Math::Angles::degrees(m_fp_stage2.lon) << "," << m_fp_stage2.d << "," << m_kf_stage2.innov.norm_p(2) << "," << m_kf_stage2.PHat.norm_p(2) << "," << m_kf_stage2.PHat.trace() << "," << time_log << std::endl;

            //>> Log Stage 3: Linearized KF Estimate
            m_o3.precision(15);
            m_o3 << m_r1.unix_milisecond_time << "," << m_r2.unix_milisecond_time << "," << m_r3.unix_milisecond_time << "," << DUNE::Math::Angles::degrees(m_fp_stage3.lat) << "," << DUNE::Math::Angles::degrees(m_fp_stage3.lon) << "," << m_fp_stage3.d << "," << m_kf_stage3.innov.norm_p(2) << "," << m_kf_stage3.PHat.norm_p(2) << "," << m_kf_stage3.PHat.trace() << "," << time_log << std::endl;
            // //>> dispatch fish position
            // IMC::RemoteSensorInfo tagPosition;
            // Matrix FishPosition = m_kf_stage2.xHat;
            //
            // // double x1_lat = r1.lat, x1_lon = r1.lon, x1_depth;
            // double x1_lat = m_refCoord.lat; double x1_lon = m_refCoord.lon; double x1_depth = m_refCoord.hae;
            // WGS84::displace(FishPosition(0,0), FishPosition(1,0), FishPosition(2,0), &(x1_lat), &(x1_lon), &(x1_depth));
            //
            // tagPosition.id = "Fish_position_est";
            // tagPosition.lat = x1_lat;
            // tagPosition.lon = x1_lon;
            // tagPosition.alt = -x1_depth;
            //
            // tagPosition.data = r1.msg+ "#" +r2.msg+ "#" + r3.msg;
            // dispatch(tagPosition);
            // debug(DTR("Fish lat: %f long: %f Depth: %f"), x1_lat, x1_lon, -x1_depth);
            //
            // m_norm_P = m_kf_stage2.PHat.norm_p(2);
            // m_o.precision(15);
            //
            // m_o << m_r1.unix_milisecond_time << "," << m_r2.unix_milisecond_time << "," << m_r3.unix_milisecond_time << "," << x1_lat << "," << x1_lon << "," << "," << -x1_depth << "," << FishPosition(2,0) << "," << m_d1 << "," << m_d2 << "," << m_d3 << "," << m_norm_innov << "," << m_norm_P << "," << m_kf_stage2.PHat.trace() << std::endl;
          }

          void
          findFishPosition(ReceiverData r1, ReceiverData r2, ReceiverData r3)
          {
            if(!m_args.isXKFInitialized)
            {
              // Initialize Kalman Filter
              initializeXKF();
              m_args.isXKFInitialized = 1;
            }

            //>> Step 1: Predict
            xkf_predict();

            //>> Step 2: If measurements available

            uint8_t m12 = ((m_r1.isMeasure || m_r2.isMeasure) && (timeShiftCorrect(&m_r1, &m_r2)));
            uint8_t m23 = ((m_r2.isMeasure || m_r3.isMeasure) && (timeShiftCorrect(&m_r2, &m_r3)));
            uint8_t m31 = ((m_r3.isMeasure || m_r1.isMeasure) && (timeShiftCorrect(&m_r3, &m_r1)));

            if((m12 == 1) || (m23 == 1) || (m31 == 1))
            {
              //>> Step 2.1: Update System Matrices and Covariance
              createNEDframe(&r1, &r2, &r3);
              // if((m12 == 1) && (m23 == 1) && (m31 == 1))
              if((m23 == 1) && (m31 == 1))
              {
                bool isSuccess = LeastSquaresEstimate(&r1, &r2, &r3);
                if(isSuccess && ~m_kf_stage2.isEKFActive)
                {
                    m_kf_stage3.isEKFActive = 1;
                    m_kf_stage2.isEKFActive = 1;
                }
              }
              //>> Step 2.2: Measurement Update - stage - 2
              // if(m_kf_stage2.isEKFActive && ((m23 == 1) || (m31 == 1))) // conditions for stage 2 update
              if(m_kf_stage2.isEKFActive && ((m23 == 1) && (m31 == 1))) // conditions for stage 2 update
              {
                xkf_constructCandR_stage2(&r1, &r2, &r3);
                xkf_update_stage2();
              }
              if(m_kf_stage3.isEKFActive)
              {
                xkf_constructCandR_stage3(&r1, &r2, &r3);
                xkf_update_stage3();
              }
            }
            //m_o3 << m_r1.unix_milisecond_time << "," << m_r2.unix_milisecond_time << "," << m_r3.unix_milisecond_time << "," << DUNE::Math::Angles::degrees(m_fp_stage3.lat) << "," << DUNE::Math::Angles::degrees(m_fp_stage3.lon) << "," << m_fp_stage3.d << "," << m_kf_stage3.innov.norm_p(2) << "," << m_kf_stage3.PHat.norm_p(2) << "," << m_kf_stage3.PHat.trace() << "," << time_log << std::endl;
              IMC::RemoteSensorInfo tagPosition;
              tagPosition.lat = m_fp_stage3.lat;
              tagPosition.lon = m_fp_stage3.lon;
              tagPosition.alt = -m_fp_stage3.d;
              tagPosition.data = std::to_string(m_r1.unix_milisecond_time) + "," + std::to_string(m_r2.unix_milisecond_time) + "," + std::to_string(m_r3.unix_milisecond_time);// + "," + m_kf_stage3.innov.norm_p(2) + "," << m_kf_stage3.PHat.norm_p(2) + "," << m_kf_stage3.PHat.trace();
              dispatch(tagPosition);
            logFishPosition();
          } //End function findFishPosition

          //! Main loop.
          void
          task(void)
          {
            if(m_args.isXKFActive)
            {
              findFishPosition(m_r1, m_r2, m_r3);
            } // End of if
          } // End of task function
        };
      }
    }
  }
}

DUNE_TASK
