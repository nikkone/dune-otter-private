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
      namespace FishTagEKFTDoA
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

          //! Structure to hold the EKF data
          EKFilter m_ekf;

          //! Output file only for use in simulations
          std::ofstream m_o;

          //! Reference coordinate system
          RefCoordinate m_refCoord;

          //! Additional variables used to log the EKF data
          double m_d1, m_d2, m_d3, m_norm_innov, m_norm_P, m_trace_P;

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

            param("num_states", m_ekf.nx)
            .description("Number of States")
            .defaultValue("4");

            param("num_outputs", m_ekf.ny)
            .description("Number of Outputs")
            .defaultValue("4");

            param("num_inputs", m_ekf.nu)
            .description("Number of Inputs")
            .defaultValue("0");

            param("FishPos Cov", m_ekf.qq_cov)
            .description("Fish Position Covariance")
            .defaultValue("0");

            param("ToA Cov", m_ekf.rr_cov)
            .description("Time of Arrival Covariance")
            .defaultValue("0");

            param("Depth Cov", m_ekf.rz_cov)
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
            bind<IMC::EstimatedState>(this);
          }

          //! Update internal state with new parameter values.
          void
          onUpdateParameters(void)
          {
            m_args.dt = 1/getFrequency();
            //
            m_refCoord.lat = Math::Angles::radians(m_args.reference[0]);
            m_refCoord.lon = Math::Angles::radians(m_args.reference[1]);
            m_refCoord.hae = 0.0;

            // m_refCoord.lat = m_args.reference[0];
            // m_refCoord.lon = m_args.reference[1];
            // m_refCoord.hae = 0.0;

            debug(DTR("Reference Cooridnate: %f %f"), m_args.reference[0], m_args.reference[1]);
            debug(DTR("Fish Position: %f %f %f"), m_args.position[0], m_args.position[1], m_args.position[2]);

            // m_r1_name = resolveSystemName(m_args.rec_1_vehicle_name);
            // m_r2_name = resolveSystemName(m_args.rec_2_vehicle_name);
            // m_r3_name = resolveSystemName(m_args.rec_3_vehicle_name);
            // inf(DTR("Carrier of receiver 1 name is %s, with ID %d"), m_args.rec_1_vehicle_name.c_str(), m_r1_name);
            // inf(DTR("Carrier of receiver 2 name is %s, with ID %d"), m_args.rec_2_vehicle_name.c_str(), m_r2_name);
            // inf(DTR("Carrier of receiver 3 name is %s, with ID %d"), m_args.rec_3_vehicle_name.c_str(), m_r3_name);

            //>> Used for simulation purposes only
            // m_r1_name = 42;
            // m_r2_name = 22;
            // m_r3_name = 45;

            m_r1_name = 632; // Duckling 1
            m_r2_name = 634; // Duckling 2
            m_r3_name = 631; // Duckling 3

            // m_refCoord.lat = 1.11039166;
            // m_refCoord.lon = 0.16761261;
            // m_refCoord.hae = 0.0;
            m_o.open("/home/praveen/Desktop/ExperimentDataSet/Set3/EKFTDoAoutput.txt", std::ios::out | std::ios::trunc);
            if(m_o.is_open())
            {
              debug(DTR("EKFTDoAoutput open!"));
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

          // //Function: Parse the Fish tag data
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

            // r1->N = 0.0;
            // r1->E = 0.0;
            // r1->D = 0.0;
            //
            // trace(DTR("createNEDframe: Radians R1: %f, %f,"), m_r1.lat, m_r1.lon);
            // trace(DTR("createNEDframe: Radians R2: %f, %f,"), m_r2.lat, m_r2.lon);
            // trace(DTR("createNEDframe: Radians R3: %f, %f,"), m_r3.lat, m_r3.lon);
            //
            // WGS84::displacement(r1->lat, r1->lon, 0.0, r2->lat, r2->lon, 0.0, &(r2->N), &(r2->E), &(r2->D));
            // WGS84::displacement(r1->lat, r1->lon, 0.0, r3->lat, r3->lon, 0.0, &(r3->N), &(r3->E), &(r3->D));
            //
            // Matrix FishPosition = m_ekf.xHat;
            //
            // double x1_lat = r1->lat, x1_lon = r1->lon, x1_depth = 0.0;
            // WGS84::displace(FishPosition(0,0), FishPosition(1,0), 0.0, &(x1_lat), &(x1_lon), &(x1_depth));
            // double fn,fe,fd;
            // WGS84::displacement(r1->lat, r1->lon, 0.0, x1_lat, x1_lon, 0.0, &(fn), &(fe), &(fd));
            //
            // m_ekf.xHat(0,0) = fn; m_ekf.xHat(1,0) = fe;

            // r2->unix_milisecond_time = (r2->unix_milisecond_time - r1->unix_milisecond_time)/1000.0 + 1.0;
            // r3->unix_milisecond_time = (r3->unix_milisecond_time - r1->unix_milisecond_time)/1000.0 + 1.0;
            // r1->unix_milisecond_time = (r1->unix_milisecond_time - r1->unix_milisecond_time)/1000.0 + 1.0;

            // trace(DTR("createNEDframe: NED R1: %f, %f, %f; t=%f"), r1->N, r1->E, r1->D, r1->unix_milisecond_time);
            // trace(DTR("createNEDframe: NED R2: %f, %f, %f; t=%f"), r2->N, r2->E, r2->D, r2->unix_milisecond_time);
            // trace(DTR("createNEDframe: NED R3: %f, %f, %f; t=%f"), r3->N, r3->E, r3->D, r3->unix_milisecond_time);
          }

          void
          consume(const IMC::EstimatedState* msg)
          {
            inf("Estimated arrived from %d", msg->getSource());
          }

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
                  m_ekf.isEKFActive = 1;
                }
                else if (m_r2_name == msg->getSource())
                {
                  debug(DTR("consume: Message from R2 arrived"));
                  m_r2 = tmp;
                  m_r2.isMeasure = 1;
                  tagPosition.id = "Receiver_2";
                  dispatch(tagPosition);
                  m_ekf.isEKFActive = 1;
                }
                else if (m_r3_name == msg->getSource())
                {
                  debug(DTR("consume: Message from R3 arrived"));
                  m_r3 = tmp;
                  m_r3.isMeasure = 1;
                  tagPosition.id = "Receiver_3";
                  dispatch(tagPosition);
                  m_ekf.isEKFActive = 1;
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
          initializeEKF(void)
          {
            // Set the size of Kalman filter matrices
            m_ekf.nyk = m_ekf.ny; //TODO: Do something with this
            m_ekf.nd = 3; //TODO: remove this hardcoding

            m_ekf.A.resizeAndFill(m_ekf.nx, m_ekf.nx,0.0);
            m_ekf.C.resizeAndFill(m_ekf.ny, m_ekf.nx,0.0);
            m_ekf.Ck.resizeAndFill(m_ekf.nyk, m_ekf.nx,0.0);
            m_ekf.PHat.resizeAndFill(m_ekf.nx, m_ekf.nx,0.0);
            m_ekf.Q.resizeAndFill(m_ekf.nx, m_ekf.nx,0.0);
            m_ekf.R.resizeAndFill(m_ekf.ny, m_ekf.ny,0.0);
            m_ekf.Rk.resizeAndFill(m_ekf.nyk, m_ekf.nyk,0.0);
            m_ekf.xHat.resizeAndFill(m_ekf.nx,1,0.0);
            m_ekf.innov.resizeAndFill(m_ekf.nyk,1,0.0);
            m_ekf.yk.resizeAndFill(m_ekf.nyk,1,0.0);
            m_ekf.ykest.resizeAndFill(m_ekf.nyk,1,0.0);

            if(m_ekf.nu != 0)
            {
              m_ekf.B.resizeAndFill(m_ekf.nx, m_ekf.nu, 0.0);
            }
            if(m_ekf.nd != 0)
            {
              m_ekf.D.resizeAndFill(m_ekf.nx, m_ekf.nd, 0.0);
            }
            // Set the state transition matrix (does not change)
            double A_v[] = {1, 0, 0,
                            0, 1, 0,
                            0, 0, 1};

            double Q_v[] = {m_ekf.qq_cov, 0, 0,
                            0, m_ekf.qq_cov, 0,
                            0, 0, m_ekf.qq_cov/10.0};

            // Initial position of Fish tag in [m]
            // double qn = 150.00;
            // double qe = 150.00;
            // double qd = 10.0;
            //
            // double flat = m_refCoord.lat; double flon = m_refCoord.lon; double fhae = m_refCoord.hae;
            // WGS84::displace(qn, qe, qd, &(flat), &(flon), &(fhae));
            // debug(DTR("onUpdateParameters:: FT lat: %f, lon: %f, hae: %f"), flat, flon, fhae);
            //
            // double fish_lat = 1.11041514936877;
            // double fish_lon = 0.167665401259442;
            // debug(DTR("qn: %f qe: %f fish_lat: %f fish_lon: %f"),qn, qe, fish_lat, fish_lon);
            //
            // WGS84::displacement(m_refCoord.lat, m_refCoord.lon, 0.0, fish_lat, fish_lon, 0.0, &(qn), &(qe), &(qd));
            // debug(DTR("Intial Fish position: %f, %f, %f"), qn, qe, qd);
            // double x0_v[] = {100, 100, 7.0};
            // double x0_v[] = {10, 10, 1.0};

            double x0_v[] = {m_args.position[0], m_args.position[1], m_args.position[2]};
            double P0_v[] = {100, 0, 0,
                             0, 100, 0,
                             0, 0, 10};

            double D_v[] = {m_args.dt*1, 0, 0,
                            0, m_args.dt*1, 0,
                            0, 0, m_args.dt*1};

            m_ekf.A.fill(m_ekf.nx, m_ekf.nx, A_v);
            m_ekf.PHat.fill(m_ekf.nx, m_ekf.nx, P0_v);
            m_ekf.Q.fill(m_ekf.nx, m_ekf.nx, Q_v);
            m_ekf.xHat.fill(m_ekf.nx, 1, x0_v);
            m_ekf.D.fill(m_ekf.nx, m_ekf.nd, D_v);
            // m_ekf.B.fill(m_ekf.nx, m_ekf.nu, B_v);
          }

          void
          ekf_constructCandR(ReceiverData* r1, ReceiverData* r2, ReceiverData* r3)
          {
            //>> Get the state estimate
            Matrix x_est = m_ekf.xHat;

            //>> Compute C and R matrices
            double qN = x_est.element(0,0);
            double qE = x_est.element(1,0);
            double qD = x_est.element(2,0);

            double d1 = std::sqrt(std::pow((r1->N - qN),2) + std::pow((r1->E - qE),2) + std::pow((r1->D - qD),2));
            double d2 = std::sqrt(std::pow((r2->N - qN),2) + std::pow((r2->E - qE),2) + std::pow((r2->D - qD),2));
            double d3 = std::sqrt(std::pow((r3->N - qN),2) + std::pow((r3->E - qE),2) + std::pow((r3->D - qD),2));

            m_d1 = d1;
            m_d2 = d2;
            m_d3 = d3;

            debug(DTR("ekf_constructCandR: range d1: %f, d2: %f, d3: %f"), d1, d2, d3);
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
                            c31, c32, c33,
                            c41, c42, c43};

            double R_v[] = {2*m_ekf.rr_cov, 0, 0, 0,
                            0, 2*m_ekf.rr_cov, 0, 0,
                            0, 0, 2*m_ekf.rr_cov, 0,
                            0, 0, 0, m_ekf.rz_cov};

            m_ekf.C.fill(m_ekf.ny, m_ekf.nx, C_v);
            m_ekf.R.fill(m_ekf.ny, m_ekf.ny, R_v);

            // Compute number of measurements. Plus 1 is for depth measurement.
            // This can be done this function is called only when atleast one measurement is available.

            uint8_t m12 = ((m_r1.isMeasure || m_r2.isMeasure) && (timeShiftCorrect(&m_r1, &m_r2)));
            uint8_t m23 = ((m_r2.isMeasure || m_r3.isMeasure) && (timeShiftCorrect(&m_r2, &m_r3)));
            uint8_t m31 = ((m_r3.isMeasure || m_r1.isMeasure) && (timeShiftCorrect(&m_r3, &m_r1)));

            m_ekf.nyk = m12 + m23 + m31 + 1;
            debug(DTR("Number of measurements: %d"), m_ekf.nyk);

            // Resize matrices
            m_ekf.Ck.resizeAndFill(m_ekf.nyk, m_ekf.nx, 0.0);
            m_ekf.Rk.resizeAndFill(m_ekf.nyk, m_ekf.nyk, m_ekf.rr_cov);
            m_ekf.yk.resizeAndFill(m_ekf.nyk,1,0.0);
            m_ekf.ykest.resizeAndFill(m_ekf.nyk,1,0.0);
            m_ekf.innov.resizeAndFill(m_ekf.nyk, 1, 0.0);

            debug(DTR("resizeAndFill successful"));
            uint8_t temp_ind = 0;
            // Construct Ck, Rk and yk
            if(m12 == 1)
            {
              double rangeDiff1 = m_args.speed_of_sound_in_water*(r1->unix_milisecond_time - r2->unix_milisecond_time)/1000;
              double rangeEst1 = d1 - d2;
              debug(DTR("R12: Range difference computed %f"), rangeDiff1);
              Matrix temp_C1 = m_ekf.C.row(0);
              debug(DTR("R12: Extracted first row from C"));
              m_ekf.Ck.put(temp_ind, 0, temp_C1);
              debug(DTR("R12: Put first row into Ck"));
              m_ekf.Rk(temp_ind, temp_ind) = 2*m_ekf.rr_cov;
              debug(DTR("R12: Update Rk"));
              m_ekf.yk(temp_ind, 0) = rangeDiff1;
              m_ekf.ykest(temp_ind, 0) = rangeEst1;
              temp_ind++;
              debug(DTR("Measurement update from R12"));
              m_r1.isMeasure = 0; m_r2.isMeasure = 0;
            }

            if(m23 == 1)
            {
              double rangeDiff2 = m_args.speed_of_sound_in_water*(r2->unix_milisecond_time - r3->unix_milisecond_time)/1000;
              double rangeEst2 = d2 - d3;
              debug(DTR("R23: Range difference computed %f"), rangeDiff2);
              m_ekf.Ck.put(temp_ind, 0, m_ekf.C.row(1));
              m_ekf.Rk(temp_ind, temp_ind) = 2*m_ekf.rr_cov;
              m_ekf.yk(temp_ind, 0) = rangeDiff2;
              m_ekf.ykest(temp_ind, 0) = rangeEst2;
              temp_ind++;
              debug(DTR("Measurement update from R23"));
              m_r2.isMeasure = 0; m_r3.isMeasure = 0;
            }

            if(m31 == 1)
            {
              double rangeDiff3 = m_args.speed_of_sound_in_water*(r3->unix_milisecond_time - r1->unix_milisecond_time)/1000;
              double rangeEst3 = d3 - d1;
              debug(DTR("R31: Range difference computed %f"), rangeDiff3);
              m_ekf.Ck.put(temp_ind, 0, m_ekf.C.row(2));
              m_ekf.Rk(temp_ind, temp_ind) = 2*m_ekf.rr_cov;
              m_ekf.yk(temp_ind, 0) = rangeDiff3;
              m_ekf.ykest(temp_ind, 0) = rangeEst3;
              temp_ind++;
              debug(DTR("Measurement update from R31"));
              m_r3.isMeasure = 0; m_r1.isMeasure = 0;
            }

            double depth = m12*(r1->sensor_data + r2->sensor_data)/2 +
                           m23*(r2->sensor_data + r3->sensor_data)/2 +
                           m31*(r3->sensor_data + r1->sensor_data)/2;
            debug(DTR("ekf_constructCandR:: Depth before %f"), depth);
            depth = depth/(m12 + m23 + m31);
            debug(DTR("ekf_constructCandR:: Depth after %f by %d"), depth, (m12 + m23 + m31));

            Matrix R_temp = 0.0*m_ekf.Rk;
            m_ekf.Rk.put(temp_ind, 0, R_temp.row(temp_ind));
            m_ekf.Rk.put(0, temp_ind, R_temp.column(temp_ind));

            m_ekf.Ck.put(temp_ind, 0, m_ekf.C.row(3));
            m_ekf.Rk(temp_ind, temp_ind) = m_ekf.rz_cov;
            m_ekf.yk(temp_ind, 0) = depth;
            m_ekf.ykest(temp_ind, 0) = qD;
            debug(DTR("Update of C and R matrices successful"));
          } // End ekf_constructCandR()

          void
          ekf_predict()
          {
            //TODO: Manage situations when inputs are present
            m_ekf.xHat = m_ekf.A*m_ekf.xHat;
            m_ekf.PHat = m_ekf.A*m_ekf.PHat*transpose(m_ekf.A) + m_ekf.D*m_ekf.Q*transpose(m_ekf.D);

            // debug(DTR("Predicted state: %f %f %f"), m_ekf.xHat.element(0,0), m_ekf.xHat.element(1,0), m_ekf.xHat.element(2,0));
            // debug(DTR("Predict successful"));
          }

          void
          ekf_update()
          {
            Matrix K = m_ekf.PHat*transpose(m_ekf.Ck)*inverse(m_ekf.Ck*m_ekf.PHat*transpose(m_ekf.Ck) + m_ekf.Rk);
            Matrix I(m_ekf.nx); // create identity matrix

            m_ekf.innov = m_ekf.yk - m_ekf.ykest;
            m_ekf.xHat = m_ekf.xHat + K*m_ekf.innov;
            m_ekf.PHat = (I - K*m_ekf.Ck)*m_ekf.PHat;
            m_norm_innov = m_ekf.innov.norm_p(2);

            // Make covariance matrix symmetric
            // m_ekf.PHat = 0.5*(m_ekf.PHat + transpose(m_ekf.PHat));
            debug(DTR("Estimated state: %f %f %f"), m_ekf.xHat.element(0,0), m_ekf.xHat.element(1,0), m_ekf.xHat.element(2,0));
            // debug(DTR("Innovation: %f %f %f %f"), m_ekf.innov(0,0), m_ekf.innov(1,0), m_ekf.innov(2,0), m_ekf.innov(3,0));
            debug(DTR("Update successful"));
          }

          void
          findFishPosition(ReceiverData r1, ReceiverData r2, ReceiverData r3)
          {
            if(!m_ekf.isEKFInitialized)
            {
              // Initialize Kalman Filter
              initializeEKF();
              m_ekf.isEKFInitialized = 1;
            }

            //>> Step 1: Predict
            ekf_predict();

            //>> Step 2: If measurements available

            uint8_t m12 = ((m_r1.isMeasure || m_r2.isMeasure) && (timeShiftCorrect(&m_r1, &m_r2)));
            uint8_t m23 = ((m_r2.isMeasure || m_r3.isMeasure) && (timeShiftCorrect(&m_r2, &m_r3)));
            uint8_t m31 = ((m_r3.isMeasure || m_r1.isMeasure) && (timeShiftCorrect(&m_r3, &m_r1)));

            if((m12 == 1) || (m23 == 1) || (m31 == 1))
            {
              //>> Step 2.1: Update System Matrices and Covariance
              createNEDframe(&r1, &r2, &r3);
              ekf_constructCandR(&r1, &r2, &r3);
              //>> Step 2.2: Measurement Update
              ekf_update();
            }

            //>> dispatch fish position
            IMC::RemoteSensorInfo tagPosition;
            Matrix FishPosition = m_ekf.xHat;

            // double x1_lat = r1.lat, x1_lon = r1.lon, x1_depth;
            double x1_lat = m_refCoord.lat; double x1_lon = m_refCoord.lon; double x1_depth = m_refCoord.hae;
            WGS84::displace(FishPosition(0,0), FishPosition(1,0), FishPosition(2,0), &(x1_lat), &(x1_lon), &(x1_depth));

            tagPosition.id = "Fish_position_est";
            tagPosition.lat = x1_lat;
            tagPosition.lon = x1_lon;
            tagPosition.alt = -x1_depth;

            tagPosition.data = r1.msg+ "#" +r2.msg+ "#" + r3.msg;
            dispatch(tagPosition);
            debug(DTR("Fish lat: %f long: %f Depth: %f"), x1_lat, x1_lon, -x1_depth);
            double time_log = Clock::getSinceEpochMsec();
            m_norm_P = m_ekf.PHat.norm_p(2);
            m_o.precision(15);

            m_o << m_r1.unix_milisecond_time << "," << m_r2.unix_milisecond_time << "," << m_r3.unix_milisecond_time << "," << x1_lat << "," << x1_lon << "," << "," << -x1_depth << "," << FishPosition(2,0) << "," << m_d1 << "," << m_d2 << "," << m_d3 << "," << m_norm_innov << "," << m_norm_P << "," << m_ekf.PHat.trace() << "," << time_log << std::endl;
          } //End function findFishPosition

          //! Main loop.
          void
          task(void)
          {
            if(m_ekf.isEKFActive)
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
