//***************************************************************************
// Copyright 2007-2017 Universidade do Porto - Faculdade de Engenharia      *
// Laboratório de Sistemas e Tecnologia Subaquática (LSTS)                  *
//***************************************************************************
// This file is part of DUNE: Unified Navigation Environment.               *
//                                                                          *
// Commercial Licence Usage                                                 *
// Licencees holding valid commercial DUNE licences may use this file in    *
// accordance with the commercial licence agreement provided with the       *
// Software or, alternatively, in accordance with the terms contained in a  *
// written agreement between you and Faculdade de Engenharia da             *
// Universidade do Porto. For licensing terms, conditions, and further      *
// information contact lsts@fe.up.pt.                                       *
//                                                                          *
// Modified European Union Public Licence - EUPL v.1.1 Usage                *
// Alternatively, this file may be used under the terms of the Modified     *
// EUPL, Version 1.1 only (the "Licence"), appearing in the file LICENCE.md *
// included in the packaging of this file. You may not use this work        *
// except in compliance with the Licence. Unless required by applicable     *
// law or agreed to in writing, software distributed under the Licence is   *
// distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF     *
// ANY KIND, either express or implied. See the Licence for the specific    *
// language governing permissions and limitations at                        *
// https://github.com/LSTS/dune/blob/master/LICENCE.md and                  *
// http://ec.europa.eu/idabc/eupl.html.                                     *
//***************************************************************************
// Author: Praveen Jain                                                     *
//***************************************************************************

// DUNE headers.
#include <DUNE/DUNE.hpp>
#include <fstream>

namespace
{
  namespace Navigation
  {
    namespace FishEstimation
    {
      namespace TBRTest
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
          double sensor_data; // CHANGED: Was int, so results in paper is probably worse than needs be.
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

        struct Task: public DUNE::Tasks::Periodic
        {
          double m_time_ms;
          bool isInit = 0;

          std::ifstream m_v1;
          std::ifstream m_v2;
          std::ifstream m_v3;

          std::ofstream m_o;

          ReceiverData m_r1, m_r2, m_r3;
          IMC::TBRFishTag tag1, tag2, tag3;
          //! Constructor.
          //! @param[in] name task name.
          //! @param[in] ctx context.
          Task(const std::string& name, Tasks::Context& ctx):
          DUNE::Tasks::Periodic(name, ctx)
          {

          }

          //! Update internal state with new parameter values.
          void
          onUpdateParameters(void)
          {
            inf("onUpdateParameters");
            //>>--------------- When using three raw data files
            m_v1.open("/home/nikolai/ExperimentDataSet/Set4/Receiver_1.txt", std::ios::in);
            if(m_v1.is_open())
            {
              debug(DTR("Receiver-1 open!"));
            }

            m_v2.open( "/home/nikolai/ExperimentDataSet/Set4/Receiver_2.txt", std::ios::in);
            if(m_v2.is_open())
            {
              debug(DTR("Receiver-2 open!"));
            }

            m_v3.open("/home/nikolai/ExperimentDataSet/Set4/Receiver_3.txt", std::ios::in);
            if(m_v3.is_open())
            {
              debug(DTR("Receiver-3 open!"));
            }

            //>>--------------- When using combined good data files

            // m_v1.open("/home/praveen/Desktop/TBR/TBRData.txt", std::ios::in);
            // if(m_v1.is_open())
            // {
            //   debug(DTR("TBRData open!"));
            // }

            /*m_o.open("/home/nikolai/ExperimentDataSet/TBRoutput.txt", std::ios::out | std::ios::trunc);
            if(m_o.is_open())
            {
              debug(DTR("TBRoutput open!"));
            }*/
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

          IMC::TBRFishTag toIMCTag(ReceiverData input) {
            IMC::TBRFishTag output;
            output.lat = input.lat;
            output.lon = input.lon;
            output.serial_no = input.serial_no;
            output.unix_timestamp = static_cast<uint32_t>(input.unix_timestamp);
            output.millis = static_cast<uint16_t>(input.milisecond_timestamp);
            //inf("Timecompare: %f, %f",input.unix_milisecond_time, static_cast<double>((double)output.unix_timestamp*1000+output.millis));
            output.trans_id = input.tag_id;
            output.trans_data = uint16_t(input.sensor_data);
            //inf("Depthcompare: %f, %d",input.sensor_data, output.trans_data);
            //output.trans_freq = input.
            //output.trans_protocol = input.code_type;
            output.snr = input.snr;
            output.recv_mem_addr = input.memory_number;
            return output;
          }

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
              data.sensor_data = std::atof(data_fields[6].c_str()); // assumed depth
              debug(DTR("Depth %f"), data.sensor_data);

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

            /*IMC::TBRFishTag output;
            output.lat = data.lat;
            output.lon = data.lon;
            output.serial_no = data.serial_no;
            output.unix_timestamp = data.unix_timestamp;
            output.millis = data.milisecond_timestamp;
            inf("Timecompare: %f, %f",data.unix_milisecond_time, static_cast<double>((double)data.unix_timestamp*1000+data.milisecond_timestamp));
            output.trans_id = data.tag_id;
            output.trans_data = uint16_t(data.sensor_data/0.392);
            inf("Depthcompare: %f, %d",data.sensor_data, output.trans_data);
            //output.trans_freq = input.
            //output.trans_protocol = input.code_type;
            output.snr = data.snr;
            output.recv_mem_addr = data.memory_number;*/
            }
            else
            {
              inf(DTR("Wrong number of fields in the message, received %d"), (int)data_fields.size());
            }

            return data;
          } // End Function: parseFishTagData



      //! Main loop.
      void
      task(void)
      {
        //>>--------------- When using three raw data files

        std::string temp_msg1, temp_msg2, temp_msg3;
        if(!isInit)
        {
          //>> When three separate raw data files
          do{
            std::getline(m_v1, temp_msg1);
            m_r1 = parseFishTagData(temp_msg1);
            
            // debug(DTR("Receiver_1 message %s"), temp_msg1.c_str());
            temp_msg1.clear();
          }while(!m_r1.valid);

          do{
            std::getline(m_v2, temp_msg2);
            m_r2 = parseFishTagData(temp_msg2);
                       
            // debug(DTR("Receiver_2 message %s"), temp_msg2.c_str());
            temp_msg2.clear();
          }while(!m_r2.valid);

          do{
            std::getline(m_v3, temp_msg3);
            m_r3 = parseFishTagData(temp_msg3);
                        
            // debug(DTR("Receiver_3 message %s"), temp_msg3.c_str());
            temp_msg3.clear();
          }while(!m_r3.valid);

          isInit = 1;

          double tms[] = {m_r1.unix_milisecond_time, m_r2.unix_milisecond_time, m_r3.unix_milisecond_time};
          m_time_ms = min(Matrix(tms,3,1)) - 1000;
          // m_time_ms = 1464349787*1000;
        }

        if(m_time_ms >= m_r1.unix_milisecond_time)
        {
          IMC::DevDataText imc_msg;
          imc_msg.setSource(m_r1.serial_no);
          imc_msg.value = m_r1.msg;
          m_o << m_r1.msg << std::endl;

          do{
            if(std::getline(m_v1, temp_msg1))
            {
                IMC::RemoteSensorInfo tagPosition;
                tagPosition.lat = m_r1.lat;
                tagPosition.lon = m_r1.lon;
                tagPosition.alt = m_r1.sensor_data*0.392;
                //tagPosition.data = std::to_string(m_r1.snr) + std::to_string(m_r1.code_type) + "," + std::to_string(m_pf->x(2));
                tagPosition.id = "Receiver" + std::to_string(m_r1.serial_no);
                dispatch(tagPosition);

              dispatch(imc_msg);
              if(m_r1.valid) {
                tag1 = toIMCTag(m_r1);
                tag1.setSource(m_r1.serial_no);
                dispatch(tag1);
                //inf("sent1");
              }
              // m_o << m_r1.msg << std::endl;

              m_r1 = parseFishTagData(temp_msg1);
              //m_r1.sensor_data = m_r1.sensor_data*0.392;

              
              debug(DTR("dispatch Receiver_1 data at %f"), m_time_ms);
              temp_msg1.clear();
            }
            else
            {
              inf(DTR("Reached End of File!"));
              break;
            }
          }while(!m_r1.valid);
        }

        if(m_time_ms >= m_r2.unix_milisecond_time)
        {
          IMC::RemoteSensorInfo tagPosition;
          tagPosition.lat = m_r2.lat;
          tagPosition.lon = m_r2.lon;
          tagPosition.alt = m_r2.sensor_data*0.392;
          //tagPosition.data = std::to_string(m_r1.snr) + std::to_string(m_r1.code_type) + "," + std::to_string(m_pf->x(2));
          tagPosition.id = "Receiver" + std::to_string(m_r2.serial_no);
          dispatch(tagPosition);

          IMC::DevDataText imc_msg;
          imc_msg.setSource(m_r2.serial_no);
          imc_msg.value = m_r2.msg;
          m_o << m_r2.msg << std::endl;

          do{
            if(std::getline(m_v2, temp_msg2))
            {
              dispatch(imc_msg);
              if(m_r2.valid) {
                tag2 = toIMCTag(m_r2);
                tag2.setSource(m_r2.serial_no);
                dispatch(tag2);
                //inf("sent2");
              }
              // m_o << m_r2.msg << std::endl;

              m_r2 = parseFishTagData(temp_msg2);;
              //m_r2.sensor_data = m_r2.sensor_data*0.392;
              
              debug(DTR("dispatch Receiver_2 data at %f"), m_time_ms);
              temp_msg2.clear();
            }
            else
            {
              inf(DTR("Reached End of File!"));
              break;
            }
          }while(!m_r2.valid);
        }

        if(m_time_ms >= m_r3.unix_milisecond_time)
        {
          IMC::RemoteSensorInfo tagPosition;
          tagPosition.lat = m_r3.lat;
          tagPosition.lon = m_r3.lon;
          tagPosition.alt = m_r3.sensor_data*0.392;
          //tagPosition.data = std::to_string(m_r1.snr) + std::to_string(m_r1.code_type) + "," + std::to_string(m_pf->x(2));
          tagPosition.id = "Receiver" + std::to_string(m_r3.serial_no);
          dispatch(tagPosition);

          IMC::DevDataText imc_msg;
          imc_msg.setSource(m_r3.serial_no);
          imc_msg.value = m_r3.msg;
          m_o << m_r3.msg << std::endl;
          do{
            if(std::getline(m_v3, temp_msg3))
            {
              dispatch(imc_msg);
              debug("%d", m_r3.valid);
              if(m_r3.valid) {
                tag3 = toIMCTag(m_r3);
                tag3.setSource(m_r3.serial_no);
                dispatch(tag3);
                //inf("sent3");
              }
              // m_o << m_r3.msg << std::endl;

              m_r3 = parseFishTagData(temp_msg3);
              //m_r3.sensor_data = m_r3.sensor_data*0.392;
              
              debug(DTR("dispatch Receiver_3 data at %f"), m_time_ms);
              temp_msg3.clear();
            }
            else
            {
              inf(DTR("Reached End of File!"));
              break;
            }
          }while(!m_r3.valid);
        }

        //>>--------------- When using combined good data files

        // std::string temp_msg;
        // if(!isInit)
        // {
        //   //>> When three separate raw data files
        //   do{
        //     std::getline(m_v1, temp_msg);
        //     m_r1 = parseFishTagData(temp_msg);
        //     temp_msg.clear();
        //   }while(!m_r1.valid);
        //
        //   isInit = 1;
        //
        //   m_time_ms = m_r1.unix_milisecond_time - 1000;
        // }
        //
        // if(m_time_ms >= m_r1.unix_milisecond_time)
        // {
        //   IMC::DevDataText imc_msg;
        //   imc_msg.setSource(m_r1.serial_no);
        //   imc_msg.value = m_r1.msg;
        //
        //   do{
        //     if(std::getline(m_v1, temp_msg))
        //     {
        //       dispatch(imc_msg);
        //       debug(DTR("dispatch %d data at %f"), m_r1.serial_no, m_time_ms);
        //       m_o << m_r1.msg << std::endl;
        //
        //       m_r1 = parseFishTagData(temp_msg);
        //       temp_msg.clear();
        //     }
        //     else
        //     {
        //       inf(DTR("Reached End of File!"));
        //       break;
        //     }
        //   }while(!m_r1.valid);
        // }

        // std::getline(m_v3, temp_msg1);
        // m_r1 = parseFishTagData(temp_msg1);
        // temp_msg1.clear();
        m_time_ms = m_time_ms + (1/getFrequency())*1000;

        //>>--------------- When using fake data
        // IMC::DevDataText imc_msg1, imc_msg2, imc_msg3;
        // std::string data1 = "$TBR,000042,1473077326,347,S256,1,6,52,11420,6337.25512,N,00936.21981,E,2,12,1.0,*0A";
        // std::string data2 = "$TBR,000022,1473077326,345,S256,1,6,45,8692,6337.25168,N,00936.22455,E,2,12,0.6,*32";
        // std::string data3 = "$TBR,000045,1473077326,343,S256,1,6,55,6611,6337.24978,N,00936.22341,E,1,12,0.7,*39";
        // imc_msg1.value = data1;
        // imc_msg1.setSource(000042);
        // imc_msg2.value = data2;
        // imc_msg2.setSource(000022);
        // imc_msg3.value = data3;
        // imc_msg3.setSource(000045);
        // dispatch(imc_msg1);
        // dispatch(imc_msg2);
        // dispatch(imc_msg3);
        // inf("Sent fake message!");

      }// End of task Function
    };// End of struct Task
  }
}
}
}

DUNE_TASK
