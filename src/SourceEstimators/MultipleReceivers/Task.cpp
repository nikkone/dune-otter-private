//***************************************************************************
// Copyright 2007-2020 Universidade do Porto - Faculdade de Engenharia      *
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
// Author: Nikolai Lauvås based on works by Praveen Jain and Artur Zolich   *
//***************************************************************************

// DUNE headers.
#include <DUNE/DUNE.hpp>
#include "ExtendedKalmanFilter.hpp"
namespace SourceEstimators
{
  //! Insert short task description here.
  //!
  //! Insert explanation on task behaviour here.
  //! @author Nikolai Lauvås
  namespace MultipleReceivers
  {
    using DUNE_NAMESPACES;

    struct FishPosition
    {
      double lat;
      double lon;
      double hae;
      double n;
      double e;
      double d;
    };

    struct RefCoordinate
    {
      double lat;
      double lon;
      double hae;
    };
    //! %Task arguments.
    struct Arguments
    {
      //! Tag ID
      uint32_t tag_id;
      //! Receiver 1 carrier
      uint32_t rec_1_serial;
      //! Receiver 2 carrier
      uint32_t rec_2_serial;
      //! Receiver 3 carrier
      uint32_t rec_3_serial;

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

      //! Structure that hold the parsed fish tag data
      ReceiverData m_r1, m_r2, m_r3;
      IMC::TBRFishTag tbr_r1, tbr_r2, tbr_r3;

      //! Output file only for use in simulations
      // std::ofstream m_o;

      //! Reference coordinate system
      RefCoordinate m_refCoord;

      //! Additional variables used to log the EKF data
      //double m_norm_P, m_trace_P;

      //! Lat
      float m_lat;

      //! Lon
      float m_lon;

      //! Structure to hold the EKF data
      EKFilter *m_ekf;

      //! Constructor.
      //! @param[in] name task name.
      //! @param[in] ctx context.
      Task(const std::string& name, Tasks::Context& ctx):
        DUNE::Tasks::Periodic(name, ctx)
      {
        param("Receiver 1 Serial Number", m_args.rec_1_serial)
        .description("Receiver 1 Serial Number")
        .defaultValue("47");

        param("Receiver 2 Serial Number", m_args.rec_2_serial)
        .description("Receiver 2 Serial Number")
        .defaultValue("48");

        param("Receiver 3 Serial Number", m_args.rec_3_serial)
        .description("Receiver 3 Serial Number")
        .defaultValue("49");

        param("Tag ID", m_args.tag_id)
        .description("The ID of the tracked fish tag")
        .defaultValue("40");

        param("Reference Coordinate", m_args.reference)
        .units(Units::Degree)
        .size(2)
        .description("Origin of the reference coordinate system");

        param("Initial Position", m_args.position)
        .size(3)
        .description("Initial Fish tag position in NED from reference coordinate frame");
        bind<IMC::TBRFishTag>(this);
      }

      //! Update internal state with new parameter values.
      void
      onUpdateParameters(void)
      {
        m_args.dt = 1/getFrequency();

        m_refCoord.lat = Math::Angles::radians(m_args.reference[0]);
        m_refCoord.lon = Math::Angles::radians(m_args.reference[1]);
        m_refCoord.hae = 0.0; // TOdo: Fix
        m_lat = m_refCoord.lat;
        m_lon = m_refCoord.lon;
        debug(DTR("Reference Cooridnate: %f %f"), m_args.reference[0], m_args.reference[1]);
        debug(DTR("Initial Fish Position: %f %f %f"), m_args.position[0], m_args.position[1], m_args.position[2]);

        /*m_r1_name = resolveSystemName(m_args.rec_1_vehicle_name);
        m_r2_name = resolveSystemName(m_args.rec_2_vehicle_name);
        m_r3_name = resolveSystemName(m_args.rec_3_vehicle_name);
        inf(DTR("Carrier of receiver 1 name is %s, with ID %d"), m_args.rec_1_vehicle_name.c_str(), m_r1_name);
        inf(DTR("Carrier of receiver 2 name is %s, with ID %d"), m_args.rec_2_vehicle_name.c_str(), m_r2_name);
        inf(DTR("Carrier of receiver 3 name is %s, with ID %d"), m_args.rec_3_vehicle_name.c_str(), m_r3_name);*/
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
        m_ekf = new EKFilter(3,4,4,0,3, m_args.dt);
        setEntityState(IMC::EntityState::ESTA_NORMAL, Status::CODE_ACTIVE);
        inf("Started, initialized EKFilter");
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


      void
      createNEDframe(ReceiverData *r1, ReceiverData *r2, ReceiverData *r3)
      {
        WGS84::displacement(m_refCoord.lat, m_refCoord.lon, 0.0, r1->lat, r1->lon, 0.0, &(r1->N), &(r1->E), &(r1->D));
        WGS84::displacement(m_refCoord.lat, m_refCoord.lon, 0.0, r2->lat, r2->lon, 0.0, &(r2->N), &(r2->E), &(r2->D));
        WGS84::displacement(m_refCoord.lat, m_refCoord.lon, 0.0, r3->lat, r3->lon, 0.0, &(r3->N), &(r3->E), &(r3->D));
      }
      void
      consume(const IMC::TBRFishTag* msg)
      {
        if(msg->trans_id == m_args.tag_id) {
          //tagBuffer->push_back(*msg);
          m_ekf->setActive(true);
          //updateFilter();

          IMC::RemoteSensorInfo tagPosition;
          tagPosition.lat = msg->lat;
          tagPosition.lon = msg->lon;
          tagPosition.alt = msg->trans_data;

          if(msg->serial_no == m_args.rec_1_serial)
          {
            debug(DTR("consume: Message from R1 arrived"));
            tbr_r1 = *msg;
            inf("Received tag: %u", tbr_r1.trans_id);
            m_r1.isMeasure = 1;
            m_r1.lat = Math::Angles::degrees(msg->lat);
            m_r1.lon = Math::Angles::degrees(msg->lon);
            m_r1.unix_milisecond_time = msg->unix_timestamp*1000+msg->millis;
            m_r1.sensor_data = 2;//msg->trans_data * 0.392;

            tagPosition.id = "Receiver_1";
            dispatch(tagPosition);
          }
          else if(msg->serial_no == m_args.rec_2_serial)
          {
            debug(DTR("consume: Message from R2 arrived"));
            tbr_r2 = *msg;
            m_r2.isMeasure = 1;
            m_r2.lat = Math::Angles::degrees(msg->lat);
            m_r2.lon = Math::Angles::degrees(msg->lon);
            m_r2.unix_milisecond_time = msg->unix_timestamp*1000+msg->millis;
            m_r2.sensor_data = 2;//msg->trans_data * 0.392;
            tagPosition.id = "Receiver_2";
            dispatch(tagPosition);
          }
          else if(msg->serial_no == m_args.rec_3_serial)
          {
            debug(DTR("consume: Message from R3 arrived"));
            tbr_r3 = *msg;
            m_r3.isMeasure = 1;
            m_r3.lat = Math::Angles::degrees(msg->lat);
            m_r3.lon = Math::Angles::degrees(msg->lon);
            m_r3.unix_milisecond_time = msg->unix_timestamp*1000+msg->millis;
            m_r3.sensor_data = 2;//msg->trans_data * 0.392;
            tagPosition.id = "Receiver_3";
            dispatch(tagPosition);
          }
          else
          {
            inf(DTR("consume: Message from non registered receiver with ID %d"), msg->serial_no);
          }
        } else {
          inf(DTR("consume: Message from non registered tag with ID %d"), msg->trans_id);
        }
      }
      
      void
      consume(const IMC::EstimatedState* msg)
      {
        setEntityState(IMC::EntityState::ESTA_NORMAL, Status::CODE_ACTIVE);
                spew("Estimated State arrived from %d", msg->getSource());

                if (getSystemId() == msg->getSource())
                {
                  m_lat = msg->lat;
                  m_lon = msg->lon;

                }

      }

      void
      findFishPosition(ReceiverData &r1, ReceiverData &r2, ReceiverData &r3)
      {
        if(!m_ekf->isInitialized())
        {
          cri("ERROR");
          return;
        }
        //>> Step 1: Predict
        m_ekf->predict();

        //>> Step 2: If measurements available

        uint8_t m12 = ((r1.isMeasure || r2.isMeasure) && (m_ekf->timeShiftCorrect(&r1, &r2)));
        uint8_t m23 = ((r2.isMeasure || r3.isMeasure) && (m_ekf->timeShiftCorrect(&r2, &r3)));
        uint8_t m31 = ((r3.isMeasure || r1.isMeasure) && (m_ekf->timeShiftCorrect(&r3, &r1)));

        if((m12 == 1) || (m23 == 1) || (m31 == 1))
        {
          //>> Step 2.1: Update System Matrices and Covariance
          createNEDframe(&r1, &r2, &r3);
          /*inf("R1: %f %f %f %f %f %f", r1.lat, r1.lon, r1.N, r1.E, r1.D, r1.unix_milisecond_time);
          inf("R2: %f %f %f %f %f %f", r2.lat, r2.lon, r2.N, r2.E, r2.D, r2.unix_milisecond_time);
          inf("R3: %f %f %f %f %f %f", r3.lat, r3.lon, r3.N, r3.E, r3.D, r2.unix_milisecond_time);*/
                  
          m_ekf->constructCmatrix(&r1, &r2, &r3);
          //>> Step 2.2: Measurement Update
          m_ekf->update();
        }

        //>> dispatch fish position
        IMC::RemoteSensorInfo tagPosition;
        Matrix FishPosition = m_ekf->getEstimate();

        // double x1_lat = r1.lat, x1_lon = r1.lon, x1_depth;
        double x1_lat = m_refCoord.lat; double x1_lon = m_refCoord.lon; double x1_depth = m_refCoord.hae;
        debug("%f %f %f", FishPosition(0,0), FishPosition(1,0), FishPosition(2,0));
        WGS84::displace(FishPosition(0,0), FishPosition(1,0), FishPosition(2,0), &(x1_lat), &(x1_lon), &(x1_depth));

        // Estimated position
        tagPosition.id = "Fish_position_est";
        tagPosition.lat = x1_lat;
        tagPosition.lon = x1_lon;
        tagPosition.alt = -x1_depth;

        /*
        // Compute the estimation error - ground truth - estimated state
        double eN, eE, eD;
        inf("m_lat: %f m_lon%f", m_lat, m_lon);
        WGS84::displacement(m_lat, m_lon, 0.0, x1_lat, x1_lon, 0.0, &(eN), &(eE), &(eD));

        // Additional filter signals to log
        m_norm_P = m_ekf->PHat.norm_p(2);
        m_trace_P = m_ekf->PHat.trace();

        std::stringstream ss;
        double err = std::sqrt(std::pow(eN,2) + std::pow(eE,2) + std::pow(eD,2));
        ss << m_lat << "," << m_lon << "," << err << "," << m_d1 << "," << m_d2 << "," << m_d3 << "," << m_norm_innov << "," << m_norm_P << "," << m_trace_P;*/


        //tagPosition.data = r1.msg + "," +r2.msg+ "," + r3.msg + "," + ss.str();
        dispatch(tagPosition);
        debug(DTR("Fish lat: %f long: %f Depth: %f"), x1_lat, x1_lon, -x1_depth);
            std::ofstream logOutStream;
    logOutStream.open(("log/estimates.log"), std::fstream::app);
    if (logOutStream.good()) {
        logOutStream << Math::Angles::degrees(x1_lat) << "," << Math::Angles::degrees(x1_lon) << "," << -x1_depth;
      
        logOutStream << std::endl;
        logOutStream.close();
    }  
        debug(DTR("Fish lat: %f long: %f Depth: %f"), Math::Angles::degrees(x1_lat), Math::Angles::degrees(x1_lon), -x1_depth);
        // m_P = m_ekf.PHat.norm_p(2);
        // m_o.precision(15);
        //
        // m_o << m_r1.unix_milisecond_time << "," << m_r2.unix_milisecond_time << "," << m_r3.unix_milisecond_time << "," << x1_lat << "," << x1_lon << "," << "," << -x1_depth << "," << FishPosition(2,0) << "," << m_d1 << "," << m_d2 << "," << m_d3 << "," << m_norm_innov << "," << m_P << "," << m_ekf.PHat.trace() << std::endl;
      } //End function findFishPosition


      //! Main repeated periodically.
      void
      task(void)
      {
        consumeMessages();
        if(m_ekf->isActive())
        {
          findFishPosition(m_r1, m_r2, m_r3);
        } // End of if 
      } // End of task function
    };
  }
}

DUNE_TASK
