//***************************************************************************
// Copyright 2013-2021 Norwegian University of Science and Technology (NTNU)*
// Department of Engineering Cybernetics (ITK)                              *
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
// Author: Nikolai Lauvås                                                   *
//***************************************************************************

// DUNE headers.
#include <DUNE/DUNE.hpp>


#include "ExtendedKalmanFilter.h"
#include "ParticleFilter.h"
#include "AlgebraicSolution.hpp"

#include <boost/circular_buffer.hpp>
namespace SourceEstimators
{
  //! TODO: Implement altitude in TBR and use this. More todos in code.
  //! TODO: Implement use preassure/depth from tag. 
  //! Insert explanation on task behaviour here.
  //! @author Nikolai Lauvås
  namespace SingleReceiverDepth
  {
    using DUNE_NAMESPACES;

    static const unsigned c_buffer_size = 5;

    struct Arguments
    {
      //! Fitting of SNR to range
      std::vector<double> ranging_snr_fit;
      //! Jitter upper limit (seconds)
      float max_jitter;
      //! Expected tag Period
      float tag_period;
      //! Period between when filter is run
      float filter_timestep;
      //! Tag ID
      unsigned tag_id;
      //! Time to wait in while. In practice, this controls how regular the filter timing is
      float message_wait_time;
      //! How deep the receiver is mounted in altitude
      float receiver_depth;

// Initial Parameters for calculating Speed of Sound
      //! Initial Speed of Sound in water
      float init_c_sound;
      //! Initial temperature used for calculating Speed of Sound in water
      float init_temperature;
      //! Initial salinity used for calculating Speed of Sound in water
      float init_salinity;
      //! Initial depth used for calculating Speed of Sound in water
      float init_depth;
// Parameters for updating Speed of Sound
      //! Should the Speed of Sound in water be updated with water temperature
      bool update_c_sound_temp;
      //! Entity delivering water temperatures for updating the Speed of Sound in water
      std::string entity_c_sound_temp;
      //! Should the Speed of Sound in water be updated from measurement
      bool update_c_sound;
      //! Entity providing the Speed of Sound in water
      std::string entity_c_sound;
      //! Should the Speed of Sound in water be updated with water salinity measurement
      bool update_c_sound_salinity;
      //! Entity delivering water salinity for updating the Speed of Sound in water
      std::string entity_c_sound_salinity;

      //! How far back into the buffer to attempt period matching.
      int max_correction_attempts;
// Kalman Filter
      //! Extended Kalman filter - Qm
      std::vector<double> ekf_Qm;      
      //! Extended Kalman filter - Rm
      std::vector<double> ekf_Rm; 
      //! Extended Kalman filter - P0
      std::vector<double> ekf_P0;
      std::vector<double> ekf_x0;
// Particle Filter
      //! Particle Filter Range
      double pf_range;
      //! Particle Filter width
      int pf_particle_width;
      //! Particle Filter sigma_r squared
      double pf_sigma_r;
      //! Particle Filter sigma_rd squared
      double pf_sigma_rd;
// Location settings

      uint32_t receiver_serial;
      //! Reference coordinate position (degrees)
      std::vector<double> reference;
    };
    struct Task: public DUNE::Tasks::Task
    {
      //! Task arguments.
      Arguments m_args;
      //! Buffer holding received tag detections.
      boost::circular_buffer<IMC::TBRFishTag> *tagBuffer;
      //! Timer responsible for running filter timestep
      Time::Counter<float> m_filter_timer;
      //! Current Speed of sound in water
      float m_c_speed;
      //! Current temperature used for calculating Speed of sound in water
      float m_c_speed_temp;
      //! Current salinity used for calculating Speed of sound in water
      float m_c_speed_salinity;
      //! Current depth used for calculating Speed of sound in water
      float m_c_speed_depth;      
      //! Extended Kalman filter used for estimating transmitter state (position)
      ExtendedKalmanFilter m_ekf;
      //! Particle filter used for estimating transmitter state (position)
      ParticleFilter m_pf;
      //! Particle filter used for estimating transmitter state (position)
      AlgebraicSolver m_aslv;
      //! Temperature entity label.
      int m_temp_eid;
      //! Speed of sound provider entity label.
      int m_c_sound_eid;
      //! Salinity provider entity label.
      int m_salinity_eid;

      //! How far back into the buffer to attempt period matching.
      int m_max_correction_attempts;
      //! Reference coordinate used to calculate NED frame
      double m_refCoord[3];
      //! Constructor.
      //! @param[in] name task name.
      //! @param[in] ctx context.
      Task(const std::string& name, Tasks::Context& ctx):
        DUNE::Tasks::Task(name, ctx)
      {
        param("Ranging - SNR Fit", m_args.ranging_snr_fit)
        .size(2)
        .description("Linear fit of SNR to distance (a, b, forms r=ax+b)");

        param("Tag Max Jitter", m_args.max_jitter)
        .description("The maximum jitter allowed between measurements")
        .units(Units::Second)
        .defaultValue("0.01");

        param("Tag Transmitt Period", m_args.tag_period)
        .description("The expected period between tag registration")
        .units(Units::Second)
        .defaultValue("7.0");

        param("Receiver Serial Number", m_args.receiver_serial)
        .description("The serial number of the receiver to accept tag registrations from.")
        .defaultValue("634");

        param("Receiver Depth", m_args.receiver_depth)
        .description("The depth of the receiver providing sensor messages")
        .units(Units::Meter)
        .defaultValue("3.0");

        param("Message Wait Time", m_args.message_wait_time)
        .description("The time to wait for new messages in the while loop between checking timer.")
        .units(Units::Second)
        .defaultValue("0.01");

        param("Filter Timestep", m_args.filter_timestep)
        .description("The timestep of the filter")
        .units(Units::Second)
        .defaultValue("7.0");

        param("Max Correction Attempts", m_args.max_correction_attempts)
        .description("How far back into the buffer to attempt period matching. -1 gives max allowed in used buffer")
        .defaultValue("-1");

        param("Tag ID", m_args.tag_id)
        .description("The ID of the tracked fish tag")
        .defaultValue("40");
// Initial Parameters for calculating Speed of Sound
        param("Initial Speed Of Sound", m_args.init_c_sound)
        .units(Units::MeterPerSecond)
        .description("The ID of the tracked fish tag")
        .defaultValue("1485.0");

        param("Initial Temperature", m_args.init_temperature)
        .units(Units::DegreeCelsius)
        .description("The ID of the tracked fish tag")
        .defaultValue("10.0");

        param("Initial Salinity", m_args.init_salinity)
        .description("Initial value of salinity in PPM")
        .defaultValue("25.0");

        param("Initial Depth", m_args.init_depth)
        .units(Units::Meter) 
        .description("Initial depth used to calculate Speed of Sound in water")
        .defaultValue("1.0");  
// Parameters for updating Speed of Sound
        param("Use Temperature For Speed Of Sound", m_args.update_c_sound_temp)
        .description("If any temperature measurements are to be used for updating the speed of sound.")
        .defaultValue("false");

         param("Temperature - Entity", m_args.entity_c_sound_temp)
        .description("The entity delivering water temperatures for updating the Speed of Sound in water")
        .defaultValue("Hydrophone");

        param("Use Salinity For Speed Of Sound", m_args.update_c_sound_salinity)
        .description("If any salinity measurements are to be used for updating the speed of sound.")
        .defaultValue("false");

         param("Salinity - Entity", m_args.entity_c_sound_salinity)
        .description("The entity delivering salinity for updating the Speed of Sound in water")
        .defaultValue("CTD");

        param("Use Speed Of Sound Measurement", m_args.update_c_sound)
        .description("If speed of sound is provided through IMC messages.")
        .defaultValue("false");

         param("Speed Of Sound - Entity", m_args.entity_c_sound)
        .units(Units::MeterPerSecond)
        .description("The entity delivering the Speed of Sound in water")
        .defaultValue("CTD");        


// Kalman Filter Parameters
        param("EKF - x0", m_args.ekf_x0)
        .size(3)
        .description("Initial X value for the extended Kalman filter")
        .defaultValue("0.0, 0.0, 0.0");

        param("EKF - P0", m_args.ekf_P0)
        .size(9)
        .description("Initial P matrix value for the extended Kalman filter, first row")
        .defaultValue("66458, -26820, 0, -26820, 12116, 0, 0, 0, 0");

        param("EKF - Rm", m_args.ekf_Rm)
        .description("Initial X value for the extended Kalman filter")
        .defaultValue("1.5*3.8706, 0, 0, 2.1638e6");

        param("EKF - Qm", m_args.ekf_Qm)
        .size(9)
        .description("Initial X value for the extended Kalman filter")
        .defaultValue("0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00001");
// Particle Filter Arguments
        param("PF - Range", m_args.pf_range)
        .description("Radius of square containing particles")
        .units(Units::Meter)
        .defaultValue("400.0");

        param("PF - Particles Per Direction", m_args.pf_particle_width)
        .description("The total amount of particles equals n^2")
        .defaultValue("50");

        param("PF - Range Sigma Squared SNR", m_args.pf_sigma_r)
        .description("The expected period between tag registration")
        .defaultValue("1470.986063836");

        param("PF - Range Sigma Squared TDOA", m_args.pf_sigma_rd)
        .description("The expected period between tag registration")
        .defaultValue("3.8706");
// Others
        param("Reference Coordinate", m_args.reference)
        .units(Units::Degree)
        .size(2)
        .description("Origin of the reference coordinate system");

        bind<IMC::TBRFishTag>(this);
        bind<IMC::Temperature>(this);
        bind<IMC::Salinity>(this);
        bind<IMC::SoundSpeed>(this);
      }

      //! Update internal state with new parameter values.
      void
      onUpdateParameters(void)
      {
        m_refCoord[0] = Math::Angles::radians(m_args.reference[0]);
        m_refCoord[1] = Math::Angles::radians(m_args.reference[1]);
        m_refCoord[2] = 0.0;

        m_c_speed=m_args.init_c_sound;
        m_c_speed_salinity = m_args.init_salinity;
        m_c_speed_temp = m_args.init_temperature;
        m_c_speed_depth = m_args.init_depth;

        if(!m_args.update_c_sound) {
          if(m_args.update_c_sound_salinity || m_args.update_c_sound_temp) {
            m_c_speed=calculateSpeedOfSound(m_c_speed_temp, m_c_speed_salinity, m_c_speed_depth);
            spew("Calculated initial c_speed: %f", m_c_speed);
          }
        }

        if(m_args.max_correction_attempts < 0)
          m_max_correction_attempts = c_buffer_size-2;
        else {
          m_max_correction_attempts = m_args.max_correction_attempts;
        }

        if(paramChanged(m_args.filter_timestep))
          m_filter_timer.setTop(m_args.filter_timestep);
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
        try
          {
            m_c_sound_eid = resolveEntity(m_args.entity_c_sound);
          }
          catch (...)
          {
            if(m_args.update_c_sound) {
              err("Could not find entity: %s", m_args.entity_c_sound.c_str());
            }
            m_c_sound_eid = 0;
          }
        try
          {
            m_temp_eid = resolveEntity(m_args.entity_c_sound_temp);
          }
          catch (...)
          {
            if(m_args.update_c_sound_temp) {            
              err("Could not find entity: %s", m_args.entity_c_sound_temp.c_str());
            }
            m_temp_eid = 0;

          }
        try
          {
            m_salinity_eid = resolveEntity(m_args.entity_c_sound_salinity);
          }
          catch (...)
          {
            if(m_args.update_c_sound_salinity) {            
              err("Could not find entity: %s", m_args.entity_c_sound_salinity.c_str());
            }
            m_salinity_eid = 0;
          }
      }

      //! Acquire resources.
      void
      onResourceAcquisition(void)
      {
        tagBuffer = new boost::circular_buffer<IMC::TBRFishTag>(c_buffer_size);

        // Create Extended Kalman Filter Instance
        arma::vec x0(m_args.ekf_x0);
        arma::mat P0( &(m_args.ekf_P0).front(), 3, 3);
        arma::mat Qm( &(m_args.ekf_Qm).front(), 3, 3);
        arma::mat Rm( &(m_args.ekf_Rm).front(), std::sqrt(m_args.ekf_Rm.size()), std::sqrt(m_args.ekf_Rm.size()));

        m_ekf = ExtendedKalmanFilter(x0, P0, Qm, Rm);
        // Create ParticleFilter instance
        m_pf = ParticleFilter(m_args.pf_particle_width, m_args.pf_range, m_args.pf_sigma_rd, m_args.pf_sigma_r);
        //std::cout << P0 << std::endl << P02 << std::endl;
      }

      void
      consume(const IMC::TBRFishTag* msg)
      {
      std::ofstream logOutStream;
      logOutStream.open(("log/tag.log"), std::fstream::app);
      if (logOutStream.good()) {
        logOutStream.precision(15);
          logOutStream << DUNE::Math::Angles::degrees(msg->lat) << "," << DUNE::Math::Angles::degrees(msg->lon) << "," << msg->unix_timestamp << "," << msg->millis<< std::endl;
          logOutStream.close();
      }   
        if (m_args.receiver_serial == msg->serial_no) {
          if(msg->trans_id == m_args.tag_id) {
            tagBuffer->push_back(*msg);
            updateFilter();
          }
          // Ignore other tags
        }
      }

      void
      consume(const IMC::SoundSpeed* msg)
      {
        if(msg->getSourceEntity() == m_c_sound_eid) {
          if(m_args.update_c_sound) {
            m_c_speed = msg->value;
            spew("Setting c_sound to: %f", msg->value);
          }
        }
      }

      void
      consume(const IMC::Temperature* msg)
      {
        if(msg->getSourceEntity() == m_temp_eid) {
          if(m_args.update_c_sound_temp) {
            m_c_speed_temp = msg->value;
            m_c_speed = calculateSpeedOfSound(m_c_speed_temp, m_c_speed_salinity, m_c_speed_depth);
            spew("Recalculating c_sound with temperature: %f, Result: %f", msg->value, m_c_speed);
          }
        }
      }

      void
      consume(const IMC::Salinity* msg)
      {
        if(msg->getSourceEntity() == m_salinity_eid) {
          if(m_args.update_c_sound_salinity) {
            m_c_speed_salinity = msg->value;
            m_c_speed = calculateSpeedOfSound(m_c_speed_temp, m_c_speed_salinity, m_c_speed_depth);
            spew("Recalculating c_sound with salinity: %f, Result: %f", msg->value, m_c_speed);
          }
        }
      }

      //! Leroys formula for calculating the speed of sound in water
      //! Temperature given in Celsius
      //! Salinity given in PPM
      //! Depth given in Meters
      double calculateSpeedOfSound(double temperature, double salinity, double depth) {
        double c_speed = 1492.9;
        c_speed += 3*(temperature -10.0);
        c_speed -= (6*pow(temperature-10.0, 2.0))/1000;
        c_speed -= (4*pow(temperature-18.0, 2.0))/100;
        c_speed += 1.2*(salinity-35);
        c_speed -= ((temperature-18-0)*(salinity-35))/100;
        c_speed += depth/61;
        return c_speed;
      }
      //! Initialize resources.
      void
      onResourceInitialization(void)
      {
        //! Set timer for periodic check of surroundings.
        m_filter_timer.setTop(m_args.filter_timestep);
      }

      //! Release resources.
      void
      onResourceRelease(void)
      {
        Memory::clear(tagBuffer);
        //Memory::clear(m_ekf);
        //Memory::clear(m_pf);
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

      //! 
      void updateFilter(void) {
        //inf("Update %ld", c_buffer_size - tagBuffer->size());
        double measurement_millis = tagBuffer->rbegin()->unix_timestamp + (double)tagBuffer->rbegin()->millis/1000;
        // TODO: Stop at m_args.max_correction_attempts
        for(boost::circular_buffer<DUNE::IMC::TBRFishTag>::reverse_iterator i=tagBuffer->rbegin()+1; i != tagBuffer->rend();i++) {
          //inf("%d - %d", tagBuffer->rbegin()->unix_timestamp, i->unix_timestamp);
          double td = measurement_millis - i->unix_timestamp - (double)i->millis/1000;
          double closestMultipleOfPeriod = m_args.tag_period*std::round(td/m_args.tag_period);
          double tdoa = td - closestMultipleOfPeriod;
          inf("delta %f %f", td, closestMultipleOfPeriod);

          if(abs(tdoa) < m_args.max_jitter || td > 60.0) {
              double altitude=m_args.receiver_depth;
              double P[2] = {m_args.ranging_snr_fit[0], m_args.ranging_snr_fit[1]}; // Found on page 54 of master
              double rdoa = m_c_speed*tdoa; // Range difference
              double rangeSNR = (tagBuffer->rbegin()->snr - P[1])/P[0];
              double depth = i->trans_data*0.392;
              arma::vec ranging = {rdoa,rangeSNR};
              arma::vec measurements = {rdoa,rangeSNR, depth};

              double NED1[3];
              double NED2[3];

              toNEDframe(*tagBuffer->rbegin(), m_refCoord, NED1);
              arma::vec position_current(NED1,3);
              toNEDframe(*i, m_refCoord, NED2);
              arma::vec position_previous(NED2,3);
              position_current(2) = altitude;
              position_previous(2) = altitude;

              
              int pf_k = m_pf.update(ranging,position_previous,position_current);

              //if(pf_k>3) {
                if (m_aslv.addMeasurement(measurements,position_previous,position_current)) {
                  double result[3] = {m_aslv.x(0), m_aslv.x(1), m_aslv.x(2)};
                  arma::vec X0={m_aslv.x(0), m_aslv.x(1), m_aslv.x(2)};
                  if(!m_ekf.isInitialized()) m_ekf.initialize(X0);
                  double latLon[3];
                  fromNEDframe(result, m_refCoord, latLon);
                  std::ofstream logOutStream;
                  logOutStream.open(("log/aslv.log"), std::fstream::app);
                  if (logOutStream.good()) {
                    logOutStream.precision(15);
                    logOutStream << m_aslv.x(0) << "," << m_aslv.x(1) << "," << m_aslv.x(2) << "," << DUNE::Math::Angles::degrees(latLon[0]) << "," << DUNE::Math::Angles::degrees(latLon[1]) << std::endl;
                    logOutStream.close();
                  }   
                }
                if(m_args.ekf_Rm.size()==1) {
                  arma::vec rdoavec = {rdoa};
                  m_ekf.measurementStep(rdoavec,position_previous,position_current);
                } else if(m_args.ekf_Rm.size()>4) {
                  m_ekf.measurementStep(measurements,position_previous,position_current);
                } else {
                  m_ekf.measurementStep(ranging,position_previous,position_current);
                }
              //}
              
              if (pf_k > 0) { // Checks if an estimate has been made
                double result[3] = {m_pf.x(0), m_pf.x(1), m_pf.x(2)};
                double latLon[3];
                fromNEDframe(result, m_refCoord, latLon);

                IMC::RemoteSensorInfo tagPosition;
                tagPosition.lat = latLon[0];
                tagPosition.lon = latLon[1];
                tagPosition.alt = -m_pf.x(2);
                tagPosition.data = std::to_string(m_pf.x(0)) + std::to_string(m_pf.x(1)) + "," + std::to_string(m_pf.x(2));
                tagPosition.id = "DepthPF" + std::to_string(m_args.receiver_serial);
                dispatch(tagPosition);
                std::ofstream logOutStream;
                  logOutStream.open(("log/pf.log"), std::fstream::app);
                  if (logOutStream.good()) {
                    logOutStream.precision(15);
                    logOutStream << m_pf.x(0) << "," << m_pf.x(1) << "," << m_pf.x(2) << "," << DUNE::Math::Angles::degrees(latLon[0]) << "," << DUNE::Math::Angles::degrees(latLon[1]) << std::endl;
                    logOutStream.close();
                  }
                spew("New PF Estimate: (N,E,D,La,Lo)= %.15f,%.15f,%.15f,%.15f, %.15f", m_pf.x(0), m_pf.x(1), m_pf.x(2),DUNE::Math::Angles::degrees(latLon[0]),DUNE::Math::Angles::degrees(latLon[1]));
                return;
            }
          }
        }
        war("No good TDOA value found");
      }

      //! Main loop.
      void
      onMain(void)
      {

        while (!stopping())
        {
          if(m_filter_timer.overflow()) {
            m_filter_timer.reset();
            if(m_ekf.isInitialized()) {
              int ekf_k = m_ekf.predictionStep(); // Filter time update
              if (ekf_k > 0 ) { // Checks if an estimate has been made
                  double lati,longi;
                  double result[3] = {m_ekf.x(0),m_ekf.x(1),m_ekf.x(2)};
                double latLon[3];
                fromNEDframe(result, m_refCoord, latLon);
                lati=latLon[0], longi=latLon[1];
                  // Send output to Neptus/DUNE log
                  IMC::RemoteSensorInfo tagPosition;
                  tagPosition.lat = lati;
                  tagPosition.lon = longi;
                  tagPosition.alt = -m_ekf.x(2);
                  tagPosition.data = std::to_string(m_ekf.x(0)) + std::to_string(m_ekf.x(1)) + "," + std::to_string(m_ekf.x(2));
                  tagPosition.id = "DepthEKF" + std::to_string(m_args.receiver_serial);
                  dispatch(tagPosition);
                  std::ofstream logOutStream;
                  logOutStream.open(("log/ekf.log"), std::fstream::app);
                  if (logOutStream.good()) {
                    logOutStream.precision(15);
                      logOutStream << m_ekf.x(0) << "," << m_ekf.x(1) << "," << m_ekf.x(2) << "," << DUNE::Math::Angles::degrees(lati) << "," << DUNE::Math::Angles::degrees(longi) << std::endl;
                      logOutStream.close();
                  }   
                  spew("New Kalman Estimate: (N,E,D,La,Lo)= %.15f,%.15f,%.15f,%.15f, %.15f", m_ekf.x(0), m_ekf.x(1), m_ekf.x(2),DUNE::Math::Angles::degrees(lati),DUNE::Math::Angles::degrees(longi));
              }
            }
          }
          waitForMessages(m_args.message_wait_time);
        }
      }
    };
  }
}

DUNE_TASK
