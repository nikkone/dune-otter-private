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

#include <boost/circular_buffer.hpp>
namespace SourceEstimators
{
  //! TODO: Implement altitude in TBR and use this. More todos in code.
  //! TODO: Implement use preassure/depth from tag. 
  //! Insert explanation on task behaviour here.
  //! @author Nikolai Lauvås
  namespace SingleReceiver
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
      bool location_northern_hemisphere;
      uint8_t locaton_utm_zone;
      uint32_t receiver_serial;
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
      //! Temperature entity label.
      int m_temp_eid;
      //! Speed of sound provider entity label.
      int m_c_sound_eid;
      //! Salinity provider entity label.
      int m_salinity_eid;

      //! How far back into the buffer to attempt period matching.
      int m_max_correction_attempts;

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
        .defaultValue("0.5");

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
        .size(4)
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

        param("UTM Zone", m_args.locaton_utm_zone)
        .description("The expected period between tag registration")
        .defaultValue("32");

        param("Northern Hemisphere", m_args.location_northern_hemisphere)
        .description("The expected period between tag registration")
        .defaultValue("true");

        bind<IMC::TBRFishTag>(this);
        bind<IMC::Temperature>(this);
        bind<IMC::Salinity>(this);
        bind<IMC::SoundSpeed>(this);
      }

      //! Update internal state with new parameter values.
      void
      onUpdateParameters(void)
      {
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
        arma::mat Rm( &(m_args.ekf_Rm).front(), 2, 2);

        m_ekf = ExtendedKalmanFilter(x0, P0, Qm, Rm);
        // Create ParticleFilter instance
        m_pf = ParticleFilter(m_args.pf_particle_width, m_args.pf_range, m_args.pf_sigma_rd, m_args.pf_sigma_r);
        //std::cout << P0 << std::endl << P02 << std::endl;
      }

      void
      consume(const IMC::TBRFishTag* msg)
      {
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

      void updateFilter(void) {
        //TODO: Add discarding older then t tag detections
        //TODO: For completeness, add offsets. Not used on the Otter, because only GPS, and no dead reckoning
        //TODO: Altitude in TBRTag IMC
        //TODO: Check UTM zone is equal, fix zone?
        double altitude=m_args.receiver_depth;
        double P[2] = {m_args.ranging_snr_fit[0], m_args.ranging_snr_fit[1]}; // Found on page 54 of master

        if(tagBuffer->size() == c_buffer_size) {
          double td = (*tagBuffer)[c_buffer_size-1].unix_timestamp + (double)(*tagBuffer)[c_buffer_size-1].millis/1000 - (*tagBuffer)[c_buffer_size-2].unix_timestamp - (double)(*tagBuffer)[c_buffer_size-2].millis/1000;
          //spew("Diff %f", td);
          double closestMultipleOfPeriod = m_args.tag_period*std::round(td/m_args.tag_period);
          double tdoa = td - closestMultipleOfPeriod;

          bool goodSample = false;
          int correction_attempts = 0;
          while(correction_attempts < m_max_correction_attempts) {
            if(abs(tdoa) < m_args.max_jitter) {
                goodSample = true;
                break;
            }
            //inf("Attempt %i %i", c_buffer_size-2-correction_attempts-1, c_buffer_size-2);
            //Bad tdoa value, try another.
            //std::cout << "Time Difference discrepancy. Correcting with older data." << std::endl;
            //std::cout << "Time Difference of Arrival discrepancy." << std::endl;
            //std::cout << "Period time was          : " << td << " seconds" << std::endl; 
            //std::cout << "Expected some multiple of: ~" << m_args.tag_period << " seconds"  << std::endl;
            //std::cout << "Fetching older sample in attempt to correct." << std::endl;

            td = (*tagBuffer)[c_buffer_size-1].unix_timestamp + (double)(*tagBuffer)[c_buffer_size-1].millis/1000 - (*tagBuffer)[c_buffer_size-2-correction_attempts-1].unix_timestamp - (double)(*tagBuffer)[c_buffer_size-2-correction_attempts-1].millis/1000;
            closestMultipleOfPeriod = m_args.tag_period*std::round(td/m_args.tag_period);
            tdoa = td - closestMultipleOfPeriod;
            correction_attempts++;
          }
                    

          if (goodSample) {
              //spew("Tdiff: %f, TDOA: %f", td,tdoa);
              double rdoa = m_c_speed*tdoa; // Range difference
              double rangeSNR = ((*tagBuffer)[c_buffer_size-1].snr - P[1])/P[0];
              arma::vec ranging = {rdoa,rangeSNR};
              double N,E;
              bool north;
              int zone;
              DUNE::Coordinates::UTM::fromWGS84((*tagBuffer)[c_buffer_size-1].lat,(*tagBuffer)[c_buffer_size-1].lon,&N,&E,&zone,&north);
              arma::vec position_current = {E, N, altitude};

              DUNE::Coordinates::UTM::fromWGS84((*tagBuffer)[c_buffer_size-2].lat,(*tagBuffer)[c_buffer_size-2].lon,&N,&E,&zone,&north);

              arma::vec position_previous = {E, N, altitude};
              double lati,longi;

              m_ekf.measurementStep(ranging,position_previous,position_current,rdoa);
              int pf_k = m_pf.update(ranging,position_previous,position_current);
              if (pf_k > 0) { // Checks if an estimate has been made
                DUNE::Coordinates::UTM::toWGS84(m_pf.x(1),m_pf.x(0),zone,north,&lati,&longi);
                IMC::RemoteSensorInfo tagPosition;
                tagPosition.lat = lati;
                tagPosition.lon = longi;
                tagPosition.alt = -m_pf.x(2);
                tagPosition.data = std::to_string(m_pf.x(0)) + std::to_string(m_pf.x(1)) + "," + std::to_string(m_pf.x(2));
                tagPosition.id = "TagPositionPF" + std::to_string(m_args.receiver_serial);
                dispatch(tagPosition);

                  /*std::ofstream logOutStream;
                  logOutStream.open(("log/pf2.log"), std::fstream::app);
                  if (logOutStream.good()) {
                    logOutStream.precision(15);
                      logOutStream << m_pf.x(0) << "," << m_pf.x(1) << "," << m_pf.x(2) << "," << DUNE::Math::Angles::degrees(lati) << "," << DUNE::Math::Angles::degrees(longi) << std::endl;
                      logOutStream.close();
                  } */

                spew("New PF Estimate: (E,N,A,La,Lo)= %.15f,%.15f,%.15f,%.15f, %.15f", m_pf.x(0), m_pf.x(1), m_pf.x(2),DUNE::Math::Angles::degrees(lati),DUNE::Math::Angles::degrees(longi));
                if(zone != m_args.locaton_utm_zone) {
                  war("Transformation with different UTM zone than specified: %d", zone);
                }
                if(north != m_args.location_northern_hemisphere) {
                  war("Not in spesifiec hemisphere. Coordinate transformation may produce wrong results");
                }
              }
            
          } else {
              err("Filter measurement could not be updated. No good TDOA value available.");
          }
        } else {
          spew("Waiting for measurements");
        }
      }

      //! Main loop.
      void
      onMain(void)
      {

        while (!stopping())
        {
          if(m_filter_timer.overflow()) {
            //inf("%f", m_filter_timer.getElapsed()); // To check for filter jitter
            m_filter_timer.reset();
            if(m_ekf.isInitialized()) {
              int ekf_k = m_ekf.predictionStep(); // Filter time update
              if (ekf_k > 0 ) { // Checks if an estimate has been made
                  double lati,longi;
                  DUNE::Coordinates::UTM::toWGS84(m_ekf.x(1),m_ekf.x(0),m_args.locaton_utm_zone,m_args.location_northern_hemisphere,&lati,&longi);
                  std::ofstream logOutStream;
                  logOutStream.open(("log/ekf2.log"), std::fstream::app);
                  if (logOutStream.good()) {
                    logOutStream.precision(15);
                      logOutStream << m_ekf.x(0) << "," << m_ekf.x(1) << "," << m_ekf.x(2) << "," << DUNE::Math::Angles::degrees(lati) << "," << DUNE::Math::Angles::degrees(longi) << std::endl;
                      logOutStream.close();
                  } 
                  // Send output to Neptus/DUNE log
                  IMC::RemoteSensorInfo tagPosition;
                  tagPosition.lat = lati;
                  tagPosition.lon = longi;
                  tagPosition.alt = -m_ekf.x(2);
                  tagPosition.data = std::to_string(m_ekf.x(0)) + std::to_string(m_ekf.x(1)) + "," + std::to_string(m_ekf.x(2));
                  tagPosition.id = "TagPositionEKF" + std::to_string(m_args.receiver_serial);
                  dispatch(tagPosition);

                  spew("New Kalman Estimate: (E,N,A,La,Lo)= %.15f,%.15f,%.15f,%.15f, %.15f", m_ekf.x(0), m_ekf.x(1), m_ekf.x(2),DUNE::Math::Angles::degrees(lati),DUNE::Math::Angles::degrees(longi));
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
