/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using namespace std;
void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  if(is_initialized)
    return;
  num_particles = 100;  // TODO: Set the number of particles
  // initialize all particles to first position
  //Particle particles;
  default_random_engine generator;
  normal_distribution<double> dist_x(x,std[0]);
  normal_distribution<double> dist_y(y,std[1]);
  normal_distribution<double> dist_theta(theta,std[2]);
  // addd noise to each particle
  for (int i=0;i<num_particles;i++)
  {
    Particle particle;
    particle.id = i;
    particle.x = dist_x(generator);
    particle.y = dist_y(generator);
    particle.theta = dist_theta(generator);
    particle.weight = 1.0;
    particles.push_back(particle);
  }
  is_initialized=true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  default_random_engine generator;
  normal_distribution<double> dist_x(0,std_pos[0]);
  normal_distribution<double> dist_y(0,std_pos[1]);
  normal_distribution<double> dist_theta(0,std_pos[2]);
  for(int i=0;i<num_particles;i++)
  {
    if(fabs(yaw_rate) >= 0.00001)
    {
      particles[i].x = particles[i].x + (velocity/yaw_rate*(sin(particles[i].theta + yaw_rate*delta_t)-sin(particles[i].theta)));
      particles[i].y = particles[i].y + (velocity/yaw_rate*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t)));
      particles[i].theta = particles[i].theta + yaw_rate*delta_t;
      // adding noise
    particles[i].x = particles[i].x + dist_x(generator);
    particles[i].y = particles[i].y + dist_y(generator);
    particles[i].theta = particles[i].theta + dist_theta(generator);
    }
    else
    {
      particles[i].x = particles[i].x + (velocity * delta_t * cos(particles[i].theta));
      particles[i].y = particles[i].y + (velocity * delta_t * sin(particles[i].theta));
    }
  }
}


void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
   for(unsigned int i = 0;i<observations.size();i++)
   {
     LandmarkObs obs = observations[i];
     double min_dist = numeric_limits<double>::max();
     int min_id = -9999999;
     
     for(unsigned int j = 0;j<predicted.size();j++)
     {
       LandmarkObs pred = predicted[j];
       double e_dist = dist(obs.x,obs.y,pred.x,pred.y);
       if(e_dist<min_dist)
       {
         min_dist = e_dist;
         min_id = pred.id;
       }
     }
     observations[i].id = min_id;
   } 
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  // values for gaussian equation for updating the weights
  double cnst = 1/(2*M_PI*std_landmark[0]*std_landmark[1]);
  double den1 = 2*std_landmark[0]*std_landmark[0];
  double den2 = 2*std_landmark[1]*std_landmark[1];
  
  for(int i=0; i<num_particles; i++)
    {
    vector<LandmarkObs> predicted;
        for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {
            double distance_x = map_landmarks.landmark_list[j].x_f - particles[i].x;
            double distance_y = map_landmarks.landmark_list[j].y_f - particles[i].y;
            if(fabs(distance_x)  <= sensor_range && fabs(distance_y) <=sensor_range)
            {
              predicted.push_back(LandmarkObs{map_landmarks.landmark_list[j].id_i, map_landmarks.landmark_list[j].x_f,
                                                       map_landmarks.landmark_list[j].y_f});
            }
        }
  
  // Transformation
  vector<LandmarkObs> transformed;
  for(unsigned int i=0;i<observations.size();i++)
  {
    double trans_x = cos(particles[i].theta)*observations[i].x - sin(particles[i].theta)*observations[i].y + particles[i].x;
    double trans_y = sin(particles[i].theta)*observations[i].x + cos(particles[i].theta)*observations[i].y + particles[i].y;    
    transformed.push_back(LandmarkObs{observations[i].id,trans_x,trans_y});
  }
  
   dataAssociation(predicted,transformed);
    particles[i].weight = 1.0;
     for(unsigned int i=0;i<transformed.size();i++)
     {
       double x_obs = transformed[i].x;
       double y_obs = transformed[i].y;
       double mu_x,mu_y;
       for(unsigned int j=0;j<predicted.size();j++)
       {
         if(predicted[j].id == transformed[i].id)
         {
           mu_x = predicted[j].x;
           mu_y = predicted[j].y;
         }        
       }
       double exponent =  (pow(mu_x - x_obs, 2) / den1) + (pow(mu_y - y_obs, 2) / den2);
       double localWeight = cnst * exp(-exponent);
       particles[i].weight = particles[i].weight * localWeight ;     
     }
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  vector<double> weights;
  double max_weight = numeric_limits<double>::min();
   for(int i=0;i<num_particles;i++)
   {
     weights.push_back(particles[i].weight);
     if(particles[i].weight>max_weight)
       max_weight = particles[i].weight;
   }
	
  vector<Particle> resampled;
  default_random_engine generator;
  uniform_real_distribution<double> dist_real(0.0,max_weight);
  uniform_int_distribution<int> dist_int(0,num_particles-1);
  int index = dist_int(generator);
  double beta =0.0;
  for(int i=0;i<num_particles;i++)
  {
    beta = beta + dist_real(generator)*2.0;
    while(weights[index] < beta)
    {
      beta = beta -weights[index];
      index = (index+1) % num_particles;
    }
    resampled.push_back(particles[index]);
  }
  particles = resampled;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}