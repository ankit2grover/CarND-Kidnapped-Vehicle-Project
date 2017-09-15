/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // Initialize number of particles
  num_particles = 21;
  // Generator to generate random values
  default_random_engine gen;
  double std_x = std[0];
  double std_y = std[1];
  double std_theta = std[2];
  
  // This line creates a normal (Gaussian) distribution
  normal_distribution<double> dist_x(0, std_x);
  normal_distribution<double> dist_y(0, std_y);
  normal_distribution<double> dist_theta(0, std_theta);
  particles.resize(num_particles);
  weights.resize(num_particles);
  // Intialize particles values
  for (unsigned int i = 0; i < num_particles; i++) {
    particles[i].id = i;
    particles[i].x = x;
    particles[i].y = y;
    particles[i].theta = theta;
    // Adding GPS noise on every particle
    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);
    particles[i].weight = 1.0;
  }
  is_initialized = true;
  
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  // Adding prediction noise on every particles using guassian distribution
  default_random_engine gen;
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);
  for (unsigned int i = 0; i < num_particles; i++) {
    // If yaw rate is 0 then use below equation
    if (abs(yaw_rate) == 0) {
      particles[i].x += velocity * delta_t * cos(particles[i].theta);
      particles[i].y += velocity * delta_t * sin(particles[i].theta);
    } else {
      // Predict particle positions if yaw rate is not 0
      particles[i].x +=  (velocity / yaw_rate)*
        (sin(particles[i].theta + (yaw_rate * delta_t))- sin(particles[i].theta));
      
      particles[i].y +=   (velocity / yaw_rate)*
        (cos(particles[i].theta) - cos(particles[i].theta + (yaw_rate * delta_t)));
      
      particles[i].theta +=  (yaw_rate * delta_t);
    }
    
    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
  // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
  //   observed measurement to this particular landmark.
  // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
  //   implement this method and use it as a helper during the updateWeights phase.
  /** Didn't find that method much useful, leaving it empty **/
}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   std::vector<LandmarkObs> observations, Map map_landmarks) {
  // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
  //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
  // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
  //   according to the MAP'S coordinate system. You will need to transform between the two systems.
  //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
  //   The following is a good resource for the theory:
  //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
  //   and the following is a good resource for the actual equation to implement (look at equation
  //   3.33
  //   http://planning.cs.uiuc.edu/node99.html
  
  // Update weight constants
  const double denominator0 = (2 * std_landmark[0] * std_landmark[0]);
  const double denominator1 = (2 * std_landmark[1] * std_landmark[1]);
  const double guass_norm = 1/(2 * M_PI * std_landmark[0] * std_landmark[1]);
  // Iterate over every particle to associate it with correct landmark observations.
  for (unsigned int i = 0; i < num_particles; i++) {
    // Iterate over every particle
    double theta = particles[i].theta;
    double probability = 1.0;
    
    for (unsigned int j = 0; j < observations.size(); j++) {
      // Stores vector of landmarks distance from the observation
      std::vector<double> landmark_obs_dist(map_landmarks.landmark_list.size());
      // Transform landmark observations from vehicle coordinate system
      // to map coordinate system with respect to every particle.
      double x_trans_obs = particles[i].x + (cos(theta) * observations[j].x) - (sin(theta) * observations[j].y);
      double y_trans_obs = particles[i].y + (sin(theta) * observations[j].x) + (cos(theta) * observations[j].y);
      // Iterate over every landmark and stores its distance in vector
      for (unsigned int k = 0; k < map_landmarks.landmark_list.size(); k++) {
        double landmark_x = map_landmarks.landmark_list[k].x_f;
        double landmark_y = map_landmarks.landmark_list[k].y_f;
        double particle_landmark_dist = dist(particles[i].x, particles[i].y, landmark_x, landmark_y);
        if (particle_landmark_dist <= sensor_range) {
          double calc_dist = dist(x_trans_obs, y_trans_obs, landmark_x, landmark_y);
          landmark_obs_dist[k] = calc_dist;
        } else {
          // Set value as highest to make sure that 0 values are not stored.
          landmark_obs_dist[k] = 999999.0;
        }
      }
      // Find min index of the landmark that is associated with observation using neighrest neighbour algorithm.
      int min_pos = distance(landmark_obs_dist.begin(), min_element(landmark_obs_dist.begin(), landmark_obs_dist.end()));
      // Calculate exponent term
      double exponent =
        (pow((x_trans_obs - map_landmarks.landmark_list[min_pos].x_f), 2) / denominator0) +
        (pow((y_trans_obs - map_landmarks.landmark_list[min_pos].y_f), 2) / denominator1);
      probability = probability * guass_norm * exp(-exponent);
      
    }
    // Update the weights probability.
    particles[i].weight = probability;
    weights[i] = particles[i].weight;
  }
}

void ParticleFilter::resample() {
  // Resample particles with replacement with probability proportional to their weight.
  // Resample the particles using resampling wheel algorithm explained by Sebastian.
  default_random_engine gen;
  std::vector<Particle> new_particles(num_particles);
  uniform_int_distribution<int> samples_distr(0, num_particles -1);
  int index = samples_distr(gen);
  double beta = 0.0;
  double max_weight = *max_element(std::begin(weights), std::end(weights));
  uniform_real_distribution<double> weights_distr(0, max_weight);
  for (unsigned int i= 0; i < num_particles; i++) {
    beta += 2.0 * weights_distr(gen);
    while(beta > weights[index]) {
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    new_particles[i] = particles[index];
  }
  particles = new_particles;
  
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
  //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  
  //Clear the previous associations
  particle.associations.clear();
  particle.sense_x.clear();
  particle.sense_y.clear();
  
  particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;
  
 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
  vector<int> v = best.associations;
  stringstream ss;
  copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
string ParticleFilter::getSenseX(Particle best)
{
  vector<double> v = best.sense_x;
  stringstream ss;
  copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
string ParticleFilter::getSenseY(Particle best)
{
  vector<double> v = best.sense_y;
  stringstream ss;
  copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
