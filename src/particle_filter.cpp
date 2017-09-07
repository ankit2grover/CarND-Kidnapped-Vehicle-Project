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
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 50;
	default_random_engine gen;
	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];

	// This line creates a normal (Gaussian) distribution for GPS x
	normal_distribution<double> dist_x(0, std_x);
	normal_distribution<double> dist_y(0, std_y);
	normal_distribution<double> dist_theta(0, std_theta);
  cout << "Initial values: x: " << x << " y: " << y << endl;
	for (unsigned int i = 0; i < num_particles; i++) {
		Particle particle;
		particle.id = i;
    particle.x = x;
    particle.y = y;
    particle.theta = theta;
		particle.x += dist_x(gen);
		particle.y += dist_y(gen);
		particle.theta += dist_theta(gen);
		particles.push_back(particle);
	}
  is_initialized = true;
	
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  default_random_engine gen;
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);
	for (unsigned int i = 0; i < num_particles; i++) {
		Particle particle = particles[i];
		if (fabs(yaw_rate) < 0.001) {
      particle.x += velocity * delta_t * cos(particle.theta);
      particle.y += velocity * delta_t * sin(particle.theta);
    } else {
      particle.x = particle.x + (velocity / yaw_rate)*
          (sin(particle.theta + (yaw_rate * delta_t))- sin(particle.theta));
		
      particle.y = particle.y + (velocity / yaw_rate)*
          (cos(particle.theta) - cos(particle.theta + (yaw_rate * delta_t)));

      particle.theta = particle.theta + (yaw_rate * delta_t);
    }
    
    particle.x += dist_x(gen);
    particle.y += dist_y(gen);
    particle.theta += dist_theta(gen);
    cout << "Particle i: " << i << "x: "<< particle.x << endl;

	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
  for (unsigned int i =0; i < observations.size(); i++) {
    // Associate observation by default with 0 index landmark
    observations[i].id = predicted[0].id;
    double observation_x = observations[i].x;
    double observation_y = observations[i].y;
    double dist_min = dist(observation_x, observation_y, predicted[0].x, predicted[0].y);
    for (unsigned int j = 1; j < predicted.size(); j++) {
      double landmark_x = predicted[j].x;
      double landmark_y = predicted[j].y;
      double calc_dist = dist(observation_x, observation_y, landmark_x, landmark_y);
      if (calc_dist < dist_min) {
        dist_min = calc_dist;
        observations[i].id = predicted[j].id;
      }
    }
  }
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

	// Iterate over every particle to associate it with correct landmark observations.
	for (unsigned int i = 0; i < num_particles; i++) {
		// Iterate over every particle
		Particle particle = particles[i];
		double theta = particle.theta;
    
    // Transformed observations from vehicle coordinates to map coordinates.
    std::vector<LandmarkObs> transformed_Obs;
    
    // Map Landmarks vector that will be with in sensor range.
		std::vector<LandmarkObs> predictions;
    
    for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {
      double landmark_x = map_landmarks.landmark_list[j].x_f;
      double landmark_y = map_landmarks.landmark_list[j].y_f;
      int landmark_id = map_landmarks.landmark_list[j].id_i;
      double diff_x= fabs(landmark_x - particle.x);
      double diff_y = fabs(landmark_y - particle.y);
      if (diff_x <= sensor_range && diff_y <= sensor_range) {
        LandmarkObs map_Obs;
        map_Obs.x = landmark_x;
        map_Obs.y = landmark_y;
        map_Obs.id = landmark_id;
        predictions.push_back(map_Obs);
      }
      
    }

    // Transform landmark observations from vehicle coordinate system
    // to map coordinate system with respect to every particle.
		for (unsigned int j = 0; j < observations.size(); j++) {
			LandmarkObs transformed;
			double x_trans_obs = particle.x + (cos(theta) * observations[j].x) - (sin(theta) * observations[j].y);
			double y_trans_obs = particle.y + (sin(theta) * observations[j].x) + (cos(theta) * observations[j].y);
			transformed.id = observations[j].id;
			transformed.x = x_trans_obs;
			transformed.y = y_trans_obs;
			transformed_Obs.push_back(transformed);
		}
    // Particle observations are associated with landmarks using nearest neighbors algorithm.
		this->ParticleFilter::dataAssociation(predictions, transformed_Obs);
    
    
    // Update the weight of each particle
    double associate_landmark_x = 0.0;
    double associate_landmark_y = 0.0;
    double guass_norm = 1/(2 * M_PI * std_landmark[0] * std_landmark[1]);
    if (guass_norm < 0.001) {
      guass_norm = 0.001;
    }
    double particle_weight = 1.0;
    std::vector<int> associations;
    std::vector<double> sense_x;
    std::vector<double> sense_y;
    
    for (unsigned int j=0; j < transformed_Obs.size(); j++) {
      int associate_landmark_id = transformed_Obs[j].id;
      for (unsigned int k = 0; k < predictions.size(); k++) {
        if (associate_landmark_id == predictions[k].id) {
          associations.push_back(associate_landmark_id);
          associate_landmark_x = predictions[k].x;
          associate_landmark_y = predictions[k].y;
          sense_x.push_back(transformed_Obs[j].x);
          sense_y.push_back(transformed_Obs[j].y);
          break;
        }
      }
      
      double exponent =
        (pow((transformed_Obs[j].x - associate_landmark_x), 2) / (2 * std_landmark[0] * std_landmark[0])) +
        (pow((transformed_Obs[j].y - associate_landmark_y), 2) / (2 * std_landmark[1] * std_landmark[1]));
      particle_weight = particle_weight * guass_norm * exp(-exponent);
    }
    particle.weight = particle_weight;
		weights.push_back(particle_weight);
    this->SetAssociations(particle, associations, sense_x, sense_y);
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  default_random_engine gen;
	std::vector<Particle> new_particles;
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
		new_particles.push_back(particles[index]);
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
