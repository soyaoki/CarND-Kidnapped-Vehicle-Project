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
/**
     * init Initializes particle filter by initializing particles to Gaussian
     *   distribution around first position and all the weights to 1.
     * @param x Initial x position [m] (simulated estimate from GPS)
     * @param y Initial y position [m]
     * @param theta Initial orientation [rad]
     * @param std[] Array of dimension 3 [standard deviation of x [m], 
     *   standard deviation of y [m], standard deviation of yaw [rad]]
     */
void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 20;  // TODO: Set the number of particles
  // each particles have int id, double x, double y, double theta, double weight, 
  //                       std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y;
  Particle p;
  std::default_random_engine generator;
  std::normal_distribution<double> dist_x(0.0, std[0]);
  std::normal_distribution<double> dist_y(0.0, std[1]);
  std::normal_distribution<double> dist_theta(0.0, std[2]);
  for (int i = 0; i < num_particles; i++){
    particles.push_back(p);
    particles[i].id = i;
    particles[i].x = x + dist_x(generator);
    particles[i].y = y + dist_y(generator);
    particles[i].theta = theta + dist_theta(generator);
    particles[i].weight = 1;
    }
  std::cout << "Initialized!" << std::endl;
  is_initialized = true;
}

/**
   * prediction Predicts the state for the next time step
   *   using the process model.
   * @param delta_t Time between time step t and t+1 in measurements [s]
   * @param std_pos[] Array of dimension 3 [standard deviation of x [m], 
   *   standard deviation of y [m], standard deviation of yaw [rad]]
   * @param velocity Velocity of car from t to t+1 [m/s]
   * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
   */
void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  std::default_random_engine generator;
  for (int i = 0; i < num_particles; i++){
    std::normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
    std::normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
    std::normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);
    //particles[i].id = i;
    particles[i].x = particles[i].x + velocity / yaw_rate * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta)) + dist_x(generator);
    particles[i].y = particles[i].y + velocity / yaw_rate * (-cos(particles[i].theta + yaw_rate * delta_t) + cos(particles[i].theta)) + dist_y(generator);
    particles[i].theta = particles[i].theta + yaw_rate * delta_t + dist_theta(generator);
    }
  std::cout << "Predicted!" << std::endl;
}

/**
   * dataAssociation Finds which observations correspond to which landmarks 
   *   (likely by using a nearest-neighbors data association).
   * @param predicted Vector of predicted landmark observations
   * @param observations Vector of landmark observations
   */
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
  
  //LandmarkObs has
//  int id;     // Id of matching landmark in the map.
//  double x;   // Local (vehicle coords) x position of landmark observation [m]
//  double y;   // Local (vehicle coords) y position of landmark observation [m]
  
  vector<LandmarkObs> observed_sorted;
  //std::cout << "predicted size : " << predicted.size() << std::endl;
  //std::cout << "observed size : " << observations.size() << std::endl;
  // Each predicted distance
  for (int i = 0; i < predicted.size(); i++){
    // Each measured distance
    double distance_error_min = 50;
    LandmarkObs obs;
    for (int j = 0; j < observations.size(); j++){
      if (dist(predicted[i].x, predicted[i].y, observations[j].x, observations[j].y) < distance_error_min){
        obs.x = observations[j].x;
        obs.y = observations[j].y;
      } 
    }
    obs.id = predicted[i].id;
    observed_sorted.push_back(obs);
  }
  observations = observed_sorted;
}

/**
   * updateWeights Updates the weights for each particle based on the likelihood
   *   of the observed measurements. 
   * @param sensor_range Range [m] of sensor
   * @param std_landmark[] Array of dimension 2
   *   [Landmark measurement uncertainty [x [m], y [m]]]
   * @param observations Vector of landmark observations
   * @param map Map class containing map landmarks
   */
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

  //for (int i =0; i< map_landmarks.landmark_list.size(); i++){
  //  std::cout << map_landmarks.landmark_list[i].id_i << std::endl;
  //  std::cout << map_landmarks.landmark_list[i].x_f << std::endl;
  //  std::cout << map_landmarks.landmark_list[i].y_f << std::endl;
  //}
  
  // Each Paticles
  for (int i = 0; i < num_particles; i++){
    // calcurate distances between the particle and each landmarks.
    vector<LandmarkObs> predicted;
    for (int j = 0; j < map_landmarks.landmark_list.size(); j++){
      if (dist(particles[i].x, particles[i].y, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f) < sensor_range){
        LandmarkObs pred;
        pred.id = map_landmarks.landmark_list[j].id_i;
        pred.x = map_landmarks.landmark_list[j].x_f - particles[i].x;
        pred.y = map_landmarks.landmark_list[j].y_f - particles[i].y;
        predicted.push_back(pred);
      }
    }
    // transfar observed (measured) distace
    vector<LandmarkObs> observed;
    for (int k = 0; k < observations.size(); k++){
      LandmarkObs obs;
      obs.x= cos(-particles[i].theta) * observations[k].x - sin(-particles[i].theta) * observations[k].y - particles[i].x;
      obs.y = sin(-particles[i].theta) * observations[k].x + cos(-particles[i].theta) * observations[k].y - particles[i].y;
      observed.push_back(obs);
    }
    // wight update
    dataAssociation(predicted, observed); // resize to same size
    //std::cout << "predicted size : " << predicted.size() << std::endl;
    //std::cout << "observed size : " << observed.size() << std::endl;
    
    for (int l = 0; l < predicted.size(); l++){
      particles[i].weight = particles[i].weight * multiv_prob(std_landmark[0], std_landmark[1], observations[l].x, observations[l].y, predicted[l].x, predicted[l].y);
    }
  }
}

double ParticleFilter::multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs, double mu_x, double mu_y) {
  // calculate normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

  // calculate exponent
  double exponent;
  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
    
  // calculate weight using normalization terms and exponent
  double weight;
  std::cout << gauss_norm << "," << exponent << std::endl;
  weight = gauss_norm * exp(-exponent);
    
  return weight;
}

/**
   * resample Resamples from the updated set of particles to form
   *   the new set of particles.
   */
void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
    vector<Particle> new_p;
    int index = random() / RAND_MAX * num_particles;
    double beta = 0.0;
    double mw = 0;
  
    for (int i = 0; i < num_particles; i++){
      if (particles[i].weight > mw){
        mw = particles[i].weight;
      }
    }
  
    for (int j = 0; j < num_particles; j++){
      beta += random() / RAND_MAX * 2.0 * mw;
      while ( beta > particles[index].weight ){
        beta -= particles[index].weight;
        index = (index + 1) % num_particles;
      }
      new_p.push_back(particles[index]);
    }
    particles = new_p;
}

/**
   * Set a particles list of associations, along with the associations'
   *   calculated world x,y coordinates
   * This can be a very useful debugging tool to make sure transformations 
   *   are correct and assocations correctly connected
   */
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