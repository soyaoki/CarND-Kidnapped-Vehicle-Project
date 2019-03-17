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
    particles[i].weight = 1.0;
    weights.push_back(particles[i].weight);
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
  std::normal_distribution<double> dist_x(0.0, std_pos[0]);
  std::normal_distribution<double> dist_y(0.0, std_pos[1]);
  std::normal_distribution<double> dist_theta(0.0, std_pos[2]);
  
  for (int i = 0; i < num_particles; i++){
    
    if (fabs(yaw_rate) < 0.00001) {  
      particles[i].x += velocity * delta_t * cos(particles[i].theta) + dist_x(generator);
      particles[i].y += velocity * delta_t * sin(particles[i].theta) + dist_y(generator);
    } else {
      particles[i].x += (velocity / yaw_rate) * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta)) + dist_x(generator);
      particles[i].y += (velocity / yaw_rate) * (-cos(particles[i].theta + yaw_rate * delta_t) + cos(particles[i].theta)) + dist_y(generator);
      particles[i].theta += (yaw_rate * delta_t) + dist_theta(generator);
    }
  }
  std::cout << "Predicted!" << std::endl;
  std::cout << "Pred x :" << particles[0].x << std::endl;
  std::cout << "Pred y :" << particles[0].y << std::endl;
  std::cout << "Pred theta :" << particles[0].theta << std::endl;
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
  
  vector<LandmarkObs> observed_sorted;
  for (int ii = 0; ii < predicted.size(); ii++){
    // Each measured distance
    double distance_error_min = 100;
    LandmarkObs obs;
    for (int jj = 0; jj < observations.size(); jj++){
      double d = dist(predicted[ii].x, predicted[ii].y, observations[jj].x, observations[jj].y);
      if (d <= distance_error_min){
        obs.x = observations[jj].x;
        obs.y = observations[jj].y;
        distance_error_min = d;
      } 
    }
    obs.id = predicted[ii].id;
    observed_sorted.push_back(obs);
  }
  observations.clear();
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
  
  double weight_normalizer = 0.0;
  
  // Each Paticles
  for (int i = 0; i < num_particles; i++){
    // calcurate distances between the particle and each landmarks.
    vector<LandmarkObs> predicted;
    for (int j = 0; j < map_landmarks.landmark_list.size(); j++){
      //if (dist(particles[i].x, particles[i].y, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f) <= sensor_range){
      if ( (fabs(particles[i].x - map_landmarks.landmark_list[j].x_f) <= sensor_range ) && ( fabs(particles[i].y - map_landmarks.landmark_list[j].y_f) <= sensor_range ) ){
        LandmarkObs pred;
        pred.id = map_landmarks.landmark_list[j].id_i;
        pred.x = map_landmarks.landmark_list[j].x_f;
        pred.y = map_landmarks.landmark_list[j].y_f;
        predicted.push_back(pred);
      }
    }
    // transfar observed (measured) distace to from vehicle coordinate to global coordinate
    vector<LandmarkObs> observed;
    for (int k = 0; k < observations.size(); k++){
      LandmarkObs obs;
      obs.x = particles[i].x + ( cos(particles[i].theta) * observations[k].x ) - ( sin(particles[i].theta) * observations[k].y );
      obs.y = particles[i].y + ( sin(particles[i].theta) * observations[k].x ) + ( cos(particles[i].theta) * observations[k].y );
      obs.id =  k;
      observed.push_back(obs);
    }
    
    // wight update
    dataAssociation(predicted, observed); // resize to same size

    for (int n = 0; n < predicted.size(); n++){
      std::cout << "predicted : " << predicted[n].x << "," << predicted[n].y << ", " << predicted[n].id << std::endl;
      std::cout << "observed : " << observed[n].x << "," << observed[n].y << ", " << observed[n].id << std::endl;
    }
    
    particles[i].weight = 1.0;
    
    for (int l = 0; l < predicted.size(); l++){
      //std::cout << "weight = " << particles[i].weight << std::endl;
      particles[i].weight = particles[i].weight * multiv_prob(std_landmark[0], std_landmark[1], observed[l].x, observed[l].y, predicted[l].x, predicted[l].y);
      std::cout << "weight : " << multiv_prob(std_landmark[0], std_landmark[1], observed[l].x, observed[l].y, predicted[l].x, predicted[l].y) << std::endl;
    }
    std::cout << "No." << i << " final weight : " << particles[i].weight << std::endl;
    weight_normalizer += particles[i].weight;
  }
  
  std::cout << "weight normalizer : " << weight_normalizer << std::endl;
  
  for (int i = 0; i < particles.size(); i++) {
    //particles[i].weight = particles[i].weight / weight_normalizer;
    weights[i] = particles[i].weight;
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
  std::default_random_engine generator;
  std::uniform_int_distribution<int> particle_index(0, num_particles - 1);
  int index = particle_index(generator);
  double beta = 0.0;
  double mw = -1.0;
  
    for (int i = 0; i < num_particles; i++){
      if (weights[i] > mw){
        mw = weights[i];
      }
    }
  
    for (int j = 0; j < num_particles; j++){
      std::uniform_real_distribution<double> random_weight(0.0, 2*mw);
      beta += random_weight(generator);

      while ( beta > weights[index] ){
        beta -= weights[index];
        index = (index + 1) % num_particles;
      }
      new_p.push_back(particles[index]);
    }
    particles = new_p;
  	std::cout << "Resampled! ( N = " << particles.size() << ")" << std::endl;
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