// version 21-12-2017
#ifdef _WIN32
#include <windows.h>
#include <time.h>
#elif __APPLE__
#elif __linux__
#endif

#include "zaprandom.h"
#include <cmath>
#include <random>
#include <iostream>
#include <algorithm> // for sghuffle

std::mt19937* random::generator = nullptr;
bool random::initialized = false;

void random::initialize(unsigned int init_seed){
    if(initialized) std::cerr << "WRN: random::initialize() is being called two times or more" << endl;
    initialized = true;
    std::mt19937::result_type seed = time(0);
    if(init_seed != 0){seed = init_seed;}
    generator = new std::mt19937(seed);
    std::cout << "Seed used for optimization(std::mt19937) : \t" << seed << endl;
}

int random::uniformInteger(int a, int b){
    if(!initialized) initialize();
    std::uniform_int_distribution<int> distribution = std::uniform_int_distribution<int>(a, b);
    return distribution(*generator);
}

double random::uniformDouble(double a, double b){
    if(!initialized) initialize();
    std::uniform_real_distribution<double> distribution = std::uniform_real_distribution<double>(a, b);
    return distribution(*generator);
}

double random::normal(double mean , double SD){
    if(!initialized) initialize();
    std::normal_distribution<double> normal = std::normal_distribution<double>(mean, SD);
    return normal(*generator);
}

double random::bernouilli(double p){
    if(!initialized) initialize();
    std::bernoulli_distribution bernouilli = std::bernoulli_distribution(p);
    return bernouilli(*generator);
}

int random::binomial(int nb, double p){
    if(!initialized) initialize();
    std::binomial_distribution<int> binomial = std::binomial_distribution<int>(nb, p);
    return binomial(*generator);
}

double random::cauchy(double m, double sig){
    if(!initialized) initialize();
    std::cauchy_distribution<double> cauchy = std::cauchy_distribution<double>(m, sig);
    return cauchy(*generator);
}

double random::gamma(double shape, double scale){
    if(!initialized) initialize();
    std::gamma_distribution<double> gamma = std::gamma_distribution<double>(shape, scale);
    return gamma(*generator);
}

double random::poisson(double lambda){
    if(!initialized) initialize();
    std::poisson_distribution<int> poisson = std::poisson_distribution<int>(lambda);
    return poisson(*generator);
}

double random::geometric(double lambda){
    if(!initialized) initialize();
    std::geometric_distribution<int> geometric = std::geometric_distribution<int>(lambda);
    return geometric(*generator);
}

double random::exponential(double lambda){
    if(!initialized) initialize();
    std::exponential_distribution<double> exponential = std::exponential_distribution<double>(lambda);
    return exponential(*generator);
}

double random::logNormal(double mean, double sigma){
    if(!initialized) initialize();
    std::lognormal_distribution<double> logNormal = std::lognormal_distribution<double>(mean, sigma);
    return logNormal(*generator);
}

double random::biModal(double mn, double sig, double weight, double mn2, double sig2){
    if(uniformDouble(0,1) < weight) {   // will initialize if neccessary
        return normal(mn, sig);
    } else {
        return normal(mn2, sig2);
    }
}



