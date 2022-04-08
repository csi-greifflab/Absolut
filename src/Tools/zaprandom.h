#ifndef zap_random_h
#define zap_random_h
#include <random>
#include <algorithm> // for shuffle²
using namespace std;

struct random
{
    static std::mt19937* generator;
    static bool initialized;

    static void   initialize    (unsigned int init_seed=0);
    static int    uniformInteger(int, int);	
	static size_t uniformUInteger(size_t a, size_t b);
    static double uniformDouble (double, double);
    static double normal        (double, double);
    static double bernouilli    (double p);
    static int    binomial      (int nb, double p);
    static double cauchy        (double m, double sig);
    static double gamma         (double shape, double scale);
    static double poisson       (double lambda);
    static double geometric     (double lambda);
    static double exponential   (double lambda);
    static double logNormal     (double mean, double sigma);
    static double biModal       (double mn, double sig, double weight, double mn2, double sig2);

    template <typename T>
    static void shuffle(vector<T> &v){ // template functions should be in the .h
        if(!initialized) initialize();
        std::shuffle(v.begin(), v.end(), *generator);
    }
    static string shuffle(string str){ // template functions should be in the .h
        if(!initialized) initialize();
        std::vector<char> data(str.begin(), str.end());
        std::shuffle(data.begin(), data.end(), *generator);
        return string(data.begin(), data.end());
    }
};

#endif
