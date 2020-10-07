#ifndef i_distribution
#define i_distribution
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include "zaprandom.h"
#include <vector>
using namespace std;

// General class, you can derive more versions of distributions. The only thing that matter is that
// a sirtribution should give you a random number !
struct distribution {
    virtual double getRandValue(){cerr << "Distribution class called without daughter class" << endl; return -1;}
    virtual ~distribution(){}
};

/** @Brief Class to generate random numbers from an arbitrary distribution,
 *  defined by a liste of values (classes) and their densities
 *  THE VALUES OF THE CLASSES SHOULD BE ORDERED INCREASING */
struct probaLawFromTable : public distribution {
    int size;
    vector<double> Xs;
    vector<double> NormalizedDensities;     // to probabilities
    vector<double> cumulatedValuesHigh;      // cumulated probabilities
    double normalizingCoeff;                // just for remembering the normalisation

    bool uniformizeInsideClasses;
    vector<double> lowBoundsXs;
    vector<double> highBoundsXs;

    virtual double getRandValue();
    int getRandClass();
    probaLawFromTable(vector<double> &_Xs, vector<double> &densities, bool _uniformizeInsideClasses);
    static double IndependentRandRealGen(); // If later someone wants to separate generators for seeding.
    virtual string print();

    virtual ~probaLawFromTable(){}
};

/** @Brief Class to test/control the behavior of probaLawFromTable.
 *  by being daughter cell, reimplements the getRandomValue function and
 *  does statistics on how the mother class behaves */
struct probaLawFromTableStoringResults : public probaLawFromTable
{
    vector<int> nbOfTimesPerClass;
    int totalNbEvents;

    probaLawFromTableStoringResults(vector<double> & _Xs, vector<double> & densities, bool _uniformizeInsideClasses);
    double getRandValue();
    double getFrequency(int classIndex);
    int fitDoubleInAClass(double val);
    void clearRecord();
    string print();
    static double TestprobaLawFromTable();

    static double InteractionTimeFromData();

    // not sure whether the destructor should call the mother destructors. Here, no effect, they are all empty
    ~probaLawFromTableStoringResults(){}
};

struct histogramFromDistrib{
    vector<double> lowBoundsXs;
    vector<double> highBoundsXs;
    vector<double> averageXs;
    vector<double> densities;
    double vmin;
    double vmax;
    double maxDens;
    string print(bool showClasses = false){
        stringstream ss;
        size_t S = lowBoundsXs.size();
        if(highBoundsXs.size() != S) cerr << "Wrong sizes in histogramFromDistrib " << endl;
        if(averageXs.size() != S) cerr << "Wrong sizes in histogramFromDistrib " << endl;
        if(densities.size() != S) cerr << "Wrong sizes in histogramFromDistrib " << endl;
        for(size_t i = 0; i < S; ++i){
            if(showClasses) ss << "[" << lowBoundsXs[i] << "," << highBoundsXs[i] << "]\t";
            ss << averageXs[i] << "\t" << densities[i] << endl;
        }
        return ss.str();
    }
    histogramFromDistrib(vector<double> listValues, int nbBins)
    {
        size_t S = listValues.size();
        if(S == 0) {cerr << "Histogram: empty data < 1\n"; return;}
        if(nbBins < 1) {cerr << "Histogram: nbBins < 1\n"; return;}
        lowBoundsXs.resize(nbBins, 0.);
        highBoundsXs.resize(nbBins, 0.);
        averageXs.resize(nbBins, 0.);
        densities.resize(nbBins, 0.);

        maxDens = 0;
        vmin = listValues[0];
        vmax = listValues[0];
        for(size_t i = 0; i < S; ++i){
            vmin = min(vmin, listValues[i]);
            vmax = max(vmax, listValues[i]);
        }
        double interval = (vmax - vmin) / (static_cast<double>(nbBins) - 1.);
        if(fabs(interval) < 1e-9)
        {
            lowBoundsXs.clear();
            lowBoundsXs.resize(1, vmin);
            highBoundsXs.clear();
            highBoundsXs.resize(1, vmax);
            averageXs.clear();
            averageXs.resize(1, 0.5*vmin + 0.5 * vmax );
            densities.clear();
            densities.resize(1, 1.0);
            cerr << "WRN : Histogram: too few stddev\n"; return;
        }
        for(int i = 0; i < nbBins; ++i){
            lowBoundsXs[i] = vmin + interval * ((double) i);
            highBoundsXs[i] = vmin + interval * ((double) (i + 1));
            averageXs[i] = 0.5*lowBoundsXs[i] + 0.5*highBoundsXs[i];
        }
        for(size_t i = 0; i < S; ++i){
            double val = listValues[i];
            int bin = (val - vmin) / interval;
            if(bin < 0) bin = 0;
            if(bin >= nbBins) bin = nbBins-1;
            densities[bin] += 1.0 / (double) S;
        }
        for(size_t i = 0; i < densities.size(); ++i){
            maxDens = max(maxDens, densities[i]);
        }
    }
    histogramFromDistrib(vector<double> listValues, vector<double> groupBoundaries)
    {
        size_t S = listValues.size();
        size_t GB = groupBoundaries.size();
        if(S == 0) {cerr << "ERR: Histogram: empty data < 1\n"; return;}
        if(GB == 0) {cerr << "ERR: Histogram: empty group boundaries < 1\n"; return;}

        size_t nbBins = GB+1;
        lowBoundsXs.resize(nbBins, 0.);
        highBoundsXs.resize(nbBins, 0.);
        averageXs.resize(nbBins, 0.);
        densities.resize(nbBins, 0.);
        for(size_t i = 1; i < GB; ++i){
            lowBoundsXs[i] = groupBoundaries[i-1];
            highBoundsXs[i] = groupBoundaries[i];
            averageXs[i] = (highBoundsXs[i] + lowBoundsXs[i]) / 2.;
            if(groupBoundaries[i] < groupBoundaries[i-1] - 1e-6) cerr << "ERR: Histogram, the group boundaries should be increasing : " << groupBoundaries[i-1] << "<->" << groupBoundaries[i] << endl;
        }
        // outsider values
        lowBoundsXs[0] = -1e6;
        highBoundsXs[0] = groupBoundaries[0];
        lowBoundsXs[GB] = groupBoundaries[GB-1];
        highBoundsXs[GB] = 1e6;

        maxDens = 0;
        vmin = listValues[0];
        vmax = listValues[0];
        for(size_t i = 0; i < S; ++i)// this will be useless
        {
            vmin = min(vmin, listValues[i]);
            vmax = max(vmax, listValues[i]);
        }
        for(size_t i = 0; i < S; ++i)
        {
            double val = listValues[i];
            bool found = false;
            size_t whichClass = 0;

            while(!found && (whichClass < nbBins)){
                if(val < highBoundsXs[whichClass]){
                    found = true;
                    densities[whichClass] += 1.0 / (double) S;
                } else {
                    whichClass++;
                }
            }
        }
        for(size_t i = 0; i < densities.size(); ++i)
        {
            maxDens = max(maxDens, densities[i]);
        }
    }
};

struct ministat
{

    ministat()
    {
        clear();
    }
    void clear()
    {
        N=0;
        sum=0;
        sumsq=0;
        values.clear();
        vmin = +1e30;
        vmax = -1e30;
    }
    void add(double v)
    {
        sum += v;
        sumsq += v*v;
        N++;
        values.push_back(v);
        vmin = min(vmin, v);
        vmax = max(vmax, v);
    }
    int N;
    vector<double> values;
    double sum;
    double sumsq;
    double vmax; // new added 2019-04-08
    double vmin;
    double average()
    {
        return sum / (double) N;
    }
    double stddev()
    {
        return (double) sqrt( (double) (sumsq / (double) N) - (sum / (double) N) * (sum / (double) N));
    }
    string print()
    {
        stringstream res;
        res << "[" << N << "]";
        for(int i = 0; i < N; ++i)
        {
            res << "\t" << values[i];
        }
        return res.str();
    }
    string sumUp(){
        stringstream res;
        res << N << "\tvals, [\t" << vmin << "\t" << vmax << "\t] - avg=\t" << average() << "\tstd=\t" << stddev();
        return res.str();
    }
};

/*void example()
{
    ministat a;
a.add(15);
a.add(20);
cout << a.average() << a.stddev() << endl;
a.clear();
}*/
#endif
