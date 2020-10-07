#ifndef ENERGYFILESUPPORT_HH_
#define ENERGYFILESUPPORT_HH_

//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <sstream>
#include <limits.h>

#include <biu/Alphabet.hh>
#include <biu/DistanceEnergyFunction.hh>



/**
 * Initializes a distance interval based energy function from the
 * given stream.
 *
 * @param alph the alphabet to fill
 * @param energyMatrix the energy function contact table to fill
 * @param cAlphaDistScale a factor to scale the distances of the energy file to
 *        the corresponding distances in the used lattice model
 * @param input the input stream to read the information from
 * @return Status of reading : 0 = ok, -1 = read error
 */
int
initIntervalEnergyFunction(	biu::Alphabet *& alph,
							biu::DistanceEnergyFunction *& energyD,
							const double cAlphaDistScale,
							std::istream & input )
{
	assertbiu( input.good(), "cannot read from given input");

	// remove old energy function
	if (energyD != NULL) {
		delete energyD; energyD = NULL;
	}


	  // temp data structures
	std::string line;

	/////////   READ ALPHABET   //////////////////////////////////////////////

	if (alph != NULL) {
		delete alph; alph = NULL;
	}

	  // read first line that contains alphabet info
	  // get first line content
	getline( input, line, '\n' );

	if (line.size() == 0) {
		std::cerr <<"\n   ERROR : energy file : first line does not contain allowed alphabet !\n\n";
		return -1;
	}
	  // create new alphabet
	alph = new biu::Alphabet(line, 1);



	/////////   READ INTERVAL BOUNDS   ////////////////////////////////////////

	// holds the upper interval bounds of the energy function to create
	std::vector<double> dMax;

	  // read second line that contains interval bound info
	getline( input, line, '\n' );

	if (line.size() == 0) {
		std::cerr <<"\n   ERROR : energy file : second line does not contain interval bounds !\n\n";
		return -1;
	}

	std::istringstream dMaxLine(line);

	double lastBound = (double)-INT_MAX;
	while (dMaxLine.good()) {
		double nextBound = (double)UINT_MAX;
		  // get next interval bound
		dMaxLine >> nextBound;
		if (nextBound < lastBound) {
			std::cerr <<"\n   ERROR : energy file : the interval bounds in the second line are not increasing !\n\n";
			return -1;
		}
		if (nextBound != UINT_MAX) {
			  // store scaled upper interval bound
			dMax.push_back(nextBound*cAlphaDistScale);
			  // remember this as last bound found
			lastBound = nextBound;
		}
	}


	/////////   READ ENERGY TABLES   //////////////////////////////////////////

	// holds the energy table of the currently processed interval
	biu::EnergyMatrix eTable(	alph->getAlphabetSize(),
								alph->getAlphabetSize(),
								0.0 );

	  // initialise energy function
	energyD = new biu::IntervalEnergyFunction(alph);

	biu::IntervalEnergyFunction* energy = dynamic_cast<biu::IntervalEnergyFunction*>(energyD);

	for (size_t i=0; i<dMax.size(); i++) {
		if (input.bad()) {
			std::cerr <<"\n   ERROR : energy file : can not read " <<i+1 <<". energy matrix !\n";
			return -1;
		}
		  // read next table
		input >> eTable;
		  // push table and interval to energy function
		energy->addInterval( eTable, dMax[i]);
	}

	energyD = energy;

	return 0;
}

/**
 * Initializes the data needed for a contact based energy function from the
 * given stream.
 *
 * @param alph the alphabet to fill
 * @param energyMatrix the energy function contact table to fill
 * @param input the input stream to read the information from
 * @return Status of reading : 0 = ok, -1 = read error
 */
int
initContactEnergyFunction(	biu::Alphabet *& alph,
							biu::EnergyMatrix *& energyMatrix,
							std::istream & input )
{
	assertbiu( input.good(), "cannot read from given input");

	  // temp data structures
	std::string line;

	if (alph != NULL) {
		delete alph; alph = NULL;
	}
	if (energyMatrix != NULL) {
		delete energyMatrix; energyMatrix = NULL;
	}

	  // read first line that contains alphabet info
	  // get first line content
	getline( input, line, '\n' );

	if (line.size() == 0) {
		std::cerr <<"\n   ERROR : energy file : first line does not contain allowed alphabet !\n\n";
		return -1;
	}
	  // create new alphabet
	alph = new biu::Alphabet(line, 1);

	  // create and init new matrix
	energyMatrix = new biu::EnergyMatrix(	alph->getAlphabetSize(),
											alph->getAlphabetSize(),
											0.0 );

	  // fill matrix : read remaining istream into matrix
	input >> *energyMatrix;


	return 0;
}


//////////////////////////////////////////////////////////////////////////


#endif /*ENERGYFILESUPPORT_HH_*/
