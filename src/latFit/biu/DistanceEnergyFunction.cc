// $Id: DistanceEnergyFunction.cc,v 1.2 2016/08/08 12:41:58 mmann Exp $

#include "biu/DistanceEnergyFunction.hh"
#include <limits.h>

namespace biu
{
	
	DistanceEnergyFunction::~DistanceEnergyFunction()
	{
	}

} // biu


namespace biu
{

	ContactEnergyFunction::ContactEnergyFunction(
				const Alphabet* const _alphabet, 
				const EnergyMatrix* const _energyMat,
				const LatticeModel* const _lattice) 
	 :	biu::DistanceEnergyFunction(),
		alphabet(_alphabet), 
		energyMat(_energyMat),
		lattice(_lattice),
		firstBaseVecLength(lattice==NULL?0.0:lattice->getNeighborhood().getElementByIndex(0).distance(IntPoint(0,0,0)))
	{
		assertbiu( alphabet != NULL, "no alphabet given (NULL)");
		assertbiu( energyMat != NULL, "no energy matrix given (NULL)");
		assertbiu( lattice != NULL, "no lattice model given (NULL)");
		assertbiu(	alphabet->getAlphabetSize() == energyMat->numRows() 
				&& alphabet->getAlphabetSize() == energyMat->numColumns(),
				"energy matrix has the wrong size for this alphabet size");
	}
	
	ContactEnergyFunction::~ContactEnergyFunction()
	{
	}
	
	double 
	ContactEnergyFunction::getContactEnergy( 
							const Alphabet::AlphElem& first, 
							const Alphabet::AlphElem& second) const {
		
		return energyMat->at(	alphabet->getIndex(first),
								alphabet->getIndex(second) );
	}
	
	double 
	ContactEnergyFunction::getEnergy(	const Alphabet::AlphElem& seq_i, 
										const Alphabet::AlphElem& seq_j,
										const double & distance ) const
	{
		if (distance == firstBaseVecLength ) {
			return getContactEnergy(seq_i, seq_j);
		}
		return 0.0;
	}
								
	double 
	ContactEnergyFunction::getEnergy(	const Alphabet::AlphElem & seq_i, 
										const Alphabet::AlphElem & seq_j,
										const IntPoint & cor_i,
										const IntPoint & cor_j  ) const
	{
		if (lattice->areNeighbored( cor_i, cor_j ) ) {
			return getContactEnergy( seq_i, seq_j );
		}
		return 0.0;
	}
								
	double 
	ContactEnergyFunction::getEnergy(	const Alphabet::AlphElem & seq_i, 
										const Alphabet::AlphElem & seq_j,
										const DblPoint & cor_i,
										const DblPoint & cor_j  ) const
	{
		  // use distance evaluation function
		return this->getEnergy(seq_i,seq_j,cor_i.distance(cor_j));
	}

	bool 
	ContactEnergyFunction::operator == (const DistanceEnergyFunction& ef2) const {
		if (this == &ef2)
			return true;
		  // cast
		const ContactEnergyFunction* cef2 = dynamic_cast<const ContactEnergyFunction*> (&ef2);
		  // evaluate
		return	(alphabet == cef2->alphabet || *alphabet == *(cef2->alphabet))
				&& (energyMat == cef2->energyMat || *energyMat == *(cef2->energyMat))
				&& (lattice == cef2->lattice || *lattice == *(cef2->lattice));
	}
	
	bool 
	ContactEnergyFunction::operator != (const DistanceEnergyFunction& ef2) const {
		return	!(this->operator ==(ef2));
	}
	

} // namespace biu


#include <algorithm>

namespace biu
{

	IntervalEnergyFunction::IntervalEnergyFunction(	const Alphabet* const alph_)
	 :	biu::DistanceEnergyFunction(), 
	 	alphabet(alph_),
	 	energyMat(),
	 	intervalMax()
	{
		
	}
	
	IntervalEnergyFunction::IntervalEnergyFunction( const biu::IntervalEnergyFunction & toCopy)
	 :	biu::DistanceEnergyFunction(),
	 	alphabet(toCopy.alphabet),
	 	energyMat(),
	 	intervalMax(toCopy.intervalMax)
	{
		for (size_t i=0; i<toCopy.energyMat.size(); i++) {
			energyMat.push_back(new biu::EnergyMatrix(*(toCopy.energyMat[i])));
		}
		assertbiu(energyMat.size() == intervalMax.size(),
				"different numbers of energy tables and interval boundaries");
	}

	
	IntervalEnergyFunction::~IntervalEnergyFunction()
	{
		 // remove heap data
		for (size_t i=0; i<energyMat.size(); i++) {
			delete energyMat[i];
		}
		 // clear vectors
		energyMat.clear();
		intervalMax.clear();
	}
	
	const Alphabet* const
	IntervalEnergyFunction::getAlphabet() const 
	{
		return alphabet;
	}
	
	double 
	IntervalEnergyFunction::getEnergy(	const Alphabet::AlphElem& first, 
										const Alphabet::AlphElem& second,
										const double & distance ) const 
	{
		assertbiu( distance < *(intervalMax.rbegin()), "given distance is out of interval bounds");
	
		return (*energyMat[getInterval(distance)])[alphabet->getIndex(first)][alphabet->getIndex(second)];;
	}
	
	double 
	IntervalEnergyFunction::getEnergy(	const Alphabet::AlphElem & seq_i, 
										const Alphabet::AlphElem & seq_j,
										const IntPoint & cor_i,
										const IntPoint & cor_j  ) const
	{
		// use standard distance evaluation
		return this->getEnergy(seq_i, seq_j, cor_i.distance(cor_j));
	}
								
	double 
	IntervalEnergyFunction::getEnergy(	const Alphabet::AlphElem & seq_i, 
										const Alphabet::AlphElem & seq_j,
										const DblPoint & cor_i,
										const DblPoint & cor_j  ) const
	{
		// use standard distance evaluation
		return this->getEnergy(seq_i, seq_j, cor_i.distance(cor_j));
	}


	bool 
	IntervalEnergyFunction::operator == (const DistanceEnergyFunction& ef2) const
	{
		if (&ef2 == this)
			return true;
		
		bool retVal = false;
		
		const IntervalEnergyFunction* ief = dynamic_cast<const IntervalEnergyFunction*>(&ef2);
		if (ief != NULL) {
			retVal =	(alphabet == ief->alphabet || *alphabet == *(ief->alphabet))
						&& ief->intervalMax.size() == this->intervalMax.size()
						&& ief->energyMat.size() == this->energyMat.size();
			  // compare all intervals in worst case
			for (size_t i=0; retVal && i<intervalMax.size(); i++) {
				retVal =    ief->intervalMax[i] == this->intervalMax[i]
				         && *(ief->energyMat[i]) == *(this->energyMat[i]);
			}
		}
		
		return retVal;
	}
	
	bool 
	IntervalEnergyFunction::operator != (const DistanceEnergyFunction& ef2) const
	{
		return !(this->operator ==(ef2));
	}

	size_t
	IntervalEnergyFunction::addInterval(	const biu::EnergyMatrix& energies,
											const double upperBound )
	{
		assertbiu( intervalMax.size() == 0 || upperBound > *(intervalMax.rbegin()),
					"given upperBound is NOT greater than the last in use");
		assertbiu(	alphabet->getAlphabetSize() == energies.numRows() 
					&& alphabet->getAlphabetSize() == energies.numColumns(),
					"energy matrix has the wrong size for this alphabet size");
		
		  // add energy matrix (copy)
		energyMat.push_back(new biu::EnergyMatrix(energies));
		  // add upper bound of the interval
		intervalMax.push_back(upperBound);
		
		return intervalMax.size()-1;
	}

	size_t
	IntervalEnergyFunction::getIntervalNum(void) const {
		return intervalMax.size();
	}
	
	double
	IntervalEnergyFunction::getIntervalMax(const size_t index) const {
		assertbiu(index < intervalMax.size(), "index out of bound");
		return intervalMax[index];
	}

	const biu::EnergyMatrix* const
	IntervalEnergyFunction::getIntervalMatrix(const size_t index) const {
		assertbiu(index < intervalMax.size(), "index out of bound");
		return energyMat[index];
	}

	size_t
	IntervalEnergyFunction::getInterval( double distance) const {
		if (intervalMax.size() == 0 || distance > *(intervalMax.rbegin())) {
			return UINT_MAX;
		}
		  // find the corresponding index
		size_t index = 0;
		while (distance > intervalMax[index]) {
			index++;
		}
		return index;
	}


} // biu
