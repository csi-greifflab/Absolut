// $Id: LatticeProtein_Ipnt.cc,v 1.2 2016/08/08 12:41:59 mmann Exp $


#include "biu/LatticeProtein_Ipnt.hh"
#include <limits.h>


namespace biu
{

	// initialization of the indifferent value of double variables
	const double LatticeProtein_Ipnt::NAN_DOUBLE = (double)INT_MAX;

	// construction
	LatticeProtein_Ipnt::LatticeProtein_Ipnt(
						const LatticeModel* lattice, 
						const DistanceEnergyFunction* energyFunc, 
						const Sequence* seq,
						const bool seqShared,
						const std::string& moveString, 
						const bool isAbsoluteMove)
	 :	LatticeProtein_I(lattice, energyFunc, seq, seqShared),
	 	points( new IPointVec() ), 
	 	energy(NAN_DOUBLE), 
	 	selfavoiding(MyNaN),
	 	connected(MyNaN) 
	{
		if (isAbsoluteMove) {
			*points = lattice->absMovesToPoints(
							lattice->parseMoveString(moveString));
		} else {
			*points = lattice->relMovesToPoints(
							lattice->parseMoveString(moveString));
		}
		assertbiu ( seq->size() == points->size(),
			"sequence and structure differ in size");
		 // inform the object that the structure has changed
		updateProperties();
	}
	
	LatticeProtein_Ipnt::LatticeProtein_Ipnt(const biu::LatticeProtein& latPr) 
	 :	LatticeProtein_I(latPr.getLatticeModel(), 
	 		latPr.getEnergyFunction(), 
	 		latPr.getSequenceRef(),	
	 		latPr.isSequenceShared()),	 
 		points( new IPointVec()),
	 	energy(NAN_DOUBLE), 
	 	selfavoiding(MyNaN),
	 	connected(MyNaN) 
	{
		assertbiu (lattice != NULL && energyFunc != NULL, 
			"no lattice model or energy function available");
			// check if data is based on points too
		const LatticeProtein_Ipnt* l2 = dynamic_cast<const LatticeProtein_Ipnt*>(&latPr);
		if (l2 != NULL) {	// direct point access possible
			points->resize(l2->getPointsRef()->size());	// resize vector
				// fill vector
			std::copy<IPointVec::const_iterator, IPointVec::iterator>(
				l2->getPointsRef()->begin(), l2->getPointsRef()->end(),
				points->begin());
				// copy maybe already calculated properties
			energy = l2->energy;
			selfavoiding = l2->selfavoiding;
			connected = l2->connected;
		} else {	// via copied objects
			*points = latPr.getPoints();
			 // inform the object that the structure has changed
			updateProperties();
		}
	}
	
	LatticeProtein_Ipnt::LatticeProtein_Ipnt(const biu::LatticeProtein_Ipnt& latPr) 
	 :	LatticeProtein_I(latPr), 
	 	points( new IPointVec(latPr.points->size())),
	 	energy(latPr.energy), 
	 	selfavoiding(latPr.selfavoiding),
	 	connected(latPr.connected) 
	{
		assertbiu (lattice != NULL && energyFunc != NULL, 
			"no lattice model or energy function available");
			// copy structure
		std::copy<IPointVec::const_iterator, IPointVec::iterator>(
				latPr.points->begin(), latPr.points->end(),
				points->begin());
	}
	
	LatticeProtein_I* LatticeProtein_Ipnt::clone() const {
		return new LatticeProtein_Ipnt(*this);
	}
	
	LatticeProtein_I* LatticeProtein_Ipnt::fromString(const std::string& stringRep) const {
		size_t pos = stringRep.find_first_of("(");
		std::string absMoveStr = stringRep.substr(0, pos);
		std::string seqStr = stringRep.substr(pos+1, stringRep.size()-pos-2);
		biu::Sequence seq = energyFunc->getAlphabet()->getSequence(seqStr);
		const bool seqShared = false;
		const bool isAbsMove = true;
		return new LatticeProtein_Ipnt(lattice, energyFunc, &seq, seqShared, absMoveStr, isAbsMove);
	}
	
	LatticeProtein_Ipnt::~LatticeProtein_Ipnt()
	{
		if (points != NULL) {
			delete points;		points = NULL;
		}
	}
	
	LatticeProtein_I&
	LatticeProtein_Ipnt::operator =(const LatticeProtein_I& latProt2) {
		const LatticeProtein_Ipnt* l2 = dynamic_cast<const LatticeProtein_Ipnt*>(&latProt2);
		if (l2 != NULL) {	// direct point access possible
			return this->operator=(*l2);
		} else
			return LatticeProtein_I::operator=(latProt2);
	}
	
	LatticeProtein_Ipnt&
	LatticeProtein_Ipnt::operator =(const LatticeProtein_Ipnt& latProt2) {
		if (this != &latProt2 && *this != latProt2) {
			assertbiu(sequence != NULL && points != NULL,
				"sequence and structure not initialized");
				
				// copy data
			lattice = latProt2.getLatticeModel();
			energyFunc = latProt2.getEnergyFunction();
			if (!seqShared) { // delete sequence if not shared
				delete sequence;
			}
			seqShared = latProt2.isSequenceShared();
			  // copy sequence
			if (seqShared) {
				sequence = latProt2.getSequenceRef();
			} else {
				sequence = new Sequence(*(latProt2.getSequenceRef()));
			}
			points->resize(latProt2.points->size());
			std::copy<IPointVec::const_iterator, IPointVec::iterator>(
					latProt2.points->begin(), latProt2.points->end(),
				points->begin());
				// copy maybe already calculated properties
			energy = latProt2.energy;
			selfavoiding = latProt2.selfavoiding;
			connected = latProt2.connected;
			 // inform the object that the structure has changed
			updateProperties();
		}
		return *this;
	}
	
	bool
	LatticeProtein_Ipnt::operator ==(const LatticeProtein& latProt2) const {
		assertbiu (lattice != NULL, "no lattice available");
		assertbiu (energyFunc != NULL, "no energy function available");
		assertbiu (sequence != NULL, "no sequence available");
		assertbiu (points != NULL, "no structure available");
		bool equal = 
				(lattice == latProt2.getLatticeModel() 
						|| *lattice == *(latProt2.getLatticeModel()))
				&& (energyFunc == latProt2.getEnergyFunction()
						|| *energyFunc == *(latProt2.getEnergyFunction()))
				&& (seqShared == latProt2.isSequenceShared())
				&& ((seqShared && sequence == latProt2.getSequenceRef())
						|| *sequence == *(latProt2.getSequenceRef()));
		if (equal) {
				// check if data is based on points too
			const LatticeProtein_Ipnt* l2 = dynamic_cast<const LatticeProtein_Ipnt*>(&latProt2);
			if (l2 != NULL) {	// direct point access possible
				equal &= *points == *l2->points;
			} else {	// comparison based on relative move strings
				equal &= lattice->pointsToRelMoves(*points) == latProt2.getMoveSeqRel();
			}
		}
		return equal;
	}
	
	bool
	LatticeProtein_Ipnt::operator !=(const LatticeProtein& latProt2) const {
		return	! operator == (latProt2);
	}
	
	
	bool 
	LatticeProtein_Ipnt::isSelfAvoiding() const {
		assertbiu(points != NULL, "no structure available");
		
			// if selfavoidingness already tested
		if (selfavoiding != MyNaN) {
			return ( selfavoiding == MyTrue ? true : false );
		} 
		
		// else selfavoidingness has to be tested
		
			// check by set insertion to the cost of memory allocation
//		IPointSet tmp;
			// pruefe via einfuegen, ob alle punkte unique sind
//		for (IPointVec::const_iterator it = pData->begin();  
//				it != pData->end(); it++) {
//			if (tmp.insert(*it).second == false) // already inside (not inserted)
//				return false;	// not selfavoiding
//		}
		
			// checking via (n^2)/2 comparisons  ==  inplace  
			// (faster than insertion for small structures)
		for (size_t k=0; k<points->size(); k++) {
			for (size_t l=k+1; l<points->size(); l++) {
				if (points->at(k) == points->at(l)) {
					selfavoiding = MyFalse;
					return false;
				}
			}
		}
		selfavoiding = MyTrue;
		return true;
	}

	// abstract functions 
	
	DPointVec	
	LatticeProtein_Ipnt::get3Ddata() const  {
		assertbiu (points != NULL, "no structure available");
			// generate return vector
		DPointVec dVec(points->size());
			// fill return vector via copy and type cast
		std::copy<IPointVec::const_iterator, DPointVec::iterator> (
				points->begin(), points->end(),
				dVec.begin() );
		return dVec;
	}
	
	double 
	LatticeProtein_Ipnt::getDRMSD(const BackboneStructure3D& other) const {
		DPointVec this3Ddata = get3Ddata(), other3Ddata = other.get3Ddata();

		assertbiu(this3Ddata.size() == other3Ddata.size(),
		 "Both lattice proteins have to have the same length.");

		double rmsd = 0.0;

		for(size_t i=0; i < this3Ddata.size()-1; i++) {
			for(size_t j=i+1; j < this3Ddata.size(); j++) {
				rmsd += std::pow( this3Ddata[i].distance(this3Ddata[j])
						 - other3Ddata[i].distance(other3Ddata[j]), 2 );
			}
		}

        rmsd /= (this3Ddata.size() * (this3Ddata.size()-1) /2);

        return std::sqrt(rmsd);

	}

	double	
	LatticeProtein_Ipnt::getEnergy() const {
		if (energy == NAN_DOUBLE) {
			assertbiu(lattice != NULL, "no lattice model available");
			assertbiu(energyFunc != NULL, "no energy function available");
			energy = 0.0;
			IPointVec::size_type i=0,j=0;
			for (i=0; i < points->size(); i++) {
				for (j = i+2; j < points->size(); j++) {
					  // distance based energy evaluation
					energy += energyFunc->getEnergy(	sequence->at(i),
														sequence->at(j),
														points->at(i), 
														points->at(j));
				}
			}
		}
		return energy;
	}
	
		//! Returns the relative move representation of the lattice
		//! protein structure. 
	Structure	
	LatticeProtein_Ipnt::getStructure() const {
		assertbiu(lattice != NULL, "no lattice model available");
		return lattice->pointsToRelMoves(*points);
	}
	
		/*! Returns a combination of absolute move string representation 
		 * and sequence information of this lattice protein.
		 * 
		 * ABSOLUTEMOVESTRING(SEQUENCE) */
	std::string	
	LatticeProtein_Ipnt::getStringRepresentation() const {
		assertbiu(energyFunc != NULL, "no energy function available");
		
		return	getMoveStrAbs()
				+ std::string("(")
				+ energyFunc->getAlphabet()->getString(*sequence)
				+ std::string(")");
	}
	
	
	IPointVec 
	LatticeProtein_Ipnt::getPoints() const {
		assertbiu(points != NULL, "no structure available");
		return *points;
	}
	
	const IPointVec* const
	LatticeProtein_Ipnt::getPointsRef() const {
		return points;
	}
	
	IPointVec*
	LatticeProtein_Ipnt::getPointsRef() {
		updateProperties();
		return points;
	}
	

	bool		
	LatticeProtein_Ipnt::isValid() const {
		return isConnected() && isSelfAvoiding();
	}

		//! Returns the absolute move string representation of the lattice
		//! protein structure. 
	std::string 
	LatticeProtein_Ipnt::getMoveStrAbs() const {
		assertbiu(lattice != NULL, "no lattice model available");
		return lattice->getString(lattice->pointsToAbsMoves(*points));
	}
	
		//! Sets a new absolute move string.
	void
	LatticeProtein_Ipnt::setMoveStrAbs(const std::string& moveString) {
		*points = lattice->absMovesToPoints(
							lattice->parseMoveString(moveString));
		assertbiu ( sequence->size() == points->size(),
			"sequence and structure differ in size");
		 // inform the object that the structure has changed
		updateProperties();
	}
	
		//! Returns the relative move string representation of the lattice
		//! protein structure. 
	std::string 
	LatticeProtein_Ipnt::getMoveStrRel() const {
		assertbiu(lattice != NULL, "no lattice model available");
		return lattice->getString(lattice->pointsToRelMoves(*points));
	}
	
		//! Returns the absolute move sequence of the lattice
		//! protein structure. 
	MoveSequence 
	LatticeProtein_Ipnt::getMoveSeqAbs() const {
		assertbiu(lattice != NULL, "no lattice model available");
		return lattice->pointsToAbsMoves(*points);
	}
	
		//! Returns the relative move sequence of the lattice
		//! protein structure. 
	MoveSequence 
	LatticeProtein_Ipnt::getMoveSeqRel() const {
		assertbiu(lattice != NULL, "no lattice model available");
		return lattice->pointsToRelMoves(*points);
	}
	
		//! Tests whether or not consecutive structure elements 
		//! are neighbored in the lattice. 
	bool	
	LatticeProtein_Ipnt::isConnected() const {
		assertbiu(points != NULL, "no structure available");
		assertbiu(lattice != NULL, "no lattice model available");
		
		if (connected != MyNaN ) {
			return ( connected == MyTrue ? true : false );
		}
			// check if there are some unneighbored successive points
		for(size_t i=1; i<points->size(); i++) {
			if (! lattice->areNeighbored(points->at(i-1),points->at(i)) ) { 
				connected = MyFalse;
				return false;
			}
		}
			// everything all right .. ;)
		connected = MyTrue;
		return true;
	}
	
	
	
		// additional functions (LatticeProtein_Ipnt)
	
	// Sets all calculated properties (energy, selfavoidingness and 
	// connectedness to an indifferent state.
	void
	LatticeProtein_Ipnt::updateProperties() {
			// set all properties that are calculated on demand to an 
			// indifferent state
		energy = NAN_DOUBLE;
		selfavoiding = MyNaN;
		connected = MyNaN;
	}

}
