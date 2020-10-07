// $Id: LatticeProtein_I.cc,v 1.2 2016/08/08 12:41:56 mmann Exp $


#include "biu/LatticeProtein_I.hh"
#include "biu/assertbiu.hh"

namespace biu {
	
	LatticeProtein_I::LatticeProtein_I(const LatticeModel* lattice_, 
							const DistanceEnergyFunction* energyFunc_, 
							const Sequence* seq_,
							const bool seqShared_)
	 :	lattice(lattice_), energyFunc(energyFunc_), sequence(NULL), 
	 	seqShared(seqShared_)
	{
		assertbiu (lattice != NULL && energyFunc != NULL, 
			"no lattice model or energy function available");
		assertbiu (seq_ != NULL,
			"no sequence available");
		assertbiu(energyFunc->getAlphabet()->isAlphabetSequence(*seq_),
			"sequence is not valid for the given alphabet");
		// initialize sequence
		if (seqShared) { // use given pointer
			sequence = &(*seq_);
		} else { // copy sequence
			sequence = new Sequence(*seq_);
		}
	}
	
	LatticeProtein_I::LatticeProtein_I(const LatticeProtein_I& l2)
	  : lattice(l2.lattice), energyFunc(l2.energyFunc), sequence(NULL), 
	 	seqShared(l2.seqShared)
	{
		if (seqShared)
			sequence = l2.sequence;
		else
			sequence = new Sequence(*(l2.sequence));
	}
	
	LatticeProtein_I::~LatticeProtein_I()
	{
		 // delete sequence if neccessary
		if (!seqShared && sequence != NULL) {
			delete sequence;
			sequence = NULL;
		}
	}

// abstract function implementation (LatticeProtein)

		//! Returns the lattice model this lattice protein is basing on
	const biu::LatticeModel* 
	LatticeProtein_I::getLatticeModel() const {
		return lattice;
	}
	
		//! Returns the energy function this lattice protein is basing on
	const biu::DistanceEnergyFunction* 
	LatticeProtein_I::getEnergyFunction() const {
		return energyFunc;
	}
	
		//! Returns whether or not the sequence is shared among several 
		//! LatticeProtein objects.
	bool 
	LatticeProtein_I::isSequenceShared() const {
		return seqShared;
	}

		//! Returns the pointer to the sequence 
	const Sequence* 
	LatticeProtein_I::getSequenceRef() const {
		return sequence;
	}
	

// abstract function implementation (BioMolecule)
		
	biu::Sequence	
	LatticeProtein_I::getSequence() const
	{
		assertbiu( sequence != NULL, "no sequence available");
		return *sequence;
	}
		
		//! Returns the length of the Biomolecule, i.e. the number of 
		//! monomers.
	size_t		
	LatticeProtein_I::getLength() const {
		assertbiu( sequence != NULL, "no sequence available");
		return sequence->size();
	}
	
	LatticeProtein_I&
	LatticeProtein_I::operator =(const LatticeProtein_I& latProt2) {
		if (this != &latProt2 && *this != latProt2) {
			assertbiu(sequence != NULL,
				"sequence not initialized");
				
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
		}
		return *this;
	}
	
	
}
