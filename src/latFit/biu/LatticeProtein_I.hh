// $Id: LatticeProtein_I.hh,v 1.2 2016/08/08 12:42:00 mmann Exp $
#ifndef BIU_LATTICEPROTEIN_I_HH_
#define BIU_LATTICEPROTEIN_I_HH_

#include "biu/LatticeProtein.hh"

namespace biu
{
	/*!
	 * A partial LatticeProtein implementation providing storage and access to
	 * the underlying lattice and energy function. Further the protein objects 
	 * sequence, encoded using the alphabet of the energy function, is handled.
	 */
	class LatticeProtein_I : public biu::LatticeProtein
	{

	protected:
			//! the lattice model this structure belongs to
		LatticeModel const *	lattice;		
		
			//! the distance based energy function 
		DistanceEnergyFunction const *	energyFunc;	
		
			//! the indexed sequence of the lattice protein 
		Sequence const *				sequence;	
		
			//! whether or not the sequence is shared or has to be copied
		bool					seqShared;	
	
	public:
	
		/*! Constructs a new lattice protein object and initializes some
		 * internal data structures.
		 *  
		 * @param lattice		the lattice model the protein structure is 
		 * 						basing on
		 * @param energyFunc	the energy function inclusive the alphabet
		 * 						for the sequence
		 * @param seq			the sequence encoded using the alphabet of
		 * 						the energy function object
		 * @param seqShared		false if a local copy of seq should be 
		 * 						generated; true if seq is shared among
		 * 						several LatticeProtein objects using the 
		 * 						provided seq pointer
		 */
		LatticeProtein_I(const LatticeModel* lattice, 
						const DistanceEnergyFunction* energyFunc, 
						const Sequence* seq,
						const bool seqShared);

		LatticeProtein_I(const LatticeProtein_I& toCopy);

		virtual LatticeProtein_I* clone() const = 0;
		virtual LatticeProtein_I* fromString(const std::string& stringRep) const = 0;
						

		/*! destruction
		 */
		virtual ~LatticeProtein_I();
		
		// additional functions
		
		/*!
		 * Returns a pointer to the integer-point representation of the
		 * protein move sequence for read only access.
		 */
		virtual const IPointVec* const getPointsRef() const = 0;
		
		/*!
		 * Returns a pointer to the integer-point representation of the
		 * protein for read-write access. Calling this function sets all
		 * calculated properties to an indifferent state.
		 */
		virtual IPointVec* getPointsRef() = 0;
		
		/*!
		 * Sets a new absolute move string.
		 */
		virtual void setMoveStrAbs(const std::string& moveString) = 0;
		
	// abstract function implementation (LatticeProtein)

			//! Returns the lattice model this lattice protein is basing on
		virtual const LatticeModel* getLatticeModel() const;
		
			//! Returns the energy function this lattice protein is basing on
		virtual const DistanceEnergyFunction* getEnergyFunction() const;

			//! Returns whether or not the sequence is shared among several 
			//! LatticeProtein objects.
		virtual bool isSequenceShared() const;

			//! Returns the pointer to the sequence
		virtual const Sequence* getSequenceRef() const;

	// abstract function implementation (BioMolecule)
		
		Sequence	getSequence() const;

			//! Returns the length of the Biomolecule, i.e. the number of 
			//! monomers.
		size_t		getLength() const;
		
		virtual LatticeProtein_I&	operator= (const LatticeProtein_I& latProt2);		
			
	};

}

#endif /*LATTICEPROTEIN_I_HH_*/
