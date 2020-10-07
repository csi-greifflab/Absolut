// $Id: LatticeProtein.hh,v 1.2 2016/08/08 12:41:57 mmann Exp $
#ifndef BIU_LATTICEPROTEIN_H_
#define BIU_LATTICEPROTEIN_H_


#include "biu/BackboneStructure3D.hh"
#include "biu/BioMolecule.hh"
#include "biu/LatticeModel.hh"
#include "biu/DistanceEnergyFunction.hh"

namespace biu
{
		/*! The abstract LatticeProtein interface represents a backbone protein 
		 * representation in an energy landscape
		 * whereby all backbone edges are placed on integer lattice points
		 * and consecutive positions are neighbored in the lattice.
		 * 
		 * Such a structure is valid, if and only if it is selfavoiding walk in
		 * the lattice.
		 * 
		 * The energy calulation is based on a DistanceEnergyFunction object. 
		 * A contact is formed, if two structure positions are neighbored and 
		 * not consecutive in sequence.
		 * 
		 * @author Martin Mann, Sebastian Will, Andreas Richter
		 */
	class LatticeProtein : public BioMolecule, public BackboneStructure3D
	{
	public:
			//! empty constructor 
		LatticeProtein() 
		{}
						
			//! empty desctructor
		virtual ~LatticeProtein() 
		{} 
		
	// additional abstract functions in LatticeProtein
	
		virtual bool	 		operator== (const LatticeProtein& latProt2) const = 0;
		virtual bool	 		operator!= (const LatticeProtein& latProt2) const = 0;

	
			//! Returns the absolute move string representation of the lattice
			//! protein structure. 
		virtual std::string getMoveStrAbs() const = 0;
		
			//! Returns the relative move string representation of the lattice
			//! protein structure. 
		virtual std::string getMoveStrRel() const = 0;
		
			//! Returns the absolute move sequence of the lattice
			//! protein structure. 
		virtual MoveSequence getMoveSeqAbs() const = 0;
		
			//! Returns the relative move sequence of the lattice
			//! protein structure. 
		virtual MoveSequence getMoveSeqRel() const = 0;
		
			//! Returns the successive 3D coordinates of the lattice protein
			//! structure.
		virtual IPointVec getPoints() const = 0;

			//! Tests whether or not the structure is self avoiding.
		virtual bool	isSelfAvoiding() const = 0;
		
			//! Tests whether or not consecutive structure elements 
			//! are neighbored in the lattice. 
		virtual bool	isConnected() const = 0;
		
			//! Returns the lattice model this lattice protein is basing on
		virtual const LatticeModel* getLatticeModel() const = 0;
		
			//! Returns the energy function this lattice protein is basing on
		virtual const DistanceEnergyFunction* getEnergyFunction() const = 0;
		
			//! Returns whether or not the sequence is shared among several 
			//! LatticeProtein objects.
		virtual bool isSequenceShared() const = 0;
		
			//! Returns the pointer to the sequence if shared and a new
			//! sequence object otherwise.
			//! 
			//! NOTE: This object has to be deleted  by the user manually!
		virtual const Sequence* getSequenceRef() const = 0;

	};

} // namespace biu

#endif /*LATTICEPROTEIN_H_*/
