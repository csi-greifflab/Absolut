// $Id: LatticeProtein_Ipnt.hh,v 1.2 2016/08/08 12:42:00 mmann Exp $
#ifndef BIU_LATTICEPROTEIN_PNT_HH_
#define BIU_LATTICEPROTEIN_PNT_HH_

#include "biu/LatticeProtein_I.hh"

namespace biu
{
	
	/*!
	 * A LatticeProtein implementation using a vector of 3D coordinates as 
	 * internal structure representation. Move string representations are 
	 * computed on demand.
	 * 
	 */
	class LatticeProtein_Ipnt : public biu::LatticeProtein_I
	{
	protected:
	
			//! trigger value for double variables that they should be 
			//! recalculated
		static const double NAN_DOUBLE;
		
			//! new data type to allow indifferent state of boolean variables
			//! that allow their calculation on demand
		enum MyBool { MyFalse, MyTrue, MyNaN };
	


		//! the structure in 3D coordinates
		IPointVec*		points;	

			//! the energy of this structure that is calculated lazily on demand
		mutable double	energy;
		
			//! lazy calculated selfavoidingness
		mutable MyBool selfavoiding;
		
			//! lazy calculated connectedness
		mutable MyBool connected;

	
	public:
	
		/*! Constructs a new lattice protein object based on a 3D coordinate
		 * structure representation.
		 *  
		 * @param lattice		the lattice model the protein structure is 
		 * 						basing on
		 * @param energy 		the energy function including the alphabet
		 * 						for the sequence
		 * @param seq			the sequence encoded using the alphabet of
		 * 						the energy function object
		 * @param seqShared		false if a local copy of seq should be 
		 * 						generated; true if seq is shared among
		 * 						several LatticeProtein objects using the 
		 * 						provided seq pointer
		 * @param moveString 	a move string representation of the 
		 * 						structure
		 * @param isAbsoluteMove true if moveString is an absolute move 
		 * 						string or false if it is a relative 
		 * 						move string
		 */
		LatticeProtein_Ipnt(const LatticeModel* lattice, 
						const DistanceEnergyFunction* energy, 
						const Sequence* seq,
						const bool seqShared,
						const std::string& moveString, 
						const bool isAbsoluteMove);

		LatticeProtein_I* clone() const;
		LatticeProtein_I* fromString(const std::string& stringRep) const;
		
		// copy constructors
	
		LatticeProtein_Ipnt(const LatticeProtein& latProt); 
		LatticeProtein_Ipnt(const LatticeProtein_Ipnt& latProt); 
						
		/*! Destruction of the lattice protein object.
		 */
		virtual ~LatticeProtein_Ipnt();
		
	// additional functions
	
		virtual const IPointVec* const getPointsRef() const;
		virtual IPointVec* getPointsRef(); 
		virtual void setMoveStrAbs(const std::string& moveString);

	// abstract function implementation (BackboneStructure3D)
		
		virtual DPointVec	get3Ddata() const ;
		
		virtual double		getDRMSD(const BackboneStructure3D& other) const;

	// abstract function implementation (BioMolecule)
		
			//! Returns the relative move representation of the lattice
			//! protein structure. 
		virtual Structure	getStructure() const;
		
			/*! Returns the contact energy of the lattice protein structure.
			 * Consecutive positions are not taken into account. 
			 * Each contact is count only once.
			 */
		virtual double		getEnergy() const;
		
			/*! Returns whether or not the lattice protein is connected 
			 * and self avoiding.
			 */
		virtual bool		isValid() const;
		
			/*! Returns a combination of absolute move string representation 
			 * and sequence information of this lattice protein.
			 * 
			 * ABSOLUTEMOVESTRING(SEQUENCE) */
		virtual std::string	getStringRepresentation() const;
		

	// abstract function implementation (LatticeProtein)

			//! Returns the absolute move string representation of the lattice
			//! protein structure. 
		virtual std::string getMoveStrAbs() const;
		
			//! Returns the relative move string representation of the lattice
			//! protein structure. 
		virtual std::string getMoveStrRel() const;
		
			//! Returns the absolute move sequence of the lattice
			//! protein structure. 
		virtual MoveSequence getMoveSeqAbs() const;
		
			//! Returns the relative move sequence of the lattice
			//! protein structure. 
		virtual MoveSequence getMoveSeqRel() const;
		
			//! Returns the successive 3D coordinates of the lattice protein
			//! structure.
		virtual IPointVec getPoints() const;

			//! Tests whether or not the structure is self avoiding.
		virtual bool	isSelfAvoiding() const;
		
			//! Tests whether or not consecutive structure elements 
			//! are neighbored in the lattice. 
		virtual bool	isConnected() const;
			
		virtual LatticeProtein_Ipnt&	operator= (const LatticeProtein_Ipnt& latProt2);		
		virtual LatticeProtein_I&	operator= (const LatticeProtein_I& latProt2);		
		virtual bool	 		operator== (const LatticeProtein& latProt2) const;
		virtual bool	 		operator!= (const LatticeProtein& latProt2) const;

	// additional functions (LatticeProtein_Ipnt)
		
	protected:
		
			//! Is called if the structure or sequence of the protein has 
			//! changed. Either all new properties (e.g. energy etc.) are 
			//! calculated or they are triggered to be calculated on demand.
			//!
			//! Here: Sets all calculated properties (energy, selfavoidingness and 
			//! connectedness to an indifferent state.
		virtual void updateProperties();
	
	};

}

#endif /*LATTICEPROTEIN_PNT_HH_*/
