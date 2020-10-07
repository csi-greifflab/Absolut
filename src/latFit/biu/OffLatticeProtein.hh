// $Id: OffLatticeProtein.hh,v 1.2 2016/08/08 12:41:58 mmann Exp $
#ifndef BIU_OFFLATTICEPROTEIN_H_
#define BIU_OFFLATTICEPROTEIN_H_


#include "biu/LatticeProtein.hh"

namespace biu
{
		/**
		 * A representation of an off-lattice backbone protein
		 * structure.
		 * 
		 * It provides a PDB-file interface.
		 * 
		 * @author Martin Mann, Sebastian Will, Andreas Richter
		 */
	class OffLatticeProtein : public BackboneStructure3D
	{
	protected:

		DPointVec*	pData;			//!< the structure information 
		const Alphabet* alphabet;	//!< the alphabet of the sequence 
		Sequence*	sequence;		/*!< the sequence of the lattice protein 
									 * in one letter code */

	public:
			/*! Constructs a new off lattice protein.
			 *  @param coordinatesFileName Inputfile in the format "aa x y z \n"
			 *  in which aa is the amino acid in one letter code and
			 *  x, y and z specify the coordinates of the amino acid.
			 *  @param _alphabet The alphabet of the protein sequence.
			 */
		OffLatticeProtein(	const std::string& coordinatesFileName,
							const Alphabet* const _alphabet);
							
			/*! Constructs a new off lattice protein.
			 *  @param data3D 3D double coordinates of the lattice protein.
			 *  @param _alphabet The alphabet of the protein sequence.
			 *  @param seqStr The sequence of the lattice protein.
			 *  @param data3D 3D-position data of the protein.
			 */
		OffLatticeProtein(	const DPointVec& data3D, 
							const Alphabet* const _alphabet, 
							const std::string& seqStr);
				   
		OffLatticeProtein(const OffLatticeProtein& offLatPro);
		virtual ~OffLatticeProtein();

		
		OffLatticeProtein& operator= (const OffLatticeProtein& offLatPro2);

			//! Converts the Structure to PDB-output and writes to file 
		void writePDB(const std::string& pdbFileName);
		
			/*! Approximates an off lattice backbone protein to a
			 *  given lattice using a brute force build-up algorithm
			 *  described in
			 *  Park B. H., Levitt M. (1995). The complexity and accuracy
			 *  of discrete state models of protein structure. 
			 *  Journal of Molecular Biology, 249(2), 493-507.
			 * 
			 * 	@return The returned object was allocated with NEW and has to be
			 * 		destroyed by the calling function via DELETE!
			 */
		LatticeProtein* approximateToLattice( 
							const LatticeModel* const lattice, 
							const DistanceEnergyFunction* const energy) const;

	// abstract functions (BackboneStructure3D)
		
			/*! Returns the consecutive point representation of the
			 * backbone. */
		virtual DPointVec get3Ddata() const;
		
			/*! Returns the distance root mean square deviation (DRMSD) to
			 *  an other BackboneStructure3D
			 */
		virtual double getDRMSD(const BackboneStructure3D& other) const;
	};

} // namespace biu

#endif /*OFFLATTICEPROTEIN_H_*/
