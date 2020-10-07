#ifndef BIU_LATTICEPROTEINUTIL_HH_
#define BIU_LATTICEPROTEINUTIL_HH_

#include <biu/LatticeDescriptor.hh>
#include <biu/LatticeProtein.hh>

namespace biu
{
	
	class LatticeProteinUtil
	{
	public:
		
		
		
		////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////
		// CONVERSION BETWEEN STRUCTURE REPRESENTATIONS
		////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////
		
		
		  /*! Converts a move string into a move sequence.
		   * 
		   * @param moveString the string to convert
		   * @param latDescr the lattice descriptor of underlying lattice of the
		   *         given move string
		   * @param sideChain whether or not the given move string encodes a
		   *         backbone-only or side chain lattice protein model
		   * 
		   * @return the move sequence of the encoded lattice protein
		   * 
		   */
		static
		biu::MoveSequence
		toMoveSequence(	const std::string & moveString
						, const biu::LatticeDescriptor & latDescr
						, const bool sideChain );
		
		
		
		  /*! Converts a vectors of 3D coordinates into an absolute 
		   * move sequence.
		   * 
		   * @param points the coordinates to convert
		   * @param latDescr the lattice descriptor of underlying lattice of the
		   *         given move string
		   * @param sideChain whether or not the given move string encodes a
		   *         backbone-only or side chain lattice protein model
		   * 
		   * @return the move sequence of the encoded lattice protein
		   * 
		   */
		static
		biu::MoveSequence
		toMoveSequence(	const biu::IPointVec & points 
						, const biu::LatticeDescriptor & latDescr
						, const bool sideChain );
		
		
		
		  /*! Converts a move string into a move sequence.
		   * 
		   * @param moves the move sequence to convert
		   * @param latDescr the lattice descriptor of underlying lattice of the
		   *         given move sequence
		   * @param sideChain whether or not the given move sequence encodes a
		   *         backbone-only or side chain lattice protein model
		   * 
		   * @return the move string of the encoded lattice protein
		   * 
		   */
		static
		std::string
		toString(	const biu::MoveSequence & moves
					, const biu::LatticeDescriptor & latDescr
					, const bool sideChain );
		
		
		
		  /*! Converts an absolute move sequence into 3D coordinates.
		   * 
		   * @param moves the move sequence to convert
		   * @param latDescr the lattice descriptor of underlying lattice of the
		   *         given move sequence
		   * @param sideChain whether or not the given move sequence encodes a
		   *         backbone-only or side chain lattice protein model
		   * @param cAlphaDist the C_alpha distance to scale the lattice protein
		   * 
		   * @return the 3D coordinates
		   * 
		   */
		static
		biu::DPointVec
		toDblPoints(	const biu::MoveSequence & moves
						, const biu::LatticeDescriptor & latDescr
						, const bool sideChain 
						, const double cAlphaDist = 3.8
					);
		
		
		  /*! Converts an absolute move sequence into 3D coordinates.
		   * 
		   * @param moves the move sequence to convert
		   * @param latDescr the lattice descriptor of underlying lattice of the
		   *         given move sequence
		   * @param sideChain whether or not the given move sequence encodes a
		   *         backbone-only or side chain lattice protein model
		   * 
		   * @return the 3D coordinates
		   * 
		   */
		static
		biu::IPointVec
		toIntPoints(	const biu::MoveSequence & moves
						, const biu::LatticeDescriptor & latDescr
						, const bool sideChain );
		
		
		  /*! Calculates the base vector scaling for a given lattice description
		   * that scales all neighbor vectors to the given length
		   * 
		   * @param latDescr the lattice to scale
		   * @param neighVecLength the length that the neighbor vectors should
		   *         be scaled to
		   * 
		   * @return the multiplicator for the base vectors to scale the lattice
		   */
		static
		double
		getBaseScale(	const biu::LatticeDescriptor & latDescr
						, const double neighVecLength );
		
		
		
		
		////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////
		// STRUCTURAL DISTANCE CALCULATION
		////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////
		
		
		  /*!
		   * Calculates the pairwise coordinate root mean square deviation 
		   * (cRMSD) between two point vectors.
		   * 
		   * @param pos1 first position data
		   * @param pos2 second position data (equal size to pos1)
		   * 
		   * @return the cRMSD value
		   * 
		   */
		static
		double
		cRMSD(	const biu::DPointVec & pos1
				, const biu::DPointVec & pos2 );
		
		  /*!
		   * Calculates the pairwise coordinate root mean square deviation 
		   * (cRMSD) between two point vectors.
		   * 
		   * @param pos1 first position data
		   * @param pos2 second position data (equal size to pos1)
		   * 
		   * @return the cRMSD value
		   * 
		   */
		static
		double
		cRMSD(	const biu::IPointVec & pos1
				, const biu::IPointVec & pos2 );
		
		  /*!
		   * Calculates the distance root mean square deviation 
		   * (dRMSD) between two point vectors.
		   * 
		   * @param pos1 first position data
		   * @param pos2 second position data (equal size to pos1)
		   * 
		   * @return the dRMSD value
		   * 
		   */
		static
		double
		dRMSD(	const biu::DPointVec & pos1
				, const biu::DPointVec & pos2 );
		
		  /*!
		   * Calculates the distance root mean square deviation 
		   * (dRMSD) between two point vectors.
		   * 
		   * @param pos1 first position data
		   * @param pos2 second position data (equal size to pos1)
		   * 
		   * @return the dRMSD value
		   * 
		   */
		static
		double
		dRMSD(	const biu::IPointVec & pos1
				, const biu::IPointVec & pos2 );

		
		  /*!
		   * Calculates the Global-Distance-Test_Total-Score (GDT_TS)
		   * between two point vectors.
		   * 
		   * GDT_TS = (GDT_P1 + GDT_P2 + GDT_P4 + GDT_P8)/4
		   * 
		   * where GDT_Pn denotes percent of residues under distance cutoff <= n
		   * 
		   * @param pos1 first position data
		   * @param pos2 second position data (equal size to pos1)
		   * 
		   * @return the GDT_TS value
		   * 
		   */
		static
		double
		GDT_TS(	const biu::DPointVec & pos1
				, const biu::DPointVec & pos2 );
		
		
		  /*!
		   * Calculates the Global-Distance-Test_High-Accuracy score (GDT_HA)
		   * between two point vectors.
		   *
		   * GDT_HA = (GDT_P0.5 + GDT_P1 + GDT_P2 + GDT_P4)/4,
		   * 
		   * where GDT_Pn denotes percent of residues under distance cutoff <= n
		   * 
		   * @param pos1 first position data
		   * @param pos2 second position data (equal size to pos1)
		   * 
		   * @return the GDT_HA value
		   * 
		   */
		static
		double
		GDT_HA(	const biu::DPointVec & pos1
				, const biu::DPointVec & pos2 );
		
	};

} // namespace biu

#endif /*LATTICEPROTEINUTIL_HH_*/
