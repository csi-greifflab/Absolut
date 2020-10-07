// $Id: LatticeDescriptor.hh,v 1.2 2016/08/08 12:41:58 mmann Exp $
#ifndef BIU_LATTICEDESCRIPTOR_H_
#define BIU_LATTICEDESCRIPTOR_H_


#include <string>
#include <vector>
#include <set>
#include "biu/LatticeNeighborhood.hh"

namespace biu
{
	// this definition can potentially improve performance
	// but can be saved by allowing Matrix<T> * vector<T> multiplication
	// and implicit castings IntPoint <-> vector<int>.
	// anyway this should be moved to IntPoint definition
	// or inside of class LatticeModel as private Method
	// /*! matrix point multiplication without typecast */
	// IntPoint operator*(const biu::Automorphism& am, const biu::IntPoint& p);
	
		/**
		 * An LatticeDescriptor object handles the lattice
		 * specific data to construct an integer point based lattice.
		 * It controls all available automorphisms, the lattice base 
		 * vectors and the neighborhood vectors that are available 
		 * for the managed lattice.
		 * 
		 * @author Martin Mann, Sebastian Will, Andreas Richter
		 */
	class LatticeDescriptor {

	private:		

		std::string	name;			/*!< the name of the lattice model */

	protected:
		
		AutomorphismVec	automorphisms;	/*!< all available automorphisms 
										 * in this lattice e.g. rotations */
		
		IPointVec		latBase;		/*!< the base vectors of the lattice */
		
		MoveAlphabet*	moveAlphabet;	/*!< the move string alphabet */
		
			/*! The neighbor vectors that define the neighborhood in the 
			 * lattice. 
			 * 
			 * The first NeighborVector in the set represents the default 
			 * direction for relative move string calculation.
			 * */
		LatticeNeighborhood*	latNeighborhood;
		
			/*! A shuffle information that can be used to normalize move 
			 * sequences by replacing moves with symmetric ones according to
			 * a symmetry. It is initializes by #initAutomorphisms()
			 */ 
		std::vector<MoveSequence>	symMoveReplacement;
		

			//! helper struct for data of lattice neighbors
			//!
			//! defines a lattice move with vector and matrices
		struct NeighborData { 
			const char *name;       //!< name of move
			int vec[3];       //!< absolute move vector
			int mat[3][3];    //!< matrix for applying rel. move
			int invmat[3][3]; //!< inverse of matrix (we are lazy ;-))
		};
		
		struct AutomorphismData {
			int _0[3];
			int _1[3];
			int _2[3];
		};

			/*! Gets pointer to neighbor data array. Together with 
			 * getNeighborDataSize this defines the lattice in derived classes in 
			 * a way that can be used by initNeighborhood for initializing
			 * the lattice neighborhood and move alphabet. */
		virtual const NeighborData* getNeighborData() const = 0;
		
			//! get number of neighbors
		virtual unsigned int getNeighborDataSize() const = 0; 

			/*! Gets pointer to automorphism data array. */
		virtual const AutomorphismData* getAutomorphismData() const = 0;
		
			//! get number of automorphisms
		virtual unsigned int getAutomorphismDataSize() const = 0; 

		
			/*! inits the neighborhood from the data of getNeighborData()
			  and getNeighborDataSize().
			  Generic function that should be called in the constructors
			  of derived classes
			*/
		void initNeighborhood(); 
		
			//! inits the automorphisms using getAutomorphismData*
		void initAutomorphisms();
		
	public:
		
		
		  /*! Calculates the base vector scaling 
		   * that scales all neighboring vectors to the given length. 
		   * 
		   * @param neighVecLength the length that the neighbor vectors should
		   *         be scaled to
		   * 
		   * @return the multiplicator for the base vectors to scale the lattice
		   */
		virtual
		double
		getBaseScale( const double neighVecLength ) const = 0;
	
		LatticeDescriptor(const std::string& name_);
		LatticeDescriptor(const LatticeDescriptor& toCopy);
		virtual ~LatticeDescriptor();
		
			/*! Returns the move alphabet of this lattice descriptor. */
		const Alphabet* const getAlphabet() const { 
			return moveAlphabet; 
		}

		virtual std::string	 getName() const;
		
		virtual const LatticeNeighborhood&	getNeighborhood() const;
		
		  /*! Checks whether or not a given ring size is possible as a 
		   * selfavoiding walk with equal start/end position.
		   * @param ringSize the size of the ring to check for (at least 3)
		   * @return true if a SAW ring of that size is possible; 
		   *         false otherwise
		   */
		virtual 
		bool
		isPossibleRing( const size_t ringSize ) const = 0;

			//! Returns whether or not two points in the lattice are neighbored.
			//! @param first the first point
			//! @param second the second point
			//! @return true if (second-first) is element of the neighborhood,
			//!         false otherwise
		virtual bool areNeighbored( const IntPoint &first, 
									const IntPoint &second ) const;
		
			//! Checks if the given point is a valid node of the lattice or not
			//! @param p the point to check
			//! @return true if the point is reachable via a sequence of base
			//!   vectors; false otherwise
		virtual bool isLatticeNode( const IntPoint & p ) const = 0; 
		
			/*! Returns the lattice base vectors. */
		virtual const IPointVec&	getBase() const;

			/*! Returns a vector of all automorphisms in this lattice. */
		virtual const AutomorphismVec&	getAutomorphisms() const;

			/*! Converts an internal move string sequence representation into
			 * a string representation.
			 * 
			 * Returns an empty string in error case.
			 */
		virtual std::string getString(const MoveSequence& moveSeq) const;

			/*! Converts a move string into the internal move string sequence
			 * representation.
			 * 
			 * @param moveString Has to be a valid move string.
			 */
		virtual MoveSequence 
		getSequence(const std::string& moveString) const;
		
			/*! Converts an internal move string sequence representation
			 * into a normalized form, that is the same for all symmetric
			 * structures.
			 * 
			 * @param moveSeq the move string representation to normalize
			 * @return	the normalized move string representation 
			 */
		virtual MoveSequence 
		normalizeSequence(const MoveSequence& moveSeq) const;
		
			/*! Converts an internal move string sequence representation
			 * into all symmetric forms.
			 * 
			 * @param moveSeq the move string representation of interest
			 * @return	all symmetric move string representation 
			 */
		virtual std::set< MoveSequence > 
		getAllSymmetricSequences(const MoveSequence& moveSeq) const;
		
		virtual LatticeDescriptor& operator= (const LatticeDescriptor &ld2);

		bool	operator== (const LatticeDescriptor &ld2) const;
		bool	operator!= (const LatticeDescriptor &ld2) const;
		
	};
	
}

#include "LatticeDescriptor.icc"

#endif /*LATTICEDESCRIPTOR_H_*/
