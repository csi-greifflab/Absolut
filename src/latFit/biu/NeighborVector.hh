// $Id: NeighborVector.hh,v 1.2 2016/08/08 12:42:01 mmann Exp $
#ifndef BIU_NEIGHBORVECTOR_HH_
#define BIU_NEIGHBORVECTOR_HH_


#include <biu/SquareMatrix.hh>
#include "biu/Point.hh"
#include "biu/Alphabet.hh"

namespace biu
{
		//! a matrix to calculate lattice specific automorphisms like rotation
	typedef biu::SquareMatrix<int,3> Automorphism;

		//! a std::vector of Automorphism objects 
	typedef std::vector<Automorphism> AutomorphismVec;	

		//! an alphabet over Move strings 
	typedef Alphabet MoveAlphabet;

		//! an internal single move string representation 
	typedef MoveAlphabet::AlphElem Move;

		//! an internal move string sequence representation 
	typedef MoveAlphabet::Sequence MoveSequence;
	
		
		/**
		 * A NeighborVector manages lattice specific neighborhood data.  It
		 * handles the vector dependend move string and the related rotation
		 * matrix to generate relative move strings.
		 * 
		 * @author Martin Mann, Sebastian Will, Andreas Richter
		 */
	class NeighborVector : public IntPoint {
	private:
		Move	move;	//!< the vector specific internal move representation
		
			//! the rotation automorphism to convert relative to absolute moves
		Automorphism	rel2absRotation;	
		
			//! the rotation automorphism to convert absolute to relative moves
		Automorphism	abs2relRotation;	
		
		
	public :
		NeighborVector(const int x, const int y, const int z,
					   const Move _move,
					   const Automorphism& rel2absRot,
					   const Automorphism& abs2relRot) 
			:	IntPoint(x,y,z), move(_move), rel2absRotation(rel2absRot), 
				abs2relRotation(abs2relRot)
		{}
		
		NeighborVector(const IntPoint& point) : IntPoint(point), move(0) 
		{}
		
		virtual ~NeighborVector() 
		{}
		
		const Move& getMove() const { 
			return move; 
		}
		
		const Automorphism& getRel2AbsRotation() const {
			return rel2absRotation; 
		}
		const Automorphism& getAbs2RelRotation() const {
			return abs2relRotation;
		}
		
	};
	
		//! a std::set of NeighborVector objects 
	typedef std::set<NeighborVector> NeighSet;
	
	

} // namespace biu

#endif /*NEIGHBORVECTOR_HH_*/
