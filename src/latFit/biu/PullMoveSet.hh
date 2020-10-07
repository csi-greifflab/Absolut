// $Id: PullMoveSet.hh,v 1.2 2016/08/08 12:42:00 mmann Exp $
#ifndef BIU_PULLMOVESET_HH_
#define BIU_PULLMOVESET_HH_

#include "biu/LatticeMoveSet.hh"
#include "biu/LatticeProtein.hh"
#include <vector>
#include <map>

namespace biu
{

/*!
 * This class handles pull moves on LatticeProtein objects.
 * 
 * @author Daniel Maticzka
 */
class PullMoveSet : public LatticeMoveSet
{
public:
	
	  //! just a point of IntPoint objects
	typedef std::pair< IntPoint, IntPoint > IPointPair;
	
	
	/*!
	 * This class generates all pull positions for the pull move set from the 
	 * given LatticeDescritor.
	 * 
	 * NOTE : it utilizes always 2 new points (C,L) for each move, which is only
	 *  needed for rectangular lattices 
	 * 
	 * @author Daniel Maticzka
	 */
	class PullMoveDecoder 
	{
	protected:
		//! the lattice in use
		const biu::LatticeModel* const lattice;
		
		/*! true if lattice model supports triangular moves */
		bool triangularLattice;

		/*! stores all possible neighbor vectors as reference */
		std::vector< IntPoint > neighDirections;
		std::map< IntPoint, size_t > neighDirectionsAccess;
//		IPointVec neighDirections; 		

		/*! stores the C and L pull points */
		std::vector< std::vector< IPointPair > > stdPulls;
		
		/*! stores the end pull pull points */
		std::vector< IPointPair > endPulls;
		
		/*! number of possible neighbor vectors */
		size_t directionNumber;
		
		/*! number of standard pulls */
		size_t stdPullNumber;
		/*! number of end pulls */
		size_t endPullNumber;
		
	public:
		/*!
		 * Constructs a PullMoveDecoder using the given lattice.
		 */
		PullMoveDecoder(const LatticeModel* const lattice);
		
		/*!
		 * Decodes a given moveIndex and returns the values pullFront, stdPull
		 * and movePosition.
		 */
		void lookupMove(const size_t& proteinLength, 
				size_t& moveIndex, bool& pullFront, 
				bool& stdPull, size_t& movePosition) const;

		/*!
		 * Returns c and l pull points for given rootpoint and fixpoint where
		 * rootpoint is the position which will be pulled and fixpoint is the 
		 * position which will stay fixed.
		 */
		IPointPair 
		lookupStdMove(const IntPoint& rootpoint, 
				const IntPoint& fixpoint, 
				const size_t& moveIndex) const;
		
		/*!
		 * Returns the two points where the tail will be pulled to.
		 */
		IPointPair 
		lookupEndMove(const IntPoint& rootpoint,
				const size_t& moveIndex) const;
		
		/*!
		 * Returns the number of possible moves for a protein of length
		 * proteinLength.
		 */
		size_t getMoveNumber(const size_t proteinLength) const; 
		
		/*!
		 * Returns the number of possible end moves.
		 */
		size_t getEndMoveNumber() const; 
		
		/*!
		 * Returns the number of possible pull moves for any neighboring vector.
		 */
		size_t getPullMoveNumber() const; 
		
		const LatticeModel* const getLattice(void) const;

		/*!
		 * Returns true if the lattice contains triangular elements
		 */
		const bool
		isTriangular() const;

	};
	
private:
	/*!
	 * This Class implements some basic functions on int. The "direction" of
	 * those functions can be skipped for all objects of this class.
	 * 
	 * @author Daniel Maticzka
	 */
	class RelativeInt
	{
	private:
		/*!
		 * This class defines all funcitons needed by RelativeInt.
		 * 
		 * @author Daniel Maticzka
		 */
		class Operators {
		public:
			virtual ~Operators() {}
			
			/*!
			 * @param a		first argument of addition/subtraction
			 * @param b 	second argument of addition/subtraction
			 * @return		sum/subtraction of a and b depending on implementation
			 */
			virtual const int add(const int& a, const int& b) const = 0;
			
			/*!
			 * @param a		first argument of addition/subtraction
			 * @param b 	second argument of addition/subtraction
			 * @return		sum/subtraction of a and b depending on implementation
			 */
			virtual const int sub(const int& a, const int& b) const = 0;
			
			/*!
			 * @param a		first argument of comparison
			 * @param b 	second argument of comparison
			 * @return		logical value of implemented comparison
			 */
			virtual const bool largerEq(const int& a, const int& b) const = 0;
			
			/*!
			 * @param a		first argument of comparison
			 * @param b 	second argument of comparison
			 * @return		logical value of implemented comparison
			 */
			virtual const bool smallerEq(const int& a, const int& b) const = 0;
		};
		
		/*!
		 * This subclass of Operators implements all funtions with their
		 * standard meaning.
		 * 
		 * @author Daniel Maticzka
		 */
		class StdOps : public Operators {
		public:
			virtual const int add(const int& a, const int& b) const 
				{ return a + b; }
			virtual const int sub(const int& a, const int& b) const
				{ return a - b; }
			virtual const bool largerEq(const int& a, const int& b) const 
				{ return a >= b; }
			virtual const bool smallerEq(const int& a, const int& b) const 
				{ return a <= b; }
		};
		
		/*!
		 * This sublcass of Operators implements all functions with their
		 * ivnerted meaning.
		 * 
		 * @author Daniel Maticzka
		 */
		class InvOps : public Operators {
			virtual const int add(const int& a, const int& b) const 
				{ return a - b; }
			virtual const int sub(const int& a, const int& b) const 
				{ return a + b; }
			virtual const bool largerEq(const int& a, const int& b) const 
				{ return a <= b; }
			virtual const bool smallerEq(const int& a, const int& b) const 
				{ return a >= b; }
		};

		/*! computation is delegated to this object */
		Operators* ops;
		/*! the value of the RelativeInt, it's an int! */
		int value;
		
	public:
		/*!
		 * Constructs a RelativeInt with value _value and direction
		 */
		RelativeInt(int _value, bool direction) : value(_value) {
			if (direction)
				ops=new PullMoveSet::RelativeInt::StdOps();
			else
				ops=new PullMoveSet::RelativeInt::InvOps();
		}
		
		~RelativeInt() { delete ops; }
		
		const bool
		operator >= (const int& y) const { return ops->largerEq(value,y); }
		
		const bool
		operator <= (const int& y) const { return ops->smallerEq(value,y); }	
		
		const int
		operator+ (const int& y) const { return ops->add(value,y); }
		
		const int
		operator- (const int& y) const { return ops->sub(value,y); }
		
		void
		operator --(int) { value = ops->sub(value,1); }
		
		void
		operator ++(int) { value = ops->add(value,1); }
		
		const int
		getValue() const { return value; };
	};
	
protected:
	struct UndoRecord
	{
		//! pointer to the last object that was changed
		LatticeProtein*	lastChangedObject;
		//! true if protein has been changed by last move
		bool hasChanged;
		//! pull direction at last pull
		bool pullFront;
		//! position at last pull
		size_t position;
		//! last position pulled at last pull
		size_t lastPullPosition;
		//! first overwritten position at last pull
		IntPoint lostPos0;
		//! second overwritten position at last pull
		IntPoint lostPos1;
	} undoRec;
	
	  //! the current pull move decoder in use
	PullMoveSet::PullMoveDecoder* decoder;
	
	  //! whether or not the decoder is shared and has to be copied or not
	const bool decoderIsShared;
	
public:

	/*!
	 * Constructs a PullMoveSet object operating on the LatticeModel lattice.
	 * @param lattice the lattice in use
	 */
	PullMoveSet(const LatticeModel* lattice); 
	
	/*!
	 * Constructs a PullMoveSet object operating lattice model represented by
	 * the provided decoder.
	 *
	 * @param decoder the decoder to use 
	 * @param decoderIsShared whether or not the decoder is shared 
	 *         among all derived instances or has to be copied
	 */
	PullMoveSet(PullMoveDecoder* decoder, const bool decoderIsShared);
	
	PullMoveSet(const PullMoveSet& moveSet);
	PullMoveSet& operator=(const PullMoveSet& moveSet2);
	
	PullMoveSet* clone();

	virtual ~PullMoveSet();

	/*! Applies a pull move to the given LatticeProtein on the 
	 * specified position and returns a new LatticeProtein object 
	 * as result.
	 * 
	 * @param todo the initial state
	 * @param moveIndex the move index (has to be less than the return 
	 * 					value of getMoveNumber())
	 * @return a new LatticeProtein object neighbored to the todo object
	 *
	 *  In case of invalid moves, the behaviour
	 *  of the program is not specified.
	 */
	virtual LatticeProtein* applyMove(const LatticeProtein* const todo, 
			const size_t moveIndex);

	/*! Applies a pull move to the given LatticeProtein on the 
	 * specified position INPLACE and returns the modified object
	 * as result.
	 *
	 * @param todo the initial state to modify
	 * @param moveIndex the move index (has to be less than the return 
	 * 					value of getMoveNumber())
	 * @return the modified LatticeProtein or NULL if the move is not valid
	 */
	virtual LatticeProtein* applyMoveInPlace(LatticeProtein* const todo, 
			const size_t moveIndex);

	/*!
	 * Returns number of moves that can be applied to the given lattice.
	 */
	size_t getMoveNumber ( const LatticeProtein* const lp) const; 

	/*!
	 * Undo of last pull move performed on LatticeProtein toUndo. Undo
	 * is performed in place!
	 */
	virtual LatticeProtein* undoLastMove(LatticeProtein* toUndo);
	
	//! Access to the used PullMoveDecoder.
	//! @return the current decode in use
	virtual const PullMoveDecoder* const getDecoder(void) const;
};

} // end namespace biu

#endif /*PULLMOVESET_HH_*/
