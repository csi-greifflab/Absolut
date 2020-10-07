// $Id: LatticeModel.cc,v 1.2 2016/08/08 12:42:00 mmann Exp $

#include <biu/LatticeModel.hh>
#include <biu/assertbiu.hh>

namespace biu
{
	
	LatticeModel::LatticeModel(const LatticeDescriptor* const _latDescriptor) 
	 :	latDescriptor(_latDescriptor), 
		latNeighborhood(latDescriptor->getNeighborhood())
	{
		assertbiu(latDescriptor != NULL, 
			"tried to create a LatticeModel without a LatticeDescriptor");
	}
	
	LatticeModel::LatticeModel(const LatticeModel& toCopy) :
		// pointer to latDescriptor only copied, because not deleted by LatticeModel 
		latDescriptor(toCopy.latDescriptor),	
		latNeighborhood(latDescriptor->getNeighborhood())
	{
		assertbiu(latDescriptor != NULL, 
			"tried to create a LatticeModel without a LatticeDescriptor");
	}
	
	LatticeModel::~LatticeModel()
	{
	}

	// ------------------------------------------------------------
	// conversion of move sequences
	//
	
	MoveSequence
	LatticeModel::relMovesToAbsMoves( const MoveSequence &relMoves ) const {
		
		MoveSequence absMoves(relMoves.size()); // here we store the result
		
		// base = matrix of first move
		Automorphism base=latNeighborhood.getElement(0).getRel2AbsRotation();
		
		for (MoveSequence::size_type i=0; i<relMoves.size(); ++i) {
			
			base = base
				* latNeighborhood.getElement(relMoves[i]).getRel2AbsRotation();
				
			IntPoint absVec = base*latNeighborhood.getElement(0); // the absolute move vector
			
			// get abs move corresponding to absVec
			absMoves[i] = latNeighborhood.getElement(absVec).getMove();
		}
		
		return absMoves;
	}
	
	MoveSequence
	LatticeModel::absMovesToRelMoves( const MoveSequence &absMoves ) const {
		
		MoveSequence relMoves(absMoves.size()); // here we store the result
    
		// set base, s.t. first absolute move points "forward"
		Automorphism base 
			= latNeighborhood.getElement(absMoves[0]).getAbs2RelRotation();
    
    	Move relMove;
    	
		for (MoveSequence::size_type i=0; i<absMoves.size(); ++i) {

			// multiply Matrix by vector<int> (results in vector<int>) 
			const IntPoint relVec = base * latNeighborhood.getElement(absMoves[i]);
			
			// get move with vector relVec
			relMove = latNeighborhood.getElement(relVec).getMove();
			relMoves[i] = relMove;
			
			base = latNeighborhood.getElement(relMove).getAbs2RelRotation()
				* base;
			
		}
		return relMoves;
	}

	//
	// ------------------------------------------------------------

	
	// ------------------------------------------------------------
	// conversion moves <=> coordinates
	//

	IntPoint 
	LatticeModel::applyAbsMove( const IntPoint& actPoint, 
								const Move& absMove) const {
		assertbiu(latDescriptor != NULL, 
			"LatticeModel has no LatticeDescriptor");
		
		LatticeNeighborhood::const_iterator neighbor;
		LatticeNeighborhood::const_iterator end = latNeighborhood.end();
		for (neighbor = latNeighborhood.begin(); neighbor != end; neighbor++)
			if (neighbor->getMove() == absMove) {
				return actPoint + (*neighbor);
			}
		// fehlerfall .. absmove nicht gefunden!
		assertbiu(neighbor != end,
				  "cant find the absolute move in the neighborhood"
				  "of the LatticeDescriptor");
		return IntPoint(actPoint);
	}
	
	
	IPointVec 
	LatticeModel::absMovesToPoints( const MoveSequence& absMoveSeq ) const {
		assertbiu(latDescriptor != NULL, 
			"LatticeModel has no LatticeDescriptor");
		
		IPointVec points(absMoveSeq.size()+1);	// rueckgabe vector
		points[0] = biu::IntPoint(0,0,0); // init mit lattice center
		for (MoveSequence::size_type i=0; i< absMoveSeq.size(); i++){
				// berechne naechsten auf basis des letzten punktes
			points[i+1] = applyAbsMove( points[i], absMoveSeq[i] );	
		}
		return points;
	}
	
	IPointVec 
	LatticeModel::relMovesToPoints( const MoveSequence& relMoves) const {
		
		IPointVec points(relMoves.size()+1);	// rueckgabe vector
		points[0] = biu::IntPoint(0,0,0); // init mit lattice center
		return relMovesToPoints(relMoves, points);
		
	}
	
	IPointVec& 
	LatticeModel::relMovesToPoints( const MoveSequence& relMoves, 
			IPointVec& toFill) const 
	{
		assertbiu(latDescriptor != NULL, 
			"LatticeModel has no LatticeDescriptor");

		// base = matrix of first move
		Automorphism base=latNeighborhood.getElement(0).getRel2AbsRotation();
		
		toFill.resize(relMoves.size()+1);	// rueckgabe vector
		toFill[0] = biu::IntPoint(0,0,0); // init mit lattice center

		for (MoveSequence::size_type i=0; i<relMoves.size(); ++i) {
			
			base = base
				* latNeighborhood.getElement(relMoves[i]).getRel2AbsRotation();
				
				// add the absolute move vector to last position
			toFill[i+1] = toFill[i] + 
							base*latNeighborhood.getElement(relMoves[0]);
		}
		
		return toFill;
	}
	
	MoveSequence 
	LatticeModel::pointsToAbsMoves( const IPointVec& points ) const {
		assertbiu(points.size() > 0, "points vector is empty");
		MoveSequence absmoves(points.size()-1);	 // erstelle vector der entsprechenden groesse
		MoveSequence::size_type i = 0;
		biu::IPointVec::const_iterator act = points.begin(), last = act;
		for (act++; act != points.end(); last++, act++) {
				// haenge nacheinander moves an
			absmoves[i++] = getAbsMove(*last, *act);	
		}
		return absmoves;
	}
	
	IPointSet 
	LatticeModel::getAllNeighPoints(const IntPoint& center) const
	{ 
		assertbiu(latDescriptor != NULL, 
			"LatticeModel has no LatticeDescriptor");
		IPointSet retSet;
		for (LatticeNeighborhood::const_iterator it = latNeighborhood.begin(); 
					it != latNeighborhood.end(); it++) {
			retSet.insert(center + (*it));
		}
		return retSet;
	}

} // namespace biu
