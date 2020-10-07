
#include "biu/LatticeProteinUtil.hh"

#include <biu/assertbiu.hh>

#include <cmath>

namespace biu
{


	////////////////////////////////////////////////////////////////////////////


	biu::MoveSequence
	LatticeProteinUtil::
	toMoveSequence(	const std::string & moveString
					, const biu::LatticeDescriptor & latDescr
					, const bool sideChain )
	{
		biu::MoveSequence moves;
		const size_t mLength = latDescr.getAlphabet()->getElementLength();
		
		  // side chain model
		if (sideChain) {
			  // valid move string length is == ( l*(2*n -1) + 2*n )
			assertbiu( (moveString.size()+mLength) % (2*(mLength+1)) == 0, "side chain move string length is no multiple of move encoding length");
			  // initialize move string length
			moves.resize(((moveString.size()+mLength) / (mLength+1))-1,0);
			bool nextIsSideChain = true;
			  // current move string position to process
			size_t cur=0;
			  // parse move string
			for (size_t i=0; i<moveString.size(); ) {
				  // handle next monomer
				if (nextIsSideChain) {
					assertbiu(moveString.at(i)=='(', "side chain expected but not opening bracket found");
					assertbiu(moveString.at(i+1+mLength)==')', "side chain expected but not closing bracket found");
					moves[cur] = latDescr.getAlphabet()->getElement(moveString.substr(i+1,mLength));
					i += 2+mLength;
				} else {
					moves[cur] = latDescr.getAlphabet()->getElement(moveString.substr(i,mLength));
					i += mLength;
				}
				cur++;
				nextIsSideChain = !nextIsSideChain;
			}
			assertbiu( cur == (moves.size()), "filling of move sequence was not successful");
		}
		  // backbone-only model
		else {
			assertbiu( moveString.size() % mLength == 0, "move string length is no multiple of move encoding length");
			moves = latDescr.getAlphabet()->getSequence( moveString );
		}

		return moves;
	}

	////////////////////////////////////////////////////////////////////////////


	biu::MoveSequence
	LatticeProteinUtil::
	toMoveSequence(	const biu::IPointVec & p
					, const biu::LatticeDescriptor & latDescr
					, const bool sideChain )
	{
		assertbiu(p.size() > 1, "no structure given (length < 2)");
		biu::MoveSequence moves(p.size()-1,0);
		
		const biu::LatticeNeighborhood & neigh = latDescr.getNeighborhood();
		
		  // side chain model
		if (sideChain) {
			  // valid move string length is == ( l*(2*n -1) + 2*n )
			assertbiu( p.size()%2 == 0, "different number of side chain and backbone positions");
			  // append moves one after another
			for (size_t i=1; i<p.size(); i++) {
				if (i%2 == 0) { // backbone position
					moves[i-1] = neigh.getElement(p.at(i)-p.at(i-2)).getMove();
				} else { // side chain position
					moves[i-1] = neigh.getElement(p.at(i)-p.at(i-1)).getMove();
				}
			}
		}
		  // backbone-only model
		else {
			  // append moves one after another
			for (size_t i=1; i<p.size(); i++) {
				moves[i-1] = neigh.getElement(p.at(i)-p.at(i-1)).getMove();	
			}
		}
		
		return moves;
	}
	
	
	////////////////////////////////////////////////////////////////////////////
	
	std::string
	LatticeProteinUtil::
	toString(	const biu::MoveSequence & moves
				, const biu::LatticeDescriptor & latDescr
				, const bool sideChain )
	{
		std::string moveString;
		
		  // side chain model
		if (sideChain) {
			// valid move string length is == ( l*(2*n -1) + 2*n )
			assertbiu( moves.size() % 2 == 1, "different number of side chain and backbone positions");
			bool nextIsSideChain = true;
			  // convert into move string
			for (size_t i=0; i<moves.size(); i++) {
				  // handle next monomer
				if (nextIsSideChain) {
					moveString.append("(");
					moveString.append(latDescr.getAlphabet()->getString(moves.at(i)));
					moveString.append(")");
				} else {
					moveString.append(latDescr.getAlphabet()->getString(moves.at(i)));
				}
				nextIsSideChain = !nextIsSideChain;
			}
		}
		  // backbone-only model
		else {
			moveString = latDescr.getAlphabet()->getString( moves );
		}

		return moveString;
	}

	////////////////////////////////////////////////////////////////////////////

	biu::DPointVec
	LatticeProteinUtil::
	toDblPoints(	const biu::MoveSequence & moves
					, const biu::LatticeDescriptor & latDescr
					, const bool sideChain 
					, const double cAlphaDist
					)
	{
		  // set up coordinate container
		biu::DPointVec p(moves.size()+1, biu::DblPoint(0,0,0));
		
		const double baseScale = LatticeProteinUtil::getBaseScale( latDescr, cAlphaDist );
		const biu::LatticeNeighborhood& neigh = latDescr.getNeighborhood();
		
		  // side chain model
		if (sideChain) {
			  // fill scaled coordinates
			for (size_t i=0; i < moves.size(); i++){
					// get next neighbor
				const biu::NeighborVector& n = neigh.getElement( moves.at(i) );
				p[i+1].setX( p[i-(i%2)].getX() + (baseScale*double(n.getX())) );
				p[i+1].setY( p[i-(i%2)].getY() + (baseScale*double(n.getY())) );
				p[i+1].setZ( p[i-(i%2)].getZ() + (baseScale*double(n.getZ())) );
			}
			
		}
		  // backbone-only model
		else {
			  // fill scaled coordinates
			for (size_t i=0; i < moves.size(); i++){
					// get next neighbor
				const biu::NeighborVector& n = neigh.getElement( moves.at(i) );
				p[i+1].setX( p[i].getX() + (baseScale*double(n.getX())) );
				p[i+1].setY( p[i].getY() + (baseScale*double(n.getY())) );
				p[i+1].setZ( p[i].getZ() + (baseScale*double(n.getZ())) );
			}
		}
		
		return p;
	}


	////////////////////////////////////////////////////////////////////////////

	biu::IPointVec
	LatticeProteinUtil::
	toIntPoints(	const biu::MoveSequence & moves
					, const biu::LatticeDescriptor & latDescr
					, const bool sideChain )
	{
		  // set up coordinate container
		biu::IPointVec p( moves.size()+1, biu::IntPoint(0,0,0) );
		
		const biu::LatticeNeighborhood& neigh = latDescr.getNeighborhood();
		
		  // side chain model
		if (sideChain) {
			  // fill coordinates
			for (size_t i=0; i < moves.size(); i++){
					// get next neighbor
				p[i+1] = p[i-(i%2)] + neigh.getElement( moves.at(i) );
			}
			
		}
		  // backbone-only model
		else {
			  // fill coordinates
			for (size_t i=0; i < moves.size(); i++){
					// get next neighbor
				p[i+1] = p[i] + neigh.getElement( moves.at(i) );
			}
		}
		
		return p;
	}



	////////////////////////////////////////////////////////////////////////////

	double
	LatticeProteinUtil::
	getBaseScale(	const biu::LatticeDescriptor & latDescr
					, const double neighVecLength )
	{
		double baseScale = neighVecLength;
		
		if (latDescr.getName().compare("sqr") == 0) {
			baseScale = neighVecLength;
		} else if (latDescr.getName().compare("cub") == 0) {
			baseScale = neighVecLength;
		} else if (latDescr.getName().compare("fcc") == 0) {
			baseScale = sqrt( pow(neighVecLength, 2) / 2.0 );
		} else if (latDescr.getName().compare("210") == 0) {
			baseScale = sqrt( pow(neighVecLength, 2) / 5.0 );
		} else if (latDescr.getName().compare("100") == 0) {
			baseScale = neighVecLength;
		} else if (latDescr.getName().compare("110") == 0) {
			baseScale = sqrt( pow(neighVecLength, 2) / 2.0 );
		} else {
			assertbiu( false, "the given lattice descriptor is not handled");
		}
		
		return baseScale;
	}


	////////////////////////////////////////////////////////////////////////////
	

	double
	LatticeProteinUtil::
	cRMSD(	const biu::DPointVec & pos1
			, const biu::DPointVec & pos2 )
	{
		assertbiu(pos1.size() == pos2.size(), "pos1 and pos2 are of different size");
		
		  // simple check
		if (pos1.size() < 2) {
			return 0.0;
		}

		double squareDev = 0.0;
		for (size_t i=0; i<pos1.size(); i++) {
			squareDev += pow( pos1.at(i).distance( pos2.at(i) ) ,2);
		}
		
		return sqrt( squareDev / double(pos1.size()) );
	}
	

	////////////////////////////////////////////////////////////////////////////
	
		
	double
	LatticeProteinUtil::
	cRMSD(	const biu::IPointVec & pos1
			, const biu::IPointVec & pos2 )
	{
		assertbiu(pos1.size() == pos2.size(), "pos1 and pos2 are of different size");
		
		  // simple check
		if (pos1.size() < 2) {
			return 0.0;
		}

		double squareDev = 0.0;
		for (size_t i=0; i<pos1.size(); i++) {
			squareDev += pow( pos1.at(i).distance( pos2.at(i) ) ,2);
		}
		
		return sqrt( squareDev / double(pos1.size()) );
	}


	////////////////////////////////////////////////////////////////////////////
	

	double
	LatticeProteinUtil::
	dRMSD(	const biu::DPointVec & pos1
			, const biu::DPointVec & pos2 )
	{
		assertbiu(pos1.size() == pos2.size(), "pos1 and pos2 are of different size");
		
		  // simple check
		if (pos1.size() < 2) {
			return 0.0;
		}

		double squareDev = 0.0;
		for (size_t i=0; i<pos1.size(); i++) {
			for (size_t j=(i+1); j<pos1.size(); j++) {
				squareDev += pow( pos1.at(i).distance( pos1.at(j) ) 
									- pos2.at(i).distance( pos2.at(j) ) 
								, 2);
			}
		}
//std::cerr <<" dRMSD sum = " <<(double(pos1.size() * (pos1.size()-1))/2) <<"\n";
		return sqrt( 2*squareDev / double(pos1.size() * (pos1.size()-1)) );
	}

	
	////////////////////////////////////////////////////////////////////////////
	
		
	double
	LatticeProteinUtil::
	dRMSD(	const biu::IPointVec & pos1
			, const biu::IPointVec & pos2 )
	{
		assertbiu(pos1.size() == pos2.size(), "pos1 and pos2 are of different size");
		
		  // simple check
		if (pos1.size() < 2) {
			return 0.0;
		}

		double squareDev = 0.0;
		for (size_t i=0; i<pos1.size(); i++) {
			for (size_t j=(i+1); j<pos1.size(); j++) {
			squareDev += pow( pos1.at(i).distance( pos1.at(j) )
								- pos2.at(i).distance( pos2.at(j) ) 
							, 2);
			}
		}
		
		return sqrt( 2*squareDev / double(pos1.size() * (pos1.size()-1)) );
	}
	

	////////////////////////////////////////////////////////////////////////////
	

	double
	LatticeProteinUtil::
	GDT_TS(	const biu::DPointVec & pos1
			, const biu::DPointVec & pos2 )
	{
		assertbiu(pos1.size() == pos2.size(), "pos1 and pos2 are of different size");
	  /*
	   * GDT_TS = (GDT_P1 + GDT_P2 + GDT_P4 + GDT_P8)/4
	   * where GDT_Pn denotes percent of residues under distance cutoff <= n
	   */

		  // counters
		size_t GDT_P1 = 0, GDT_P2 = 0, GDT_P4 = 0, GDT_P8 = 0;
		
		for (size_t i=0; i<pos1.size(); i++) {
			  // calculate distance
			const double dist = pos1.at(i).distance(pos2.at(i));
			  // evaluate distance
			if (dist <= 1.0) {
				GDT_P1++;
				GDT_P2++;
				GDT_P4++;
				GDT_P8++;
			} else if (dist <= 2.0) {
				GDT_P2++;
				GDT_P4++;
				GDT_P8++;
			} else if (dist <= 4.0) {
				GDT_P4++;
				GDT_P8++;
			} else if (dist <= 8.0) {
				GDT_P8++;
			}
		}
		
		  // calculate the final GDT_TS
		return (	  double(GDT_P1)/double(pos1.size())
					+ double(GDT_P2)/double(pos1.size())
					+ double(GDT_P4)/double(pos1.size())
					+ double(GDT_P8)/double(pos1.size())
				) / 4.0 ;
	}


	////////////////////////////////////////////////////////////////////////////
	

	
	double
	LatticeProteinUtil::
	GDT_HA(	const biu::DPointVec & pos1
			, const biu::DPointVec & pos2 )
	{
		assertbiu(pos1.size() == pos2.size(), "pos1 and pos2 are of different size");
	  /*
	   * GDT_HA = (GDT_P0.5 + GDT_P1 + GDT_P2 + GDT_P4)/4,
	   * where GDT_Pn denotes percent of residues under distance cutoff <= n
	   */

		  // counters
		size_t GDT_P05 = 0, GDT_P1 = 0, GDT_P2 = 0, GDT_P4 = 0;
		
		for (size_t i=0; i<pos1.size(); i++) {
			  // calculate distance
			const double dist = pos1.at(i).distance(pos2.at(i));
			  // evaluate distance
			if (dist <= 0.5) {
				GDT_P05++;
				GDT_P1++;
				GDT_P2++;
				GDT_P4++;
			} else if (dist <= 1.0) {
				GDT_P1++;
				GDT_P2++;
				GDT_P4++;
			} else if (dist <= 2.0) {
				GDT_P2++;
				GDT_P4++;
			} else if (dist <= 4.0) {
				GDT_P4++;
			}
		}
		
		  // calculate the final GDT_HA
		return (	  double(GDT_P05)/double(pos1.size())
					+ double(GDT_P1)/double(pos1.size())
					+ double(GDT_P2)/double(pos1.size())
					+ double(GDT_P4)/double(pos1.size())
				) / 4.0 ;
	}


	////////////////////////////////////////////////////////////////////////////
	

}
