// $Id: LatticeDescriptor.cc,v 1.2 2016/08/08 12:41:57 mmann Exp $

#include <biu/LatticeDescriptor.hh>
#include <algorithm>

namespace biu
{

	/*
	biu::IntPoint operator*(const biu::Automorphism& am, 
							const biu::IntPoint& p)  {
		return IntPoint(
			am[0][0]*p.getX() + am[0][1]*p.getY() + am[0][2]*p.getZ(),
			am[1][0]*p.getX() + am[1][1]*p.getY() + am[1][2]*p.getZ(),
			am[2][0]*p.getX() + am[2][1]*p.getY() + am[2][2]*p.getZ()
			);
	}
	*/

	/************************************************
	 * LatticeDescriptor
	 ***********************************************/
	
	LatticeDescriptor::LatticeDescriptor(const std::string& name_)  
		: name(name_), moveAlphabet(NULL), latNeighborhood(NULL)
	{
	}
	
	LatticeDescriptor::LatticeDescriptor(const LatticeDescriptor& toCopy)
		: name(toCopy.name),
		  automorphisms(toCopy.automorphisms), 
		  latBase(toCopy.latBase),
		  moveAlphabet(new MoveAlphabet(*(toCopy.moveAlphabet))),
		  latNeighborhood(new LatticeNeighborhood(*(toCopy.latNeighborhood)))
	{
	}
	
	LatticeDescriptor&
	LatticeDescriptor::operator= (const LatticeDescriptor &ld2) {
		if (this != &ld2) {
			name			= ld2.name;
			automorphisms	= ld2.automorphisms;
			latBase 		= ld2.latBase;
			if (moveAlphabet != NULL)
				delete moveAlphabet;
			moveAlphabet	= new MoveAlphabet(*(ld2.moveAlphabet));
			if (latNeighborhood != NULL)
				delete latNeighborhood;
			latNeighborhood	= new LatticeNeighborhood(*(ld2.latNeighborhood));
		}
    	return *this;
	}
	
	
	LatticeDescriptor::~LatticeDescriptor()
	{
		if (moveAlphabet != NULL) {
			delete moveAlphabet;	moveAlphabet = NULL;
		}
		if (latNeighborhood != NULL) {
			delete latNeighborhood;	latNeighborhood = NULL;
		}
	}
	
	
	bool	
	LatticeDescriptor::operator== (const LatticeDescriptor &ld2) const {
		return	this->name == ld2.name
				&& *(this->moveAlphabet) == *(ld2.moveAlphabet)
				&& this->automorphisms == ld2.automorphisms
				&& this->latBase == ld2.latBase
				&& std::equal(	this->latNeighborhood->begin(), 
								this->latNeighborhood->end(),
								ld2.latNeighborhood->begin());
	}

	bool	
	LatticeDescriptor::operator!= (const LatticeDescriptor &ld2) const {
		return	this->name != ld2.name
				|| *(this->moveAlphabet) != *(ld2.moveAlphabet)
				|| this->automorphisms != ld2.automorphisms
				|| this->latBase != ld2.latBase
				|| !(std::equal(	this->latNeighborhood->begin(), 
									this->latNeighborhood->end(),
									ld2.latNeighborhood->begin()));
	}
	


	void 
	LatticeDescriptor::initNeighborhood() {
		unsigned int n=getNeighborDataSize();
		
		std::vector<std::string> alphabetVec;
		
		const NeighborData *data = getNeighborData();
		
		NeighSet neighbors;
		for (unsigned int i=0; i<n; i++) {
			const NeighborData &d = data[i]; // current entry in data

			alphabetVec.push_back(d.name);
				
			Automorphism mat(d.mat);
			Automorphism invmat(d.invmat);
			
			neighbors.insert(NeighborVector(d.vec[0],d.vec[1],d.vec[2],
											i,
											mat,
											invmat
								 )
				);
		}
		
		// move alphabet
		if (moveAlphabet != NULL) delete moveAlphabet;
		moveAlphabet = new MoveAlphabet(alphabetVec);
		
		if (latNeighborhood != NULL) delete latNeighborhood;
		latNeighborhood = new LatticeNeighborhood(moveAlphabet, neighbors);
		
		assertbiu(neighbors.size() == getNeighborDataSize(),
	"neighborhood size and the initialized neighbor number have to be equal");
		
	}


	void 
	LatticeDescriptor::initAutomorphisms() {
		unsigned int n=getAutomorphismDataSize();
		
		
		// initializing automorphisms
		const AutomorphismData *data = getAutomorphismData();
		
		automorphisms.clear();
		for (unsigned int i=0; i<n; i++) {
			int d[3][3];
			for(unsigned int k=0; k<3; k++) {
				d[0][k] = data[i]._0[k];
				d[1][k] = data[i]._1[k];
				d[2][k] = data[i]._2[k];
			}
			
			Automorphism a(d);
			automorphisms.push_back(a);
		}
		
		// initializing move convertion tables for normalization
		symMoveReplacement.clear();
		
		for (unsigned int i=0; i<n; i++) {
			MoveSequence symMoves;
			  // generate symmetric neighborvector
			for (unsigned int k=0; k<latNeighborhood->size(); k++) {
				IntPoint symNV = automorphisms[i] * latNeighborhood->getElement(k);
				  // find for each symmetric neighbor vector the corresponding
				  // index in the neighborhood
				for (unsigned int l=0; l<latNeighborhood->size(); l++)
					if ( symNV == latNeighborhood->getElement(l)) {
						  // store the symmetric index
						symMoves.push_back(l);
						break;
					}
			}
			  // ensure that for each move a symmetric one is stored
			assertbiu(symMoves.size() == latNeighborhood->size(),
			"calculated symmetric move vector differs in length with neighbor vectors");
			  // save symmetric change table
			symMoveReplacement.push_back(symMoves);
		}
	}

	MoveSequence 
	LatticeDescriptor::
	normalizeSequence(const MoveSequence& origMoveSeq) const {
		  // new normalized sequence to return
		MoveSequence newSeq = MoveSequence(origMoveSeq);
		MoveSequence smallestSeq = MoveSequence(origMoveSeq);

		  // find smallest move string reprentation within all symmetric move strings
		for (std::vector<MoveSequence>::const_iterator sm = symMoveReplacement.begin();
				sm != symMoveReplacement.end(); sm++)
		{
			bool smaller=false;
			for (unsigned int i=0; i<smallestSeq.size(); i++) {
				if (smaller || (*sm)[origMoveSeq[i]] <= smallestSeq[i]) {
				    smaller = smaller || (*sm)[origMoveSeq[i]] < smallestSeq[i];
				    newSeq[i] = (*sm)[origMoveSeq[i]];
				}
				else {
					smaller = false;
					break;
				}
			}
			if (smaller) {
				smallestSeq = newSeq;
			}
		}
		
		return smallestSeq;
	}

	std::set< MoveSequence >
	LatticeDescriptor::
	getAllSymmetricSequences( const MoveSequence& origMoveSeq ) const
	{
		std::set< MoveSequence > allSymms;
		
		
		  // find smallest move string reprentation within all symmetric move strings
		for (std::vector<MoveSequence>::const_iterator sm = symMoveReplacement.begin();
				sm != symMoveReplacement.end(); sm++)
		{
			  // create symmetric sequence to add
			MoveSequence newSeq = MoveSequence(origMoveSeq);
			for (size_t i=0; i<origMoveSeq.size(); i++) {
				  // replace with symmetric move
			    newSeq[i] = (*sm)[origMoveSeq[i]];
			}
			  // add to storage if not present
			allSymms.insert(newSeq);
		}
		
		return allSymms;
	}


}	// namespace biu
