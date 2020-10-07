// $Id: LatticeFrame.cc,v 1.2 2016/08/08 12:41:57 mmann Exp $

#include <biu/LatticeFrame.hh>
#include <climits>

namespace biu
{


	LatticeFrame::LatticeFrame(	const LatticeDescriptor* const latDescriptor, 
								const unsigned int frameSize_)
		:	LatticeModel(latDescriptor), 
			frameSize(frameSize_), 
			frameSizeP2(frameSize*frameSize),
			indexedNeighborhood()
	{
		assertbiu(INT_MAX >= (frameSize*frameSize*frameSize), 
			"the frame size is to big to assure a one-to-one indexing");
		this->setFrameSize( frameSize_ );
	}

	LatticeFrame::LatticeFrame(const LatticeFrame& toCopy)
		:	LatticeModel(toCopy.latDescriptor), 
			frameSize(toCopy.frameSize), 
			frameSizeP2(toCopy.frameSizeP2),
			indexedNeighborhood(toCopy.indexedNeighborhood)
	{
	}
	
	LatticeFrame::~LatticeFrame()
	{
	}
	
	void 
	LatticeFrame::setFrameSize(const unsigned int frameSize_) {
		frameSize = frameSize_; 
		frameSizeP2 = frameSize*frameSize;
		assertbiu(INT_MAX >= (frameSize*frameSize*frameSize), 
			"the frame size is to big to assure a one-to-one indexing");
		
		  // update the list of neighborhood indices
		this->indexedNeighborhood.clear();
		for(LatticeNeighborhood::const_iterator it = latDescriptor->getNeighborhood().begin();
			it != latDescriptor->getNeighborhood().end(); it++) 
		{
			this->indexedNeighborhood.insert(getIndex(*it)); 
		}

	}	
	
	bool 
	LatticeFrame::isInFrame(const IntPoint& point) const {
		return	   point.getX() >= 0 && (unsigned int)point.getX() < frameSize 
				&& point.getY() >= 0 && (unsigned int)point.getY() < frameSize
				&& point.getZ() >= 0 && (unsigned int)point.getZ() < frameSize;
	}
		
	LatticeFrame::index_type 
	LatticeFrame::getIndex(const IntPoint& point) const {
		return point.getX() + point.getY()*frameSize + point.getZ()*frameSizeP2;
	}
	
	IntPoint 
	LatticeFrame::getPoint(const index_type& index) const {
		return IntPoint(	index%frameSize, 
							(index/frameSize)%frameSize, 
							(index/frameSizeP2)%frameSize);
	}


	bool	
	LatticeFrame::operator== (const LatticeFrame &lf2) const {
		bool retVal = (latDescriptor == lf2.latDescriptor)
						&& (frameSize == lf2.frameSize);
		return retVal;
	}
	
	bool	
	LatticeFrame::operator!= (const LatticeFrame &lf2) const {
		return !(this->operator ==( lf2 ));
	}
	
	const LatticeFrame::index_set&
	LatticeFrame::getIndexedNeighborhood() const {
		return indexedNeighborhood;
	}
	
	MoveSequence 
	LatticeFrame::indicesToAbsMoves(const std::vector<index_type> indVec) const{
		IPointVec points;
		for (std::vector<index_type>::const_iterator it = indVec.begin();
				it != indVec.end(); it++) {
			points.push_back(getPoint(*it));
		}
		return this->pointsToAbsMoves(points);
	}	
	
} // namespace biu

std::ostream & operator <<(std::ostream &os, biu::LatticeFrame::index_set &x) {
	os <<"( ";
	for (biu::LatticeFrame::index_set::const_iterator it = x.begin(); it != x.end(); it++)
		os <<*it <<", ";
	os <<")";
	return os;
}
	
std::ostream & operator <<(std::ostream &os, biu::LatticeFrame::index_vec &x) {
	os <<"( ";
	for (biu::LatticeFrame::index_vec::const_iterator it = x.begin(); it != x.end(); it++)
		os <<*it <<", ";
	os <<")";
	return os;
}
		
