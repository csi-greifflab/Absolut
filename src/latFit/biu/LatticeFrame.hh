// $Id: LatticeFrame.hh,v 1.2 2016/08/08 12:41:58 mmann Exp $
#ifndef BIU_LATTICEFRAME_HH_
#define BIU_LATTICEFRAME_HH_


#include "biu/LatticeModel.hh"

namespace biu
{
	

		/**
		 * A lattice frame is a subset of a lattice in which all points
		 * coordinates XYZ are in the range ZERO <= XYZ < FRAMESIZE.
		 * 
		 * The lattice frame provides a one-to-one mapping of the managed
		 * points to an index and vice versa.
		 * 
		 * @author Martin Mann, Sebastian Will, Andreas Richter
		 */
	class LatticeFrame : public LatticeModel
	{
	public:
	
		typedef int index_type;			//!< the type of the managed indices
		typedef std::set<biu::LatticeFrame::index_type> index_set; //!< the type of a set of indices
		typedef std::vector<biu::LatticeFrame::index_type> index_vec; //!< the type of a set of indices
		
	protected:
		
		unsigned int	frameSize;		//!< the size of the frame
		unsigned int	frameSizeP2;	//!< short access for (frameSize^2)
		index_set		indexedNeighborhood; //!< the indices of the neighboring vectors
	
	public:
		
			//! construction
			//! @param _latDescriptor This lattice property handler has to
		    //! be != NULL.
		    //! @param frameSize the size of the lattice frame
		LatticeFrame(	const LatticeDescriptor* const _latDescriptor, 
						const unsigned int frameSize);
						
		LatticeFrame(const LatticeFrame& toCopy);
		virtual ~LatticeFrame();

			//! Returns the current size of the lattice frame.		
		unsigned int getFrameSize() const { return frameSize; }
		
			//! Resizes the frame.
			//!
			//! All indices, created with the old frame size will be
			//! incompatible with the new lattice frame!
			//! @param frameSize the new frame size to set
		void setFrameSize(const unsigned int frameSize);
		
			//! Returns the center of the lattice frame.
			//! @return the centered position of the frame
		IntPoint getCenter() const  { 
			unsigned int cen = frameSize/2 ; 
			return IntPoint(cen, cen, cen); 
		}
		
			//! Returns the highest index value available in this frame
		index_type getMaxIndex() const { 
			return getIndex(IntPoint(frameSize-1,frameSize-1,frameSize-1));
		}

			//! Tests whether or not a point is in the frame or not.
		bool isInFrame(const IntPoint& point) const;
		
			//! Converts a point to its index.
			//! The Point has to be inside the frame borders 
			//! (isInFrame() == true).
		index_type getIndex(const IntPoint& point) const;
		
			//! Converts an index to its point coordinates.
			//! Assure that the index was created with the same frame
			//! size. Otherwise the convertion will lead to a differnt 
			//! point.
		IntPoint getPoint(const index_type& index) const;
		
			//! Converts a sequence of indices to the corresponding
			//! absolute move sequence.
		MoveSequence 
		indicesToAbsMoves(const std::vector<index_type> indVec) const;
		
			//! Returns the indices of the neighbor vectors, this
			//! lattice frame is based on.
		const index_set& getIndexedNeighborhood() const;
		
		bool	operator== (const LatticeFrame &lf2) const;
		bool	operator!= (const LatticeFrame &lf2) const;
		
	}; // class LatticeFrame
	
}	// namespace biu

std::ostream & operator <<(std::ostream &os, biu::LatticeFrame::index_set& x);
std::ostream & operator <<(std::ostream &os, biu::LatticeFrame::index_vec& x);

#endif /*LATTICEFRAME_HH_*/
