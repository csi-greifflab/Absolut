// $Id: BackboneStructure3D.hh,v 1.2 2016/08/08 12:41:59 mmann Exp $
#ifndef BIU_BACKBONESTRUCTURE3D_H_
#define BIU_BACKBONESTRUCTURE3D_H_

#include "biu/Point.hh"

namespace biu
{

		/*! A BackboneStructure3D object represents a linear structure in 
		 * the 3-dimensional space. 
		 * 
		 * @author Martin Mann, Sebastian Will, Andreas Richter
		 */	
	class BackboneStructure3D
	{
	public:
	
		BackboneStructure3D()
		{}
		
		virtual ~BackboneStructure3D() 
		{}
		
			//! Returns the consecutive point representation of the backbone.
		virtual DPointVec get3Ddata() const = 0;
		
			/*! Returns the distance root mean square deviation (DRMSD) to an 
			 * other BackboneStructure3D
			 */
		virtual double getDRMSD(const BackboneStructure3D& other) const = 0;
	};

} // namespace biu

#endif /*BACKBONESTRUCTURE3D_H_*/
