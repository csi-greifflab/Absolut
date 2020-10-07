// $Id: LatticeDescriptorSQR.cc,v 1.2 2016/08/08 12:41:57 mmann Exp $

#include <biu/LatticeDescriptorSQR.hh>
#include <cstdlib>

namespace biu
{

	LatticeDescriptorSQR::LatticeDescriptorSQR() 
		: LatticeDescriptor("sqr")
	{
		
		// gitterbasis
		latBase.push_back(IntPoint(1,0,0));	// x
		latBase.push_back(IntPoint(0,1,0));	// y
		latBase.push_back(IntPoint(0,0,1));	// z
		
		initNeighborhood();
		initAutomorphisms();
	}
	
	LatticeDescriptorSQR::~LatticeDescriptorSQR() {
	}

	unsigned int 
	LatticeDescriptorSQR::getNeighborDataSize() const {
		return 4;
	}
	
	const LatticeDescriptor::NeighborData *
	LatticeDescriptorSQR::getNeighborData() const {
			// init
		static const NeighborData squareNeighborData[] = {
			{"F",{ 1, 0, 0},
			 {{1,0,0},{0,1,0},{0,0,1}},  {{1,0,0},{0,1,0},{0,0,1}}},
			{"R",{ 0, 1, 0},
			 {{0,-1,0},{1,0,0},{0,0,1}}, {{0,1,0},{-1,0,0},{0,0,1}}},
			{"B",{-1, 0, 0},
			 {{-1,0,0},{0,-1,0},{0,0,1}},{{-1,0,0},{0,-1,0},{0,0,1}}},
			{"L",{ 0,-1, 0},
			 {{0,1,0},{-1,0,0},{0,0,1}}, {{0,-1,0},{1,0,0},{0,0,1}}},
		};
		
		return squareNeighborData;
	}

	unsigned int LatticeDescriptorSQR::getAutomorphismDataSize() const {
		return 8;
	}

	const LatticeDescriptorSQR::AutomorphismData *
	LatticeDescriptorSQR::getAutomorphismData() const {
		static const AutomorphismData data[] = {
			/* rotations */
			{{1,0,0},{0,1,0},{0,0,1}},   // 360
			{{0,1,0},{-1,0,0},{0,0,1}},  // 90 left
			{{-1,0,0},{0,-1,0},{0,0,1}}, // 180
			{{0,-1,0},{1,0,0},{0,0,1}},  // 90 right
			/* reflections */
			{{-1,0,0},{0,1,0},{0,0,1}},  // X-axis
			{{0,-1,0},{-1,0,0},{0,0,1}}, // -XY diagonal
			{{0,1,0},{1,0,0},{0,0,1}},   // +XY diagonal
			{{1,0,0},{0,-1,0},{0,0,1}},  // Y-axis
		};
		return data;
	}

	bool 
	LatticeDescriptorSQR::areNeighbored( const IntPoint &first, 
								const IntPoint &second ) const
	{
		int sum = abs(second.getX()-first.getX());
		return	sum < 2 
				&& (sum+=abs(second.getY()-first.getY())) < 2 
				&& (sum+=abs(second.getZ()-first.getZ())) == 1;
	}
	
	double
	LatticeDescriptorSQR::
	getBaseScale( const double neighVecLength ) const
	{
		return neighVecLength;
	}
	
	bool 
	LatticeDescriptorSQR::
	isLatticeNode( const IntPoint & p ) const
	{
		  // only XY plane points are valid nodes
		return (p.getZ() == 0); 
	}

	bool
	LatticeDescriptorSQR::
	isPossibleRing( const size_t ringSize ) const
	{
		  // check if ring size is even and at least 4
		return (ringSize > 3) && (ringSize%2==0);
	}

}	// namespace biu
