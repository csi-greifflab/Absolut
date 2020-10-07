// $Id: LatticeDescriptorFCC.cc,v 1.2 2016/08/08 12:41:59 mmann Exp $

#include <biu/LatticeDescriptorFCC.hh>
#include <cstdlib>

namespace biu
{
	
	LatticeDescriptorFCC::LatticeDescriptorFCC() 
		: LatticeDescriptor("fcc")
	{
		
		// gitterbasis
		latBase.push_back(IntPoint(1,0,0));	// x
		latBase.push_back(IntPoint(0,1,0));	// y
		latBase.push_back(IntPoint(0,0,1));	// z
		
		initNeighborhood();
		initAutomorphisms();
	}
	
	
	LatticeDescriptorFCC::~LatticeDescriptorFCC() {
	}

	unsigned int 
	LatticeDescriptorFCC::getNeighborDataSize() const {
		return 12;
	}
	
	const LatticeDescriptor::NeighborData *
	LatticeDescriptorFCC::getNeighborData() const {
			// init
		static const NeighborData fccNeighborData[] = {
			{"FL",{1,1,0},
			 {{1,0,0},{0,1,0},{0,0,1}},   {{1,0,0},{0,1,0},{0,0,1}}},
			{"LU",{0,1,1},
			 {{0,0,-1},{0,1,0},{1,0,0}},  {{0,0,1},{0,1,0},{-1,0,0}}},
			{"FU",{1,0,1},
			 {{1,0,0},{0,0,-1},{0,1,0}},  {{1,0,0},{0,0,1},{0,-1,0}}},
			{"BL",{-1,1,0},
			 {{0,-1,0},{1,0,0},{0,0,1}},  {{0,1,0},{-1,0,0},{0,0,1}}},
			{"RU",{0,-1,1},
			 {{0,0,-1},{-1,0,0},{0,1,0}}, {{0,-1,0},{0,0,1},{-1,0,0}}},
			{"BU",{-1,0,1},
			 {{0,-1,0},{0,0,-1},{1,0,0}}, {{0,0,1},{-1,0,0},{0,-1,0}}},
			{"FR",{1,-1,0},
			 {{0,1,0},{-1,0,0},{0,0,1}},  {{0,-1,0},{1,0,0},{0,0,1}}},
			{"LD",{0,1,-1},
			 {{0,0,1},{0,1,0},{-1,0,0}},  {{0,0,-1},{0,1,0},{1,0,0}}},
			{"FD",{1,0,-1},
			 {{1,0,0},{0,0,1},{0,-1,0}},  {{1,0,0},{0,0,-1},{0,1,0}}},
			{"BR",{-1,-1,0},
			 {{-1,0,0},{0,-1,0},{0,0,1}}, {{-1,0,0},{0,-1,0},{0,0,1}}},
			{"RD",{0,-1,-1},
			 {{0,0,1},{-1,0,0},{0,-1,0}}, {{0,-1,0},{0,0,-1},{1,0,0}}},
			{"BD",{-1,0,-1},
			 {{0,-1,0},{0,0,1},{-1,0,0}}, {{0,0,-1},{-1,0,0},{0,1,0}}}
			};
		
		return fccNeighborData;
	}

	
	unsigned int LatticeDescriptorFCC::getAutomorphismDataSize() const {
		return 48;
	}

	const LatticeDescriptorFCC::AutomorphismData *
	LatticeDescriptorFCC::getAutomorphismData() const {
		static const AutomorphismData data[] = 
			{{{1,0,0},{0,1,0},{0,0,1}}, {{0,0,1},{0,1,0},{1,0,0}},
			 {{0,0,1},{0,1,0},{-1,0,0}}, {{0,0,1},{0,-1,0},{1,0,0}},
			 {{0,0,1},{0,-1,0},{-1,0,0}}, {{0,0,1},{1,0,0},{0,1,0}},
			 {{0,0,1},{1,0,0},{0,-1,0}}, {{0,0,1},{-1,0,0},{0,1,0}},
			 {{0,0,1},{-1,0,0},{0,-1,0}}, {{0,0,-1},{0,1,0},{1,0,0}},
			 {{0,0,-1},{0,1,0},{-1,0,0}}, {{0,0,-1},{0,-1,0},{1,0,0}},
			 {{0,0,-1},{0,-1,0},{-1,0,0}}, {{0,0,-1},{1,0,0},{0,1,0}},
			 {{0,0,-1},{1,0,0},{0,-1,0}}, {{0,0,-1},{-1,0,0},{0,1,0}},
			 {{0,0,-1},{-1,0,0},{0,-1,0}}, {{0,1,0},{0,0,1},{1,0,0}},
			 {{0,1,0},{0,0,1},{-1,0,0}}, {{0,1,0},{0,0,-1},{1,0,0}},
			 {{0,1,0},{0,0,-1},{-1,0,0}}, {{0,1,0},{1,0,0},{0,0,1}},
			 {{0,1,0},{1,0,0},{0,0,-1}}, {{0,1,0},{-1,0,0},{0,0,1}},
			 {{0,1,0},{-1,0,0},{0,0,-1}}, {{0,-1,0},{0,0,1},{1,0,0}},
			 {{0,-1,0},{0,0,1},{-1,0,0}}, {{0,-1,0},{0,0,-1},{1,0,0}},
			 {{0,-1,0},{0,0,-1},{-1,0,0}}, {{0,-1,0},{1,0,0},{0,0,1}},
			 {{0,-1,0},{1,0,0},{0,0,-1}}, {{0,-1,0},{-1,0,0},{0,0,1}},
			 {{0,-1,0},{-1,0,0},{0,0,-1}}, {{1,0,0},{0,0,1},{0,1,0}},
			 {{1,0,0},{0,0,1},{0,-1,0}}, {{1,0,0},{0,0,-1},{0,1,0}},
			 {{1,0,0},{0,0,-1},{0,-1,0}}, {{1,0,0},{0,1,0},{0,0,-1}},
			 {{1,0,0},{0,-1,0},{0,0,1}}, {{1,0,0},{0,-1,0},{0,0,-1}},
			 {{-1,0,0},{0,0,1},{0,1,0}}, {{-1,0,0},{0,0,1},{0,-1,0}},
			 {{-1,0,0},{0,0,-1},{0,1,0}}, {{-1,0,0},{0,0,-1},{0,-1,0}},
			 {{-1,0,0},{0,1,0},{0,0,1}}, {{-1,0,0},{0,1,0},{0,0,-1}},
			 {{-1,0,0},{0,-1,0},{0,0,1}}, {{-1,0,0},{0,-1,0},{0,0,-1}}};
		return data;
	}
	
	bool 
	LatticeDescriptorFCC::areNeighbored( const IntPoint &first, 
								const IntPoint &second ) const
	{
		int dx = abs(second.getX()-first.getX()), dy=0, dz=0;
		return (!(	dx > 1 
					|| (dy = abs(second.getY()-first.getY())) > 1 
					|| (dz = abs(second.getZ()-first.getZ())) > 1 
					|| dx+dy+dz!=2));
	}
	
	double
	LatticeDescriptorFCC::
	getBaseScale( const double neighVecLength ) const
	{
		return  sqrt( pow(neighVecLength, 2) / 2.0 );
	}

	bool 
	LatticeDescriptorFCC::
	isLatticeNode( const IntPoint & p ) const
	{
		  // check if the coordinate sum is even
		return (abs(p.getX()+p.getY()+p.getZ()) % 2) == 0;
	}

	bool
	LatticeDescriptorFCC::
	isPossibleRing( const size_t ringSize ) const
	{
		  // check if ring size is at least 3
		return (ringSize > 2);
	}

	
} // namespace biu

