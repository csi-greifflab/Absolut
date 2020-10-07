// $Id: LatticeDescriptorCUB.cc,v 1.2 2016/08/08 12:42:00 mmann Exp $

#include <biu/LatticeDescriptorCUB.hh>
#include <cstdlib>
#include <biu/assertbiu.hh>

namespace biu
{
	
	LatticeDescriptorCUB::LatticeDescriptorCUB() 
		: LatticeDescriptor("cub")
	{
		
		// gitterbasis
		latBase.push_back(IntPoint(1,0,0));	// x
		latBase.push_back(IntPoint(0,1,0));	// y
		latBase.push_back(IntPoint(0,0,1));	// z
		
		initNeighborhood();
		initAutomorphisms();
		
			// for faster "isElement" call use special neighborhood class
		assertbiu(latNeighborhood != NULL, "no lattice neighborhood given");
		LatticeNeighborhood * tmp = new LatticeNeighborhoodCUB(*latNeighborhood);
		delete latNeighborhood;
		latNeighborhood = tmp;
	}

	LatticeDescriptorCUB::LatticeDescriptorCUB(const LatticeDescriptorCUB& toCopy)
		: LatticeDescriptor(toCopy)
	{
			// for faster "isElement" call use special neighborhood class
		assertbiu(latNeighborhood != NULL, "no lattice neighborhood given");
		LatticeNeighborhood * tmp = new LatticeNeighborhoodCUB(*latNeighborhood);
		delete latNeighborhood;
		latNeighborhood = tmp;
	}
	
	LatticeDescriptorCUB::~LatticeDescriptorCUB() {
		// NOTE: automatically calls destructor of ancestor LatticeDescriptor 
	}
	
	unsigned int 
	LatticeDescriptorCUB::getNeighborDataSize() const {
		return 6;
	}
	
	const LatticeDescriptor::NeighborData *
	LatticeDescriptorCUB::getNeighborData() const {
		static const NeighborData cubicNeighborData[] = {
			{"F",{ 1, 0, 0},
			 {{1,0,0},{0,1,0},{0,0,1}},  {{1,0,0},{0,1,0},{0,0,1}}},
			{"R",{ 0, 1, 0},
			 {{0,-1,0},{1,0,0},{0,0,1}}, {{0,1,0},{-1,0,0},{0,0,1}}},
			{"U",{ 0, 0, 1},
			 {{0,0,-1},{0,1,0},{1,0,0}}, {{0,0,1},{0,1,0},{-1,0,0}}},
			{"B",{-1, 0, 0},
			 {{-1,0,0},{0,-1,0},{0,0,1}},{{-1,0,0},{0,-1,0},{0,0,1}}},
			{"L",{ 0,-1, 0},
			 {{0,1,0},{-1,0,0},{0,0,1}}, {{0,-1,0},{1,0,0},{0,0,1}}},
			{"D",{ 0, 0,-1},
			 {{0,0,1},{0,1,0},{-1,0,0}}, {{0,0,-1},{0,1,0},{1,0,0}}}
		};
		return cubicNeighborData;
	}


	unsigned int LatticeDescriptorCUB::getAutomorphismDataSize() const {
		return 48;
	}

	const LatticeDescriptorCUB::AutomorphismData *
	LatticeDescriptorCUB::getAutomorphismData() const {
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
	LatticeDescriptorCUB::areNeighbored( const IntPoint &first, 
								const IntPoint &second ) const
	{
		int sum = abs(second.getX()-first.getX());
		return	sum < 2 
				&& (sum+=abs(second.getY()-first.getY())) < 2 
				&& (sum+=abs(second.getZ()-first.getZ())) == 1;
	}
	
	double
	LatticeDescriptorCUB::
	getBaseScale( const double neighVecLength ) const
	{
		return  neighVecLength;
	}

	bool 
	LatticeDescriptorCUB::
	isLatticeNode( const IntPoint & p ) const
	{
		  // all IntPoints are reachable
		return true;
	}
	
	
	bool
	LatticeDescriptorCUB::
	isPossibleRing( const size_t ringSize ) const
	{
		  // check if ring size is even and at least 4
		return (ringSize > 3) && (ringSize%2==0);
	}

}
