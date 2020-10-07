/*
 *  Main authors:
 *     Martin Mann <mmann@informatik.uni-freiburg.de>
 *
 *  Contributing authors:
 *     Sebastian Will <will@informatik.uni-freiburg.de>
 *
 *  Copyright:
 *     Martin Mann, 2007
 *
 *  This file is part of the CPSP-tools package:
 *     http://www.bioinf.uni-freiburg.de/sw/cpsp/
 *
 *  See the file "LICENSE" for information on usage and
 *  redistribution of this file, and for a
 *     DISCLAIMER OF ALL WARRANTIES.
 *
 */

#ifndef HPCONVERT_HH_
#define HPCONVERT_HH_

#include <biu/OptionParser.hh>
#include <biu/LatticeModel.hh>


	/**
	 * HPconvert allows the conversion of different lattice structure
	 * descriptions into each other.
	 * Supported formats are
	 *  - absolute moves
	 *  - relative moves
	 *  - absolute points in XYZ-file format
	 *  - absolute points in CML-file format 
	 */
	 class HPconvert {
	 protected:
	 
			//! types of conversion modes supported 
		enum Mode { A2R, A2P, R2A, R2P, P2A, P2R, A2C, A2N, R2N, R2C, A2X, NONE };
	 
			//! Inits the allowed parameters for HPconvert.
		void initAllowedArguments(biu::OptionMap & allowedArgs,
											std::string &infoText ) const;
											
			//! Converts an absolute move string into a relative one and prints
			//! it to stream.									
		void abs2rel(	const std::string & absMoves, 
						const biu::LatticeModel * lattice,
						std::ostream & out = std::cout) const;
						
			//! Converts an absolute move string into absolute positions and 
			//! prints it to stream.									
		void abs2xyz(	const std::string & absMoves, 
						const biu::LatticeModel * lattice,
						std::ostream & out = std::cout) const;
						
			//! Converts an absolute move string into absolute positions and 
			//! prints it to stream in CML-file format.									
		void abs2cml(	const std::string & absMoves, 
						const biu::LatticeModel * lattice,
						std::ostream & out = std::cout) const;
						
			//! Converts an absolute move string into absolute positions and 
			//! prints it to stream in PDB-file format.									
		void abs2pdb(	const std::string & absMoves, 
						const biu::LatticeModel * lattice,
						std::ostream & out = std::cout,
						const std::string* const seq = NULL,
						const double cAlphaDist = 3.8) const;
						
			//! Converts a relative move string into an absolute one and prints
			//! it to stream.									
		void rel2abs(	const std::string & absMoves, 
						const biu::LatticeModel * lattice,
						std::ostream & out = std::cout) const;
		
			//! Converts a relative move string into absolute positions and 
			//! prints it to stream.									
		void rel2xyz(	const std::string & absMoves, 
						const biu::LatticeModel * lattice,
						std::ostream & out = std::cout) const;
						
			//! Converts an relative move string into absolute positions and 
			//! prints it to stream in CML-file format.									
		void rel2cml(	const std::string & relMoves, 
						const biu::LatticeModel * lattice,
						std::ostream & out = std::cout) const;
						
			//! Converts absolute positions given by an XYZ-file to an 
			//! absolute move string and prints it to stream.									
		void xyz2abs(	const std::string & xyzFile, 
						const biu::LatticeModel * lattice,
						std::ostream & out = std::cout) const;
						
			//! Converts absolute positions given by an XYZ-file to a 
			//! relative move string and prints it to stream.									
		void xyz2rel(	const std::string & xyzFile, 
						const biu::LatticeModel * lattice,
						std::ostream & out = std::cout) const;
						
			//! Writes a point vector to stream in the XYZ-file format.
		void writeXYZ(	const biu::IPointVec & points,
						const biu::LatticeModel * lattice,
						std::ostream & out = std::cout) const;
						
			//! Writes a point vector to stream in the CML-file format.
		void writeCML(	const biu::IPointVec & points,
						const biu::LatticeModel * lattice,
						std::ostream & out = std::cout) const;
						
			//! Normalizes a move string and prints it to stream.									
		void mov2norm(	const std::string & moves, 
						const biu::LatticeModel * lattice,
						std::ostream & out = std::cout) const;
						
			//! Reads absolute points from stream in XYZ-file format.
			//! @param in the stream to read
			//! @param points the point vector to write the positions in
			//! @return the execution status of reading
		int readXYZ( std::istream & in, biu::IPointVec & points) const;
				
	 public:
	 		//! construction
	 	HPconvert();
	 	
	 		//! destruction
	 	~HPconvert();
	 	
	 		//! does conversion according to the given program parameters
	 		//! @return the execution status of parameter parsing and conversion
	 	int convert(int argc, char** argv);
	 	
	 }; // HPconvert

#endif /*HPCONVERT_HH_*/
