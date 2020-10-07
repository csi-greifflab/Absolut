#ifndef PDBSUPPORT_HH_
#define PDBSUPPORT_HH_

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <time.h>


#include <biu/LatticeModel.hh>




std::string getTimeString() {
	  // return string
	std::string strTime = "";
	  // get current time
	time_t tim = time(NULL);
	tm *curTime = localtime(&tim);
	  // generate time string
	strTime += char(48+curTime->tm_mday/10);
	strTime += char(48+curTime->tm_mday%10);
	strTime += '-';
	switch (curTime->tm_mon) {
		case 0  : strTime += "JAN"; break;
		case 1  : strTime += "FEB"; break;
		case 2  : strTime += "MAR"; break;
		case 3  : strTime += "APR"; break;
		case 4  : strTime += "MAY"; break;
		case 5  : strTime += "JUN"; break;
		case 6  : strTime += "JUL"; break;
		case 7  : strTime += "AUG"; break;
		case 8  : strTime += "SEP"; break;
		case 9  : strTime += "OCT"; break;
		case 10 : strTime += "NOV"; break;
		case 11 : strTime += "DEC"; break;
	}
	strTime += '-';
	int year10 = curTime->tm_year%10;
	int year100 = (curTime->tm_year-year10)%100;
	strTime += char(48+year100);
	strTime += char(48+year10);
	  // return time string
	return strTime;
}


void
writePDBseq( const std::vector<std::string> & seq,
			char chainID,
			std::ostream & out) 
{
	  // print sequence information
	size_t seqResLine = 1;
	out <<"SEQRES "<<std::setw(3) <<seqResLine++ <<" " <<chainID <<" "<<std::setw(4) <<seq.size()<<" ";
	for (size_t i=1; i<=seq.size(); i++) {
		out <<" " <<seq[i-1];
		if (i%13==0 && (i+1) < seq.size()) {
			out <<"\nSEQRES "<<std::setw(3) <<seqResLine++ <<" " <<chainID <<" "<<std::setw(4) <<seq.size()<<" ";
		}
	}
	out <<"\n";
}

void
writePDBhetatm( const std::vector<std::string> & seq,
			const biu::DPointVec & points, 
			char chainID,
			bool hasSideChain,
			std::ostream & out,
			size_t firstID) 
{
	size_t atomID = 0;
	for (size_t i=0; i<seq.size(); i++) {
		
		out <<"HETATM" <<std::setw(5)<<(atomID+1+firstID) <<"  "<<"CA " <<" "<<seq[i]
			<<" " <<chainID <<std::setw(4) <<(i+1) 
			<<" " <<"   " <<std::fixed 
			<<std::setprecision(3) <<std::setw(8) <<points[atomID].getX() 
			<<std::setprecision(3) <<std::setw(8) <<points[atomID].getY() 
			<<std::setprecision(3) <<std::setw(8) <<points[atomID].getZ()
			<<std::setprecision(2) <<std::setw(6) <<1.0
			<<std::setprecision(2) <<std::setw(6) <<0.0
			<<"           C  "
			<<"\n";
		atomID++;
		if (hasSideChain) {
			out <<"HETATM" <<std::setw(5)<<(atomID+1+firstID) <<"  "<<"CB " <<" "<<seq[i]
				<<" " <<chainID <<std::setw(4) <<(i+1) 
				<<" " <<"   " <<std::fixed 
				<<std::setw(8) <<std::setprecision(3) <<points[atomID].getX() 
				<<std::setw(8) <<std::setprecision(3) <<points[atomID].getY() 
				<<std::setw(8) <<std::setprecision(3) <<points[atomID].getZ()
				<<std::setw(6) <<std::setprecision(2) <<1.0
				<<std::setw(6) <<std::setprecision(2) <<0.0
				<<"          Cl  "
				<<"\n";
			atomID++;
		}
	}
}

void
writePDBconect( const std::vector<std::string> & seq,
			const biu::DPointVec & points, 
			char chainID,
			bool hasSideChain,
			std::ostream & out,
			size_t firstID) 
{
//                  1         2         3         4         5         6         7         8         
//         12345678901234567890123456789012345678901234567890123456789012345678901234567890
	 // CA -> CB
	 
	if (hasSideChain) {
		  // give connections explicitly
		size_t i = firstID;
		out	<<"CONECT" <<std::setw(5) <<(2*i)+1 
			<<std::setw(5) <<((2*i)+1)+1 <<std::setw(5) <<(2*(i+1))+1
			<<"                                                           \n";
		out	<<"CONECT" <<std::setw(5) <<((2*i)+1)+1 
			<<std::setw(5) <<(2*i)+1
			<<"                                                                \n";
		for (i++; (i+1)< (firstID+seq.size()); i++) {
			out	<<"CONECT" <<std::setw(5) <<(2*i)+1 
				<<std::setw(5) <<(2*(i-1))+1 <<std::setw(5) <<((2*i)+1)+1 <<std::setw(5) <<(2*(i+1))+1
				<<"                                                      \n";
			out	<<"CONECT" <<std::setw(5) <<((2*i)+1)+1 
				<<std::setw(5) <<(2*i)+1
				<<"                                                                \n";
		}
		out	<<"CONECT" <<std::setw(5) <<(2*i)+1 
			<<std::setw(5) <<(2*(i-1))+1 <<std::setw(5) <<((2*i)+1)+1
			<<"                                                           \n";
		out	<<"CONECT" <<std::setw(5) <<((2*i)+1)+1 
			<<std::setw(5) <<(2*i)+1
			<<"                                                                \n";
	} else {
		  // give connections explicitly
		size_t i = firstID;
		out	<<"CONECT" <<std::setw(5) <<i+1 
			<<std::setw(5) <<(i+1)+1
			<<"                                                                \n";
		for (i++; (i+1)< (firstID+seq.size()); i++) {
			out	<<"CONECT" <<std::setw(5) <<i+1 
				<<std::setw(5) <<i <<std::setw(5) <<(i+1)+1
				<<"                                                           \n";
		}
		out	<<"CONECT" <<std::setw(5) <<i+1 
			<<std::setw(5) <<i
			<<"                                                                \n";
	}	
}
void 
writePDB(	const std::string & sourceInfo,
			const std::vector<std::string> & seq,
			const biu::DPointVec & points, 
			const biu::LatticeModel & lattice,
			std::ostream & out,
			bool hasSideChain,
			const biu::LatticeModel & sclattice)  
{
	const char chainID = 'L';
	  // get lattice name
	std::string latName = lattice.getDescriptor()->getName();
	latName.resize(58,' ');

//                  1         2         3         4         5         6         7         8         
//         12345678901234567890123456789012345678901234567890123456789012345678901234567890
	out	<<"HEADER    LATTICE PROTEIN STRUCTURE               "<<getTimeString()<<"                     \n"              
		<<"TITLE     HP LATTICE PROTEIN STRUCTURE                                          \n" 
		<<"COMPND    MOL_ID: 1;                                                            \n"
		<<"COMPND   2 MOLECULE: LATTICE PROTEIN;                                           \n"        
		<<"COMPND   3 CHAIN: "<<chainID<<";                                                            \n"
		<<"COMPND   4 ENGINEERED: YES                                                      \n"
		<<"SOURCE    MOL_ID: 1                                                             \n"
		<<"KEYWDS    LATTICE PROTEIN MODEL                                                 \n"
		<<"EXPDTA    THEORETICAL MODEL                                                     \n"
		<<"AUTHOR    CPSP-TOOLS SOFTWARE                                                   \n"               
		<<"REMARK  40                                                                      \n"
		<<"REMARK  40 GENERATED WITH LATPACK (C) MARTIN MANN 2008                          \n";
	out	<<"REMARK  40  LATTICE = "<<latName <<"\n";
	if (hasSideChain) {		
		std::string latName = sclattice.getDescriptor()->getName();
		latName.resize(47,' ');
	out	<<"REMARK  40  SIDE CHAIN MODEL                                                    \n"
		<<"REMARK  40  SIDE CHAIN LATTICE = " <<latName <<"\n";
	}
	out	<<"REMARK 220                                                                      \n"
		<<"REMARK 220 EXPERIMENTAL DETAILS                                                 \n"
		<<"REMARK 220 EXPERIMENT TYPE : THEORETICAL MODELLING                              \n"
		<<"REMARK 220                                                                      \n"
		<<"REMARK 220 PARAMETER SET APPLIED TO LATFIT :                                    \n"
		<<"REMARK 220                                                                      \n";
	  // write parameter set out of infoString
	size_t cut = 0;
	while(cut < sourceInfo.size() && sourceInfo[cut]==' ')
	{ cut++; }
	while(cut < sourceInfo.size() && cut != std::string::npos) {
		size_t start = cut;
		cut = sourceInfo.find(" ", start+1);
		std::string output = "REMARK 220  ";
		output += sourceInfo.substr(start,cut-start);
		output.resize(80,' ');
		out <<output <<"\n";
		if (cut!=std::string::npos)
			cut++;
	}
		
	out	<<"REMARK 220                                                                      \n"
		<<"REMARK 225                                                                      \n"
		<<"REMARK 225 THEORETICAL MODEL                                                    \n"
		<<"REMARK 225 THE COORDINATES IN THIS ENTRY REPRESENT A LATTICE MODEL STRUCTURE.   \n"
		<<"REMARK 225                                                                      \n"
		<<"REMARK 225 THE LATTICE PROTEIN STRUCTURE IS GIVEN AS CHAIN "<<chainID<<".                   \n"
		<<"REMARK 225 THE BACKBONE IS REPRESENTED BY A CA-ATOMS AS A 'C'-elements.         \n"
		<<(hasSideChain?"REMARK 225 THE SIDE CHAIN IS REPRESENTED BY A CB-ATOMS AS A 'Cl'-elements.      \n":"")
		<<"REMARK 225                                                                      \n"
		;

	  // print sequence information
	writePDBseq( seq, chainID, out);
	
	  // print coordinates
	writePDBhetatm( seq, points, chainID, hasSideChain, out, 0);
	  // mark end of chain
	out <<"TER                                                                             \n";
	
	  // print contacts
	writePDBconect( seq, points, chainID, hasSideChain, out, 0);
	
	  // mark EOF
	out <<"END   "
		<<std::endl;

	out.flush();	
}


#endif /*PDBSUPPORT_HH_*/
