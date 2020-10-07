// $Id: OffLatticeProtein.cc,v 1.2 2016/08/08 12:41:56 mmann Exp $

#include <biu/OffLatticeProtein.hh>
#include <fstream>
#include <stdlib.h>


namespace biu
{
	OffLatticeProtein::OffLatticeProtein(const std::string& coordinatesFileName,
	 					const Alphabet* const _alphabet)
		:	pData(new DPointVec()), 
	 		alphabet(_alphabet), 
	 		sequence(NULL) 
	 {
		std::ifstream ifs;
		ifs.open(coordinatesFileName.c_str());
		if (!ifs) {
			std::cerr<< "ERROR OffLatticeProtein::OffLatticeProtein : "+
			 coordinatesFileName+" kann nicht ge�ffnet werden.";
			delete pData;
			exit(-1);
		}

		double x, y, z;
		char aa;
		std::string remark, seqStr = "";
		
		// ignore remarks (#...)
		while(ifs.peek() == '#') {
			getline(ifs, remark);
		}
		while(ifs >> aa >> x >> y >> z) {
			pData->push_back(DblPoint(x, y, z));
			seqStr += aa;
		}
		ifs.close();

		assertbiu(alphabet->isAlphabetString(seqStr),
		 "AA sequence contains characters which are not in the alphabet.");
		sequence = new Sequence(alphabet->getSequence(seqStr));
	}

	OffLatticeProtein::OffLatticeProtein(	const DPointVec& data3D, 
											const Alphabet* const _alphabet, 
											const std::string& seqStr) 
	 : pData(new DPointVec(data3D)), alphabet(_alphabet), sequence(NULL) 
	{
		assertbiu(alphabet != NULL,  "alphabet is not allowed to be NULL.");
		sequence = new Sequence(alphabet->getSequence(seqStr));
		assertbiu(pData->size() == sequence->size(),
			"Structure and sequence have to have the same length.");
	}

	OffLatticeProtein::OffLatticeProtein( const OffLatticeProtein& offLatPro)
	 :	pData(new DPointVec(*(offLatPro.pData))), 
	 	alphabet(offLatPro.alphabet),
		sequence(new Sequence(*(offLatPro.sequence))) 
	{
		assertbiu(alphabet != NULL,  "alphabet is not allowed to be NULL.");
		assertbiu(pData->size() == sequence->size(),
		 "Structure and sequence have to have the same length.");
	}
	
	OffLatticeProtein::~OffLatticeProtein()	{
		delete pData; pData = NULL;
		delete sequence; sequence = NULL;
	}

	
	OffLatticeProtein& 
	OffLatticeProtein::operator= ( const OffLatticeProtein& offLatPro2) 
	{
		if (this != &offLatPro2) {
			delete pData;
			pData = new DPointVec(*(offLatPro2.pData));
			delete sequence;
			sequence = new Alphabet::Sequence(*(offLatPro2.sequence));
			alphabet = offLatPro2.alphabet;
		}
		return *this;
	}

	void		
	OffLatticeProtein::writePDB(const std::string& pdbFileName) {
	//TODO
	}
	
	LatticeProtein*	
	OffLatticeProtein::approximateToLattice(const LatticeModel* const lattice,	
				const DistanceEnergyFunction* const energy) const 
	{
	//TODO
		return NULL;
	}

	DPointVec	
	OffLatticeProtein::get3Ddata() const{
		return *pData;	
	}
	
	double 		
	OffLatticeProtein::getDRMSD(const BackboneStructure3D& other) const{
// TODO implementieren
	/*	DPointVec other3Ddata = other.get3Ddata();
		assertbiu(pData->size() == other3Ddata->size(),
		 "Structure and sequence have to have the same length.");
        double rmsd = 0;
        for(int k=0; k<i; k++)
                for(int l=k+1; l<=i; l++)
                        rmsd += pow( pos[k].dist(pos[l]) - (scaling * conformation[k].dist(conformation[l])), 2);
        rmsd /= ((i+1)*i /2);   // Index i => Kettenlaenge i+1
        return sqrt(rmsd);*/
		return 0.0;

	}




}
