#include "nucleotides.h"
#include <map>
#include <set>
#include <string>
#include <vector>
#include "../Tools/zaprandom.h"
using namespace std;

std::map<int,char> nuclToAA;
std::map<char, vector<string> > AAtoNucl;
bool AAtoNuclMapsLoaded = false;

int singleCode(char c){
    if(c == 'A') return 0;
    if(c == 'T') return 1;
    if(c == 'G') return 2;
    if(c == 'C') return 3;
    cerr << "Bad nucleotide " << c << endl;
    return -1000;
}

int code(string S, unsigned int startPos = 0){
    if(S.size() < startPos + 3) {cerr << "code(" << S << "," << startPos << " out of bounds" << endl; return -1000;}
    return 16*singleCode(S[startPos]) + 4*singleCode(S[startPos+1]) + singleCode(S[startPos+2]);
}

int code(char c1, char c2, char c3){
    return 16*singleCode(c1) + 4*singleCode(c2) + singleCode(c3);
}

void add(string S, char C){
    int Ncode = code(S);
    if(nuclToAA.find(Ncode) == nuclToAA.end()){ // ie not found
        nuclToAA.insert(std::pair<int,char>(Ncode, C));
    } else {
        cout << "WRN: loading multiple times the nucl to AA definition for " << S << " new=" << C << " old = " << nuclToAA.find(Ncode)->second << " code=" << Ncode << endl ;
    }
    if(AAtoNucl.find(C) == AAtoNucl.end()){ // ie not found
        AAtoNucl.insert(std::pair<char, vector<string> >(C, {S}));
    } else {
        AAtoNucl[C].push_back(S);
    }
}

void loadAAtoNuclMaps(){
    cout << "   ... loading dictionnary to convert nucleodides into AAs" << endl;
    add("TTT" ,'F'); /*  Phe */   add("TCT" ,'S');/*  Ser */  add("TAT" ,'Y');/*  Tyr */  add("TGT" ,'C');/*  Cys */
    add("TTC" ,'F'); /*  Phe */   add("TCC" ,'S');/*  Ser */  add("TAC" ,'Y');/*  Tyr */  add("TGC" ,'C');/*  Cys */
    add("TTA" ,'L'); /*  Leu */   add("TCA" ,'S');/*  Ser */  add("TAA" ,'!');/*  Stop*/  add("TGA" ,'!');/* Stop */
    add("TTG" ,'L'); /*  Leu */   add("TCG" ,'S');/*  Ser */  add("TAG" ,'!');/*  Stop*/  add("TGG" ,'W');/*  Trp */

    add("CTT" ,'L'); /*  Leu */ 	add("CCT" ,'P');/*  Pro */  add("CAT" ,'H');/*  His */  add("CGT" ,'R');/*  Arg  */
    add("CTC" ,'L'); /*  Leu */ 	add("CCC" ,'P');/*  Pro */  add("CAC" ,'H');/*  His */  add("CGC" ,'R');/*  Arg  */
    add("CTA" ,'L'); /*  Leu */ 	add("CCA" ,'P');/*  Pro */  add("CAA" ,'Q');/*  Gln */  add("CGA" ,'R');/*  Arg  */
    add("CTG" ,'L'); /*  Leu */ 	add("CCG" ,'P');/*  Pro */  add("CAG" ,'Q');/*  Gln */  add("CGG" ,'R');/*  Arg  */

    add("ATT" ,'I'); /*  Ile */ 	add("ACT" ,'T');/*  Thr */  add("AAT" ,'N');/*  Asn */  add("AGT" ,'S');/*  Ser  */
    add("ATC" ,'I'); /*  Ile */ 	add("ACC" ,'T');/*  Thr */  add("AAC" ,'N');/*  Asn */  add("AGC" ,'S');/*  Ser  */
    add("ATA" ,'I'); /*  Ile */ 	add("ACA" ,'T');/*  Thr */  add("AAA" ,'K');/*  Lys */  add("AGA" ,'R');/*  Arg  */
    add("ATG" ,'M'); /*  Met */     add("ACG" ,'T');/*  Thr */  add("AAG" ,'K');/*  Lys */  add("AGG" ,'R');/*  Arg  */

    add("GTT" ,'V'); /*  Val */ 	add("GCT" ,'A');/*  Ala */  add("GAT" ,'D');/*  Asp */  add("GGT" ,'G');/*  Gly  */
    add("GTC" ,'V'); /*  Val */ 	add("GCC" ,'A');/*  Ala */  add("GAC" ,'D');/*  Asp */  add("GGC" ,'G');/*  Gly  */
    add("GTA" ,'V'); /*  Val */ 	add("GCA" ,'A');/*  Ala */  add("GAA" ,'E');/*  GlT */  add("GGA" ,'G');/*  Gly  */
    add("GTG" ,'V'); /*  Val */ 	add("GCG" ,'A');/*  Ala */  add("GAG" ,'E');/*  GlT */  add("GGG" ,'G');/*  Gly  */
    AAtoNuclMapsLoaded = true;
}

char getAA(string nucleotideSequence, int startPos){
    if(!AAtoNuclMapsLoaded) loadAAtoNuclMaps();
    int Ncode = code(nucleotideSequence, startPos);
    std::map<int,char>::iterator it = nuclToAA.find(Ncode);
    if(it != nuclToAA.end()) return it->second;
    cerr << "ERR: " << nucleotideSequence << ", from pos " << startPos << ", code unknown" << endl;
    return '?';
}

string randomDNA(int size){
    string allNucls = "ATGC";
    string res = string(size, '?');
    for(int i = 0; i < size; ++i){
        int newNuclId = random::uniformInteger(0,3);
        //int (double (rand()) * 20 / (double (RAND_MAX) + double (1)));
        if(newNuclId == 4) cerr << "ERR of iRandom (2)" << endl;
        res[i] = allNucls[newNuclId];
    }
    return res;
}


// this function is not 100% tested
string mutateDNA(string seq, bool onlyNonSynMuts, bool mutateAgainIfStop){
    int cpt = 0;
    size_t size = seq.size();
    bool hasStop = false;
    string initialAAseq = convertToProtein(seq, hasStop);
    if(hasStop && (mutateAgainIfStop)){
        cerr << "ERR: mutateDNA(" << seq << ") requested to generate without stop codons but already contains some... Please give stop free sequences" << endl;
    }
    while(cpt < 100){
        size_t position = static_cast<size_t>(random::uniformInteger(0, size - 1));
        if ((position >= size) || (position < 0)) {
           cerr << "Random does shit, mutateDNA";
           return string("");
        }
        char test = randomNucleotide();
        if(test != seq.at(position)){
            bool seqIsOk = true;
            seq.at(position) = test;
            bool containsStop = true;
            // silent mutation
            if(!initialAAseq.compare(convertToProtein(seq, containsStop))){
                if(onlyNonSynMuts) seqIsOk = false;
            }
            if(containsStop && (!mutateAgainIfStop)){
                if(containsStop) seqIsOk = false;
            }
            if(seqIsOk) return seq; //convertToProtein(seq);
        }
        cpt++;
    }

    cerr << "ERR: mutateDNA, got trouble to mutate an AA" << endl;
    return string(""); // the sequence was not changed (for instance all conserved)
}

char randomNucleotide(){
    string allNucls = "ATGC";
    int newNuclId = random::uniformInteger(0,3);
    if(newNuclId == 4) cerr << "ERR of iRandom (3)" << endl;
    return allNucls[newNuclId];
}

string clearStops(string nucl){
    string DNA = nucl;
    size_t L = DNA.size();
    if((L%3) != 0) {cerr << "ERR: convertToProtein " << DNA << ", size not a multiple of 3" << endl; return string("");}
    for(size_t i = 0; i < L/3; ++i){
        char C = getAA(DNA, 3*i);
        int cpt = 0;
        while((C == '!') && (cpt < 100)){
            for(int k = 3*i; k < 3*i+3; ++k){
                DNA.at(k) = randomNucleotide();
            }
            C = getAA(DNA, 3*i);
            cpt++;
        }
        if(cpt > 99) cerr << "ERR: clearStops failed to remove stops ! Either problem or very bad luck!" << endl;
    }
    return DNA;
}

string convertToProtein(string nucl){
    bool useless = false;
    return convertToProtein(nucl, useless); // a bit too much string copying, could make it better
}

string convertToProtein(string nucl, bool& containsStop){
    //cout << "Input " << nucl << endl;
    containsStop = false;
    int L = nucl.size();
    if((L%3) != 0) {cerr << "ERR: convertToProtein " << nucl << ", size not a multiple of 3" << endl; return string("");}
    string res(L/3, '?');
    for(int i = 0; i < L/3; ++i){
        res.at(i) = getAA(nucl, 3*i);
        if(res.at(i) == '!') containsStop = true;
        //cout << nucl[3*i] << nucl[3*i+1] << nucl[3*i+2] << "->" << getAA(nucl, 3*i) << endl;
    }
    //cout << "Output " << res << endl;
    return res;
}

void testNucleotides(){
    string s1 = "TAATATTAGGCCGAT"; // note: stop codon inside
    bool hasStop;
    cout << s1 << "\t" << convertToProtein(s1, hasStop) << endl;
    if(hasStop) cout << "Has Stop\n"; else  cout << "Without Stop\n";
    string s4 = clearStops(s1);
    cout << s4 << "\t" << convertToProtein(s4, hasStop) << endl;
    if(hasStop) cout << "Has Stop\n"; else  cout << "Without Stop\n";
    string s2 = "TTATATTAGGCCGA"; // wrong size
    cout << s2 << "\t" << convertToProtein(s2) << endl;
    string s3 = "TTATACTACGCCGAZ"; // wrong character
    cout << s3 << "\t" << convertToProtein(s3, hasStop) << endl;
    if(hasStop) cout << "Has Stop\n"; else  cout << "Without Stop\n";
}
