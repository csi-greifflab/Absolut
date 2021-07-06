#include "generatemutants.h"




/* functions from quality.cpp
 * // mask=1 = mutable
string mutateNPositions(string content, int nMut, string mask) {
    size_t size = content.size();
    if(mask.size() != size) {
        cerr << "ERR:mutateOnePosition, wrong mask size " << endl;
        return string("Wrong Mask");
    }
   int cpt = 0;
   int success = 0;
   while(cpt < 25){
       char test = randomAA();
       size_t position = random::uniformInteger(0,size-1);
       if ((position >= size) || (position < 0)) {
          cerr << "Random does shit, sequencespace.cpp::mutateOnePosition";
          return content; // fail
       }
       if((test != content[position]) && (mask[position] == '1')){
            content[position] = test;
            success++;
            mask[position] = '0'; // do not mutate again there
            if(success == nMut){
                return content; // success
            }
       }
       cpt++;
   }
   cerr << "ERR: foldedFree::mutateOnePosition, got trouble to mutate an AA" << endl;
   return string("MutateFailed"); // silent
}


// starts from an AA sequence, test the quality, then looks as N-points mutants and looks at the quality as well,
// Mask = 0100011 with 1 are positions that can be mutated.
// note: the random sequence is NOT following the mask, it's completely independent of the original one
// use NbMuts < 0 tto completely randomize each time
void evaluateAntigen(string antigenID, int nrAGs, int nBCRseqPerAG, bool startFromRandomAAseq, int NbMuts, string maskMut, string initialSeq){
    cout << "Now real mask" << endl;

    precise_stopwatch stopwatch; // careful, doesnt work for more than 2 hours

    std::pair<superProtein*, vector<int> > AG = getAntigen(antigenID);

    string AAseq = AG.first->getAAseq();

    // if gave initial AG sequence, takes it
    if(initialSeq.size() > 0) AAseq = initialSeq;
    size_t nAAs = AAseq.size();

    // checks if mask has good size. If not put  everywhere
    if(maskMut.size() != nAAs){
        cerr << "ERR: evaluateAntigen(" << antigenID << ", mask= " << maskMut << " is not same size as AA sequence from antigen " << AAseq << endl;
        maskMut = string(nAAs, '1');
    }

    // if want random sequence, randomizes only mutable positions
    if(startFromRandomAAseq){
        string newProt = randomProt(nAAs);
        for(size_t i = 0; i < newProt.size(); ++i){
            if(maskMut[i] == '1'){
                AAseq[i] = newProt[i];
            }
        }
    }

    int sizeReceptors = 10;
    int minInteract = 11;

    // store the original AA sequence
    string originalSeq = AAseq;


    for(int KL = 0; KL < nrAGs; ++KL){
        // the first time, use the original sequence
        if(KL == 0){
            AAseq = originalSeq;
        } else {
            // each time want a new sequence, will do mutations around
            if(NbMuts > 0)
                AAseq = mutateNPositions(originalSeq, NbMuts, maskMut);
            else {
                // randomizes all mutable positions
                string newProt = randomProt(nAAs);
                for(size_t i = 0; i < newProt.size(); ++i){
                    if(maskMut[i] == '1'){
                        AAseq[i] = newProt[i];
                    }
                }
            }
        }


        AG.first->setAAs(AAseq);
        affinityOneLigand T1 = affinityOneLigand(AG.first, sizeReceptors, minInteract, -1, 1, AG.second);

        unsigned int actual_wait_time = stopwatch.elapsed_time<unsigned int, std::chrono::microseconds>();
        cerr << actual_wait_time / 1000. << "ms Elapsed to compute structures" << endl;

        cout << "Computing 1000 affinities for ligand " << antigenID << " with sequence " << AAseq << " as reference, receptorSize " << sizeReceptors << " minI=" << minInteract << endl;

        vector<double> listAff;
        double vmin = 0;
        double avg = 0;
        int cpt = 0;
        for(int i = 0; i < nBCRseqPerAG; ++i){
            string Px = randomProt(sizeReceptors+1);
            std::pair<double, double> res = T1.affinity(Px);
            vmin = min(vmin, res.first);
            avg += res.first;
            listAff.push_back(res.first);
            //cout << "affinity(" << Px << ") = \t" << res.first << "\t" << res.second << endl;
        }
        cout << "AG:" << AAseq << ", avg binding " << avg / (double) cpt << endl;

        unsigned int actual_wait_time2 = stopwatch.elapsed_time<unsigned int, std::chrono::microseconds>();
        cerr << actual_wait_time2 / 1000.<< " ms Elapsed to compute 1000 affinities" << endl;

//        vector<double> intBoundaries;
//        for(double d = -120; d < -10; ++d){
//            intBoundaries.push_back(d-0.5);
//        }
//        cout << "   -> Distributions of best energies" << endl;
//        histogramFromDistrib h1(listAff, intBoundaries);
//        cout << h1.print(true) << endl;
    }

    //    cout << "Best Aff(type Best) so far:" << vmin << endl;
    //    for(int i = 0; i < 1; ++i){
    //        string Px = randomProt(sizeReceptors+1);
//        cout << "Details of the structures and affinities" << endl; // for " << simpleAccessible << " (Ag=" << AAsimple << ", BCR=" << Px << "), receptors " << receptorSize << " minI=4" << endl;
//        std::pair<double, double> res = T1->affinity(Px, true);
//        cout << "affinity(" << Px << ") = \t" << res.first << "\t" << res.second << endl;
//    }

}
*/
generateMutants::generateMutants()
{

}
