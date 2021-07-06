#include "quality.h"
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include "antigenLib.h"
#include "../Tools/stopwatch.h"
using namespace std;

void testQuality(){
    evaluateAntigen("L10paper2AGs(V1a)", 2, 200, true, 2, "101101111101011111100000000000100000000");
    evaluateAntigen("L10paper2AGs(V1a)", 2, 200, true, -1, "101101111101011111100000000000100000000");
}

// mask=1 = mutable
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

        /*        vector<double> intBoundaries;
        for(double d = -120; d < -10; ++d){
            intBoundaries.push_back(d-0.5);
        }
        cout << "   -> Distributions of best energies" << endl;
        histogramFromDistrib h1(listAff, intBoundaries);
        cout << h1.print(true) << endl; */
    }

    //    cout << "Best Aff(type Best) so far:" << vmin << endl;
    //    for(int i = 0; i < 1; ++i){
    //        string Px = randomProt(sizeReceptors+1);
//        cout << "Details of the structures and affinities" << endl; // for " << simpleAccessible << " (Ag=" << AAsimple << ", BCR=" << Px << "), receptors " << receptorSize << " minI=4" << endl;
//        std::pair<double, double> res = T1->affinity(Px, true);
//        cout << "affinity(" << Px << ") = \t" << res.first << "\t" << res.second << endl;
//    }

}

/*
vector<double> distribEnergies(
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

*/


//string binarySequence::testAffinityFunctions(double L, double R, int maxClusters, int typeAffinityFunction) {

//    stringstream subFolder;
//    subFolder << "L" << L << "_r" << R << nameTypeAff(typeAffinityFunction) << "_C" << maxClusters;
//    string folder = string("C:/Users/Philippe/Desktop/Sequences/") + subFolder.str() + string("/");
//    createFolder(folder);

//    //#define out cout
//    stringstream out;
//    out << "Testing the properties ot the affinity function for the following parameters : \n";
//    out << "   ->     L= " << L << "\t(Size of sequences)" << endl;
//    out << "   ->     R= " << R << "\t(specificity parameter)" << endl;
//    out << "   -> maxCl= " << maxClusters << "\t(cluster size scale)" << endl;
//    switch (typeAffinityFunction) {
//    case seqAff: {out << "   -> Using standard affinity (Saham's)\n"; break;}
//    case seqAffNorm: {out << "   -> Using standard affinity normalized by maxCl^r\n"; break;}
//    case seqAffWindow: {out << "   -> Using the maximum affinity of a sliding window\n"; break;}
//    }
//    out << "==== Part 1 : enumerates all (if possible), or a lot of sequences and sort them by affinity to get the best ones : ====" << endl;

//#define resolutiondistrib 100
//#define maxSequencesToEnumeate 500000

//    // The reference antigen is 00000...
//    binarySequence * ref = new binarySequence(L);

//    // will store a large list of random sequences, with their affinity to ref
//    vector<pair<double, binarySequence*> > store;


//    out << "1 ----------- Distribution of affinities ----------------------" << endl;
//    vector<double> distribution;
//    vector<double> logDistrib;
//    distribution.resize(resolutiondistrib + 1);

//    // in case L is big, will only sample maxSequencesToEnumeate sequences.
//    int total = 0;
//    int maxim = pow(2, L);
//    if (L > 26) {
//        maxim = maxSequencesToEnumeate + 1;           // to avoid it to become negative ...
//    }

//    // if possible to enumerate all, do it one by one.
//    bool enumerateAll = (maxim < maxSequencesToEnumeate);
//    if (enumerateAll) {
//        for (int i = 0; i < maxim; ++i) {
//            binarySequence * a = new binarySequence(L, (long) i);
//            double affi = binarySequence::affinity(a, ref, R, maxClusters, typeAffinityFunction);
//            distribution[(int) (((double) resolutiondistrib) * affi)] += 1.0; // put into the histogram
//            store.push_back(pair<double, binarySequence*> (affi,a));
//            logDistrib.push_back(log10(affi+1e-6));
//            total++;
//        }
//    } else {
//        // if not, than sample randomly
//        for (int i = 0; i < maxSequencesToEnumeate; ++i) {
//            if(((i%100000) == 0) || (i == maxSequencesToEnumeate-1)) cout << i << "/" << maxSequencesToEnumeate << endl;
//            binarySequence * a = new binarySequence(L);
//            a->randomize();
//            double affi = binarySequence::affinity(a, ref, R, maxClusters, typeAffinityFunction);
//            distribution[(int) (((double) resolutiondistrib) * affi)] += 1.0;
//            store.push_back(pair<double, binarySequence*> (affi,a));
//            logDistrib.push_back(log10(affi+1e-6));
//            total++;
//        }
//    }

//    out << "Distribution of affinities\n";
//    for (int i = 0; i < (int) distribution.size(); ++i) {
//        distribution[i] /= (double) total;
//        out << i << "\t"  << double (i) * (1.0 / (double) resolutiondistrib)
//            << "\t" << double (i + 1) * (1.0 / (double) resolutiondistrib)
//            << "\t" << distribution[i] <<  endl;
//    }
//    out << "Distribution of affinities in LOG scale\n";
//    histogramFromDistrib v(logDistrib, 100);
//    out << v.print() << endl;

//    out << "\nSequences and affinity, "
//        << ((enumerateAll) ? " in the order of ID\n" : " randomly generated\n");
//    for (int i = 0; i < 200; ++i) {
//        out << i << "\t" << store[i].second->print() << "\t" << store[i].first << "\n";
//    }

//    out << "\nSequences, sorted from the best, out of the "
//        << min(maxim,(int) maxSequencesToEnumeate) << " evaluated sequences\n";
//    std::sort(store.begin(), store.end(), compSequences);
//    for (int i = 0; i < 200; ++i) {
//        out << i << "\t" << store[i].second->print() << "\t" << store[i].first << "\n";
//    }

//    if (enumerateAll) {
//        out << "\nAffinity of sequences taken randomly\n";
//        for (int i = 0; i < 100; ++i) {
//            binarySequence * seqtmp = new binarySequence(L);
//            seqtmp->randomize();
//            out << i << "\t" << binarySequence::affinity(seqtmp,ref,R,maxClusters,typeAffinityFunction) << "\t" << seqtmp->print() << "\n";
//        }
//    }

//    out << "2 ------------------------- Mutation histograms -----------------------" << endl;

//    // Cuts the distribution of affinities in blocs of 5% of the sequences
//#define sizeClasses 0.05
//    out << "\nEqual classes of affinities (best percentiles)" << endl;
//    vector<double> classes;
//    vector<double> valuesClasses;

//    std::reverse(store.begin(), store.end()); // now to put it increasing

//    // every 5% of sequences browsed (fraction), store the value of affinity
//    double fraction = 0.00;
//    for(unsigned int i = 0; i < store.size(); ++i){
//        //out << store[i].second->print() << " Aff " << store[i].first << endl;
//        if((((double) i / (double) (store.size())) >= fraction - 1e-6) || (i == store.size() - 1)){
//            //out << "]" << fraction - sizeClasses << "," << fraction << "]\t" << store[i].first << endl;
//            valuesClasses.push_back(store[i].first);
//            classes.push_back(fraction);
//            fraction += sizeClasses;
//        }
//    }
//    out << "Affinity classes: cell at percent XX has the affinity value of: (last line=log, 0->-10)" << endl;
//    out << printVec(classes) << endl;
//    out << printVec(valuesClasses) << endl;
//    for(unsigned int i = 0; i < valuesClasses.size(); ++i){
//        out << "\t" << log10(valuesClasses[i] + 1e-10);
//    } out << endl;

//    out << " Now, for each class of sequences, makes the distribution of affinities." << endl;

//    vector<double> newAffinitiesByMutation;      // collecting new affinities only inside one class
//    vector<double> foldAffinitiesByMutation;
//    vector< vector<double>> tableAbsMutations;   // to store the histograms inside each class.
//    vector< vector<double>> tableFoldMutations;
//    vector<double> TOTALnewAffinitiesByMutation; // collecting new affinities for all sequences.
//    vector<double> TOTALfoldAffinitiesByMutation;
//    // For making histograms in term of fold induction, the following groups/classes will be used to make histograms
//    vector<double> classesFoldInd = {0,0.1,0.2,0.4,0.6,0.8,0.9,0.99,1.01,1.1,1.2,1.4,2,2.5,3.25,5,7.5,10};
//    // For making histograms in term of affinity, the classes inside valueclasses will be used.

//    classes.push_back(1.0); // to avoid seg fault

//    // Now will browse again, class by class. Values of sequences inside valuesClasses[currentClass-1] and valuesClasses[currentClass]
//    fraction = sizeClasses;
//    int currentClass = 1;

//    for(unsigned int i = 0; i < store.size(); ++i){
//        binarySequence* thisSeq = store[i].second;
//        double oldAff = binarySequence::affinity(thisSeq, ref, R, maxClusters, typeAffinityFunction);
//        if((oldAff > valuesClasses[currentClass] + 1e-6) || (oldAff < valuesClasses[currentClass-1] - 1e-6) ) cerr << "Class problems, sequence outside its class" << endl;

//        // within one class, add the new possible affinities
//        for(int j = 0; j < thisSeq->size * 2; ++j){ // might get sequence of different lengths !
//            binarySequence seq(thisSeq);
//            seq.mutateOnePosition();
//            double newAff = binarySequence::affinity(&seq, ref, R, maxClusters, typeAffinityFunction);
//            newAffinitiesByMutation.push_back(newAff);
//            if(oldAff > 0) foldAffinitiesByMutation.push_back(newAff / oldAff);
//            TOTALnewAffinitiesByMutation.push_back(newAff);
//            if(oldAff > 0) TOTALfoldAffinitiesByMutation.push_back(newAff / oldAff);
//            //out << thisSeq->print() << " Aff " << oldAff << "->" << seq.print() << " Newaff " << newAff << " Fold " << newAff / oldAff << endl;
//        }
//        if((((double) i / (double) (store.size())) >= fraction - 1e-6) || (i == store.size() - 1)){
//            //out << "Histograms for mutation for the class [" << valuesClasses[currentClass-1] << "," << valuesClasses[currentClass] << "]" << endl;
//            histogramFromDistrib resAbsForThisClass(newAffinitiesByMutation,valuesClasses);
//            tableAbsMutations.push_back(resAbsForThisClass.densities);
//            //out << "New affinities" << endl;
//            //out << resAbsForThisClass.print(true);
//            histogramFromDistrib resFoldForThisClass(foldAffinitiesByMutation,classesFoldInd);
//            tableFoldMutations.push_back(resFoldForThisClass.densities);
//            //out << "Fold increase in the affinity" << endl;
//            //out << resFoldForThisClass.print(true);
//            //out << endl << endl << endl;
//            newAffinitiesByMutation.clear();
//            foldAffinitiesByMutation.clear();
//            currentClass++;
//            fraction += sizeClasses;
//        }

//    }
//    out << "Outputing the affinity changes from all sequences" << endl;
//    histogramFromDistrib TOTALresAbsForThisClass(TOTALnewAffinitiesByMutation,valuesClasses);
//    out << "New affinities" << endl;
//    out << TOTALresAbsForThisClass.print(true);
//    histogramFromDistrib TOTALresFoldForThisClass(TOTALfoldAffinitiesByMutation,classesFoldInd);
//    out << "Fold increase in the affinity" << endl;
//    out << TOTALresFoldForThisClass.print(true);

//    out << "Outputing the affinity changes (new affinity value), from each class of affinity" << endl;
//    out << "\t";
//    for(unsigned int j = 0; j < tableAbsMutations.size(); ++j){
//        out << "[" << valuesClasses[j] << "," << valuesClasses[j+1] << "]\t";
//    }
//    out << "\n";
//    for(unsigned int i = 0; i < TOTALresAbsForThisClass.densities.size(); ++i){
//        out << "[" << TOTALresAbsForThisClass.lowBoundsXs[i] << "," << TOTALresAbsForThisClass.highBoundsXs[i] << "]";
//        for(unsigned int j = 0; j < tableAbsMutations.size(); ++j){
//            out << "\t" << tableAbsMutations[j][i];
//        }
//        out << endl;
//    }

//    out << "Outputing the affinity changes (new affinity value), from each class of affinity" << endl;
//    out << "\t";
//    for(unsigned int j = 0; j < tableFoldMutations.size(); ++j){
//        out << "[" << valuesClasses[j] << "," << valuesClasses[j+1] << "]\t";
//    }
//    out << "\n";
//    for(unsigned int i = 0; i < TOTALresFoldForThisClass.densities.size(); ++i){
//        out << "[" << TOTALresFoldForThisClass.lowBoundsXs[i] << "," << TOTALresFoldForThisClass.highBoundsXs[i] << "]";
//        for(unsigned int j = 0; j < tableFoldMutations.size(); ++j){
//            out << "\t" << tableFoldMutations[j][i];
//        }
//        out << endl;
//    }




//    out << "==== Part 2 : Evaluating cross-reactivity in the system : ====" << endl;

//    for(int k = 0; k < 10; ++k){

//        binarySequence ref2(ref);
//        ref2.mutateOnePosition();
//        ref2.mutateOnePosition();

//        stringstream fname;
//        fname << folder << "DotPlot2CloseAG-L" << L << "_r" << R << nameTypeAff(typeAffinityFunction) << "_C" << maxClusters << "rep" << k << "_" << ref->print() << "to" << ref2.print() << ".txt";
//        fstream of(fname.str().c_str(), ios::out);

//        vector<double> affRef;
//        vector<double> affRef2;

//        for(int i = 0; i < 20000; ++i){
//            binarySequence test(ref2);
//            test.randomize();
//            double aff1 = binarySequence::affinity(&test, ref, R, maxClusters, typeAffinityFunction);
//            double aff2 = binarySequence::affinity(&test, &ref2, R, maxClusters, typeAffinityFunction);
//            affRef.push_back(aff1);
//            affRef2.push_back(aff2);
//            of << aff1 << "\t" << aff2 << endl;
//        }
//        of.close();
//    }

//    for(int k = 0; k < 10; ++k){

//        binarySequence ref2(ref);
//        ref2.mutateOnePosition();
//        ref2.mutateOnePosition();
//        ref2.mutateOnePosition();
//        ref2.mutateOnePosition();
//        ref2.mutateOnePosition();
//        ref2.mutateOnePosition();
//        ref2.mutateOnePosition();

//        stringstream fname;
//        fname << folder << "DotPlot6CloseAG-L" << L << "_r" << R << nameTypeAff(typeAffinityFunction) << "_C" << maxClusters << "rep" << k << "_" << ref->print() << "to" << ref2.print() << ".txt";
//        fstream of(fname.str().c_str(), ios::out);

//        vector<double> affRef;
//        vector<double> affRef2;

//        for(int i = 0; i < 20000; ++i){
//            binarySequence test(ref2);
//            test.randomize();
//            double aff1 = binarySequence::affinity(&test, ref, R, maxClusters, typeAffinityFunction);
//            double aff2 = binarySequence::affinity(&test, &ref2, R, maxClusters, typeAffinityFunction);
//            affRef.push_back(aff1);
//            affRef2.push_back(aff2);
//            of << aff1 << "\t" << aff2 << endl;
//        }
//        of.close();
//    }

//    for(int k = 0; k < 10; ++k){
//        binarySequence ref2(ref);
//        ref2.randomize();

//        stringstream fname;
//        fname << folder << "DotPlotRandCloseAG-L" << L << "_r" << R << nameTypeAff(typeAffinityFunction) << "_C" << maxClusters << "rep" << k << "_" << ref->print() << "to" << ref2.print() << ".txt";
//        fstream of(fname.str().c_str(), ios::out);


//        vector<double> affRef;
//        vector<double> affRef2;

//        for(int i = 0; i < 20000; ++i){
//            binarySequence test(ref2);
//            test.randomize();
//            double aff1 = binarySequence::affinity(&test, ref, R, maxClusters, typeAffinityFunction);
//            double aff2 = binarySequence::affinity(&test, &ref2, R, maxClusters, typeAffinityFunction);
//            affRef.push_back(aff1);
//            affRef2.push_back(aff2);
//            of << aff1 << "\t" << aff2 << endl;
//        }
//        of.close();
//    }

//    ofstream ffin(folder + string("Output.txt"));
//    ffin << out.str();
//    ffin.close();
//    return out.str();


//    int nbAntigens = 10;
//    out << "Generating randomly " << nbAntigens << " antigens " << endl;
//    vector<binarySequence*> ags;
//    for (int i = 0; i < nbAntigens; ++i) {
//        binarySequence * seq = new binarySequence(L);
//        seq->randomize();
//        ags.push_back(seq);
//        out << "\tAg nr " << i << "\t" << seq->print() << endl;
//    }
//    out << "\nNumber of antigens recognized by randomly generated sequences, based on threshold\n";

//    out << "  -> (for the first 100 sequences : ) In the case of random sequences" << endl;
//    total = 0;
//#define thresholdRecoAg 0.1
//    int nbDiscardedSeq = 0;  // sequences that don't recognize anything
//    int countprint = 0;
//    for (int k = 0; k < min(maxim, (int) maxSequencesToEnumeate); ++k) {
//        if (k == 100) {
//            out
//                    <<
//                       "  -> (for the remaining sequences) for sequences recognizing at least an antigen with affinity 0.1"
//                    << endl;
//        }
//        total++;

//        // for each sequence,
//        bool recoAtLeastOne = false;
//        vector<double> nbRecoDepThresh(10, 0.0);
//        vector<double> affinityEach(nbAntigens, 0.0);
//        binarySequence * seqtmp = new binarySequence(L);
//        seqtmp->randomize();
//        for (int j = 0; j < nbAntigens; ++j) {
//            double thisAff = binarySequence::affinity(seqtmp,
//                                                      ags[j],
//                                                      R,
//                                                      maxClusters,
//                                                      typeAffinityFunction);
//            if ((thisAff > thresholdRecoAg) || (k < 100)) {
//                recoAtLeastOne = true;
//            } else { nbDiscardedSeq++; }
//            affinityEach[j] = thisAff;
//            for (int i = 0; i <= (int) (9.99 * thisAff); ++i) {
//                if (i < 10) { nbRecoDepThresh[i]++; }
//            }
//        }
//        if (recoAtLeastOne && (countprint < 5000)) {
//            countprint++;
//            out << "RandSeq " << k << ", " << seqtmp->print() << " ";
//            out << "nbAgPerThreshold:";
//            for (int i = 0; i < 10; ++i) {
//                out << "\t" << nbRecoDepThresh[i];
//            }
//            out << "\taffPerAg:";
//            for (int i = 0; i < nbAntigens; ++i) {
//                out << "\t" << affinityEach[i];
//            }
//            out << endl;
//        }
//        delete seqtmp;
//    }
//    out << "   ... Nb of sequences analyzed: " << total << endl;
//    out << "   ... Nb of sequences discarded: " << nbDiscardedSeq
//        << "(except the 100 first ones, i.e. among the :" << total - 100 << " remaining)" << endl;

//    out << "==== Part 3 : Evaluating the effect of mutations : ====" << endl;

//    binarySequence * start = new binarySequence(L);    // starting by '0000' : the best sequence
//    out << "NbMut\tsequence\taffinity\n";
//    for (int i = 0; i < 2 * L; ++i) {
//        out << i << "\t" << start->print() << "\t" << binarySequence::affinity(start,
//                                                                               ref,
//                                                                               R,
//                                                                               maxClusters,
//                                                                               typeAffinityFunction)
//            << endl;
//        start->mutateOnePosition();
//    }

//    out << "\tReaching a good affinity\t" << endl;
//    start->randomize();
//    double prevaff = binarySequence::affinity(start, ref, R, maxClusters, typeAffinityFunction);

//    bool stop = false;
//    for (int i = 0; (i < L) && (!stop); ++i) {
//        out << "sequence : " << start->print() << "\tAff:\t" << prevaff << "\t";
//        out << "PossibleMut:";
//        vector<int> posGoodMutations;
//        for (int i = 0; i < L; ++i) {
//            binarySequence stmp = binarySequence(start);
//            stmp.content[i] = !stmp.content[i];
//            double newaff
//                    = binarySequence::affinity(&stmp, ref, R, maxClusters, typeAffinityFunction);
//            out << "\t" << newaff;
//            if (newaff > prevaff) { posGoodMutations.push_back(i); }
//        }
//        out << endl;
//        if (posGoodMutations.size() > 0) {
//            int nextmut = random::uniformInteger(0,posGoodMutations.size() - 1);
//            start->content[posGoodMutations[nextmut]] = !start->content[posGoodMutations[nextmut]];
//            prevaff = binarySequence::affinity(start, ref, R, maxClusters, typeAffinityFunction);
//        } else {
//            stop = true;
//        }
//    }


//    for (int i = 0; i < (int) store.size(); ++i) {
//        delete store[i].second;
//    }


//    out << "\tReaching a good affinity\t" << endl;
//    start->randomize();
//    /*double prevaff = binarySequence::seq_affinity(start, ref, R, maxClusters, typeAffinityFunction);

//       bool stop = false;
//       for (int i = 0; (i < L) && (!stop); ++i) {
//          out << "sequence : " << start->print() << "\tAff:\t" << prevaff << "\t";
//          out << "PossibleMut:";
//          vector<int> posGoodMutations;
//          for (int i = 0; i < L; ++i) {
//             binarySequence stmp = binarySequence(start);
//             stmp.content[i] = !stmp.content[i];
//             double newaff
//                = binarySequence::seq_affinity(&stmp, ref, R, maxClusters, typeAffinityFunction);
//             out << "\t" << newaff;
//             if (newaff > prevaff) { posGoodMutations.push_back(i); }
//          }
//          out << endl;
//          if (posGoodMutations.size() > 0) {
//             int nextmut = irandom(posGoodMutations.size() - 1);
//             start->content[posGoodMutations[nextmut]] = !start->content[posGoodMutations[nextmut]];
//             prevaff = binarySequence::seq_affinity(start, ref, R, maxClusters, typeAffinityFunction);
//          } else {
//             stop = true;
//          }
//       }*/



//    out << "2D plot of affinities to 2 antigens" << endl;
//    out << "Case 1 : two close antigens, one mutation away" << endl;

//    return out.str();
//}







//#include "quality.h"
//#include <vector>
//#include <string>
//#include <iostream>
//#include <sstream>
//#include <fstream>
//#include <antigenLib.h>
//#include "../Tools/stopwatch.h"
//using namespace std;

//void testQuality(){
//    evaluateAntigen("L10paper2AGs(V1a)", true, "101101111101011111100000000000100000000", 2);
//}

//// starts from an AA sequence, test the quality, then looks as N-points mutants and looks at the quality as well,
//// Mask = 0100011 with 1 are positions that can be mutated.
//// note: the random sequence is NOT following the mask, it's completely independent of the original one
//void evaluateAntigen(string antigenID, bool startFromRandomAAseq, string maskMut, int NbMuts){
//    cout << "Now real mask" << endl;
//    for(int KL = 0; KL < 100; ++KL){
//        precise_stopwatch stopwatch; // careful, doesnt work for more than 2 hours

//        std::pair<superProtein*, vector<int> > AG = getAntigen(antigenID);

//        string AAseq = AG.first->getAAseq();
//        size_t nAAs = AAseq.size();
//        if(maskMut.size() != nAAs){
//            cerr << "ERR: evaluateAntigen(" << antigenID << ", mask= " << maskMut << " is not same size as AA sequence from antigen " << AAseq << endl;
//        }
//        if(startFromRandomAAseq){
//            string newProt = randomProt(nAAs);
//            for(size_t i = 0; i < newProt.size(); ++i){
//                if(maskMut[i] == '1'){
//                    AAseq[i] = newProt[i];
//                }
//            }
//            AG.first->setAAs(AAseq);
//        }

//        int sizeReceptors = 8;
//        int minInteract = 4;
//        affinityOneLigand T1 = affinityOneLigand(AG.first, sizeReceptors, minInteract, -1, 1, AG.second);

//        unsigned int actual_wait_time = stopwatch.elapsed_time<unsigned int, std::chrono::microseconds>();
//        cerr << actual_wait_time / 1000. << "ms Elapsed to compute structures" << endl;

//        cout << "Computing 1000 affinities for ligand " << antigenID << " with sequence " << AAseq << " as reference, receptorSize " << sizeReceptors << " minI=" << minInteract << endl;

//        vector<double> listAff;
//        double vmin = 0;
//        double avg = 0;
//        int cpt = 0;
//        for(int i = 0; i < 200; ++i){
//            string Px = randomProt(sizeReceptors+1);
//            std::pair<double, double> res = T1.affinity(Px);
//            vmin = min(vmin, res.first);
//            avg += res.first;
//            listAff.push_back(res.first);
//            cout << "affinity(" << Px << ") = \t" << res.first << "\t" << res.second << endl;
//        }
//        cout << "AG:" << AAseq << ", avg binding " << avg / (double) cpt << endl;

//        unsigned int actual_wait_time2 = stopwatch.elapsed_time<unsigned int, std::chrono::microseconds>();
//        cerr << actual_wait_time2 / 1000.<< " ms Elapsed to compute 1000 affinities" << endl;

///*        vector<double> intBoundaries;
//        for(double d = -120; d < -10; ++d){
//            intBoundaries.push_back(d-0.5);
//        }
//        cout << "   -> Distributions of best energies" << endl;
//        histogramFromDistrib h1(listAff, intBoundaries);
//        cout << h1.print(true) << endl; */
//    }

////    cout << "Best Aff(type Best) so far:" << vmin << endl;
////    for(int i = 0; i < 1; ++i){
////        string Px = randomProt(sizeReceptors+1);
////        cout << "Details of the structures and affinities" << endl; // for " << simpleAccessible << " (Ag=" << AAsimple << ", BCR=" << Px << "), receptors " << receptorSize << " minI=4" << endl;
////        std::pair<double, double> res = T1->affinity(Px, true);
////        cout << "affinity(" << Px << ") = \t" << res.first << "\t" << res.second << endl;
////    }

//}




////string binarySequence::testAffinityFunctions(double L, double R, int maxClusters, int typeAffinityFunction) {

////    stringstream subFolder;
////    subFolder << "L" << L << "_r" << R << nameTypeAff(typeAffinityFunction) << "_C" << maxClusters;
////    string folder = string("C:/Users/Philippe/Desktop/Sequences/") + subFolder.str() + string("/");
////    createFolder(folder);

////    //#define out cout
////    stringstream out;
////    out << "Testing the properties ot the affinity function for the following parameters : \n";
////    out << "   ->     L= " << L << "\t(Size of sequences)" << endl;
////    out << "   ->     R= " << R << "\t(specificity parameter)" << endl;
////    out << "   -> maxCl= " << maxClusters << "\t(cluster size scale)" << endl;
////    switch (typeAffinityFunction) {
////    case seqAff: {out << "   -> Using standard affinity (Saham's)\n"; break;}
////    case seqAffNorm: {out << "   -> Using standard affinity normalized by maxCl^r\n"; break;}
////    case seqAffWindow: {out << "   -> Using the maximum affinity of a sliding window\n"; break;}
////    }
////    out << "==== Part 1 : enumerates all (if possible), or a lot of sequences and sort them by affinity to get the best ones : ====" << endl;

////#define resolutiondistrib 100
////#define maxSequencesToEnumeate 500000

////    // The reference antigen is 00000...
////    binarySequence * ref = new binarySequence(L);

////    // will store a large list of random sequences, with their affinity to ref
////    vector<pair<double, binarySequence*> > store;


////    out << "1 ----------- Distribution of affinities ----------------------" << endl;
////    vector<double> distribution;
////    vector<double> logDistrib;
////    distribution.resize(resolutiondistrib + 1);

////    // in case L is big, will only sample maxSequencesToEnumeate sequences.
////    int total = 0;
////    int maxim = pow(2, L);
////    if (L > 26) {
////        maxim = maxSequencesToEnumeate + 1;           // to avoid it to become negative ...
////    }

////    // if possible to enumerate all, do it one by one.
////    bool enumerateAll = (maxim < maxSequencesToEnumeate);
////    if (enumerateAll) {
////        for (int i = 0; i < maxim; ++i) {
////            binarySequence * a = new binarySequence(L, (long) i);
////            double affi = binarySequence::affinity(a, ref, R, maxClusters, typeAffinityFunction);
////            distribution[(int) (((double) resolutiondistrib) * affi)] += 1.0; // put into the histogram
////            store.push_back(pair<double, binarySequence*> (affi,a));
////            logDistrib.push_back(log10(affi+1e-6));
////            total++;
////        }
////    } else {
////        // if not, than sample randomly
////        for (int i = 0; i < maxSequencesToEnumeate; ++i) {
////            if(((i%100000) == 0) || (i == maxSequencesToEnumeate-1)) cout << i << "/" << maxSequencesToEnumeate << endl;
////            binarySequence * a = new binarySequence(L);
////            a->randomize();
////            double affi = binarySequence::affinity(a, ref, R, maxClusters, typeAffinityFunction);
////            distribution[(int) (((double) resolutiondistrib) * affi)] += 1.0;
////            store.push_back(pair<double, binarySequence*> (affi,a));
////            logDistrib.push_back(log10(affi+1e-6));
////            total++;
////        }
////    }

////    out << "Distribution of affinities\n";
////    for (int i = 0; i < (int) distribution.size(); ++i) {
////        distribution[i] /= (double) total;
////        out << i << "\t"  << double (i) * (1.0 / (double) resolutiondistrib)
////            << "\t" << double (i + 1) * (1.0 / (double) resolutiondistrib)
////            << "\t" << distribution[i] <<  endl;
////    }
////    out << "Distribution of affinities in LOG scale\n";
////    histogramFromDistrib v(logDistrib, 100);
////    out << v.print() << endl;

////    out << "\nSequences and affinity, "
////        << ((enumerateAll) ? " in the order of ID\n" : " randomly generated\n");
////    for (int i = 0; i < 200; ++i) {
////        out << i << "\t" << store[i].second->print() << "\t" << store[i].first << "\n";
////    }

////    out << "\nSequences, sorted from the best, out of the "
////        << min(maxim,(int) maxSequencesToEnumeate) << " evaluated sequences\n";
////    std::sort(store.begin(), store.end(), compSequences);
////    for (int i = 0; i < 200; ++i) {
////        out << i << "\t" << store[i].second->print() << "\t" << store[i].first << "\n";
////    }

////    if (enumerateAll) {
////        out << "\nAffinity of sequences taken randomly\n";
////        for (int i = 0; i < 100; ++i) {
////            binarySequence * seqtmp = new binarySequence(L);
////            seqtmp->randomize();
////            out << i << "\t" << binarySequence::affinity(seqtmp,ref,R,maxClusters,typeAffinityFunction) << "\t" << seqtmp->print() << "\n";
////        }
////    }

////    out << "2 ------------------------- Mutation histograms -----------------------" << endl;

////    // Cuts the distribution of affinities in blocs of 5% of the sequences
////#define sizeClasses 0.05
////    out << "\nEqual classes of affinities (best percentiles)" << endl;
////    vector<double> classes;
////    vector<double> valuesClasses;

////    std::reverse(store.begin(), store.end()); // now to put it increasing

////    // every 5% of sequences browsed (fraction), store the value of affinity
////    double fraction = 0.00;
////    for(unsigned int i = 0; i < store.size(); ++i){
////        //out << store[i].second->print() << " Aff " << store[i].first << endl;
////        if((((double) i / (double) (store.size())) >= fraction - 1e-6) || (i == store.size() - 1)){
////            //out << "]" << fraction - sizeClasses << "," << fraction << "]\t" << store[i].first << endl;
////            valuesClasses.push_back(store[i].first);
////            classes.push_back(fraction);
////            fraction += sizeClasses;
////        }
////    }
////    out << "Affinity classes: cell at percent XX has the affinity value of: (last line=log, 0->-10)" << endl;
////    out << printVec(classes) << endl;
////    out << printVec(valuesClasses) << endl;
////    for(unsigned int i = 0; i < valuesClasses.size(); ++i){
////        out << "\t" << log10(valuesClasses[i] + 1e-10);
////    } out << endl;

////    out << " Now, for each class of sequences, makes the distribution of affinities." << endl;

////    vector<double> newAffinitiesByMutation;      // collecting new affinities only inside one class
////    vector<double> foldAffinitiesByMutation;
////    vector< vector<double>> tableAbsMutations;   // to store the histograms inside each class.
////    vector< vector<double>> tableFoldMutations;
////    vector<double> TOTALnewAffinitiesByMutation; // collecting new affinities for all sequences.
////    vector<double> TOTALfoldAffinitiesByMutation;
////    // For making histograms in term of fold induction, the following groups/classes will be used to make histograms
////    vector<double> classesFoldInd = {0,0.1,0.2,0.4,0.6,0.8,0.9,0.99,1.01,1.1,1.2,1.4,2,2.5,3.25,5,7.5,10};
////    // For making histograms in term of affinity, the classes inside valueclasses will be used.

////    classes.push_back(1.0); // to avoid seg fault

////    // Now will browse again, class by class. Values of sequences inside valuesClasses[currentClass-1] and valuesClasses[currentClass]
////    fraction = sizeClasses;
////    int currentClass = 1;

////    for(unsigned int i = 0; i < store.size(); ++i){
////        binarySequence* thisSeq = store[i].second;
////        double oldAff = binarySequence::affinity(thisSeq, ref, R, maxClusters, typeAffinityFunction);
////        if((oldAff > valuesClasses[currentClass] + 1e-6) || (oldAff < valuesClasses[currentClass-1] - 1e-6) ) cerr << "Class problems, sequence outside its class" << endl;

////        // within one class, add the new possible affinities
////        for(int j = 0; j < thisSeq->size * 2; ++j){ // might get sequence of different lengths !
////            binarySequence seq(thisSeq);
////            seq.mutateOnePosition();
////            double newAff = binarySequence::affinity(&seq, ref, R, maxClusters, typeAffinityFunction);
////            newAffinitiesByMutation.push_back(newAff);
////            if(oldAff > 0) foldAffinitiesByMutation.push_back(newAff / oldAff);
////            TOTALnewAffinitiesByMutation.push_back(newAff);
////            if(oldAff > 0) TOTALfoldAffinitiesByMutation.push_back(newAff / oldAff);
////            //out << thisSeq->print() << " Aff " << oldAff << "->" << seq.print() << " Newaff " << newAff << " Fold " << newAff / oldAff << endl;
////        }
////        if((((double) i / (double) (store.size())) >= fraction - 1e-6) || (i == store.size() - 1)){
////            //out << "Histograms for mutation for the class [" << valuesClasses[currentClass-1] << "," << valuesClasses[currentClass] << "]" << endl;
////            histogramFromDistrib resAbsForThisClass(newAffinitiesByMutation,valuesClasses);
////            tableAbsMutations.push_back(resAbsForThisClass.densities);
////            //out << "New affinities" << endl;
////            //out << resAbsForThisClass.print(true);
////            histogramFromDistrib resFoldForThisClass(foldAffinitiesByMutation,classesFoldInd);
////            tableFoldMutations.push_back(resFoldForThisClass.densities);
////            //out << "Fold increase in the affinity" << endl;
////            //out << resFoldForThisClass.print(true);
////            //out << endl << endl << endl;
////            newAffinitiesByMutation.clear();
////            foldAffinitiesByMutation.clear();
////            currentClass++;
////            fraction += sizeClasses;
////        }

////    }
////    out << "Outputing the affinity changes from all sequences" << endl;
////    histogramFromDistrib TOTALresAbsForThisClass(TOTALnewAffinitiesByMutation,valuesClasses);
////    out << "New affinities" << endl;
////    out << TOTALresAbsForThisClass.print(true);
////    histogramFromDistrib TOTALresFoldForThisClass(TOTALfoldAffinitiesByMutation,classesFoldInd);
////    out << "Fold increase in the affinity" << endl;
////    out << TOTALresFoldForThisClass.print(true);

////    out << "Outputing the affinity changes (new affinity value), from each class of affinity" << endl;
////    out << "\t";
////    for(unsigned int j = 0; j < tableAbsMutations.size(); ++j){
////        out << "[" << valuesClasses[j] << "," << valuesClasses[j+1] << "]\t";
////    }
////    out << "\n";
////    for(unsigned int i = 0; i < TOTALresAbsForThisClass.densities.size(); ++i){
////        out << "[" << TOTALresAbsForThisClass.lowBoundsXs[i] << "," << TOTALresAbsForThisClass.highBoundsXs[i] << "]";
////        for(unsigned int j = 0; j < tableAbsMutations.size(); ++j){
////            out << "\t" << tableAbsMutations[j][i];
////        }
////        out << endl;
////    }

////    out << "Outputing the affinity changes (new affinity value), from each class of affinity" << endl;
////    out << "\t";
////    for(unsigned int j = 0; j < tableFoldMutations.size(); ++j){
////        out << "[" << valuesClasses[j] << "," << valuesClasses[j+1] << "]\t";
////    }
////    out << "\n";
////    for(unsigned int i = 0; i < TOTALresFoldForThisClass.densities.size(); ++i){
////        out << "[" << TOTALresFoldForThisClass.lowBoundsXs[i] << "," << TOTALresFoldForThisClass.highBoundsXs[i] << "]";
////        for(unsigned int j = 0; j < tableFoldMutations.size(); ++j){
////            out << "\t" << tableFoldMutations[j][i];
////        }
////        out << endl;
////    }




////    out << "==== Part 2 : Evaluating cross-reactivity in the system : ====" << endl;

////    for(int k = 0; k < 10; ++k){

////        binarySequence ref2(ref);
////        ref2.mutateOnePosition();
////        ref2.mutateOnePosition();

////        stringstream fname;
////        fname << folder << "DotPlot2CloseAG-L" << L << "_r" << R << nameTypeAff(typeAffinityFunction) << "_C" << maxClusters << "rep" << k << "_" << ref->print() << "to" << ref2.print() << ".txt";
////        fstream of(fname.str().c_str(), ios::out);

////        vector<double> affRef;
////        vector<double> affRef2;

////        for(int i = 0; i < 20000; ++i){
////            binarySequence test(ref2);
////            test.randomize();
////            double aff1 = binarySequence::affinity(&test, ref, R, maxClusters, typeAffinityFunction);
////            double aff2 = binarySequence::affinity(&test, &ref2, R, maxClusters, typeAffinityFunction);
////            affRef.push_back(aff1);
////            affRef2.push_back(aff2);
////            of << aff1 << "\t" << aff2 << endl;
////        }
////        of.close();
////    }

////    for(int k = 0; k < 10; ++k){

////        binarySequence ref2(ref);
////        ref2.mutateOnePosition();
////        ref2.mutateOnePosition();
////        ref2.mutateOnePosition();
////        ref2.mutateOnePosition();
////        ref2.mutateOnePosition();
////        ref2.mutateOnePosition();
////        ref2.mutateOnePosition();

////        stringstream fname;
////        fname << folder << "DotPlot6CloseAG-L" << L << "_r" << R << nameTypeAff(typeAffinityFunction) << "_C" << maxClusters << "rep" << k << "_" << ref->print() << "to" << ref2.print() << ".txt";
////        fstream of(fname.str().c_str(), ios::out);

////        vector<double> affRef;
////        vector<double> affRef2;

////        for(int i = 0; i < 20000; ++i){
////            binarySequence test(ref2);
////            test.randomize();
////            double aff1 = binarySequence::affinity(&test, ref, R, maxClusters, typeAffinityFunction);
////            double aff2 = binarySequence::affinity(&test, &ref2, R, maxClusters, typeAffinityFunction);
////            affRef.push_back(aff1);
////            affRef2.push_back(aff2);
////            of << aff1 << "\t" << aff2 << endl;
////        }
////        of.close();
////    }

////    for(int k = 0; k < 10; ++k){
////        binarySequence ref2(ref);
////        ref2.randomize();

////        stringstream fname;
////        fname << folder << "DotPlotRandCloseAG-L" << L << "_r" << R << nameTypeAff(typeAffinityFunction) << "_C" << maxClusters << "rep" << k << "_" << ref->print() << "to" << ref2.print() << ".txt";
////        fstream of(fname.str().c_str(), ios::out);


////        vector<double> affRef;
////        vector<double> affRef2;

////        for(int i = 0; i < 20000; ++i){
////            binarySequence test(ref2);
////            test.randomize();
////            double aff1 = binarySequence::affinity(&test, ref, R, maxClusters, typeAffinityFunction);
////            double aff2 = binarySequence::affinity(&test, &ref2, R, maxClusters, typeAffinityFunction);
////            affRef.push_back(aff1);
////            affRef2.push_back(aff2);
////            of << aff1 << "\t" << aff2 << endl;
////        }
////        of.close();
////    }

////    ofstream ffin(folder + string("Output.txt"));
////    ffin << out.str();
////    ffin.close();
////    return out.str();


////    int nbAntigens = 10;
////    out << "Generating randomly " << nbAntigens << " antigens " << endl;
////    vector<binarySequence*> ags;
////    for (int i = 0; i < nbAntigens; ++i) {
////        binarySequence * seq = new binarySequence(L);
////        seq->randomize();
////        ags.push_back(seq);
////        out << "\tAg nr " << i << "\t" << seq->print() << endl;
////    }
////    out << "\nNumber of antigens recognized by randomly generated sequences, based on threshold\n";

////    out << "  -> (for the first 100 sequences : ) In the case of random sequences" << endl;
////    total = 0;
////#define thresholdRecoAg 0.1
////    int nbDiscardedSeq = 0;  // sequences that don't recognize anything
////    int countprint = 0;
////    for (int k = 0; k < min(maxim, (int) maxSequencesToEnumeate); ++k) {
////        if (k == 100) {
////            out
////                    <<
////                       "  -> (for the remaining sequences) for sequences recognizing at least an antigen with affinity 0.1"
////                    << endl;
////        }
////        total++;

////        // for each sequence,
////        bool recoAtLeastOne = false;
////        vector<double> nbRecoDepThresh(10, 0.0);
////        vector<double> affinityEach(nbAntigens, 0.0);
////        binarySequence * seqtmp = new binarySequence(L);
////        seqtmp->randomize();
////        for (int j = 0; j < nbAntigens; ++j) {
////            double thisAff = binarySequence::affinity(seqtmp,
////                                                      ags[j],
////                                                      R,
////                                                      maxClusters,
////                                                      typeAffinityFunction);
////            if ((thisAff > thresholdRecoAg) || (k < 100)) {
////                recoAtLeastOne = true;
////            } else { nbDiscardedSeq++; }
////            affinityEach[j] = thisAff;
////            for (int i = 0; i <= (int) (9.99 * thisAff); ++i) {
////                if (i < 10) { nbRecoDepThresh[i]++; }
////            }
////        }
////        if (recoAtLeastOne && (countprint < 5000)) {
////            countprint++;
////            out << "RandSeq " << k << ", " << seqtmp->print() << " ";
////            out << "nbAgPerThreshold:";
////            for (int i = 0; i < 10; ++i) {
////                out << "\t" << nbRecoDepThresh[i];
////            }
////            out << "\taffPerAg:";
////            for (int i = 0; i < nbAntigens; ++i) {
////                out << "\t" << affinityEach[i];
////            }
////            out << endl;
////        }
////        delete seqtmp;
////    }
////    out << "   ... Nb of sequences analyzed: " << total << endl;
////    out << "   ... Nb of sequences discarded: " << nbDiscardedSeq
////        << "(except the 100 first ones, i.e. among the :" << total - 100 << " remaining)" << endl;

////    out << "==== Part 3 : Evaluating the effect of mutations : ====" << endl;

////    binarySequence * start = new binarySequence(L);    // starting by '0000' : the best sequence
////    out << "NbMut\tsequence\taffinity\n";
////    for (int i = 0; i < 2 * L; ++i) {
////        out << i << "\t" << start->print() << "\t" << binarySequence::affinity(start,
////                                                                               ref,
////                                                                               R,
////                                                                               maxClusters,
////                                                                               typeAffinityFunction)
////            << endl;
////        start->mutateOnePosition();
////    }

////    out << "\tReaching a good affinity\t" << endl;
////    start->randomize();
////    double prevaff = binarySequence::affinity(start, ref, R, maxClusters, typeAffinityFunction);

////    bool stop = false;
////    for (int i = 0; (i < L) && (!stop); ++i) {
////        out << "sequence : " << start->print() << "\tAff:\t" << prevaff << "\t";
////        out << "PossibleMut:";
////        vector<int> posGoodMutations;
////        for (int i = 0; i < L; ++i) {
////            binarySequence stmp = binarySequence(start);
////            stmp.content[i] = !stmp.content[i];
////            double newaff
////                    = binarySequence::affinity(&stmp, ref, R, maxClusters, typeAffinityFunction);
////            out << "\t" << newaff;
////            if (newaff > prevaff) { posGoodMutations.push_back(i); }
////        }
////        out << endl;
////        if (posGoodMutations.size() > 0) {
////            int nextmut = random::uniformInteger(0,posGoodMutations.size() - 1);
////            start->content[posGoodMutations[nextmut]] = !start->content[posGoodMutations[nextmut]];
////            prevaff = binarySequence::affinity(start, ref, R, maxClusters, typeAffinityFunction);
////        } else {
////            stop = true;
////        }
////    }


////    for (int i = 0; i < (int) store.size(); ++i) {
////        delete store[i].second;
////    }


////    out << "\tReaching a good affinity\t" << endl;
////    start->randomize();
////    /*double prevaff = binarySequence::seq_affinity(start, ref, R, maxClusters, typeAffinityFunction);

////       bool stop = false;
////       for (int i = 0; (i < L) && (!stop); ++i) {
////          out << "sequence : " << start->print() << "\tAff:\t" << prevaff << "\t";
////          out << "PossibleMut:";
////          vector<int> posGoodMutations;
////          for (int i = 0; i < L; ++i) {
////             binarySequence stmp = binarySequence(start);
////             stmp.content[i] = !stmp.content[i];
////             double newaff
////                = binarySequence::seq_affinity(&stmp, ref, R, maxClusters, typeAffinityFunction);
////             out << "\t" << newaff;
////             if (newaff > prevaff) { posGoodMutations.push_back(i); }
////          }
////          out << endl;
////          if (posGoodMutations.size() > 0) {
////             int nextmut = irandom(posGoodMutations.size() - 1);
////             start->content[posGoodMutations[nextmut]] = !start->content[posGoodMutations[nextmut]];
////             prevaff = binarySequence::seq_affinity(start, ref, R, maxClusters, typeAffinityFunction);
////          } else {
////             stop = true;
////          }
////       }*/



////    out << "2D plot of affinities to 2 antigens" << endl;
////    out << "Case 1 : two close antigens, one mutation away" << endl;

////    return out.str();
////}


