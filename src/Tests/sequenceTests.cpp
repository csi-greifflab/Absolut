#include <iostream>
#include <vector>
#include <string>
using namespace std;

#include "../binaryLineages.h"
#include "../binarySequences.h"
#include "../foldedFree.h"
#include "../foldedCubic.h"
#include "../probabilisticArup.h"
#include "../perelsonSequence.h"
#include "../Tools/distribution.h"
#include "../perelsonGeneralized.h"
#include <algorithm>
#include "../Tools/zaprandom.h"
#include "common.h"
#include <fstream>
void testAllSequences();
string testeAffinityFunctions(double L, double R, int maxClusters, int typeAffinityFunction);

string printVec(vector<double> v){
    stringstream res;
    for(unsigned int i = 0; i < v.size(); ++i){
        res << "\t" << v[i];
    }
    return res.str();
}

//template <class T>
//void GenericTestSequence

/* For all sequence representations, the functions are :
   void randomize();
   void mutateOnePosition();
   static double seq_affinity( *** parameters ***);
   static double hamming(sequence * x, sequence * y);

   Different hypotheses for representation :

   BinarySequences:
   static double seq_affinity(sequence * x, sequence * y, double r, int maxSizeClusters, int type_affinity_function);
   1/ model original with clusters, seqAff

        static double seq_affinity(sequence * x, sequence * y, double r, XXX, seqAff)

   2/ seqAffNorm


   3/ seqAffNorm

   2/ seqAffWindow



   What we want to know:

   1/ for a target epitope sequence,
        what is the best antibody / how many are best
        distribution of affinities of different antibodies
        =>  how many antibodies with initial recognition
            or,what is the threshold for a good basal affinity

   2/ Key mutations:
        when do they happen ?
        distance when doing mutations
        probabilities of affinity change (distribution)

            Show trees ?

   How do we define high affinity
   3/ polyreactivity:
        two sequence (unrelated)
            => distance is more than a few mutations
        Can two antibodies recognize with high affinity ?
            In which cases it can happen ?
        And: is there a continuous path between two antibodies ?
        Can we find a settings of epitopes where it exists ?

   4/ conservation:
        set of a lot of sequences, with a conserved region,
        => how to model it ?
        affinity of antibodies recognizing conserved versus variable
            => what is best to variable,
            => what is best to conserved,
                => what affinities ?

   5/ promiscuous
        a lot of similar sequences
        => can one antibody recognize them all ? or how many ?

   6/ 2 antigens representation:
        => dot plot random populations (log scale ?)
        => dot plots what happens by mutation ?
   heatmaps ...

   */

/*
int main3(void){
    cout << testBinaryLineage(14) << endl;
    //binaryLineage::test();
    //perelsonSequence::test();
    //perelsonGeneralized::testUpgrade();
    //testAllSequences();

    //testeAffinityFunctions(double L, double R, int maxClusters, int typeAffinityFunction);
    //cout << testeAffinityFunctions(30,2,15, seqAffSlidingWindow);
    //        return 0;

    cout << "Starting" << endl;



}*/

//seqAffNorm, seqAffWindow, seqAffSlidingWindow, seqBiggestSubsequence, seqMultiEpitopeAffinity, seqNumberTypesAffFunctions
void testAllSequences(){

}


bool compSequences(pair<double, binarySequence*> a, pair<double, binarySequence*> b) {
   return a.first > b.first;
}

string binarySequence::testAffinityFunctions(double L, double R, int maxClusters, int typeAffinityFunction) {

    stringstream subFolder;
    subFolder << "L" << L << "_r" << R << nameTypeAff(typeAffinityFunction) << "_C" << maxClusters;
    string folder = string("C:/Users/Philippe/Desktop/Sequences/") + subFolder.str() + string("/");
    createFolder(folder);

    //#define out cout
    stringstream out;
    out << "Testing the properties ot the affinity function for the following parameters : \n";
    out << "   ->     L= " << L << "\t(Size of sequences)" << endl;
    out << "   ->     R= " << R << "\t(specificity parameter)" << endl;
    out << "   -> maxCl= " << maxClusters << "\t(cluster size scale)" << endl;
    switch (typeAffinityFunction) {
    case seqAff: {out << "   -> Using standard affinity (Saham's)\n"; break;}
    case seqAffNorm: {out << "   -> Using standard affinity normalized by maxCl^r\n"; break;}
    case seqAffWindow: {out << "   -> Using the maximum affinity of a sliding window\n"; break;}
    }
    out << "==== Part 1 : enumerates all (if possible), or a lot of sequences and sort them by affinity to get the best ones : ====" << endl;

#define resolutiondistrib 100
#define maxSequencesToEnumeate 500000

    // The reference antigen is 00000...
    binarySequence * ref = new binarySequence(L);

    // will store a large list of random sequences, with their affinity to ref
    vector<pair<double, binarySequence*> > store;


    out << "1 ----------- Distribution of affinities ----------------------" << endl;
    vector<double> distribution;
    vector<double> logDistrib;
    distribution.resize(resolutiondistrib + 1);

    // in case L is big, will only sample maxSequencesToEnumeate sequences.
    int total = 0;
    int maxim = pow(2, L);
    if (L > 26) {
        maxim = maxSequencesToEnumeate + 1;           // to avoid it to become negative ...
    }

    // if possible to enumerate all, do it one by one.
    bool enumerateAll = (maxim < maxSequencesToEnumeate);
    if (enumerateAll) {
        for (int i = 0; i < maxim; ++i) {
            binarySequence * a = new binarySequence(L, (long) i);
            double affi = binarySequence::affinity(a, ref, R, maxClusters, typeAffinityFunction);
            distribution[(int) (((double) resolutiondistrib) * affi)] += 1.0; // put into the histogram
            store.push_back(pair<double, binarySequence*> (affi,a));
            logDistrib.push_back(log10(affi+1e-6));
            total++;
        }
    } else {
        // if not, than sample randomly
        for (int i = 0; i < maxSequencesToEnumeate; ++i) {
            if(((i%100000) == 0) || (i == maxSequencesToEnumeate-1)) cout << i << "/" << maxSequencesToEnumeate << endl;
            binarySequence * a = new binarySequence(L);
            a->randomize();
            double affi = binarySequence::affinity(a, ref, R, maxClusters, typeAffinityFunction);
            distribution[(int) (((double) resolutiondistrib) * affi)] += 1.0;
            store.push_back(pair<double, binarySequence*> (affi,a));
            logDistrib.push_back(log10(affi+1e-6));
            total++;
        }
    }

    out << "Distribution of affinities\n";
    for (int i = 0; i < (int) distribution.size(); ++i) {
        distribution[i] /= (double) total;
        out << i << "\t"  << double (i) * (1.0 / (double) resolutiondistrib)
            << "\t" << double (i + 1) * (1.0 / (double) resolutiondistrib)
            << "\t" << distribution[i] <<  endl;
    }
    out << "Distribution of affinities in LOG scale\n";
    histogramFromDistrib v(logDistrib, 100);
    out << v.print() << endl;

    out << "\nSequences and affinity, "
        << ((enumerateAll) ? " in the order of ID\n" : " randomly generated\n");
    for (int i = 0; i < 200; ++i) {
        out << i << "\t" << store[i].second->print() << "\t" << store[i].first << "\n";
    }

    out << "\nSequences, sorted from the best, out of the "
        << min(maxim,(int) maxSequencesToEnumeate) << " evaluated sequences\n";
    std::sort(store.begin(), store.end(), compSequences);
    for (int i = 0; i < 200; ++i) {
        out << i << "\t" << store[i].second->print() << "\t" << store[i].first << "\n";
    }

    if (enumerateAll) {
        out << "\nAffinity of sequences taken randomly\n";
        for (int i = 0; i < 100; ++i) {
            binarySequence * seqtmp = new binarySequence(L);
            seqtmp->randomize();
            out << i << "\t" << binarySequence::affinity(seqtmp,ref,R,maxClusters,typeAffinityFunction) << "\t" << seqtmp->print() << "\n";
        }
    }

    out << "2 ------------------------- Mutation histograms -----------------------" << endl;

    // Cuts the distribution of affinities in blocs of 5% of the sequences
#define sizeClasses 0.05
    out << "\nEqual classes of affinities (best percentiles)" << endl;
    vector<double> classes;
    vector<double> valuesClasses;

    std::reverse(store.begin(), store.end()); // now to put it increasing

    // every 5% of sequences browsed (fraction), store the value of affinity
    double fraction = 0.00;
    for(unsigned int i = 0; i < store.size(); ++i){
        //out << store[i].second->print() << " Aff " << store[i].first << endl;
        if((((double) i / (double) (store.size())) >= fraction - 1e-6) || (i == store.size() - 1)){
            //out << "]" << fraction - sizeClasses << "," << fraction << "]\t" << store[i].first << endl;
            valuesClasses.push_back(store[i].first);
            classes.push_back(fraction);
            fraction += sizeClasses;
        }
    }
    out << "Affinity classes: cell at percent XX has the affinity value of: (last line=log, 0->-10)" << endl;
    out << printVec(classes) << endl;
    out << printVec(valuesClasses) << endl;
    for(unsigned int i = 0; i < valuesClasses.size(); ++i){
        out << "\t" << log10(valuesClasses[i] + 1e-10);
    } out << endl;

    out << " Now, for each class of sequences, makes the distribution of affinities." << endl;

    vector<double> newAffinitiesByMutation;      // collecting new affinities only inside one class
    vector<double> foldAffinitiesByMutation;
    vector< vector<double>> tableAbsMutations;   // to store the histograms inside each class.
    vector< vector<double>> tableFoldMutations;
    vector<double> TOTALnewAffinitiesByMutation; // collecting new affinities for all sequences.
    vector<double> TOTALfoldAffinitiesByMutation;
    // For making histograms in term of fold induction, the following groups/classes will be used to make histograms
    vector<double> classesFoldInd = {0,0.1,0.2,0.4,0.6,0.8,0.9,0.99,1.01,1.1,1.2,1.4,2,2.5,3.25,5,7.5,10};
    // For making histograms in term of affinity, the classes inside valueclasses will be used.

    classes.push_back(1.0); // to avoid seg fault

    // Now will browse again, class by class. Values of sequences inside valuesClasses[currentClass-1] and valuesClasses[currentClass]
    fraction = sizeClasses;
    int currentClass = 1;

    for(unsigned int i = 0; i < store.size(); ++i){
        binarySequence* thisSeq = store[i].second;
        double oldAff = binarySequence::affinity(thisSeq, ref, R, maxClusters, typeAffinityFunction);
        if((oldAff > valuesClasses[currentClass] + 1e-6) || (oldAff < valuesClasses[currentClass-1] - 1e-6) ) cerr << "Class problems, sequence outside its class" << endl;

        // within one class, add the new possible affinities
        for(int j = 0; j < thisSeq->size * 2; ++j){ // might get sequence of different lengths !
            binarySequence seq(thisSeq);
            seq.mutateOnePosition();
            double newAff = binarySequence::affinity(&seq, ref, R, maxClusters, typeAffinityFunction);
            newAffinitiesByMutation.push_back(newAff);
            if(oldAff > 0) foldAffinitiesByMutation.push_back(newAff / oldAff);
            TOTALnewAffinitiesByMutation.push_back(newAff);
            if(oldAff > 0) TOTALfoldAffinitiesByMutation.push_back(newAff / oldAff);
            //out << thisSeq->print() << " Aff " << oldAff << "->" << seq.print() << " Newaff " << newAff << " Fold " << newAff / oldAff << endl;
        }
        if((((double) i / (double) (store.size())) >= fraction - 1e-6) || (i == store.size() - 1)){
            //out << "Histograms for mutation for the class [" << valuesClasses[currentClass-1] << "," << valuesClasses[currentClass] << "]" << endl;
            histogramFromDistrib resAbsForThisClass(newAffinitiesByMutation,valuesClasses);
            tableAbsMutations.push_back(resAbsForThisClass.densities);
            //out << "New affinities" << endl;
            //out << resAbsForThisClass.print(true);
            histogramFromDistrib resFoldForThisClass(foldAffinitiesByMutation,classesFoldInd);
            tableFoldMutations.push_back(resFoldForThisClass.densities);
            //out << "Fold increase in the affinity" << endl;
            //out << resFoldForThisClass.print(true);
            //out << endl << endl << endl;
            newAffinitiesByMutation.clear();
            foldAffinitiesByMutation.clear();
            currentClass++;
            fraction += sizeClasses;
        }

    }
    out << "Outputing the affinity changes from all sequences" << endl;
    histogramFromDistrib TOTALresAbsForThisClass(TOTALnewAffinitiesByMutation,valuesClasses);
    out << "New affinities" << endl;
    out << TOTALresAbsForThisClass.print(true);
    histogramFromDistrib TOTALresFoldForThisClass(TOTALfoldAffinitiesByMutation,classesFoldInd);
    out << "Fold increase in the affinity" << endl;
    out << TOTALresFoldForThisClass.print(true);

    out << "Outputing the affinity changes (new affinity value), from each class of affinity" << endl;
    out << "\t";
    for(unsigned int j = 0; j < tableAbsMutations.size(); ++j){
        out << "[" << valuesClasses[j] << "," << valuesClasses[j+1] << "]\t";
    }
    out << "\n";
    for(unsigned int i = 0; i < TOTALresAbsForThisClass.densities.size(); ++i){
        out << "[" << TOTALresAbsForThisClass.lowBoundsXs[i] << "," << TOTALresAbsForThisClass.highBoundsXs[i] << "]";
        for(unsigned int j = 0; j < tableAbsMutations.size(); ++j){
            out << "\t" << tableAbsMutations[j][i];
        }
        out << endl;
    }

    out << "Outputing the affinity changes (new affinity value), from each class of affinity" << endl;
    out << "\t";
    for(unsigned int j = 0; j < tableFoldMutations.size(); ++j){
        out << "[" << valuesClasses[j] << "," << valuesClasses[j+1] << "]\t";
    }
    out << "\n";
    for(unsigned int i = 0; i < TOTALresFoldForThisClass.densities.size(); ++i){
        out << "[" << TOTALresFoldForThisClass.lowBoundsXs[i] << "," << TOTALresFoldForThisClass.highBoundsXs[i] << "]";
        for(unsigned int j = 0; j < tableFoldMutations.size(); ++j){
            out << "\t" << tableFoldMutations[j][i];
        }
        out << endl;
    }




    out << "==== Part 2 : Evaluating cross-reactivity in the system : ====" << endl;

    for(int k = 0; k < 10; ++k){

        binarySequence ref2(ref);
        ref2.mutateOnePosition();
        ref2.mutateOnePosition();

        stringstream fname;
        fname << folder << "DotPlot2CloseAG-L" << L << "_r" << R << nameTypeAff(typeAffinityFunction) << "_C" << maxClusters << "rep" << k << "_" << ref->print() << "to" << ref2.print() << ".txt";
        fstream of(fname.str().c_str(), ios::out);

        vector<double> affRef;
        vector<double> affRef2;

        for(int i = 0; i < 20000; ++i){
            binarySequence test(ref2);
            test.randomize();
            double aff1 = binarySequence::affinity(&test, ref, R, maxClusters, typeAffinityFunction);
            double aff2 = binarySequence::affinity(&test, &ref2, R, maxClusters, typeAffinityFunction);
            affRef.push_back(aff1);
            affRef2.push_back(aff2);
            of << aff1 << "\t" << aff2 << endl;
        }
        of.close();
    }

    for(int k = 0; k < 10; ++k){

        binarySequence ref2(ref);
        ref2.mutateOnePosition();
        ref2.mutateOnePosition();
        ref2.mutateOnePosition();
        ref2.mutateOnePosition();
        ref2.mutateOnePosition();
        ref2.mutateOnePosition();
        ref2.mutateOnePosition();

        stringstream fname;
        fname << folder << "DotPlot6CloseAG-L" << L << "_r" << R << nameTypeAff(typeAffinityFunction) << "_C" << maxClusters << "rep" << k << "_" << ref->print() << "to" << ref2.print() << ".txt";
        fstream of(fname.str().c_str(), ios::out);

        vector<double> affRef;
        vector<double> affRef2;

        for(int i = 0; i < 20000; ++i){
            binarySequence test(ref2);
            test.randomize();
            double aff1 = binarySequence::affinity(&test, ref, R, maxClusters, typeAffinityFunction);
            double aff2 = binarySequence::affinity(&test, &ref2, R, maxClusters, typeAffinityFunction);
            affRef.push_back(aff1);
            affRef2.push_back(aff2);
            of << aff1 << "\t" << aff2 << endl;
        }
        of.close();
    }

    for(int k = 0; k < 10; ++k){
        binarySequence ref2(ref);
        ref2.randomize();

        stringstream fname;
        fname << folder << "DotPlotRandCloseAG-L" << L << "_r" << R << nameTypeAff(typeAffinityFunction) << "_C" << maxClusters << "rep" << k << "_" << ref->print() << "to" << ref2.print() << ".txt";
        fstream of(fname.str().c_str(), ios::out);


        vector<double> affRef;
        vector<double> affRef2;

        for(int i = 0; i < 20000; ++i){
            binarySequence test(ref2);
            test.randomize();
            double aff1 = binarySequence::affinity(&test, ref, R, maxClusters, typeAffinityFunction);
            double aff2 = binarySequence::affinity(&test, &ref2, R, maxClusters, typeAffinityFunction);
            affRef.push_back(aff1);
            affRef2.push_back(aff2);
            of << aff1 << "\t" << aff2 << endl;
        }
        of.close();
    }

    ofstream ffin(folder + string("Output.txt"));
    ffin << out.str();
    ffin.close();
    return out.str();


    int nbAntigens = 10;
    out << "Generating randomly " << nbAntigens << " antigens " << endl;
    vector<binarySequence*> ags;
    for (int i = 0; i < nbAntigens; ++i) {
        binarySequence * seq = new binarySequence(L);
        seq->randomize();
        ags.push_back(seq);
        out << "\tAg nr " << i << "\t" << seq->print() << endl;
    }
    out << "\nNumber of antigens recognized by randomly generated sequences, based on threshold\n";

    out << "  -> (for the first 100 sequences : ) In the case of random sequences" << endl;
    total = 0;
#define thresholdRecoAg 0.1
    int nbDiscardedSeq = 0;  // sequences that don't recognize anything
    int countprint = 0;
    for (int k = 0; k < min(maxim, (int) maxSequencesToEnumeate); ++k) {
        if (k == 100) {
            out
                    <<
                       "  -> (for the remaining sequences) for sequences recognizing at least an antigen with affinity 0.1"
                    << endl;
        }
        total++;

        // for each sequence,
        bool recoAtLeastOne = false;
        vector<double> nbRecoDepThresh(10, 0.0);
        vector<double> affinityEach(nbAntigens, 0.0);
        binarySequence * seqtmp = new binarySequence(L);
        seqtmp->randomize();
        for (int j = 0; j < nbAntigens; ++j) {
            double thisAff = binarySequence::affinity(seqtmp,
                                                      ags[j],
                                                      R,
                                                      maxClusters,
                                                      typeAffinityFunction);
            if ((thisAff > thresholdRecoAg) || (k < 100)) {
                recoAtLeastOne = true;
            } else { nbDiscardedSeq++; }
            affinityEach[j] = thisAff;
            for (int i = 0; i <= (int) (9.99 * thisAff); ++i) {
                if (i < 10) { nbRecoDepThresh[i]++; }
            }
        }
        if (recoAtLeastOne && (countprint < 5000)) {
            countprint++;
            out << "RandSeq " << k << ", " << seqtmp->print() << " ";
            out << "nbAgPerThreshold:";
            for (int i = 0; i < 10; ++i) {
                out << "\t" << nbRecoDepThresh[i];
            }
            out << "\taffPerAg:";
            for (int i = 0; i < nbAntigens; ++i) {
                out << "\t" << affinityEach[i];
            }
            out << endl;
        }
        delete seqtmp;
    }
    out << "   ... Nb of sequences analyzed: " << total << endl;
    out << "   ... Nb of sequences discarded: " << nbDiscardedSeq
        << "(except the 100 first ones, i.e. among the :" << total - 100 << " remaining)" << endl;

    out << "==== Part 3 : Evaluating the effect of mutations : ====" << endl;

    binarySequence * start = new binarySequence(L);    // starting by '0000' : the best sequence
    out << "NbMut\tsequence\taffinity\n";
    for (int i = 0; i < 2 * L; ++i) {
        out << i << "\t" << start->print() << "\t" << binarySequence::affinity(start,
                                                                               ref,
                                                                               R,
                                                                               maxClusters,
                                                                               typeAffinityFunction)
            << endl;
        start->mutateOnePosition();
    }

    out << "\tReaching a good affinity\t" << endl;
    start->randomize();
    double prevaff = binarySequence::affinity(start, ref, R, maxClusters, typeAffinityFunction);

    bool stop = false;
    for (int i = 0; (i < L) && (!stop); ++i) {
        out << "sequence : " << start->print() << "\tAff:\t" << prevaff << "\t";
        out << "PossibleMut:";
        vector<int> posGoodMutations;
        for (int i = 0; i < L; ++i) {
            binarySequence stmp = binarySequence(start);
            stmp.content[i] = !stmp.content[i];
            double newaff
                    = binarySequence::affinity(&stmp, ref, R, maxClusters, typeAffinityFunction);
            out << "\t" << newaff;
            if (newaff > prevaff) { posGoodMutations.push_back(i); }
        }
        out << endl;
        if (posGoodMutations.size() > 0) {
            int nextmut = random::uniformInteger(0,posGoodMutations.size() - 1);
            start->content[posGoodMutations[nextmut]] = !start->content[posGoodMutations[nextmut]];
            prevaff = binarySequence::affinity(start, ref, R, maxClusters, typeAffinityFunction);
        } else {
            stop = true;
        }
    }


    for (int i = 0; i < (int) store.size(); ++i) {
        delete store[i].second;
    }


    out << "\tReaching a good affinity\t" << endl;
    start->randomize();
    /*double prevaff = binarySequence::seq_affinity(start, ref, R, maxClusters, typeAffinityFunction);

       bool stop = false;
       for (int i = 0; (i < L) && (!stop); ++i) {
          out << "sequence : " << start->print() << "\tAff:\t" << prevaff << "\t";
          out << "PossibleMut:";
          vector<int> posGoodMutations;
          for (int i = 0; i < L; ++i) {
             binarySequence stmp = binarySequence(start);
             stmp.content[i] = !stmp.content[i];
             double newaff
                = binarySequence::seq_affinity(&stmp, ref, R, maxClusters, typeAffinityFunction);
             out << "\t" << newaff;
             if (newaff > prevaff) { posGoodMutations.push_back(i); }
          }
          out << endl;
          if (posGoodMutations.size() > 0) {
             int nextmut = irandom(posGoodMutations.size() - 1);
             start->content[posGoodMutations[nextmut]] = !start->content[posGoodMutations[nextmut]];
             prevaff = binarySequence::seq_affinity(start, ref, R, maxClusters, typeAffinityFunction);
          } else {
             stop = true;
          }
       }*/



    out << "2D plot of affinities to 2 antigens" << endl;
    out << "Case 1 : two close antigens, one mutation away" << endl;

    return out.str();
}





/*void sequenceSpace::testSequenceSpace() {
       // use cerr to catch errors and avoid delays between cout/cerr
        #define OutputTest cerr
       OutputTest << "============ Testing sequences ... ============\n" << endl;
       sequence * a = new sequence(10);
       OutputTest << "Empty sequence  a :" << a->print() << endl;
       a->mutateOnePosition();
       OutputTest << "One Mutation    a :" << a->print() << endl;
       a->mutateOnePosition();
       OutputTest << "Another one     a :" << a->print() << endl;
       a->randomize();
       OutputTest << "Randomized :    a :" << a->print() << endl;
       sequence * b = new sequence(string("01111111111110"));
       OutputTest << "New sequence    b :" << b->print() << endl;
       // b->mutate(0.5);
       OutputTest << "Mutate 50%/base b:" << b->print() << endl;
       sequence * c = new sequence(b);
       OutputTest << "New sequence  c=b :" << c->print() << endl;
       OutputTest << ((sequence::compSeq(b,c)) ? "c equals b" : "problem : c != b") << endl;
       OutputTest << "affinity b-c (r=2): " << sequence::seq_affinity(b, c, 2.0, -1, seqAff) << endl;
       sequence * d = new sequence(string("1111100000"));
       sequence * e = new sequence(string("1111011111"));
       OutputTest << "affinity d-e (r=3): " << sequence::seq_affinity(d, e, 3.0, -1, seqAff)
                  << " between " << d->print() << "\t" << e->print() << endl;
       OutputTest << "hamming(d,e) = " << sequence::hamming(d,e) << endl;

       // Antigen::number_of_bins = 10;
       OutputTest << "Getting type of a sequence : in enum, " << e->getType() << " and as string: "
                  << e->typeInString() << "\t" << e->print() << endl;
       BCR * s1 = new BCR(10);
       OutputTest << "Getting type of a BCR      : in enum, " << s1->getType() << " and as string: "
                  << s1->typeInString() << "\t" << s1->print() << endl;
       TCR * s2 = new TCR(10);
       OutputTest << "Getting type of a TCR      : in enum, " << s2->getType() << " and as string: "
                  << s2->typeInString() << "\t" << s2->print() << endl;
       Antigen * s3 = new Antigen(10);
       OutputTest << "Getting type of an Antigen : in enum, " << s3->getType() << " and as string: "
                  << s3->typeInString() << "\t" << s3->print() << endl;
       s1->add_producer(0.15, 4);
        * s1->add_producer(0.35, 5);
        * s1->rem_producer(0.31, 2);
        * s1->rem_producer(0.0, 1);
        * s1->add_producer(0.0, 2);
        * s1->printNbCells();
        * s1->print();

       OutputTest << "============ Testing sequenceSpace ... ============\n" << endl;

       ofstream ana("testFile.out");
       Parameter par;
       par.Value.size_sequences = 10;
       par.Value.init_antigen_sequences = 8;
       par.Value.initAntigenSeqs.resize(100, string(""));
       par.Value.initAntigenSeqs[0] = string("0000000001");
       par.Value.initAntigenSeqs[0] = string("-1");
       par.Value.max_hamming_antigens = 2;
       par.Value.min_hamming_antigens = 1;
       par.Value.totalBss = 5;  // nb seeder cells
       par.Value.initBCRSeqs.resize(100, string(""));
       par.Value.initBCRSeqs[0] = string("1111111110");
       par.Value.initBCRSeqs[0] = string("-1");
       par.Value.max_hamming_BCRs = 2;
       par.Value.min_initial_affinity_BCRs = 0.1;
       par.Value.max_initial_affinity_BCRs = 1.0;
       par.Value.initTCRSeqs.resize(100, string(""));
       par.Value.initTCRSeqs[0] = string("1111111110");
       par.Value.initTCRSeqs[0] = string("-1");
       par.Value.max_hamming_TCRs = 2;
       par.Value.min_initial_affinity_TCRs = 0.1;
       par.Value.max_initial_affinity_TCRs = 1.0;
       par.Value.R_affinity = 2.0;
       par.Value.pm_differentiation_time = 0.5;
       ana << "hi" << endl;
       sequenceSpace sp(par, ana);

       OutputTest << "Sequence space at initialisation : " << endl;
       OutputTest << sp.printSequences() << endl;

       sequence * newseq = new sequence("0001111110");
       BCR * aBCR = new BCR(newseq);
       BCR * bBCR = new BCR("1100111000");

       OutputTest
                  <<
          "when trying to add a sequence (without type) to the sequencespace, should raise an error. "
                  << endl;
       sp.index_adding_sequence(newseq);   // should raise an error

       OutputTest << endl;
       long newIda = sp.index_adding_sequence(aBCR);   // should raise an error
       long newIdb = sp.index_adding_sequence(bBCR);   // should raise an error
       OutputTest
          <<
          "inserting two BCRs to the space, and got the IDs, and best affinity to antigens (automatically updated :"
          << endl;
       OutputTest << sp.getSequence(newIda)->print() << " with ID " << newIda << "\t"
                  << sp.best_affinity(newIda) << endl;
       OutputTest << sp.getSequence(newIdb)->print() << " with ID " << newIdb << "\t"
                  << sp.best_affinity(newIdb) << endl;

       sp.add_cell(soutext, newIda);
       sp.add_cell(soutext, newIda);
       sp.add_cell(soutext, newIda);
       OutputTest
                  <<
          "After adding three cell to the first BCR sequence, now list of cells (n_cell) for this sequence:"
                  << endl;
       OutputTest << sp.getSequence(newIda)->printNbCells() << endl;
       cout << "differentiating the three cells for dt = 0.5 with pm_differentiation_time = 0.5\n";
       sp.PM_differentiate(soutext, soutextproduce, 0.5);
       OutputTest << "new state of the sequence space : " << endl;
       OutputTest << sp.printSequences(true) << endl;

       Antigen * newA = new Antigen(string("01000001010"));
       Antigen * newB = new Antigen(string("01111001010"));
       long idNewA = sp.add_Antigen(newA);
       long idNewB = sp.add_Antigen(newB);
       OutputTest << "Now, adding two antigens : " << newA->print() << " ID=" << idNewA << "\t"
                  << newB->print() << " ID=" << idNewB << endl;
       OutputTest << sp.printSequences(true) << endl;
    }
    string sequence::typeCellInString(cells index)


            string arupProtein::testeAffinityFunctions() {

             *   stringstream out;
             *   out << "Testing the properties ot the affinity function for the following parameters : \n";
             *   out << "   ->     L= " << L << "\t(Size of arupProteins)" << endl;
             *   out << "   ->     R= " << R << "\t(specificity parameter)" << endl;
             *   out << "   -> maxCl= " << maxClusters << "\t(cluster size scale)" << endl;
             *   switch(typeAffinityFunction){
             *       case seqAff:{out << "   -> Using standard affinity (Saham's)\n"; break;}
             *       case seqAffNorm:{out << "   -> Using standard affinity normalized by maxCl^r\n"; break;}
             *       case seqAffWindow:{out << "   -> Using the maximum affinity of a sliding window\n"; break;}
             *   }
             *
             *   out << "==== Part 1 : enumerates all (if possible), or a lot of arupProteins and sort them by
             * affinity to get the best ones : ====" << endl;
             *
             #define resolutiondistrib 100
             #define maxarupProteinsToEnumeate 1100000
             *
             *   vector< pair<double, arupProtein*> > store;
             *   arupProtein* ref = new arupProtein(L);
             *   vector<double> distribution;
             *   distribution.resize(resolutiondistrib + 1);
             *
             *   int total = 0;
             *   int maxim = pow(2, L);
             *   if(L > 26) maxim = maxarupProteinsToEnumeate + 1;  // to avoid it to become negative ...
             *   bool enumerateAll = (maxim < maxarupProteinsToEnumeate);
             *   if(enumerateAll){
             *       for(int i = 0; i < maxim; ++i){
             *           arupProtein* a = new arupProtein(L, (long) i);
             *           double affi = arupProtein::seq_affinity(a, ref, R, maxClusters, typeAffinityFunction);
             *           distribution[(int) (((double) resolutiondistrib)*affi)] += 1.0;     // put into the
             * histogram classes for affinity
             *           store.push_back( pair<double, arupProtein*>(affi,a));
             *           total++;
             *       }
             *   } else {
             *       for(int i = 0; i < maxarupProteinsToEnumeate; ++i){
             *           arupProtein* a = new arupProtein(L);
             *           a->randomize();
             *           double affi = arupProtein::seq_affinity(a, ref, R, maxClusters, typeAffinityFunction);
             *           distribution[(int) (((double) resolutiondistrib)*affi)] += 1.0;
             *           store.push_back( pair<double, arupProtein*>(affi,a));
             *           total++;
             *       }
             *   }
             *
             *   out << "Distribution of affinities\n";
             *   for(int i = 0; i < (int) distribution.size(); ++i){
             *       distribution[i] /= (double) total;
             *       out << i << "\t" << distribution[i] << "\t" << double(i) * (1.0 / (double)
             * resolutiondistrib) << "\t" << double(i + 1) * (1.0 / (double) resolutiondistrib) << endl;
             *   }
             *
             *   out << "\narupProteins and affinity, " << ((enumerateAll) ? " in the order of ID\n" : "
             * randomly generated\n");
             *   for(int i = 0; i < 200; ++i){
             *       out << i << "\t" << store[i].second->print() << "\t" << store[i].first << "\n";
             *   }
             *
             *   out << "\narupProteins, sorted from the best, out of the " << min(maxim, (int)
             * maxarupProteinsToEnumeate) << " evaluated arupProteins\n";
             *   std::sort(store.begin(), store.end(), comparupProteins);
             *   for(int i = 0; i < 200; ++i){
             *       out << i << "\t" << store[i].second->print() << "\t" << store[i].first << "\n";
             *   }
             *   if(enumerateAll){
             *   out << "\nAffinity of arupProteins taken randomly\n";
             *       for(int i = 0; i < 100; ++i){
             *           arupProtein* seqtmp = new arupProtein(L);
             *           seqtmp->randomize();
             *           out << i << "\t" << arupProtein::seq_affinity(seqtmp, ref, R, maxClusters,
             * typeAffinityFunction) << "\t" << seqtmp->print() << "\n";
             *       }
             *   }
             *   for(int i = 0; i < (int) store.size(); ++i)
             *       delete store[i].second;
             *
             *
             *
             *   out << "==== Part 2 : Evaluating cross-reactivity in the system : ====" << endl;
             *
             *   int nbAntigens = 10;
             *   out << "Generating randomly " << nbAntigens << " antigens " << endl;
             *   vector<arupProtein*> ags;
             *   for(int i = 0; i < nbAntigens; ++i){
             *       arupProtein* seq = new arupProtein(L);
             *       seq->randomize();
             *       ags.push_back(seq);
             *       out << "\tAg nr " << i << "\t" << seq->print() << endl;
             *   }
             *   out << "\nNumber of antigens recognized by randomly generated arupProteins, based on
             * threshold\n";
             *
             *   out << "  -> (for the first 100 arupProteins : ) In the case of random arupProteins" << endl;
             *   total = 0;
             #define thresholdRecoAg 0.1
             *   int nbDiscardedSeq = 0; // arupProteins that don't recognize anything
             *   int countprint = 0;
             *   for(int k = 0; k < min(maxim, (int) maxarupProteinsToEnumeate); ++k){
             *       if(k == 100) out << "  -> (for the remaining arupProteins) for arupProteins recognizing at
             * least an antigen with affinity 0.1" << endl;
             *       total++;
             *
             *       // for each arupProtein,
             *       bool recoAtLeastOne = false;
             *       vector<double> nbRecoDepThresh(10, 0.0);
             *       vector<double> affinityEach(nbAntigens, 0.0);
             *       arupProtein* seqtmp = new arupProtein(L);
             *       seqtmp->randomize();
             *       for(int j = 0; j < nbAntigens; ++j){
             *           double thisAff = arupProtein::seq_affinity(seqtmp, ags[j], R, maxClusters,
             * typeAffinityFunction);
             *           if((thisAff > thresholdRecoAg) || (k < 100)) recoAtLeastOne = true; else nbDiscardedSeq
             * ++;
             *           affinityEach[j] = thisAff;
             *           for(int i = 0; i <= (int) (9.99 * thisAff); ++i){
             *               if(i < 10) nbRecoDepThresh[i] ++;
             *           }
             *       }
             *       if(recoAtLeastOne && (countprint < 5000)){
             *           countprint++;
             *           out << "RandSeq " << k << ", " << seqtmp->print() << " ";
             *           out << "nbAgPerThreshold:";
             *           for(int i = 0; i < 10; ++i){
             *               out << "\t" << nbRecoDepThresh[i];
             *           }
             *           out << "\taffPerAg:";
             *           for(int i = 0; i < nbAntigens; ++i){
             *               out << "\t" << affinityEach[i];
             *           }
             *           out << endl;
             *       }
             *       delete seqtmp;
             *   }
             *   out << "   ... Nb of arupProteins analyzed: " << total << endl;
             *   out << "   ... Nb of arupProteins discarded: " << nbDiscardedSeq << "(except the 100 first
             * ones, i.e. among the :" << total - 100 << " remaining)"<< endl;
             *
             *
             *
             *   out << "==== Part 3 : Evaluating the effect of mutations : ====" << endl;
             *
             *   arupProtein* start = new arupProtein(L);   // starting by '0000' : the best arupProtein
             *   out << "NbMut\tarupProtein\taffinity\n";
             *   for(int i = 0; i < 2*L; ++i){
             *       out << i << "\t" << start->print() << "\t" << arupProtein::seq_affinity(start, ref, R,
             * maxClusters, typeAffinityFunction) << endl;
             *       start->mutateOnePosition();
             *   }
             *
             *   out << "\tReaching a good affinity\t" << endl;
             *   start->randomize();
             *   double prevaff = arupProtein::seq_affinity(start, ref, R, maxClusters, typeAffinityFunction);
             *
             *   bool stop = false;
             *   for(int i = 0; (i < L) && (!stop); ++i){
             *       out << "arupProtein : " << start->print() << "\tAff:\t" << prevaff << "\t";
             *       out << "PossibleMut:";
             *       vector<int> posGoodMutations;
             *       for(int i = 0; i < L; ++i){
             *           arupProtein stmp = arupProtein(start);
             *           stmp.content[i] = !stmp.content[i];
             *           double newaff = arupProtein::seq_affinity(&stmp, ref, R, maxClusters,
             * typeAffinityFunction);
             *           out << "\t" << newaff;
             *           if(newaff > prevaff) posGoodMutations.push_back(i);
             *       }
             *       out << endl;
             *       if(posGoodMutations.size() > 0){
             *           int nextmut = irandom(posGoodMutations.size()-1);
             *           start->content[posGoodMutations[nextmut]] = !
             * start->content[posGoodMutations[nextmut]];
             *           prevaff = arupProtein::seq_affinity(start, ref, R, maxClusters, typeAffinityFunction);
             *       } else {
             *           stop = true;
             *       }
             *   }
             *
             *   return out.str();*/
/*               return string("");
            }
            void arupSpace::testarupSpace() { }
             *
             *   //use cerr to catch errors and avoid delays between cout/cerr
             #define OutputTest cerr
             *   OutputTest << "============ Testing arupProteins ... ============\n" << endl;
             *   arupProtein* a = new arupProtein(10);
             *   OutputTest << "Empty arupProtein  a :" << a->print() << endl;
             *   a->mutateOnePosition();
             *   OutputTest << "One Mutation    a :" << a->print() << endl;
             *   a->mutateOnePosition();
             *   OutputTest << "Another one     a :" << a->print() << endl;
             *   a->randomize();
             *   OutputTest << "Randomized :    a :" << a->print() << endl;
             *   arupProtein* b = new arupProtein(string("01111111111110"));
             *   OutputTest << "New arupProtein    b :" << b->print() << endl;
             *   b->mutate(0.5);
             *   OutputTest << "Mutate 50%/base b:" << b->print() << endl;
             *   arupProtein* c = new arupProtein(b);
             *   OutputTest << "New arupProtein  c=b :" << c->print() << endl;
             *   OutputTest << ((arupProtein::compSeq(b,c)) ? "c equals b" : "problem : c != b") << endl;
             *   OutputTest << "affinity b-c (r=2): " << arupProtein::seq_affinity(b, c, 2.0, -1, seqAff) <<
             * endl;
             *   arupProtein* d = new arupProtein(string("1111100000"));
             *   arupProtein* e = new arupProtein(string("1111011111"));
             *   OutputTest << "affinity d-e (r=3): " << arupProtein::seq_affinity(d, e, 3.0, -1, seqAff) << "
             * between " << d->print() << "\t" << e->print() << endl;
             *   OutputTest << "hamming(d,e) = " << arupProtein::hamming(d,e) << endl;
             *
             *   //Antigen::number_of_bins = 10;
             *   OutputTest << "Getting type of a arupProtein : in enum, " << e->getType() << " and as string: "
             * << e->typeInString() << "\t" << e->print() << endl;
             *   BCR* s1 = new BCR(10);
             *   OutputTest << "Getting type of a BCR      : in enum, " << s1->getType() << " and as string: "
             * << s1->typeInString() << "\t" << s1->print() << endl;
             *   TCR* s2 = new TCR(10);
             *   OutputTest << "Getting type of a TCR      : in enum, " << s2->getType() << " and as string: "
             * << s2->typeInString() << "\t" << s2->print() << endl;
             *   Antigen* s3 = new Antigen(10);
             *   OutputTest << "Getting type of an Antigen : in enum, " << s3->getType() << " and as string: "
             * << s3->typeInString() << "\t" << s3->print() << endl;
             *
             *   s1->add_producer(0.15, 4);
             *   s1->add_producer(0.35, 5);
             *   s1->rem_producer(0.31, 2);
             *   s1->rem_producer(0.0, 1);
             *   s1->add_producer(0.0, 2);
             *   s1->printNbCells();
             *   s1->print();*/

            /*
             *
             *   OutputTest << "============ Testing arupSpace ... ============\n" << endl;
             *
             *   ofstream ana("testFile.out");
             *   Parameter par;
             *   par.Value.size_arupProteins = 10;
             *   par.Value.init_antigen_Sequences = 8;
             *   par.Value.initAntigenSeqs.resize(100, string(""));
             *   par.Value.initAntigenSeqs[0] = string("0000000001");
             *   par.Value.initAntigenSeqs[0] = string("-1");
             *   par.Value.max_hamming_antigens = 2;
             *   par.Value.min_hamming_antigens = 1;
             *   par.Value.totalBss = 5; // nb seeder cells
             *   par.Value.initBCRSeqs.resize(100, string(""));
             *   par.Value.initBCRSeqs[0] = string("1111111110");
             *   par.Value.initBCRSeqs[0] = string("-1");
             *   par.Value.max_hamming_BCRs = 2;
             *   par.Value.min_initial_affinity_BCRs = 0.1;
             *   par.Value.max_initial_affinity_BCRs = 1.0;
             *   par.Value.initTCRSeqs.resize(100, string(""));
             *   par.Value.initTCRSeqs[0] = string("1111111110");
             *   par.Value.initTCRSeqs[0] = string("-1");
             *   par.Value.max_hamming_TCRs = 2;
             *   par.Value.min_initial_affinity_TCRs = 0.1;
             *   par.Value.max_initial_affinity_TCRs = 1.0;
             *   par.Value.R_affinity = 2.0;
             *   par.Value.pm_differentiation_time = 0.5;
             *   ana << "hi" << endl;
             *   arupSpace sp(par, ana);
             *
             *
             *   OutputTest << "arupProtein space at initialisation : "<< endl;
             *   OutputTest << sp.printarupProteins() << endl;
             *
             *   arupProtein* newseq = new arupProtein("0001111110");
             *   BCR* aBCR = new BCR(newseq);
             *   BCR* bBCR = new BCR("1100111000");
             *
             *   OutputTest << "when trying to add a arupProtein (without type) to the arupSpace, should raise
             * an error. " << endl;
             *   sp.index_adding_arupProtein(newseq);  // should raise an error
             *
             *   OutputTest << endl;
             *   long newIda = sp.index_adding_arupProtein(aBCR);  // should raise an error
             *   long newIdb = sp.index_adding_arupProtein(bBCR);  // should raise an error
             *   OutputTest << "inserting two BCRs to the space, and got the IDs, and best affinity to antigens
             * (automatically updated :" << endl;
             *   OutputTest << sp.getSequence(newIda)->print() << " with ID " << newIda << "\t" <<
             * sp.best_affinity(newIda) << endl;
             *   OutputTest << sp.getSequence(newIdb)->print() << " with ID " << newIdb << "\t" <<
             * sp.best_affinity(newIdb) << endl;
             *
             *   sp.add_cell(soutext, newIda);
             *   sp.add_cell(soutext, newIda);
             *   sp.add_cell(soutext, newIda);
             *   OutputTest << "After adding three cell to the first BCR arupProtein, now list of cells (n_cell)
             * for this arupProtein:" << endl;
             *   OutputTest << sp.getSequence(newIda)->printNbCells() << endl;
             *   cout << "differentiating the three cells for dt = 0.5 with pm_differentiation_time = 0.5\n";
             *   sp.PM_differentiate(soutext, soutextproduce, 0.5);
             *   OutputTest << "new state of the arupProtein space : " << endl;
             *   OutputTest << sp.printarupProteins(true) << endl;
             *
             *   Antigen* newA = new Antigen(string("01000001010"));
             *   Antigen* newB = new Antigen(string("01111001010"));
             *   long idNewA = sp.add_Antigen(newA);
             *   long idNewB = sp.add_Antigen(newB);
             *   OutputTest << "Now, adding two antigens : " << newA->print() << " ID=" << idNewA << "\t" <<
             * newB->print() << " ID=" << idNewB << endl;
             *   OutputTest << sp.printarupProteins(true) << endl;
             *
             *  }*/






