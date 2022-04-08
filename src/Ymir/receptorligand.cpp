#include "receptorligand.h"
#include "plot3d.h"
#include <fstream>
#include <sstream>
#include <set>
#include <map>
using namespace std;

#include "../Tools/md5.h"



// very important function :
// for each interaction between the structure s2 and the ligand,
// returns : <posS2, -100+AAligand>
// for inside the structure s2
// returns : <posS2, posS2>. Doublets not allowed: only print it with first < second.
//vector<std::pair<int, int> > interactions(protein& ligand, protein & s2){
//    vector<std::pair<int, int>> res;
//    if(ligand.sequence.size() == 0) return res;
//    if(s2.sequence.size() == 0) return res;

//    int SL = ligand.points.size();
//    int S2L = s2.points.size();
//    for(int i = 0; i < SL; ++i){
//        for(int j = 0; j < S2L; ++j){
//            if(lattice::areNeighbors(ligand.points[i].IDposition, s2.points[j].IDposition)){
//                res.push_back(std::pair<int,int>(j,-NB_AAs-1+(int) ligand.points[i].TypeResidue));
//            }
//        }
//    }
//    for(int i = 0; i < S2L; ++i){
//        for(int j = i+2; j < S2L; ++j){ // j=i1 are neighbors by definition !
//            if(lattice::areNeighbors(s2.points[i].IDposition, s2.points[j].IDposition)){
//                res.push_back(std::pair<int,int>(i,j));
//            }
//        }
//    }
//    return res;
//}

//    stringstream codeInters;
//    vector<std::pair<int, int> > ints = interactions(ligand, possibleReceptor); // for linux g++ compiler, keep a space in > >
//    int nInts = ints.size();
//    for(int j = 0; j < nInts; ++j){
//        codeInters << (char) ('a' + (char) ints[j].first);
//        if(ints[j].second > 0)
//            codeInters << (char) ('a' + (char) ints[j].second);
//        else codeInters << AAname((AA) ((int) NB_AAs + 1 + ints[j].second));
//    }



string interactions(superProtein& ligand, superProtein & s2){
    stringstream codeInters;
    //if(ligand.structure->sequence.size() == 0) return codeInters.str();
    //if(s2.sequence.size() == 0) return codeInters.str();
    if((!ligand.structure) || (!s2.structure) || (ligand.points.size() == 0) || (s2.points.size() == 0)) return codeInters.str();
    vector<std::pair<int, int>> res;
    size_t SL = ligand.points.size();
    size_t S2L = s2.points.size();
    for(size_t i = 0; i < SL; ++i){
        for(size_t j = 0; j < S2L; ++j){
            if(lattice::areNeighbors(ligand.points[i].IDposition, s2.points[j].IDposition)){
                res.push_back(std::pair<int,int>(j,-NB_AAs-1+(int) ligand.points[i].TypeResidue));
                codeInters << (char) ('a' + (char) j) << AAname(ligand.points[i].TypeResidue);
            }
        }
    }
    for(size_t i = 0; i < S2L; ++i){
        for(size_t j = i+2; j < S2L; ++j){ // j=i1 are neighbors by definition !
            if(lattice::areNeighbors(s2.points[i].IDposition, s2.points[j].IDposition)){
                res.push_back(std::pair<int,int>(i,j));
                codeInters << (char) ('a' + (char) i) << (char) ('a' + (char) j);
            }
        }
    }
    return codeInters.str();
}

string interactions(struct3D& ligand, struct3D & s2, string AAsequenceLigand){
    //££tag
    superProtein ProtLigand = superProtein(ligand.sequence, ligand.startingPosition);

    // TODO much faster!
    if(ProtLigand.points.size() != AAsequenceLigand.size()){
        cerr << "ERR: putSequenceLigand: size error, should not happen" << endl;
    }
    for(unsigned int i = 0; i < ProtLigand.points.size(); ++i){
        ProtLigand.points[i].TypeResidue = AA_ID(AAsequenceLigand[i]);
    }
    //££tag
    superProtein ProtRecept = superProtein(s2.sequence, s2.startingPosition);
    return interactions(ProtLigand, ProtRecept);
}


//vector<std::pair<int, int> > interactions(struct3D& ligand, struct3D & s2, string AAsequenceLigand){
//    protein ProtLigand = protein(ligand.sequence, ligand.startingPosition);

//    // TODO much faster!
//    if(ProtLigand.points.size() != AAsequenceLigand.size()){
//        cerr << "ERR: putSequenceLigand: size error, should not happen" << endl;
//    }
//    for(unsigned int i = 0; i < ProtLigand.points.size(); ++i){
//        ProtLigand.points[i].TypeResidue = AA_ID(AAsequenceLigand[i]);
//    }

//    protein ProtRecept = protein(s2.sequence, s2.startingPosition);
//    return interactions(ProtLigand, ProtRecept);
//}

vector<std::pair<int, int>> selfInteractions(superProtein & s2);
vector<std::pair<int, int>> selfInteractions(struct3D & s2){
    //££tag
    superProtein newS = superProtein(s2.sequence, s2.startingPosition);
    return selfInteractions(newS);
}

vector<std::pair<int, int>> selfInteractions(superProtein & s2){
    vector<std::pair<int, int>> res;
    //if(s2.structure->sequence.size() == 0) return res;
    if(s2.points.size() == 0) return res;
    size_t S2L = s2.points.size();
    for(size_t i = 0; i < S2L; ++i){
        for(size_t j = i+2; j < S2L; ++j){ // j=i1 are neighbors by definition !
            if(lattice::areNeighbors(s2.points[i].IDposition, s2.points[j].IDposition)){
                res.push_back(std::pair<int,int>(i,j));
            }
        }
    }
    return res;
}

string codeSelfInteractions(struct3D& s){
    vector<std::pair<int, int> > ints = selfInteractions(s); // for linux g++ compiler, keep a space in > >
    stringstream codeInters;
    size_t nInts = ints.size();
    for(size_t j = 0; j < nInts; ++j){
        codeInters << (char) ('a' + (char) ints[j].first);
        codeInters << (char) ('a' + (char) ints[j].second);
    }
    return codeInters.str();
}





void displayLigand(superProtein* P, vector<int> listForbiddenPositions, bool alwaysVisible, bool doNotLoop){
    #ifdef ALLOW_GRAPHICS
    glDisplay();
    addToDisplay(new superProtein(*P), alwaysVisible);
    set<int>* s = new set<int>(listForbiddenPositions.begin(), listForbiddenPositions.end());
    addToDisplay(s);
    if(!doNotLoop) glutMainLoop();
    #endif
}
void displayLigand(string structureSeq, int startPosition, vector<int> listForbiddenPositions){

#ifdef ALLOW_GRAPHICS
    char* a[] = {(char*) "banana", (char*) "apple"}; // needs an arbitrary argc and argv
    glDisplay(1,a);

    set<int>* merged = generateForbidden(listForbiddenPositions);


    addToDisplay(merged);

    //struct3D* ligand = new struct3D(string("ULSRDSLRDUS"));
    //struct3D* ligand = new struct3D(structureSeq,UnDefined,startPosition);

    //££tag
    superProtein* ligandProt = new superProtein(structureSeq,startPosition);
    cout << "Displaying protein " << print(*ligandProt) << endl;
    //struct3D* ligand = new struct3D(string("ULSRDSLRDUS"));
    //set<int>* ligand2 = new set<int>(ligand->occupiedPositions);

    addToDisplay(ligandProt);

    glutMainLoop();
#endif

}



// Use:
//
set<int>* generateForbidden(vector<int> listForbiddenPositions){
    set<int>* merged = new set<int>();
    if((listForbiddenPositions.size() == 0) || (listForbiddenPositions[0] == POSITIONS_SQUARE_BLOCKED)){
        // by default, we always assume this same forbidden area, so can read again from files...
        // defines a forbidden area (by making a flat protein)
        struct3D* lowB1 = new struct3D(string("LSSSSRRSSSSSSSSSSRRSSSSS"),UnDefined, lattice::idFromPosisition(25,32,32));
        struct3D* lowB2 = new struct3D(string("LSSSSRRSSSSSSSSSSRRSSSSS"),UnDefined, lattice::idFromPosisition(27,32,32));
        struct3D* lowB3 = new struct3D(string("LSSSSRRSSSSSSSSSSRRSSSSS"),UnDefined, lattice::idFromPosisition(29,32,32));
        struct3D* lowB4 = new struct3D(string("LSSSSRRSSSSSSSSSSRRSSSSS"),UnDefined, lattice::idFromPosisition(31,32,32));
        struct3D* lowB5 = new struct3D(string("LSSSSRRSSSSSSSSSSRRSSSSS"),UnDefined, lattice::idFromPosisition(33,32,32));
        struct3D* lowB6 = new struct3D(string("LSSSSRRSSSSSSSSSSRRSSSSS"),UnDefined, lattice::idFromPosisition(35,32,32));
        struct3D* lowB7 = new struct3D(string("LSSSSRRSSSSSSSSSSRRSSSSS"),UnDefined, lattice::idFromPosisition(37,32,32));
        struct3D* lowB8 = new struct3D(string("LSSSSRRSSSSSSSSSSRRSSSSS"),UnDefined, lattice::idFromPosisition(39,32,32));
        set<int> s1 = lowB1->occupiedPositions;
        set<int> s2 = union_sets(lowB2->occupiedPositions, s1);
        set<int> s3 = union_sets(lowB3->occupiedPositions, s2);
        set<int> s4 = union_sets(lowB4->occupiedPositions, s3);
        set<int> s5 = union_sets(lowB5->occupiedPositions, s4);
        set<int> s6 = union_sets(lowB6->occupiedPositions, s5);
        set<int> s7 = union_sets(lowB7->occupiedPositions, s6);
        set<int> s8 = union_sets(lowB8->occupiedPositions, s7);
        delete merged;
        merged = new set<int>(s8);
        delete lowB1;
        delete lowB2;
        delete lowB3;
        delete lowB4;
        delete lowB5;
        delete lowB6;
        delete lowB7;
        delete lowB8;
    }
    if(listForbiddenPositions.size() > 0){
        if(listForbiddenPositions[0] != POSITIONS_ALL_FREE){
            for(size_t i = 0; i < listForbiddenPositions.size(); ++i){
               //cerr << "insert " << listForbiddenPositions[i] << endl;
               if(listForbiddenPositions[i] >= 0) merged->insert(listForbiddenPositions[i]);
            }
        }
    }
    return merged;
}

//note: this function seems redundant with lattice::getIDfromMove, or so.,
int position(int startPos, moveDirection dir){
    vector<int> startPosVec = lattice::positionFromID(startPos);
    vector<int> newPos (3,0);
    vector<int> dirVec = moveVector(dir);
    for(size_t i = 0; i < 3; ++i){
        newPos[i] = startPosVec[i] + dirVec[i];
    }
    return lattice::idFromPosisition(newPos);
}




// Idea : wants all the proteins starting from a point and direction + should not touch ligand + optional forbidden zone
// the initial direction SHOULD ALREADY HAVE BEEN CHECKED FOR TOUCHING !!
// the number of touching points does not include the initial point, but includes the first move point.

// CAREFUL, we mix absolute and relative directions here ...
// TODO make an error in the copy constructor of complex structures, to make sure there is no copy !!
// TODO : replace vector<struct3D*> by vector<string>
// LENGTH = NB of MOVES !!!
#define DBGgenProts 0
vector<struct3D*> generateProteins(int length, int IDposStart, moveDirection absStartDir, int minNbInteract, bool shouldNotTouchProt, struct3D & prot, set<int> & forbiddenVolume){

    // note: the existing struct should finish at IDposStart, and the absStartDir should already have been checked
    // both initial position and the one pointed by the direction should not collide with the prot nor be in forbidden area
    int newPos = lattice::getIdFromAbsoluteMove(IDposStart, absStartDir);
    if(DBGgenProts){ // test input of the function. Don't test touching because this point might be chosen to touch on purpose.
        if((contains(prot.occupiedPositions, IDposStart)) || (contains(forbiddenVolume, IDposStart))){
            cout << "ERR: generateProteins, the starting position " << IDposStart << "(" << printVector(lattice::positionFromID(IDposStart)) << " is inside protein or forbidden places" << endl;
        }
        if((contains(prot.occupiedPositions, newPos)) || (contains(forbiddenVolume, newPos))){
            cout << "ERR: generateProteins, the given direction (reaches pos " << newPos << "(" << printVector(lattice::positionFromID(newPos)) << ") is inside protein or forbidden places" << endl;
        }
        if(shouldNotTouchProt && touch(prot, newPos)){
            cerr << "ERR: generateProteins, the given direction (reaches pos " << newPos << "(" << printVector(lattice::positionFromID(newPos)) << " that touches the protein, while asked not to" << endl;
        }
        if(shouldNotTouchProt && (minNbInteract > 0)){
            cerr <<  "ERR: generateProteins, inconsistency in the options : the protein should not touch, but should at least interact " << minNbInteract << " times" << endl;
            return vector<struct3D*>();
        }
        /*if(existingStruct.endingPosition != IDposStart){
            cerr << "ERR: generateProteins, the existingStructure should end at the 'IDposStart'" << endl;
        }*/
        //if(contains(existingStruct.occupiedPositions, neighbor1)){
        //    cerr << "ERR: generate proteins, the given direction collides with the pre-existing protein part at position (" << neighbor1 << "(" << printVector(lattice::positionFromID(neighbor1)) << " is inside protein or forbidden places " << endl;
        //}
    }

    int nbPoints = nbTouchPoints(prot, newPos);
    if(nbPoints < 0) cerr << "generateSelfFoldings should be called from non-empty start" << endl;
    if(DBGgenProts) cout << " recursive L= " << length << " Position " << printVector(lattice::positionFromID(IDposStart)) << " with dir " << intToMoveChar(absStartDir) << ", touches " << nbPoints << endl;
    vector<struct3D*> res;

    if(length == 0) return res;

    if(length == 1) {
        struct3D* single = new struct3D(string(1,intToMoveChar(absStartDir)), UnDefined, IDposStart); // to make string from char, string(1,c)
        if(nbPoints >= minNbInteract){      // note: the first move is already checked for colliding and touching => just take it
            //if(contains(forbiddenVolume, ))
            res.push_back(single);
        } else {
            delete single;
        }
        if(DBGgenProts) cout << " recursive L= " << length << " Position " << printVector(lattice::positionFromID(IDposStart)) << " with dir " << intToMoveChar(absStartDir) << " should touch " << minNbInteract << " New pt touches " << nbPoints << " , gave " << res.size() << " sequences " << endl;
        return res;
    }

    // for all next move, checks whether it is consistent with 1/ not colliding 2/ not touching if asked


    for(int d1 = 0; d1 < Nb_Moves_Relative; ++d1){ // It is relative because it starts from a Ox and Oy observer coordinates, so only enumerates relative moves (not colliding) and transforms into absolute.

        // there is a simpler function for that, no ?
        moveDirection nextAbsDir = nextAbsoluteMove(absStartDir, initialYaxis(absStartDir), (moveDirection) d1).first;
        if(nextAbsDir == UnDefined) {
            cerr << "ERR: inside generateProteins, couldn't find the nextAbsoluteMove(Ox=" << intToMoveChar(absStartDir) << ", Oy=" << intToMoveChar(initialYaxis(absStartDir)) << ", nextRel=" << intToMoveChar(d1) << ")" << endl;
            exit(-1);
        }
        int neighbor1 = lattice::getIdFromAbsoluteMove(newPos, nextAbsDir);

        if((!contains(prot.occupiedPositions, neighbor1)) && (!contains(forbiddenVolume, neighbor1))){
            if(! (shouldNotTouchProt && (touch(prot, neighbor1)))){
                // here, means that the point is compatible with the touching and not colliding.

                // now, gathers all the proteins possible, from this point,
                int nbPoints = nbTouchPoints(prot, newPos);
                if(DBGgenProts) cout << " Calling recursive " << length - 1 << ", pos =" << printVector(lattice::positionFromID(neighbor1)) << ", dir " << intToMoveChar(nextAbsDir) << endl;
                vector<struct3D*> recursive = generateProteins(length - 1, newPos, nextAbsDir, max(0,minNbInteract - nbPoints) /*,struct3D & existingStruct*/,  shouldNotTouchProt, prot, forbiddenVolume);
                size_t nbS = recursive.size();
                for(size_t i = 0; i < nbS; ++i){
                    struct3D* currentS = recursive[i];
                    if(DBGgenProts) cout << " Receiving from recursive " << length - 1 << ", pos =" << printVector(lattice::positionFromID(neighbor1)) << ", dir " << intToMoveChar(nextAbsDir) << endl;
                    if(DBGgenProts) cout << " -> received : " << currentS->sequence;
                    if(DBGgenProts) if(currentS->sequence.size() > 0) cout << ", start at " << printVector(lattice::positionFromID(currentS->startingPosition)) << endl;
                    if((int) currentS->sequence.size() != length - 1){
                        cerr << "ERR: generateProteins(), recursive call gave prots of size " << currentS->sequence.size() << " instead of length - 1 = " << length - 1 << endl;
                    }
                    // All the proteins are supposed to be 1/ not colliding with ligand; 2/ of good length and 3/ with proper min nr of interactions
                    // needs to be checked : - that they don't auto-fold, i.e. they don't touch the current point. Upstream resursive calls will solve previous positions
                    if((!contains(currentS->occupiedPositions, IDposStart))){
                        string newS = fuse(string(1,intToMoveChar(absStartDir)), currentS->sequence);

                        /* test fuse */
                        if(currentS->sequence.size() > 0){
                            struct3D st1(newS, UnDefined, IDposStart);
                            struct3D stBase(string(1,intToMoveChar(absStartDir)), UnDefined, IDposStart);
                            struct3D stExtension(currentS->sequence, UnDefined, currentS->startingPosition);
                            set<int> merge = union_sets(stBase.occupiedPositions, stExtension.occupiedPositions);
                            std::set<int>::iterator it;
                            bool foundErr = false;
                            if(merge.size() != st1.occupiedPositions.size()){
                                foundErr = true;
                            }
                            for(it = merge.begin(); it != merge.end(); ++it){
                                if(!contains(st1.occupiedPositions, *it)){
                                    foundErr = true;
                                }
                            }
                            if(foundErr){
                                cout << "ERR on fuse (" << string(1,intToMoveChar(absStartDir)) << ", " << currentS->sequence << " gave " << newS << endl;
                                cout << "BASE " << print(stBase) << endl;
                                cout << "EXT " << print(stExtension) << endl;
                                cout << "GIVEN BY FUSE " << print(st1) << endl;
                            }
                        }
                        res.push_back(new struct3D(newS, UnDefined, IDposStart));
                    }
                    delete currentS;
                }
            }
        }
    }
    if(DBGgenProts) cout << " recursive L= " << length << " Position " << printVector(lattice::positionFromID(IDposStart)) << " with dir " << intToMoveChar(absStartDir) << " should touch " << minNbInteract << " New pt touches " << nbPoints << " , gave " << res.size() << " sequences " << endl;
    return res;
}




void receptorLigand::generateReceptors(){
    //cerr << "Start Generate" << endl;
    set<int> startingPositions = neighborPositions(ligand, forbiddenVolume);
    int nR = 0;

    if(DBGgenRecept) cout << "for each starting position (" << startingPositions.size() << ")" << endl;
    std::set<int>::iterator it;
    size_t counter = 0;
    int cptPos = 0;
    cout << "   -> Start computing structures from each bordering position. Nr of bordering/Starting positions: " << startingPositions.size() << endl;
    for (it = startingPositions.begin(); it != startingPositions.end(); ++it)
    {
        if(DBGgenRecept) cout << " ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----" << endl;
        cptPos++;
        int pos = *it;
        if(!contains(forbiddenVolume, pos)){

            // either the protein starts by an interaction, so starts here. Only one direction
            for(int d1 = 0; d1 < Nb_Moves_Absolute; ++d1){
                int neighbor1 = position(pos, (moveDirection) d1);
                if((!contains(ligand.occupiedPositions, neighbor1)) && (!contains(forbiddenVolume, neighbor1))){

                    vector<struct3D*> toAdd = generateProteins(sizeReceptor, pos, (moveDirection) d1, minimalNInteract - 1, false, ligand, forbiddenVolume);
                    size_t NTA = toAdd.size();
                    for(size_t i = 0; i < NTA; ++i){
                        possibleReceptors.push_back(toAdd[i]);
                        nR++;
                    }
                    if(DBGgenRecept) cout << nR << "\tPosition " << printVector(lattice::positionFromID(pos));
                    if(DBGgenRecept) cout << " Dir " << intToMoveChar(d1) << ", " << NTA << " sequences starting from an interaction there" << endl;
                }
            }

            // for each couple of possible directions ...
            vector<std::pair<moveDirection, moveDirection> > possibleDirections;
            for(int d1 = 0; d1 < Nb_Moves_Absolute; ++d1){
                for(int d2 = d1+1; d2 < Nb_Moves_Absolute; ++d2){// the two moves should be different.

                    // conditions on the two directions :
                    // direction 1 should not touch again the protein,
                    //           and the following structure should not touch the protein either
                    // direction 2 can touch again.
                    // if dir1 and dir2 don't touch, explore in both orders.
                    int neighbor1 = position(pos, (moveDirection) d1);
                    int neighbor2 = position(pos, (moveDirection) d2);

                    // if both directions are free
                    if((!contains(ligand.occupiedPositions, neighbor1)) &&
                            (!contains(ligand.occupiedPositions, neighbor2)) && (!contains(forbiddenVolume, neighbor1)) && (!contains(forbiddenVolume, neighbor2))){
                        // condition : first direction should not touch. If both do so, both order possible.
                        if(!touch(ligand, neighbor1)){
                            possibleDirections.push_back(std::pair<moveDirection, moveDirection>((moveDirection) d1, (moveDirection) d2));
                        }
                        if(!touch(ligand, neighbor2)){
                            possibleDirections.push_back(std::pair<moveDirection, moveDirection>((moveDirection) d2, (moveDirection) d1));
                        }
                    }
                }
            }
            size_t nDir = possibleDirections.size();
            for(size_t i = 0; i < nDir; ++i){
                for(int subLength = 1; subLength < sizeReceptor; subLength++){
                    // note: the function 'proteins' should be dynamically programmed to memorize.
                   // if want to display left and right  vector<struct3D*> possibleStructLeft = ProteinsNoTouch(pos, possibleDirections[i].first, subLength-1);
                   // if want to display left and right  vector<struct3D*> possibleStructRight = ProteinsCanTouch(pos, possibleDirections[i].first, sizeReceptor-subLength-1);
                    // could save time by computing the smallest set first, and not computing the next one if empty
                    vector<struct3D*> structRight = generateProteins(sizeReceptor - subLength, pos, possibleDirections[i].second, minimalNInteract-1 /* because already touches ! */, false, ligand, forbiddenVolume);
                    size_t NSR = structRight.size();
                    if(NSR > 0){
                        vector<struct3D*> structLeft = generateProteins(subLength, pos, possibleDirections[i].first, 0, true, ligand, forbiddenVolume);
                        size_t NSL = structLeft.size();
                        if(DBGgenRecept) cout << nR << "\tPosition " << printVector(lattice::positionFromID(pos));
                        if(DBGgenRecept) cout << " Dir " << intToMoveChar(possibleDirections[i].first) << "," << intToMoveChar(possibleDirections[i].second);
                        if(DBGgenRecept) cout << " Fuse: L=" << subLength << " (nb= " << NSL  << ") x L=" << sizeReceptor - subLength << " (nb= " << NSR << ")\t->(before foldtest) " << NSL * NSR << endl;
                        if(NSL * NSR > 0){
/*                            for(int s1 = 0; s1 < NSL; ++s1){
                                addToDisplay(structLeft[s1], false);
                            }
                            for(int s2 = 0; s2 < NSR; ++s2){
                                addToDisplay(structRight[s2], false);
                            }*/
                            for(size_t s1 = 0; s1 < NSL; ++s1){
                                for(size_t s2 = 0; s2 < NSR; ++s2){
                                    //int startPosComb = structLeft[s1]->points[structLeft[s1]->points.size() - 1].IDposition;
                                    int startPosComb = structLeft[s1]->endingPosition;
                                    struct3D* test = new struct3D(fuse(revert(structLeft[s1]->sequence), structRight[s2]->sequence), UnDefined, startPosComb);
                                    if(test->properlyFolded) {
                                        possibleReceptors.push_back(test);
                                        //addToDisplay(test, false);
                                        nR++;
                                    }
                                    else delete test;
                                }
                            }
                            for(size_t s1 = 0; s1 < NSL; ++s1){
                                delete structLeft[s1];
                            }
                            for(size_t s2 = 0; s2 < NSR; ++s2){
                                delete structRight[s2];
                            }
                        }
         /*               if(nR > 1000){
                            glutMainLoop();
                        }*/
                        // merge compatible structures
                    }
                }
            }
        }
        cout << "   ... " << cptPos << "/" << startingPositions.size() << " completed, starting pos " << pos << "(" << printVector(lattice::positionFromID(pos)) << "), + " << possibleReceptors.size() - counter << " structures -> " << possibleReceptors.size() << " total " << endl;
        counter = possibleReceptors.size();
    }
    if(DBGgenRecept) cout << "FINAL number of generated receptors : " << nR << " - " << possibleReceptors.size() << endl;
}






// Two possible strategies:
// 1 Branch and bound: from the previously generated structure, add one by one by calling the remaining tail with min nb of interactions.
//    interest: eliminates the structures as fast as there is a problem.
//          problem: 1/ can not do dynamical programming because it always is a different pre-formed structyre
//                   2/ at the beginning, a very high number of recursive will be called
// 2 Divide and conquer:
//          for each direction, enumerate all right and left, (size half), and then merges. Problem: can not predict where the middle will form.
// ultilmate divide and conquer, with the existing half tail !!!
// 3

// do branch and bound.
vector<struct3D*> generateSelfFoldings(int length, int minNbInteract){
    vector<struct3D*> res;
    std::set<string> normalizedSeqStructures;
    int red = 0;
    for(int i = 1; i < length; ++i){
        string seq = "";
        for(int j = 0; j < i; ++j){
            seq += string("S");
        }

        //        seq += string("U"); => instead, the first move will be up
        // CAREFUL, WE ENUMERATE EVERYTHING TWICE, even with the receptor !!
        struct3D startingTail = struct3D(seq, UnDefined, lattice::centralPosition());
        vector<struct3D*> resSSS = generateSelfFoldings(length - i, startingTail.endingPosition, Up, minNbInteract, startingTail);
        // cout << seq <<"U*: " << resSSS.size() << " found" << endl;

        for(size_t j = 0; j < resSSS.size(); ++j){
            struct3D* fused = new struct3D(fuse(seq,resSSS[j]->sequence),UnDefined,lattice::centralPosition());
            string thisOne = normalizeAbsolute(revert(fused->sequence));
            set<string>::iterator it = normalizedSeqStructures.find(thisOne);
            if(it != normalizedSeqStructures.end()){
                //cout << fused->sequence << " was already enumerated from other side as " << thisOne << endl;
                red++;
            } else {
                normalizedSeqStructures.insert(fused->sequence); // anyway they start by
                //cout << fused->sequence << " added (rev=" << thisOne << ")" << endl;
                res.push_back(fused);
            }
            delete resSSS[j];
        }
    }     
    //cout << red << "/" << res.size() << " discarded" << endl;
    return res;
}


vector<struct3D*> generateSelfFoldings(int length, int IDposStart, moveDirection absStartDir, int minNbInteract, struct3D & alreadyFoldedTail){

    // trick: already inside interactions + the asked remaining nb interactions should always be = initial min Nb Interact
    // note: the existing struct should finish at IDposStart, and the absStartDir should already have been checked
    // both initial position and the one pointed by the direction should not collide with the prot nor be in forbidden area
    int newPos = lattice::getIdFromAbsoluteMove(IDposStart, absStartDir);   
    if((contains(alreadyFoldedTail.occupiedPositions, newPos)) /* || (contains(forbiddenVolume, newPos))*/ ){
        cout << "ERR: generateProteins, the given direction (reaches pos " << newPos << "(" << printVector(lattice::positionFromID(newPos)) << ") is inside protein or forbidden places" << endl;
    }
    // Number of self-contacts added by the new point
    int nbPoints = nbTouchPoints(alreadyFoldedTail, newPos) - 1; // to remove the covalent link
    if(DBGgenProts) cout << " BranchBoundFold L= " << length << " Position " << printVector(lattice::positionFromID(IDposStart)) << " with dir " << intToMoveChar(absStartDir) << ", touches " << nbPoints << endl;
    vector<struct3D*> res;

    //if(contains(forbiddenVolume, IDposStart) || contains(forbiddenVolume, newPos)) return res;

    if(length == 0) return res;

    if(length == 1) { // $$ could create single inside the if...
        struct3D* single = new struct3D(string(1,intToMoveChar(absStartDir)), UnDefined, IDposStart); // to make string from char, string(1,c)
        if(nbPoints >= minNbInteract){
            res.push_back(single);
        } else {
            delete single;
        }
        if(DBGgenProts) cout << " recursive L= " << length << " Position " << printVector(lattice::positionFromID(IDposStart)) << " with dir " << intToMoveChar(absStartDir) << " should touch " << minNbInteract << " New pt touches " << nbPoints << " , gave " << res.size() << " sequences " << endl;
        return res;
    }

    // for all next move, checks whether it is consistent with 1/ not colliding 2/ min nb of touching
    for(int d1 = 0; d1 < Nb_Moves_Relative; ++d1){ // It is relative because it starts from a Ox and Oy observer coordinates, so only enumerates relative moves (not colliding) and transforms into absolute.

        // there is a simpler function for that, no ?
        moveDirection nextAbsDir = nextAbsoluteMove(absStartDir, initialYaxis(absStartDir), (moveDirection) d1).first;
        if(nextAbsDir == UnDefined) {
            cerr << "ERR: inside generateSelfFoldings, couldn't find the nextAbsoluteMove(Ox=" << intToMoveChar(absStartDir) << ", Oy=" << intToMoveChar(initialYaxis(absStartDir)) << ", nextRel=" << intToMoveChar(d1) << ")" << endl;
            exit(-1);
        }
        int neighbor1 = lattice::getIdFromAbsoluteMove(newPos, nextAbsDir);

        if((!contains(alreadyFoldedTail.occupiedPositions, neighbor1)) /*&& (!contains(forbiddenVolume, neighbor1))*/){
            //if(! (shouldNotTouchProt && (touch(alreadyFoldedTail, neighbor1)))){
                // here, means that the point is compatible with the touching and not colliding.

                // now, gathers all the proteins possible, from this point,
                //int nbPoints2 = nbTouchPoints(alreadyFoldedTail, neighbor1);
                if(DBGgenProts) cout << " Calling recursive " << length - 1 << ", pos =" << printVector(lattice::positionFromID(neighbor1)) << ", dir " << intToMoveChar(nextAbsDir) << endl;

                // fuse à faire, ici !!
                struct3D expandedCurrentFolding = struct3D(alreadyFoldedTail);
                bool succeded = expandedCurrentFolding.pushBackAbsoluteMove(absStartDir); // fuse could also have worked
                if(succeded && expandedCurrentFolding.properlyFolded){ // properlyfolder should be ok because tested before
                    vector<struct3D*> recursive = generateSelfFoldings(length - 1, newPos, nextAbsDir, max(0,minNbInteract - nbPoints), expandedCurrentFolding);
                    size_t nbS = recursive.size();
                    for(size_t i = 0; i < nbS; ++i){
                        struct3D* currentS = recursive[i];
                        if(DBGgenProts) cout << " Receiving from recursive " << length - 1 << ", pos =" << printVector(lattice::positionFromID(neighbor1)) << ", dir " << intToMoveChar(nextAbsDir) << endl;
                        if(DBGgenProts) cout << " -> received : " << currentS->sequence;
                        if(DBGgenProts) if(currentS->sequence.size() > 0) cout << ", start at " << printVector(lattice::positionFromID(currentS->startingPosition)) << endl;
                        if(static_cast<int>(currentS->sequence.size()) != length - 1){
                            cerr << "ERR: generateSelfFoldings(), recursive call gave prots of size " << currentS->sequence.size() << " instead of length - 1 = " << length - 1 << endl;
                        }
                        // All the proteins are supposed to be 1/ not colliding with the pre-folded structure; 2/ of good length and 3/ with proper min nr of interactions
                        // needs to be checked : - that they don't auto-fold, i.e. they don't touch the current point. Upstream resursive calls will solve previous positions
                        if((!contains(currentS->occupiedPositions, IDposStart))){
                            string newS = fuse(string(1,intToMoveChar(absStartDir)) /*expandedCurrentFolding.sequence*/, currentS->sequence);

                            /* test fuse */
    /*                        if(currentS->sequence.size() > 0){
                                struct3D st1(newS, UnDefined, IDposStart);
                                struct3D stExtension(currentS->sequence, UnDefined, currentS->points[0].IDposition);
                                set<int> merge = union_sets(stBase.occupiedPositions, stExtension.occupiedPositions);
                                std::set<int>::iterator it;
                                bool foundErr = false;
                                if(merge.size() != st1.occupiedPositions.size()){
                                    foundErr = true;
                                }
                                for(it = merge.begin(); it != merge.end(); ++it){
                                    if(!contains(st1.occupiedPositions, *it)){
                                        foundErr = true;
                                    }
                                }
                                if(foundErr){
                                    cout << "ERR on fuse (" << string(1,intToMoveChar(absStartDir)) << ", " << currentS->sequence << " gave " << newS << endl;
                                    cout << "BASE " << print(stBase) << endl;
                                    cout << "EXT " << print(stExtension) << endl;
                                    cout << "GIVEN BY FUSE " << print(st1) << endl;
                                }
                            }*/
                            res.push_back(new struct3D(newS, UnDefined, IDposStart));
                        }
                        delete currentS;
                    }
                }
            //}
        }
    }
    if(DBGgenProts) cout << " recursive L= " << length << " Position " << printVector(lattice::positionFromID(IDposStart)) << " with dir " << intToMoveChar(absStartDir) << " should touch " << minNbInteract << " New pt touches " << nbPoints << " , gave " << res.size() << " sequences " << endl;
    return res;
}


#define DBGgenRecept 0
//    receptorLigand::receptorLigand(struct3D & _ligand, int _sizeReceptor, int _minimalNInteract, set<int> & _forbiddenVolume) :
//        ligand(_ligand), sizeReceptor(_sizeReceptor), minimalNInteract(_minimalNInteract), forbiddenVolume(_forbiddenVolume) {}

// inputs:
//struct3D ligand;
//int sizeReceptor;
//int minimalNInteract;
//set<int> forbiddenVolume;


// Want :
// table : for different ligands :
//      L1: easy, accessible,
//      L2: ligand with easy parts and hard parts
//      L3: ligands random with different sizes
//      L4: HIV: lots of forbidden places => compromise between outside and inside.
//              -> also generate possible mutants.
//
// Want :   X=size receptor,
//          Y=Min Nb interactions (to antigen)
//          Table of number of sequences and number of compressedInteract. => will help decide which constraints make sense.
//          Table of the number

// analysis things filled by printToFile ...


// make it 'append to file', to be called during generateReceptors and free memory...
/*int receptorLigand::printToFile(string fnameAllInfo, string fnameCompressed){
    if(ligand.sequence.size() == 0) return 0;
    ofstream res(fnameAllInfo.c_str()); // with linux g++ compiler, do not say ofstream res = ofstream(...) but rather use this syntax. operator = is private
    res << "Ligand=\t" << ligand.startingPosition << "\t" << ligand.sequence << endl;
    res << "ForbiddenVolume\t" << forbiddenVolume.size();
    std::set<int>::iterator it;
    for(it = forbiddenVolume.begin(); it != forbiddenVolume.end(); ++it){
        res << "\t" << *it;
    }
    res << endl;
    res << "LatticeWidthXYZ\t" << XWidth << "\t" << YWidth << "\t" << ZWidth << endl;
    int nRecept = possibleReceptors.size();
    res << "Receptors\tSize=\t" << sizeReceptor << "\tMinInteract=\t" << minimalNInteract << endl;
    res << nRecept << endl;

    set<string> compressedReceptors; // only based on the interactions. Two receptors with the same interaction profile (inside & outside) can be merged
    for(int i = 0; i < nRecept; ++i){
        res << possibleReceptors[i]->points[0].IDposition << "\t" << possibleReceptors[i]->sequence << "\t";
        stringstream codeInters;
        //codeInters << nbTouchPoints(*(possibleReceptors[i]), ligand);
        vector<std::pair<int, int> > ints = XXXXinteractions(ligand, *possibleReceptors[i]); // for linux g++ compiler, keep a space in > >
        int nInts = ints.size();
        //codeInters << nInts;
        for(int j = 0; j < nInts; ++j){
            codeInters << ints[j].first;
            if(ints[j].second > 0)
                codeInters << ints[j].second;
            else codeInters << AAname((AA) ((int) NB_AAs + 1 + ints[j].second));
        }
        compressedReceptors.insert(codeInters.str());
        res << codeInters.str() << endl;
    }
    res.close();
    ofstream res2 (fnameCompressed.c_str());
    std::set<string>::iterator it2;
    res2 << compressedReceptors.size() << "\t" << AAsequence << "\t" << sizeReceptor << endl;
    for(it2 = compressedReceptors.begin(); it2 != compressedReceptors.end(); ++it2){
        res2 << *it2 << endl;
    }
    res2.close();
    return compressedReceptors.size();
}*/



// The last argument is a small trick. I made sure only the residues of the ligand are used (and their position) throughout fastaffinities and
// the enumeration of receptors, such that unconventional proteins can be used (ony one chain). So, old way: the file contains the structure.
// New way: force the ligand, then the structure stored in the file is not used to reconstruct the ligand
int reGenerateCompressedStructures(string fnameAllStructures, string newSequenceLigand, string fnameOutAll, string fnameOutCompressed, superProtein* forceLigand){
    ifstream fstr(fnameAllStructures);
    if(!fstr) {cerr << "ERR: reGenerateCompressedStructures, file not found:" << fnameAllStructures << endl; return -1;}
    string trash;
    int posLig = -1; string strLig;
    fstr >> trash >> posLig >> strLig;

    superProtein* ligand; // still problems with deleting structures when copying a ligand => Prefer to use pointers
    if(forceLigand == nullptr){
        //struct3D ligand = struct3D(strLig, UnDefined, posLig);
        ligand = new superProtein(strLig, posLig);
        ligand->setAAs(newSequenceLigand);
    } else {
        cout << "       Re-reading file " << fnameAllStructures << "\n       But with new/special ligand " << forceLigand->getAAseq() <<  endl; //print(*forceLigand) << endl;
        //ligand = superProtein(*forceLigand);
        ligand = forceLigand;
    }

    int nRecept = -1;
    fstr >> trash >> nRecept;
    cout << "       " << nRecept << " structures to read" << endl;
    if((nRecept < 0) || (nRecept > 1e10)){
        cerr << "ERR: reGenerateCompressedStructures , " << fnameAllStructures << ", should start by an appropriate number of structures (and < 1e10). got " << nRecept << endl;
        return -1;
    }

    #define writeFAll
    #ifdef writeFAll
    ofstream res(fnameOutAll.c_str());

    // Bad way to
    res << "Ligand=\t" << posLig << "\t" << strLig << endl;
    res << "ForbiddenVolume_lookAtFileName\t" << 0;
    //std::set<int>::iterator it;
    //for(it = forbiddenVolume.begin(); it != forbiddenVolume.end(); ++it){
    //    res << "\t" << *it;
    //}
    res << endl;
    res << "LatticeWidthXYZ\t" << XWidth << "\t" << YWidth << "\t" << ZWidth << endl;
    res << "Receptors\tSize=\t" << -1 << "\tMinInteract=\t" << -1 << endl;
    res << nRecept << endl;
    #endif

    int sizeReceptor = -1;
    // this is copy paste from the next function (print to File). Maybe make a function, later
    cout << "      ";
    map<string, int> compressedReceptors;   // only based on the interactions. Two receptors with the same interaction profile (inside & outside) can be merged
    for(int i = 0; i < nRecept; ++i){
        int pos; string seqStr;
        fstr >> pos >> seqStr;

        if((i % 100000) == 0) cout << ".";
        if(sizeReceptor == -1){sizeReceptor = seqStr.size();}
        else {if((int) seqStr.size() != sizeReceptor) cerr << "ERR: reGenerateCompressedStructures, non-match of receptor sizes " << endl;}

        //struct3D possibleReceptor = struct3D(seqStr, UnDefined, pos);
        //££tag
        superProtein possibleReceptor = superProtein(seqStr, pos);
        #ifdef writeFAll
        res << pos << "\t" << seqStr << "\t";
        #endif
//        stringstream codeInters;
//        // codeInters << nbTouchPoints(*(possibleReceptors[i]), ligand);
//        // Interactions = <pos, pos> or <pos, ligandAA>. Position encoded a...z, ligandAA encoded A..Z
//        vector<std::pair<int, int> > ints = interactions(ligand, possibleReceptor); // for linux g++ compiler, keep a space in > >
//        int nInts = ints.size();
//        //codeInters << nInts;
//        for(int j = 0; j < nInts; ++j){
//            codeInters << (char) ('a' + (char) ints[j].first);
//            if(ints[j].second > 0)
//                codeInters << (char) ('a' + (char) ints[j].second);
//            else codeInters << AAname((AA) ((int) NB_AAs + 1 + ints[j].second));
//        }
        string codeInters = interactions(*ligand, possibleReceptor); // for linux g++ compiler, keep a space in > >

        if(compressedReceptors.find(codeInters) != compressedReceptors.end()){
            compressedReceptors[codeInters]++;
        } else {
            compressedReceptors[codeInters] = 1;
        }
        #ifdef writeFAll
        res << codeInters << endl;
        #endif
    }
    #ifdef writeFAll
    res.close();
    #endif

    ofstream res2 (fnameOutCompressed.c_str());
    std::map<string,int>::iterator it2;
    res2 << compressedReceptors.size() << "\t" << newSequenceLigand << "\t" << sizeReceptor << endl;
    for(it2 = compressedReceptors.begin(); it2 != compressedReceptors.end(); ++it2){
        res2 << it2->first << "\t" << it2->second << endl;
    }
    res2.close();

    fstr.close();
    return compressedReceptors.size();
}


// This is how the output files look like:

/** File1 (structures): DRRDSRSRSDSDSRSRRLLRSSRSSRSSDDSSRSSRS-7-14-40f251eb48cc63b45d63da402b6a5828Structures
Ligand=	145376	DRRDSRSRSDSDSRSRRLLRSSRSSRSSDDSSRSSRS
Receptors=	929
137054	SULSSUU
137054	SULSUDD
137054	SULUDSD
137054	SULUDDU
137054	SULUDDR
137054	ULSSDDS
...

File 2 (Compact): Depends on antigen structure + antigen AA sequence. Regenerated from the structures file when Antigen AA sequence changes.
868
HILKGALYVPPVYHILVKSPLIVIAAKLGRALYPNVGR	7
aHbIaGaLcLdYaVeVePeYgYdHhHcIbKgGfRbe	1
aHbIaGaLcLdYaVeVgVhPePeYdHcIbKfRhRbe	1
aHbIcLaGcGdAaLeLfYaVgVgPgYfHeIdLcVbKhRbebg	1
aHbIcPdVeSeLfIfIdLcAhLhPgNgGbRbebgch	1
aHbIcPdVeSeLfIfIgAgKdLcAhLhPbRbechdg	1
aHbIcPdVfIfIeAeKdLcAhLhPgNgGbRbgcheh	1
aHbIcPdVhIhIeAeKdLcAfLfPgNgGbRbgcfeh	1
aHbIcPfVgSgLhIhIeAeKfLcAdLdPbRbgcfeh	1
aHbIcPhVfIfIgAgKhLcAdLdPeNeGbRbechdg	1
aHbIcPhVgSgLfIfIeAeKhLcAdLdPbRbgcheh	1
aHbIcPhVgSgLfIfIhLcAdLdPeNeGbRbebgch	1
aHbIdVcScLfIfIeAeKdLhLhPgNgGbRbgcfeh	1
...

File 3 (Full Details): gives the interaction profile for each structure => So it's possible to retreive all structures that have a particular interaction profile
Ligand=	145376	DRRDSRSRSDSDSRSRRLLRSSRSSRSSDDSSRSSRS
ForbiddenVolume	204	132761	132762	132763	132764	132765	132766	132767	132768	132769	132770	132771	132772	132773	132774	132775	132776	132825	132826	132827	132828	132829	132830	132831	132832	132833	132834	132835	132836	132837	132838	132839	132840	132889	132890	132891	132892	132893	132894	132895	132896	132897	132898	132899	132900	132901	132902	132903	132904	132953	132954	132955	132956	132957	132958	132959	132960	132961	132962	132963	132964	132965	132966	132967	132968	133017	133018	133019	133020	133021	133022	133023	133024	133025	133026	133027	133028	133029	133030	133031	133032	133081	133082	133083	133084	133085	133086	133087	133088	133089	133090	133091	133092	133093	133094	133095	133096	133145	133146	133147	133148	133149	133150	133151	133152	133153	133154	133155	133156	133157	133158	133159	133160	133209	133210	133211	133212	133213	133214	133215	133216	133217	133218	133219	133220	133221	133222	133223	133224	133273	133274	133275	133276	133277	133278	133279	133280	133281	133282	133283	133284	133285	133286	133287	133288	133337	133338	133339	133340	133341	133342	133343	133344	133345	133346	133347	133348	133349	133350	133351	133352	133401	133402	133403	133404	133405	133406	133407	133408	133409	133410	133411	133412	133413	133414	133415	133416	133465	133466	133467	133468	133469	133470	133471	133472	133473	133474	133475	133476	133477	133478	133479	133480	137375	137376	137377	137378	137379	137380	141471	141472	141473	141474	141475	141476
LatticeWidthXYZ	64	64	64
Receptors	Size=	7	MinInteract=	14
929
137054	SULSSUU	eIdLdGhLgYfVfPfYgHhIaLbVdVeKeh
137054	SULSUDD	eIdLdGfLgYhVhPhYgHfIaLbVdVeKeh
137054	SULUDSD	dLdGeAfLgYhVhPhYgHfIaLeLbVdV
137054	SULUDDU	gIdLdGeAfLhVhPhYfIaLeLbVdVgKdg
137054	SULUDDR	hHgIdLdGhGeAfLhLhVfIaLeLbVdVgKdg
137054	ULSSDDS	gIhLhGcAdLeYfVfPfYeHdIaLcLhVgKchdg
137054	ULSSDDR	hHgIhGcAdLhLeYfVhVfPfYeHdIaLcLgKdg
137054	ULDUSSL	eIdLdGcAfVfPfYhYaLcLdVeKhGgR
137054	ULDUSUU	eIdLdGcAhLgYfVfPfYgHhIaLcLdVeKcheh
137054	ULDUUDD	eIdLdGcAfLgYhVhPhYgHfIaLcLdVeKcfeh
...
**/

void receptorLigand::printToFile(string fnameOnlyStruct){
    //////if(ligand.sequence.size() == 0) return 0;
    ofstream resStruct(fnameOnlyStruct);

    size_t nRecept = possibleReceptors.size();
    resStruct << "Ligand=\t" << ligand.startingPosition << "\t" << ligand.sequence << endl;
    resStruct << "Receptors=\t" << nRecept << endl;
    for(size_t i = 0; i < nRecept; ++i){
        resStruct << possibleReceptors[i]->startingPosition << "\t" << possibleReceptors[i]->sequence << "\n";
    }
    resStruct.close();
}

// old function, don't use it, works only for continuous structures
int receptorLigand::printToFile(string fnameOnlyStruct, string fnameAllInfo, string fnameCompressed){
    printToFile(fnameOnlyStruct);

    size_t nRecept = possibleReceptors.size();
    ofstream res(fnameAllInfo.c_str()); // with linux g++ compiler, do not say ofstream res = ofstream(...) but rather use this syntax. operator = is private
    res << "Ligand=\t" << ligand.startingPosition << "\t" << ligand.sequence << endl;
    res << "ForbiddenVolume\t" << forbiddenVolume.size();
    std::set<int>::iterator it;
    for(it = forbiddenVolume.begin(); it != forbiddenVolume.end(); ++it){
        res << "\t" << *it;
    }
    res << endl;
    res << "LatticeWidthXYZ\t" << XWidth << "\t" << YWidth << "\t" << ZWidth << endl;
    res << "Receptors\tSize=\t" << sizeReceptor << "\tMinInteract=\t" << minimalNInteract << endl;
    res << nRecept << endl;


    map<string, int> compressedReceptors; // only based on the interactions. Two receptors with the same interaction profile (inside & outside) can be merged
    for(int i = 0; i < nRecept; ++i){
        res << possibleReceptors[i]->startingPosition << "\t" << possibleReceptors[i]->sequence << "\n";
        // resStruct << possibleReceptors[i]->startingPosition << "\t" << possibleReceptors[i]->sequence << "\n";
//        stringstream codeInters;
//        //codeInters << nbTouchPoints(*(possibleReceptors[i]), ligand);
//        // Interactions = <pos, pos> or <pos, ligandAA>. Position encoded a...z, ligandAA encoded A..Z
//        vector<std::pair<int, int> > ints = interactions(ligand, *possibleReceptors[i], AAsequence); // for linux g++ compiler, keep a space in > >
//        int nInts = ints.size();
//        //codeInters << nInts;
//        for(int j = 0; j < nInts; ++j){
//            codeInters << (char) ('a' + (char) ints[j].first);
//            if(ints[j].second > 0)
//                codeInters << (char) ('a' + (char) ints[j].second);
//            else codeInters << AAname((AA) ((int) NB_AAs + 1 + ints[j].second));
//        }
        string codeInters = interactions(ligand, *possibleReceptors[i], AAsequence);
        if(compressedReceptors.find(codeInters) != compressedReceptors.end()){
            compressedReceptors[codeInters]++;
        } else {
            compressedReceptors[codeInters] = 1;
        }
        res << codeInters << endl;
    }
    res.close();
    //resStruct.close();
    ofstream res2 (fnameCompressed.c_str());
    std::map<string,int>::iterator it2;
    res2 << compressedReceptors.size() << "\t" << AAsequence << "\t" << sizeReceptor << endl;
    for(it2 = compressedReceptors.begin(); it2 != compressedReceptors.end(); ++it2){
        res2 << it2->first << "\t" << it2->second << endl;
    }
    res2.close();
    return compressedReceptors.size();
}

//vector<struct3D*> possibleReceptors;

//string AAsequence;
void receptorLigand::putSequenceLigand(string _AAsequence){
    AAsequence = _AAsequence;
    if(ligand.sequence.size() != (AAsequence.size() - 1)){
        cerr << "ERR: putSequenceLigand(" << _AAsequence << " of size " << _AAsequence.size() <<"), is not compatible with the structure of the ligand : " << ligand.sequence << ", of size " << ligand.sequence.size() << " moves. And there should be one more AA than moves." << endl;
        return;
    }
    /*if(ligand.points.size() != AAsequence.size()){
        cerr << "ERR: putSequenceLigand: size error, should not happen" << endl;
    }
    for(unsigned int i = 0; i < ligand.points.size(); ++i){
        ligand.points[i].TypeResidue = AA_ID(AAsequence[i]);
    }*/
}

void receptorLigand::setForbiddenVolume(vector<int> positions){
    forbiddenVolume.clear();
    for(unsigned int i = 0; i < positions.size(); ++i){
        forbiddenVolume.insert(positions[i]);
    }
}

void receptorLigand::setForbiddenVolume(set<int> positions){
    forbiddenVolume.clear();
    for(set<int>::iterator it = positions.begin(); it != positions.end(); ++it){
        forbiddenVolume.insert(*it);
    }
}


// the IDs start at 0
void showStructures(string fileCompact, bool containsThirdColumn, int fromNr, int toNr){


    cout << "Opening " << fileCompact << endl;
    ifstream f(fileCompact);
    if(!f) {cerr << "File not found " << fileCompact << endl; return;}

    #ifdef ALLOW_GRAPHICS
    char *c[] = {(char*)"Hello",nullptr};
    glDisplay(0,c);
    #else
    cerr << "WRN: You are calling showStructures when ALLOW_GRAPHICS is not #defined, inside plot3d.h -> will just print to screen" << endl;
    #endif

    string structure;
    string interactionInfos;
    int startPos;
    int cpt = 0;
    while((f >> startPos)){
        f >> structure;
        if(containsThirdColumn) f >> interactionInfos;

        if((cpt >= fromNr) && (cpt <= toNr)){
            superProtein* test = new superProtein(structure);

            #ifdef ALLOW_GRAPHICS
            addToDisplay(test, false);
            #else
            cout << print(*test) << "\n";
            #endif
        }
        cpt++;
    }
    #ifdef ALLOW_GRAPHICS
    glutMainLoop();
    #endif
}

void testRecepLigands(){
    string compactS1 = "SSUSSULRDRDUS";
    string absS1 = "BSSULRDRDUS";
    string absS2 = "LSSULRDRDUS";
    string absS3 = "SSSULRDRDUS";

    // test conversion to absolute and back:
    cout << "Compact  " << compactS1 << " ID=" << relativeToInt(compactS1) << " - " << intToRelative(relativeToInt(compactS1)) << endl;
    cout << "Absolute " << absS1     << " ID=" << absoluteToInt(absS1)   << " - " << intToAbsolute(absoluteToInt(absS1)) << " => encoded as " << intToRelative(absoluteToInt(absS1)) << endl;
    cout << "Absolute " << absS2     << " ID=" << absoluteToInt(absS2)   << " - " << intToAbsolute(absoluteToInt(absS2)) << " => encoded as " << intToRelative(absoluteToInt(absS2)) << endl;
    cout << "Absolute " << absS3     << " ID=" << absoluteToInt(absS3)   << " - " << intToAbsolute(absoluteToInt(absS3)) << " => encoded as " << intToRelative(absoluteToInt(absS3)) << endl;

    testRotate();

    cout << "test compact " << endl;
    cout << string("BUSLUSD") << "->" << intToRelative(compacte(absoluteToInt(string("BUSLUSD")))) << endl;
    cout << string("SUSLUSD") << "->" << intToRelative(compacte(absoluteToInt(string("SUSLUSD")))) << endl;
    cout << string("BSSSDULUSD") << "->" << intToRelative(compacte(absoluteToInt(string("BSSSDULUSD")))) << endl;
    cout << string("SSSUUSD") << "->" << intToRelative(compacte(absoluteToInt(string("SSSUUSD")))) << endl;
    cout << string("DSSSSSS") << "->" << intToRelative(compacte(absoluteToInt(string("DSSSSSS")))) << endl;
    cout << string("BDULUSD") << "->" << intToRelative(compacte(absoluteToInt(string("BDULUSD")))) << endl;
    // note: at some point, make a function to check absolute sequences.

    cout << "testing position (start, move) (note: absolute moves without memory)" << endl;
    int startPos = lattice::centralPosition();
    int nextPos = 0;
    string move = string("BDUUSLLUSSDDD");
    for(unsigned int i = 0; i < move.size(); ++i){
        nextPos = position(startPos,(moveDirection) charMoveToInt(move[i]));
        cout << printVector(lattice::positionFromID(startPos)) << "(" << startPos << ") -> Move:" << move[i] << " -> " << printVector(lattice::positionFromID(nextPos)) << " (" << nextPos << ")" << endl;
        startPos = nextPos;
    }

    struct3D s1 = struct3D(string("BUULDSSS"));
    set<int> emptyset = set<int>();
    receptorLigand rl = receptorLigand(s1, 6, 3, emptyset);
    rl.minimalNInteract = 2;
    rl.sizeReceptor = 5;
    rl.generateReceptors();

/*
    struct3D ligand;
    set<int> forbiddenVolume;

    vector<struct3D*> possibleReceptors;
    int sizeReceptor;
    int minimalNInteract;

    void generateReceptors(){
    }*/
}

//vector<struct3D*> generateProteins(int length, int IDposStart, moveDirection absStartDir, int minNbInteract, bool shouldNotTouchProt, struct3D & prot, set<int> & forbiddenVolume){

void testGenerateProteins(){
    // defines a forbidden area (by making a flat protein)
    struct3D* lowB1 = new struct3D(string("LSSSSRRSSSSSSSSSSRRSSSSS"),UnDefined, lattice::idFromPosisition(25,32,32));
    struct3D* lowB2 = new struct3D(string("LSSSSRRSSSSSSSSSSRRSSSSS"),UnDefined, lattice::idFromPosisition(27,32,32));
    struct3D* lowB3 = new struct3D(string("LSSSSRRSSSSSSSSSSRRSSSSS"),UnDefined, lattice::idFromPosisition(29,32,32));
    struct3D* lowB4 = new struct3D(string("LSSSSRRSSSSSSSSSSRRSSSSS"),UnDefined, lattice::idFromPosisition(31,32,32));
    struct3D* lowB5 = new struct3D(string("LSSSSRRSSSSSSSSSSRRSSSSS"),UnDefined, lattice::idFromPosisition(33,32,32));
    struct3D* lowB6 = new struct3D(string("LSSSSRRSSSSSSSSSSRRSSSSS"),UnDefined, lattice::idFromPosisition(35,32,32));
    struct3D* lowB7 = new struct3D(string("LSSSSRRSSSSSSSSSSRRSSSSS"),UnDefined, lattice::idFromPosisition(37,32,32));
    struct3D* lowB8 = new struct3D(string("LSSSSRRSSSSSSSSSSRRSSSSS"),UnDefined, lattice::idFromPosisition(39,32,32));
    set<int> s1 = lowB1->occupiedPositions;
    set<int> s2 = union_sets(lowB2->occupiedPositions, s1);
    set<int> s3 = union_sets(lowB3->occupiedPositions, s2);
    set<int> s4 = union_sets(lowB4->occupiedPositions, s3);
    set<int> s5 = union_sets(lowB5->occupiedPositions, s4);
    set<int> s6 = union_sets(lowB6->occupiedPositions, s5);
    set<int> s7 = union_sets(lowB7->occupiedPositions, s6);
    set<int> s8 = union_sets(lowB8->occupiedPositions, s7);
    set<int>* merged = new set<int>(s8);

    //set<int>* s1 = new set<int>(neighborPositions(*lowB1));

    #ifdef ALLOW_GRAPHICS
    char *c[] = {(char*)"Hello",nullptr};
    glDisplay(0,c);
    /*addToDisplay(lowB1);
    addToDisplay(lowB2);
    addToDisplay(lowB3);
    addToDisplay(lowB4);
    addToDisplay(lowB5);
    addToDisplay(lowB6);
    addToDisplay(lowB7);
    addToDisplay(lowB8);*/
    addToDisplay(merged);
    //std::set_union()
    #endif

    superProtein* ligand = new superProtein(string("ULSRDSLRDUS"));
    //set<int>* ligand2 = new set<int>(ligand->occupiedPositions);

    #ifdef ALLOW_GRAPHICS
    addToDisplay(ligand);
    #endif

    /*
    vector<struct3D*> res = generateProteins(2, lattice::idFromPosisition(33,32,34), Left, 1, false, *ligand, *merged);
    cout << "Now testing proteins from position " << lattice::idFromPosisition(33,32,34) << " that should interact with the protein at least once" << endl;
    int Nres = res.size();
    cout << "Nr of solutions found : " << Nres << endl;
    for(int i = 0; i < Nres; ++i){
        addToDisplay(res[i]);
        cout << print(* res[i]) << endl;
    }
    */

    vector<struct3D*> res = generateProteins(05, lattice::idFromPosisition(33,32,35), Left, 0, true, *ligand->structure, *merged);
    cout << "Now testing proteins from position " << lattice::idFromPosisition(33,32,35) << " that should interact with the protein at least once" << endl;
    int Nres = res.size();
    cout << "Nr of solutions found : " << Nres << endl;
    for(int i = 0; i < Nres; ++i){
        #ifdef ALLOW_GRAPHICS
        superProtein* protForDisplay = new superProtein(*(res[i]));
        addToDisplay(protForDisplay, false);
        #endif
        //cout << print(* res[i]) << endl;
    }

    #ifdef ALLOW_GRAPHICS
    glutMainLoop();
    #endif

}

string codePos(vector<int> pos){
    stringstream res;
    for(unsigned int i = 0; (i < 2) && (i < pos.size()); ++i){
        res << pos[i] << "-";
    }
    res << pos.size();
    if(res.str().size() > 5) return md5(res.str());
    return res.str();
}

string shrt(string longString){
    if(longString.size() > 50) return longString.substr(0,5) + md5(longString);
    return longString;
}


string hashCode(superProtein &p){
    // string that concatenates ID residue, AA and position for each residue, then makes a hash
    stringstream cct;
    int NP = p.points.size();
    cct<<NP;
    for(int i = 0; i < NP; ++i){
        cct << p[i].IDresidue << AAname(p[i].TypeResidue) << p[i].IDposition;
    }
    return md5(cct.str());
}

string hashWithoutAAs(superProtein &p){
    // string that concatenates ID residue, AA and position for each residue, then makes a hash
    stringstream cct;
    int NP = p.points.size();
    cct<<NP;
    for(int i = 0; i < NP; ++i){
        cct << p[i].IDresidue << p[i].IDposition;
    }
    return md5(cct.str());
}


// careful, threshold is in energy, should be lower
string fnameLibrary(string ligandStructSeq, string agSeq, int sizeReceptors, int minimalNInteract, double threshold, vector<int> forbiddenPos){
    stringstream fname;
    fname << "Library" << shrt(ligandStructSeq) << "+" << shrt(agSeq) << "-" << sizeReceptors << "-" << threshold << "-" << minimalNInteract << "-" << codePos(forbiddenPos) << ".txt";
    return fname.str();
}
string fnameStructures(superProtein* ligand, int sizeReceptors, int minimalNInteract, vector<int> forbiddenPos){
    stringstream fname;
    if(!ligand->structure) return string("Ligand has empty structure!");
    if(ligand->contiguous()){
        string ligandStructSeq = ligand->structure->sequence;
        fname << shrt(ligandStructSeq) << "-" << sizeReceptors << "-" << minimalNInteract << "-"  << codePos(forbiddenPos) << "Structures.txt";
    } else {
        fname << hashWithoutAAs(*ligand) << "-" << sizeReceptors << "-" << minimalNInteract << "-" << codePos(forbiddenPos) << "Structures.txt";
    }
    return fname.str();
}

string fnameStructuresAndCompactForAASeqLigand(superProtein* ligand, int sizeReceptors, int minimalNInteract, vector<int> forbiddenPos){
    stringstream fname;
    if(!ligand->structure) return string("Ligand has empty structure!");
    string ligandAAseq = ligand->getAAseq();
    if(ligand->contiguous()){
        string ligandStructSeq = ligand->structure->sequence;
        fname << manualStructureFilesLocation << shrt(ligandStructSeq) << "+" << shrt(ligandAAseq) << "-" << sizeReceptors << "-" << minimalNInteract << "-" << codePos(forbiddenPos) << ".txt";
    } else {
        fname << manualStructureFilesLocation << hashWithoutAAs(*ligand) << "+" << shrt(ligandAAseq) << "-" << sizeReceptors << "-" << minimalNInteract << "-" << codePos(forbiddenPos) << ".txt";
    }
    return fname.str();
}

string fileNameCompactForAASeqLigand(superProtein* ligand, int sizeReceptors, int minimalNInteract, vector<int> forbiddenPos){
    stringstream fname;
    if(!ligand->structure) return string("Ligand has empty structure!");
    string ligandAAseq = ligand->getAAseq();
    if(ligand->contiguous()){
        string ligandStructSeq = ligand->structure->sequence;
        fname << manualStructureFilesLocation << shrt(ligandStructSeq) << "+" << shrt(ligandAAseq) << "-" << sizeReceptors << "-" << minimalNInteract << "-"  << codePos(forbiddenPos) << "Compact.txt";
    } else {
        fname << manualStructureFilesLocation << hashWithoutAAs(*ligand) << "+" << shrt(ligandAAseq) << "-" << sizeReceptors << "-" << minimalNInteract << "-" << codePos(forbiddenPos) << "Compact.txt";
    }
    return fname.str();
}



void testBigReceptorLigand(int sizeReceptor, int minimalNInteract, string structSeq, int startPos, vector<int> listForbidden){


    // defines a forbidden area (by making a flat protein)
//    struct3D* lowB1 = new struct3D(string("LSSSSRRSSSSSSSSSSRRSSSSS"),UnDefined, lattice::idFromPosisition(25,32,32));
//    struct3D* lowB2 = new struct3D(string("LSSSSRRSSSSSSSSSSRRSSSSS"),UnDefined, lattice::idFromPosisition(27,32,32));
//    struct3D* lowB3 = new struct3D(string("LSSSSRRSSSSSSSSSSRRSSSSS"),UnDefined, lattice::idFromPosisition(29,32,32));
//    struct3D* lowB4 = new struct3D(string("LSSSSRRSSSSSSSSSSRRSSSSS"),UnDefined, lattice::idFromPosisition(31,32,32));
//    struct3D* lowB5 = new struct3D(string("LSSSSRRSSSSSSSSSSRRSSSSS"),UnDefined, lattice::idFromPosisition(33,32,32));
//    struct3D* lowB6 = new struct3D(string("LSSSSRRSSSSSSSSSSRRSSSSS"),UnDefined, lattice::idFromPosisition(35,32,32));
//    struct3D* lowB7 = new struct3D(string("LSSSSRRSSSSSSSSSSRRSSSSS"),UnDefined, lattice::idFromPosisition(37,32,32));
//    struct3D* lowB8 = new struct3D(string("LSSSSRRSSSSSSSSSSRRSSSSS"),UnDefined, lattice::idFromPosisition(39,32,32));
//    set<int> s1 = lowB1->occupiedPositions;
//    set<int> s2 = union_sets(lowB2->occupiedPositions, s1);
//    set<int> s3 = union_sets(lowB3->occupiedPositions, s2);
//    set<int> s4 = union_sets(lowB4->occupiedPositions, s3);
//    set<int> s5 = union_sets(lowB5->occupiedPositions, s4);
//    set<int> s6 = union_sets(lowB6->occupiedPositions, s5);
//    set<int> s7 = union_sets(lowB7->occupiedPositions, s6);
//    set<int> s8 = union_sets(lowB8->occupiedPositions, s7);
    set<int>* merged = generateForbidden(listForbidden);// new set<int>(s8);

    //struct3D* ligand = new struct3D(string("ULSRDSLRDUS"));
    struct3D ligand1 = struct3D(structSeq, UnDefined, startPos);
    superProtein* ligand = new superProtein(ligand1);
    //set<int>* s1 = new set<int>(neighborPositions(*lowB1));

    //int sizeReceptor = 6;
    ///for(sizeReceptor = 10; sizeReceptor < 11; ++sizeReceptor){

        //for(int minimalNInteract = 2; minimalNInteract < 12 ; ++minimalNInteract){
            //int minimalNInteract = 4;
            cout << "Performing testBigReceptorLigand with size receptor " << sizeReceptor << " and " << minimalNInteract << " Minimal number of interactions " << endl;
            receptorLigand* a = new receptorLigand(*ligand->structure, sizeReceptor, minimalNInteract, *merged);
            string AAseq = randomProt(structSeq.size()+1);// "AGNTGYMPARNW";

            //cout << print(a->ligand) << endl;
            a->generateReceptors();
            cout << " Total Nb : " << sizeReceptor << "\t" << minimalNInteract << "\t" << a->possibleReceptors.size()  <<  endl;
            //stringstream fname; fname << ligand->sequence << sizeReceptor << "-" << minimalNInteract;

            string fOnlyStruct  = fnameStructures(ligand, sizeReceptor, minimalNInteract, listForbidden);
            string fDetail  = fnameStructuresAndCompactForAASeqLigand(ligand,  sizeReceptor, minimalNInteract, listForbidden);
            string fCompact = fileNameCompactForAASeqLigand(ligand, sizeReceptor, minimalNInteract, listForbidden);

            a->putSequenceLigand(AAseq);
            cout << "After Reduc " << a->printToFile(fOnlyStruct, fDetail, fCompact) << "\n";
            int NRs = a->possibleReceptors.size();
            vector<int> classesNbInteract(50, 0);
            for(int i = 0; i < NRs; ++i){
                classesNbInteract[nbTouchPoints(*(a->possibleReceptors[i]), *ligand->structure)]++;
                //delete a->possibleReceptors[i];
            }
            int sum = 0;
            for(int i = 0; i < (int) classesNbInteract.size(); ++i){
                if( classesNbInteract[i] > 0){
                    cout << "   " << i << " Interact, \t" << classesNbInteract[i] << endl;
                    sum += classesNbInteract[i];
                }
            }
            if(sum != NRs) cerr << "ERROR, size inconsistency " << endl;
         //}

    #ifdef ALLOW_GRAPHICS
    if(true){
        cerr << "Now displaying the found structures " << endl;
        char *c[] = {(char*)"Hello",nullptr};
        glDisplay(0,c);
        addToDisplay(merged);
        addToDisplay(ligand, true);
        size_t NRs = a->possibleReceptors.size();
        for(size_t i = 0; i < min(static_cast<size_t>(100000),NRs); ++i){
            if(!a->possibleReceptors[i]) {cerr << "WTF for recept " << i << endl; }
            else addToDisplay(a->possibleReceptors[i], false);
        }
        glutMainLoop();
    }
    #endif
    ///}

}

string addBeforeFileType(string fname, string toAdd){
    int L = fname.size();
    if(L == 0) return toAdd + string(".txt");
    for(int i = L-1; i >= 0 ; --i){
        //cout << fname[i];
        if(fname[i] == '.'){
            return fname.substr(0, i) + toAdd + fname.substr(i, L-i);
        }
    }
    return fname + toAdd + string(".txt");
}


int exportSelfInteractions(vector<struct3D*> &structures, string fname){

    ofstream res(fname.c_str());
    if(!res) cerr << "ERR: exportSelfInteractions, could not write into " << fname << endl;

    map<string, int> compressedReceptors; // only based on the interactions. Two receptors with the same interaction profile (inside & outside) can be merged

    int nRecept = structures.size();
    for(int i = 0; i < nRecept; ++i){
        string codeInters = codeSelfInteractions(*(structures[i]));
        // hope the structure is not empty ...
        res << structures[i]->startingPosition << "\t" << structures[i]->sequence << "\t" << codeInters << endl;
        if(compressedReceptors.find(codeInters) != compressedReceptors.end()){
            compressedReceptors[codeInters]++;
        } else {
            compressedReceptors[codeInters] = 1;
        }
    }
    res.close();

    string compressed = addBeforeFileType(fname, string("Compact"));
    ofstream res2 (compressed.c_str());
    if(!res2) cerr << "ERR: exportSelfInteractions, could not write into " << compressed << endl;
    std::map<string,int>::iterator it2;
    res2 << compressedReceptors.size() << endl;
    for(it2 = compressedReceptors.begin(); it2 != compressedReceptors.end(); ++it2){
        res2 << it2->first << "\t" << it2->second << endl;
    }
    res2.close();
    return compressedReceptors.size();
}

// not very important function,
// this will take care of the files behind
vector<struct3D*> loadOrGenerateSelfFoldings(int sizeReceptor, int minInteractForSelfFold, string folderToLookFor){
    vector<struct3D*> selfFoldingsForThisSize;

    if(minInteractForSelfFold < 1) cerr << "ERR: forbidden to call loadOrGenerateSelfFoldings with less than one interaction" << endl;
    stringstream fname0;
    fname0 << folderToLookFor << "/selfFoldingsL=" << sizeReceptor << "minI=" << minInteractForSelfFold << ".txt";
    string fname = fname0.str();
    ifstream f(fname.c_str());
    if(!f){
        cout << "The self-foldings for L=" << sizeReceptor << " and minI=" << minInteractForSelfFold << " are being generated " << endl;
        selfFoldingsForThisSize = generateSelfFoldings(sizeReceptor,minInteractForSelfFold);
        cout << " saving self-foldings in " << fname << endl;
        int newNS = exportSelfInteractions(selfFoldingsForThisSize, fname);
        cout << "saved " << newNS << " non-redundant structures " << endl;
    }
    else {
        cout << "loading self-foldings already pre-computed in " << fname << endl;
        int startPos;
        string structure;
        string interactions;
        while((f >> startPos)){
            f >> structure;
            f >> interactions;
            selfFoldingsForThisSize.push_back(new struct3D(structure, UnDefined, startPos)); // I don't do .reserve() because anyway this function is called only once. no need for optimizations here
        }
        f.close();
    }
    return selfFoldingsForThisSize;
}






void testFoldingStructures(){

    //struct3D* ligand = new struct3D(string("ULSRDSLRDUS"));
    int sizeReceptor = 10;
    int minimalNInteract = 5;

    vector<struct3D*> res = generateSelfFoldings(sizeReceptor,minimalNInteract);
    cout << "Performing testFoldingStructures with size receptor " << sizeReceptor << " and " << minimalNInteract << " Minimal number of interactions " << endl;

    cout << " Total Nb : " << res.size()  <<  endl;
     /*       stringstream fname; fname << ligand->sequence << sizeReceptor << "-" << minimalNInteract;
            cout << "After Reduc " << a->printToFile((fname.str() + string(".txt")).c_str(), (fname.str() + string("Compact.txt")).c_str()) << "\n";
            int NRs = a->possibleReceptors.size();
            vector<int> classesNbInteract(50, 0);
            for(int i = 0; i < NRs; ++i){
                classesNbInteract[nbTouchPoints(*(a->possibleReceptors[i]), *ligand)]++;
                //delete a->possibleReceptors[i];
            }
            int sum = 0;
            for(int i = 0; i < classesNbInteract.size(); ++i){
                if( classesNbInteract[i] > 0){
                    cout << "   " << i << " Interact, \t" << classesNbInteract[i] << endl;
                    sum += classesNbInteract[i];
                }
            }
            if(sum != NRs) cerr << "ERROR, size inconsistency " << endl;
         //}
    ///}*/

    #ifdef ALLOW_GRAPHICS
    if(true){
        cerr << "Now displaying the found structures " << endl;
        char *c[] = {(char*)"Hello",nullptr};
        glDisplay(0,c);
        //addToDisplay(merged);
        //addToDisplay(ligand, true);
        size_t NRs = res.size();
        for(size_t i = 0; i < min(static_cast<size_t>(100000),NRs); ++i){
            if(!res[i]) {cerr << "WTF for recept " << i << endl; }
            else addToDisplay(res[i], false);
        }
        glutMainLoop();
    }
    #endif
}




/// Different types of representation :
///     CompactStructure (unrooted in space)  ex: SSULDRS
///     SpaceStructure : a position, an absolute direction for first move
///             ex: (1,3,5) - USLUDRS   (the first sign replaces the first S, the first non S gives possibles rotations)
///         Compact to Space : position + 1'S becomes anything + 1'U becomes anything. 25 possible rotations
///         Space to Compact : replace first by S. The first non S becomes U.
///     Volume = set of occupied positions in the lattice.
///         Space to volume : follow the movement.
///         Volume to space (init point) : more tricky, because multiple choices. Don't do it here.
///
/// To fuse two CompactStructures- fuses the last segment with the first segment of the second protein,
/// (warning: removes 2 AA in the way, because 1 link is fused)
/// ex: SSULDRL + SSURD => the first S disappears (merged with the last L). Choses the intial direction
///    12345678  123456    the first U is changed by whatever direction
///                     => possible fusions : SSULDRL-SRRD
///                                           SSULDRL-SLRD
///                                           SSULDRL-SURD
///                                           SSULDRL-SDRD
///                                           SSULDRL-SSRD
///                                          12345678 9012
///
/// The ligand is (0,0,0) + CompactStructure (starts S(S)*U(?)* ).
///     + volume (list of 'surface points' decided inaccessible)
///     note: would also be possible to center on the receptor (AB), but I think it requires more computation (for each structure see if the ligand fit)
///           whereas here, start from ligand, and construct all possible receptor sequences
/// The receptor,
///     has a first point of interaction in the ligand (residue X)
///     has two directions of leaving this point
///         constraints : both go outside the protein (don't enter the core nor inaccessible space).
///         enumeration : either both directions make proteins that don't touch the ligand again
///                       or one direction touches again, not the other one (double amount of cases : touch - nontouch and nontouch-touch)
///
/// Ligand-receptor =
///     CombiStruct
///     (0,0,0) + Compact structure (ligand)
///     + volume (inaccessible space)
///     + first point(residue) of impact on the ligand
///     + two absolute structures  // compact structures won't work because multiple rotations are possible.
///
///
///     Other representations of L-R:
///        BiProtStruct:
///            (0,0,0) - Compact
///            (x,y,z) - Space Structure
///        Match
///            senseMatch: AAi(ligand) => AAj(receptor)
///            antisenseMatch: AAi(ligand) => AA(n-j) receptor
///
///
/// Note: for filling with AAs, each structures will need to be filled in both directions.



