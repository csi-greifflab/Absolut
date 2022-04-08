#include "compact.h"
#include "lattice.h"
#include "../Tools/zaprandom.h"
#include "plot3d.h"
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <set>
#include <map>
#include <algorithm>
#include <cstring>    // for strcpy
//#include <random>     // necessary in Linux but complains about C++11 - should not need it anymore
using namespace std;






// =============================================== Part A: Defining movements. ============================================

// Default initial y axis. Defined as if the observer was following x with the y axis.
moveDirection initialYaxis(moveDirection dir){
    switch(dir){
        case Right: return Straight;
        case Left: return Backwards;
        case Up: return Left;
        case Down: return Left;     //case Down: return moveID(vectorialProduct(moveVector(Straight), moveVector(dir)));
        case Straight: return Left; // default
        case UnDefined: return UnDefined; // should raise error later maybe
        case Backwards: return Right; // default
    }
    return UnDefined;
}




moveDirection randRelMove(){
    return static_cast<moveDirection>(random::uniformInteger(0,(int) Nb_Moves_Relative-1));
}
moveDirection randAbsMove(){
    return static_cast<moveDirection>(random::uniformInteger(0,(int) Nb_Moves_Absolute-1));
}


char intToMoveChar(int m){
    switch((moveDirection) m){
    case Right: return 'R';
    case Left: return 'L';
    case Up: return 'U';
    case Down: return 'D';
    case Straight: return 'S';
    case UnDefined: return '-';
    case Backwards: return 'B';
    }
    return '?';
}

char intToMoveChar(moveDirection m){
    switch(m){
    case Right: return 'R';
    case Left: return 'L';
    case Up: return 'U';
    case Down: return 'D';
    case Straight: return 'S';
    case UnDefined: return '-';
    case Backwards: return 'B';
    }
    return '?';
}

int charMoveToInt(char c){
    switch(c){
    case 'R': {return Right;}
    case 'L': {return Left;}
    case 'U': {return Up;}
    case 'D': {return Down;}
    case 'S': {return Straight;}
    case '-': {return UnDefined;}  // Actually, calling charMoveToInt('-') doesn't make much sense.
    case 'B': {return Backwards;}  // only for the first move of an absolute structure
    default:  {
        cerr << "ERR: Proteins Structures are supposed to only contain -,R,L,U,D and S, and not " << c << endl;
        return UnDefined;
    }
    }
}


vector<int> moveVector(moveDirection dir){
    vector<int> res(3,0);
    switch(dir){
        case Right:     {res[1] = -1; break;}
        case Left:      {res[1] = +1; break;}
        case Up:        {res[2] = +1; break;}
        case Down:      {res[2] = -1; break;}
        case Straight:  {res[0] = +1; break;}
        case Backwards: {res[0] = -1; break;}
        case UnDefined:  {/*cerr << "ERR: moveVector(UnDefined) is not an allowed direction" << endl This error always happens during init */ break;}
    }
    return res;
}

moveDirection moveID(vector<int> v){
    if(v.size() != 3) cerr << "ERR: moveID, vector with wrong length" << endl;
    if((abs(v[0]) + abs(v[1]) + abs(v[2])) != 1) cerr << "ERR: moveID, vector is not a unit direction " << v[0] << "," << v[1] << "," << v[2] << endl;
    if(v[1] == -1) return Right;
    if(v[1] == +1) return Left;
    if(v[2] == +1) return Up;
    if(v[2] == -1) return Down;
    if(v[0] == +1) return Straight;
    if(v[0] == -1) return Backwards;
    cerr << "ERR: moveID, move unknown : " << v[0] << "," << v[1] << "," << v[2] << endl;
    return UnDefined;
}

vector<int> vectorialProduct(vector<int> x, vector<int> y){
    if((x.size() != 3) || (y.size() != 3)) cerr << "ERR: Vectorial product with bad sized vectors " << endl;
    vector<int> z(3,0);
    z[0] = x[1] * y[2] - y[1] * x[2];
    z[1] = - x[0] * y[2] + y[0] * x[2];
    z[2] = x[0] * y[1] - y[0] * x[1];
    return z;
}

moveDirection vectorialProduct(moveDirection x, moveDirection y){
    // Pre-stores the result for all combinations of x and y into a map (static) when the function is called first.
    // Then, each further call is just an access to the map.
    static map<int, moveDirection> products;
    static bool loaded = false;
    if(!loaded){
        //cerr << "Loading vectorial product combinations" << endl;
        for(int i = 0; i < Nb_Moves_Absolute + 1; ++i){
            for(int j = 0; j < Nb_Moves_Absolute + 1; ++j){
                int comb = i + 100 * j;
                if((i == UnDefined) || (j == UnDefined)) products[comb] = UnDefined;
                else {
                    vector<int> z = vectorialProduct(moveVector(static_cast<moveDirection> (i)), moveVector(static_cast<moveDirection>(j)));
                    if(norm2(z) == 0) products[comb] = UnDefined;
                    else {
                        moveDirection mz = moveID(z);
                        products[comb] = mz;
                    }
                    //cerr << intToMoveChar(i) << " ^ " << intToMoveChar(j) << " = " << intToMoveChar(products[comb]) << endl;
                }
            }
        }
        loaded = true;
    }
    return products[static_cast<int>(x) + 100 * static_cast<int>(y)];
}

void testVectorsDirections(){
    // do everything in CERR such that errors occur at the good place
    cerr << "============= testing moveVector() adn moveID(). Expected errors (normal)!" << endl;
    for(int i = 0; i < Nb_Moves_Absolute+1; ++i){ // goes one step too much (+1) to reach 'Undefined' or raise an error
        cerr << "Direction : " << i << " (" << intToMoveChar(i) << ") " << printVector(moveVector(static_cast<moveDirection>(i))) << " IDmove" << moveID(moveVector(static_cast<moveDirection>(i))) << endl;
    }

    cerr << "============= testing vectorial product " << endl;
    for(int i = 0; i < 10; ++i){
        moveDirection i1 = randAbsMove();
        moveDirection i2 = randAbsMove();
        vector<int> res= vectorialProduct(moveVector(i1), moveVector(i2));
        cerr << i1 << " " << i2 << " => " << intToMoveChar(i1) << " ^ " << intToMoveChar(i2) << " = " << intToMoveChar(moveID(res)) << endl;
        cerr << "     as vectors : " << printVector(moveVector(i1)) << " ^ " << printVector(moveVector(i2)) << " = " << printVector(res) << endl;
    }

    cerr << "============= testing initialYaxis" << endl;
    for(int i = 0; i < Nb_Moves_Absolute+1; ++i){ // goes one step too much (+1) to reach 'Undefined' or raise an error
        cerr << "Direction : " << i << " (" << intToMoveChar(i) << ") " << initialYaxis(static_cast<moveDirection>(i)) << " (" << intToMoveChar(initialYaxis(static_cast<moveDirection>(i))) <<endl;
    }
}








// ------------------ tool functions necessary for the nextAbsoluteMove function (not shown in the .h): pre-computation of all possible (current X, current Y, next move) --------------------

// 1/ an ID is given to each combination (current observer direction (Ox), current observer y axis (Oy), and next relative move)
int combinedID(moveDirection previous, moveDirection yaxis, moveDirection nextToTranslate){
    return static_cast<int>(previous) * (Nb_Moves_Absolute + 1) * (Nb_Moves_Absolute + 1) +
            static_cast<int>(yaxis) * (Nb_Moves_Absolute + 1) +
            static_cast<int>(nextToTranslate);  // note: the +1's are in case someone puts UnDefined.
}

vector<int> decombine(int ID){
    vector<int> res(3,-1);
    res[0] = ID / ((Nb_Moves_Absolute + 1) * (Nb_Moves_Absolute + 1));
    res[1] = (ID - res[0] * ((Nb_Moves_Absolute + 1) * (Nb_Moves_Absolute + 1))) / (Nb_Moves_Absolute + 1);
    res[2] = ID % (Nb_Moves_Absolute + 1);
    return res;
}

//cerr << "============= testing combined ID" << endl;
//for(int i = 0; i < 10; ++i){
//    moveDirection i1 = randAbsMove();
//    moveDirection i2 = randAbsMove();
//    moveDirection i3 = randAbsMove();
//    int res = combinedID(i1, i2, i3);
//    cerr << " Combined Directions " << i1 << " " << intToMoveChar(i1) << ",";
//    cerr << " Combined Directions " << i2 << " " << intToMoveChar(i2) << ",";
//    cerr << " Combined Directions " << i3 << " " << intToMoveChar(i3) << " => " << res << " decombined:" << printVector(decombine(res)) << endl;
//}
//cerr << "============= testing movementMap" << endl;
//cerr << print(movementMap());

string showComb(int ID){
    stringstream res;
    vector<int> vec = decombine(ID);
    res << "(";
    res << intToMoveChar(vec[0]) << ",";
    res << intToMoveChar(vec[1]) << ") + ";
    res << intToMoveChar(vec[2]);
    return res.str();
}

int norm2(vector<int> toTest){
    return toTest[0] * toTest[0] + toTest[1] * toTest[1] + toTest[2] * toTest[2];
}

double norm2(vector<double> toTest){
    return toTest[0] * toTest[0] + toTest[1] * toTest[1] + toTest[2] * toTest[2];
}

int norm2(vector<int> v1, vector<int> v2){
    return norm2(vector<int>({v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]}));
}

double norm2(vector<double> v1, vector<double> v2){
    return norm2(vector<double>({v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]}));
}

// 2/ creates a 'movement map', where, for each combination, the next direction (Ox) and y axis (Oy) is pre-computed / stored
std::map<int, std::pair<moveDirection, moveDirection> > movementMap(){
    std::map<int, std::pair<moveDirection, moveDirection> > decision;

    // i is a current direction, including Backwards
    for(int i = 0; i <= Nb_Moves_Absolute; ++i){
        vector<int> currentX = moveVector( (moveDirection) i);
        //cout << " current " << intToMoveChar(i) << ":" << printVector(currentX) << endl;

        // j is a current Yaxis plane (including Backwards).
        for(int j = 0; j <= Nb_Moves_Absolute; ++j){
            vector<int> currentY = moveVector((moveDirection) j);
            vector<int> currentZ = vectorialProduct(currentX, currentY);


            // the yaxis vector has to be non-colinear !
            if(norm2(currentZ) != 0){

                //cout << "    ... current Y " << intToMoveChar(j) << ":" << printVector(currentY) << endl;
                //cout << "           => current Z " << intToMoveChar(moveID(currentZ)) << ":" << printVector(currentZ) << endl;

                // k is the next (relative) move.
                for(int k = 0; k < Nb_Moves_Relative; ++k){
                    vector<int> moveXYZ = moveVector( (moveDirection) k);
                    vector<int> moveInBase(3,0);
                    // the move 'X' means following currentX, 'Y' means follow current Y, etc
                    // new move, in XYZ basis
                    moveInBase[0] = moveXYZ[0] * currentX[0] + moveXYZ[1] * currentY[0] + moveXYZ[2] * currentZ[0];
                    moveInBase[1] = moveXYZ[0] * currentX[1] + moveXYZ[1] * currentY[1] + moveXYZ[2] * currentZ[1];
                    moveInBase[2] = moveXYZ[0] * currentX[2] + moveXYZ[1] * currentY[2] + moveXYZ[2] * currentZ[2];

                    // new Yaxis vector : if Straight, kept. if not, vectorial(prev, new)
                    vector<int> newY = currentY; // current Yaxix

                    // if left or right, keep same plane of walking, but needs to turn. => Zaxis vector kept, Y=N^newX
                    // this is a trick: instead of
                    if((k == Right) || (k == Left)){
                        newY = vectorialProduct(currentZ, moveInBase);
                    }

                    // if up or right, the plane is turned (like piloting a plane). Y is kept. The Zaxis vector is changed.
                    if((k == Up) || (k == Down)){
                        newY = currentY;
                    }

                    //cout << " Move " << intToMoveChar(k) << " results in : newDir " << intToMoveChar(moveID(moveInBase)) << " and new Yaxis: " << intToMoveChar(moveID(newY)) << endl;
                    decision[combinedID((moveDirection) i,(moveDirection) j,(moveDirection) k)] = std::pair<moveDirection, moveDirection> (moveID(moveInBase), (norm2(newY) == 0) ? UnDefined : moveID(newY) ); // moveID raises an error if gets a wrong vectror
                }
            }
        }
    }
    return decision;
}

string print( std::map<int, std::pair<moveDirection, moveDirection> > theMap){ // note: if would define the function with passing reference (&), then calling print(movementMap()) gets refused. Would need to create it outside.
    stringstream res;
    for(map<int, std::pair<moveDirection, moveDirection> >::iterator it=theMap.begin(); it!=theMap.end(); ++it) {
        res << it->first << showComb(it->first) << " => " << intToMoveChar(it->second.first) << "," << intToMoveChar(it->second.second) << endl;
    }
    return res.str();
}

// Tool function: movement map when browsing a structure from the terminal end and walking backwards.

// Interest: When fusing two proteins, one might want to fuse the last residue of the first structure
// to the beginning of next structure according to a specific rotation (of the first structure)
// in that case, the sequence to fuse needs to be transformed from its last residue to its first.
// it requires a 'reverse movement' : if the observer is (Ox, Oy), and the previous movement was M,
// what was the position of the observer before ?

std::map<int, std::pair<moveDirection, moveDirection> > reverseMovementMap(){
    std::map<int, std::pair<moveDirection, moveDirection> > decision;

    // i is a current direction, including Backwards
    for(int i = 0; i < Nb_Moves_Absolute; ++i){
        vector<int> currentX = moveVector( (moveDirection) i);
        //cout << " current " << intToMoveChar(i) << ":" << printVector(currentX) << endl;

        // j is a current Yaxis plane (including Backwards).
        for(int j = 0; j < Nb_Moves_Absolute; ++j){
            vector<int> currentY = moveVector((moveDirection) j);
            vector<int> currentZ = vectorialProduct(currentX, currentY);
            if(norm2(currentZ) != 0){

                // desired target move
                for(int k = 0; k < Nb_Moves_Absolute; ++k){

                    // now finds the relative move that makes it go to the desired direction.
                    for(int l = 0; l < Nb_Moves_Relative; ++l){
                        std::pair<moveDirection, moveDirection> possible = nextAbsoluteMove((moveDirection) i, (moveDirection) j, (moveDirection) l);
                        if(possible.first == (moveDirection) k) {
                            decision[combinedID((moveDirection) i,(moveDirection) j,(moveDirection) k)] =
                                    std::pair<moveDirection, moveDirection> ((moveDirection) l, possible.second);
                            cout << "Abs " << intToMoveChar(i) << " with Y=" << intToMoveChar(j) << " and wanted move " << intToMoveChar(k) << " needs the absolute movement " << intToMoveChar(l) << endl;
                        }
                    }
                    /*if(newRelativeMove  == Nb_Moves){
                        cerr << "ERR: fuse(" << seq1 << " -> absolutes " << s.seqAbsMoves << "," << seq2 << " -> absolutes " << s2.seqAbsMoves;
                        cerr << " at position " << i << "inside seq2 there is an opposite (absolute) direction than the last mobe. Forbidden to go backwards !" << endl;
                        return string();
                    }*/
                }
            }
        }
    }
    return decision;
}

std::pair<moveDirection, moveDirection> convertMoveWithNewYAxis(moveDirection currentDirection, moveDirection newYAxis, moveDirection desiredNewAbsDirection){
    static std::map<int, std::pair<moveDirection, moveDirection> > decision = reverseMovementMap();
    int IDmove = combinedID(currentDirection, newYAxis, desiredNewAbsDirection);
    std::map<int, std::pair<moveDirection, moveDirection> >::iterator it = (decision.find(IDmove));
    if (it == decision.end()){
        cerr << "ERR: convertMoveWithNewYAxis(" << intToMoveChar(currentDirection) << "," << intToMoveChar(newYAxis) << "," << intToMoveChar(desiredNewAbsDirection) << "), this movement is not authorized " << endl;
        return std::pair<moveDirection, moveDirection>(UnDefined, UnDefined); // todo: make a proper error enum
    }
    std::pair<moveDirection, moveDirection> res = decision[IDmove]; /// Phi to change:: do  *it.reference
    return res;
//    return *it.value??;
}

// ------------------- end of tool functions ---------------------



// Main function to move along an absolute sequence, in space.
// returns 1/ next absolute direction and 2/ next y axis
std::pair<moveDirection, moveDirection> nextAbsoluteMove(moveDirection previous, moveDirection yaxis, moveDirection nextToTranslate){
    static std::map<int, std::pair<moveDirection, moveDirection> > decision = movementMap();

    int IDmove = combinedID(previous, yaxis, nextToTranslate);
    std::map<int, std::pair<moveDirection, moveDirection> >::iterator it = (decision.find(IDmove));
    if (it == decision.end()){
        cerr << "ERR: nextAbsoluteMove(" << intToMoveChar(previous) << "," << intToMoveChar(yaxis) << "," << intToMoveChar(nextToTranslate) << "), this movement is not authorized " << endl;
        return std::pair<moveDirection, moveDirection>(UnDefined, UnDefined);
    }
    std::pair<moveDirection, moveDirection> res = decision[IDmove]; // Can be improved to *it.reference
    //cout << res.first << "," << res.second << endl;
    // Note: in non debug mode, this function could be shorter.
    return res;
}






// =========================== Part 2: Manipulating structures in 3D ===========================


struct3D::struct3D(string AbsoluteSequence, moveDirection initYAxis, int IDinitposition)
    : sequence(AbsoluteSequence), properlyFolded(true){

    //cerr << "S3D constructor" << AbsoluteSequence << " x=" << IDinitposition << " initY=" << intToMoveChar(initYAxis) << endl;

    //cout << "struct3D constructor " << AbsoluteSequence << endl;
    if(IDinitposition == -1) IDinitposition = lattice::centralPosition();
    startingPosition = IDinitposition;
    int cptIDresidue = 0;
    vector<int> currentPos = lattice::positionFromID(IDinitposition);
    if(!lattice::testPos(currentPos)) {
        cerr << "ERR: struct3D(" << sequence << ", pos=" << IDinitposition << ", protein going out of bounds at residue " << cptIDresidue  << ", pos=" << currentPos[0] << "," << currentPos[1] << "," << currentPos[2] << ". Suggestions : change starting point of the protein, or change X/Y/ZWidth" << endl;
        properlyFolded = false;
    } // leaving area is an exclusion criterion

    // even if the proteins is 1 residue (empty moves), puts the initial points in.
    //§§§§§points.push_back(residue(IDinitposition, cptIDresidue));
    occupiedPositions.insert(IDinitposition);
    if(sequence.size() == 0) return;
    cptIDresidue++;

    // First move. Note: here, take the initial axis. Would have been similar to start by straight and follow a decision, but in case of Backwards, it doesn't work.

    moveDirection currentDir = (moveDirection) charMoveToInt(sequence[0]);
    moveDirection currentYaxis = initialYaxis((moveDirection) charMoveToInt(sequence[0]));
    if(initYAxis != UnDefined){ // note:backwards also possible.
        currentYaxis = initYAxis;
        if(vectorialProduct(currentDir, currentYaxis) == UnDefined){
            cerr << "ERR: Struct3d(" << sequence << "), the initial direction : " << intToMoveChar(currentDir) << " is not compatible (not perpendicular) with the requested Y axis : " << intToMoveChar(initYAxis) << " : prod=" << intToMoveChar(vectorialProduct(currentDir, currentYaxis)) << endl;
        }
    }
    vector<int> actualMove = moveVector(currentDir);
    currentPos[0] += actualMove[0];
    currentPos[1] += actualMove[1];
    currentPos[2] += actualMove[2];
    // No need for folding checking yet.
    int IDcurrentPos = lattice::idFromPosisition(currentPos);
    //§§§§§points.push_back(residue(IDcurrentPos, cptIDresidue));
    listYAxis.push_back(intToMoveChar(currentYaxis));
    separatedSingleMoves.push_back(intToMoveChar(currentDir));
    occupiedPositions.insert(IDcurrentPos);
    cptIDresidue++;


    for(size_t i = 1; i < sequence.size(); ++i){
        std::pair<moveDirection, moveDirection> next = nextAbsoluteMove(currentDir, currentYaxis, (moveDirection) charMoveToInt(sequence[i]));
        //cout << "Current Pos " << printVector(currentPos) << "\tCurrentDir/Yaxis= " << intToMoveChar(currentDir) << "," << intToMoveChar(currentYaxis);
        //cout << " -> " << sequence[i] << " -> next direction/Yaxis = " << intToMoveChar(next.first) << "," << intToMoveChar(next.second);
        //if(next.first == UnDefined) cerr << "ERR: struct3D(), couldn't find the next move for " << intToMoveChar(currentDir) << ", " << intToMoveChar(currentYaxis) << endl;
        currentDir = next.first;
        currentYaxis = next.second;
        listYAxis.push_back(intToMoveChar(next.second));
        separatedSingleMoves.push_back(intToMoveChar(next.first));
        vector<int> actualMove = moveVector(currentDir);
        currentPos[0] += actualMove[0];
        currentPos[1] += actualMove[1];
        currentPos[2] += actualMove[2];
        //cout << "\tNewPos= " << printVector(currentPos) << endl;
        if(!lattice::testPos(currentPos)) {
            cerr << "ERR: struct3D(" << sequence << ", pos=" << IDinitposition << ", protein going out of bounds at residue " << cptIDresidue  << ", pos=" << currentPos[0] << "," << currentPos[1] << "," << currentPos[2] << ". Suggestions : change starting point of the protein, or change X/Y/ZWidth" << endl;
            properlyFolded = false;
        } // leaving area is an exclusion criterion
        IDcurrentPos = lattice::idFromPosisition(currentPos);
        //§§§§§points.push_back(residue(IDcurrentPos, cptIDresidue));
        properlyFolded = properlyFolded && ((occupiedPositions.find(IDcurrentPos)) == occupiedPositions.end());
        occupiedPositions.insert(IDcurrentPos);
        cptIDresidue++;
    }

    endingPosition = IDcurrentPos;
}



struct3D::struct3D(const struct3D &toCopy)
    // 2019-11-04 originally, this constructor would recreate the full structure, it allowed to check consistency. But this takes too much time so now just copying all the content.
    //: struct3D(toCopy.sequence, (toCopy.listYAxis.size() > 0) ? (moveDirection) charMoveToInt(toCopy.listYAxis[0]) : UnDefined, toCopy.startingPosition)
{
    startingPosition = toCopy.startingPosition;
    sequence = toCopy.sequence;
    separatedSingleMoves = toCopy.separatedSingleMoves;
    listYAxis = toCopy.listYAxis;
    endingPosition = toCopy.endingPosition;
    occupiedPositions = set<int>(toCopy.occupiedPositions);
    //cerr << &occupiedPositions << " VS " << &toCopy.occupiedPositions << endl;

    properlyFolded = toCopy.properlyFolded;
}

struct3D::struct3D(){
    // First, generates an empty structure. Note, the empty struct constructor does not exist because it would screw many functions
    // using struct3D, as they always expect a starting position, etc...
    startingPosition = 1e16;
    sequence = "";
    separatedSingleMoves = "";
    listYAxis = "";
    endingPosition = 1e16;
    occupiedPositions = set<int>();
    properlyFolded = true;
}

bool struct3D::pushBackAbsoluteMove(moveDirection d){
    moveDirection prevX = (moveDirection) charMoveToInt(separatedSingleMoves[sequence.size() - 1]);
    moveDirection prevY = (moveDirection) charMoveToInt(listYAxis[sequence.size() - 1]);
    moveDirection newRelativeMove = (moveDirection) UnDefined;
    for(int j = 0; j < Nb_Moves_Relative; ++j){
        std::pair<moveDirection, moveDirection> possible = nextAbsoluteMove(prevX,prevY, (moveDirection) j);
        if(possible.first == d){
            newRelativeMove = (moveDirection) j;
            separatedSingleMoves.push_back(intToMoveChar(possible.first));
            listYAxis.push_back(intToMoveChar(possible.second));
            sequence.push_back(intToMoveChar(newRelativeMove));
            vector<int> currentPos = lattice::positionFromID(endingPosition);
            vector<int> actualMove = moveVector(d);
            currentPos[0] += actualMove[0];
            currentPos[1] += actualMove[1];
            currentPos[2] += actualMove[2];
            if(!lattice::testPos(currentPos)) {cerr << "ERR: struct3D::pushBackAbsoluteMove(" << sequence << ", protein going out of bounds, pos=" << currentPos[0] << "," << currentPos[1] << "," << currentPos[2] << ". Suggestions : change starting point of the protein, or change X/Y/ZWidth" << endl;
                properlyFolded = false;} // leaving area is an exclusion criterion
            int IDcurrentPos = lattice::idFromPosisition(currentPos);
            //§§§§ points.push_back(residue(IDcurrentPos, points.back().IDresidue + 1));
            properlyFolded = properlyFolded && ((occupiedPositions.find(IDcurrentPos)) == occupiedPositions.end());
            //if(! properlyFolded) cerr << "WRN: you generated a self-colliding protein" << endl;
            occupiedPositions.insert(IDcurrentPos);
            endingPosition = IDcurrentPos;
            return true;
        }
    }
    //cerr << "ERR: couldn't add move" << endl;
    return false;
}


moveDirection lastAbsDirection(struct3D & s){
    return (moveDirection) charMoveToInt(s.separatedSingleMoves.back());
}

string print(set<int> &s){
    stringstream res;
    std::set<int>::iterator it;
    for(it = s.begin(); it != s.end(); ++it){
        res << *it << " ";
    }
    return res.str();
}

string print(set<string> &s){
    stringstream res;
    std::set<string>::iterator it;
    for(it = s.begin(); it != s.end(); ++it){
        res << *it << " ";
    }
    return res.str();
}


string print(struct3D& s){
    stringstream res;
    res << "struct3D(" << s.sequence << "), abs: " << s.separatedSingleMoves << " |-:" << s.listYAxis << ", " << ((s.properlyFolded) ? "Properly Folded" : "!! Auto-Colliding !!") << endl;
    res << "Positions: Start " << s.startingPosition << ", End " << s.endingPosition << ", Occupied positions : ";
    std::set<int>::iterator it;
    for(it = s.occupiedPositions.begin(); it != s.occupiedPositions.end(); ++it){
        res << *it << " ";
    }
    res << endl;
    return res.str();
}


bool contains(set<int> s, int v){
    return(!((s.find(v)) == s.end()));
}

bool collide(struct3D& s1, struct3D& s2){
    std::vector<int> common_points;
    set_intersection(s1.occupiedPositions.begin(),s1.occupiedPositions.end(),s2.occupiedPositions.begin(),s2.occupiedPositions.end(), std::back_inserter(common_points));
    return(common_points.size() > 0);
}

set<int> union_sets(set<int>& s1, set<int>& s2){
    set<int> union_points(s1);
    union_points.insert(s2.begin(), s2.end());

    //???set<int> union_points();
    //???set_union(s1.begin(),s1.end(),s2.begin(),s2.end(), std::set::insert(union_points));
    return union_points;
}

/*set<int> intersection_sets(set<int> &s1, set<int> &s2){
    // so there is no intersection sets function in C++, it takes VECTORS that need to be SORTED (hopefully, sets are)
    // see there: https://cukic.co/2018/06/03/set-intersection-in-cxx/

    vector<int> v1 = vector<int>(s1.begin(), s1.end());
    vector<int> v2 = vector<int>(s2.begin(), s2.end());

    // Precondition: make sure they're sorted
    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());

    std::vector<int> v_intersection;

    std::set_intersection(v1.begin(), v1.end(),
                          v2.begin(), v2.end(),
                          std::back_inserter(v_intersection));

    set<int> res = set<int>(v_intersection.begin(), v_intersection.end());
    return res;
}*/

//template <typename T>
//set<T> intersection_sets(set<T> &s1, set<T> &s2){
//    vector<T> v1 = vector<T>(s1.begin(), s1.end());
//    vector<T> v2 = vector<T>(s2.begin(), s2.end());

//    // Precondition: make sure they're sorted
//    std::sort(v1.begin(), v1.end());
//    std::sort(v2.begin(), v2.end());

//    std::vector<T> v_intersection;

//    std::set_intersection(v1.begin(), v1.end(),
//                          v2.begin(), v2.end(),
//                          std::back_inserter(v_intersection));

//    set<T> res = set<T>(v_intersection.begin(), v_intersection.end());
//    return res;
//}




vector<int> intersection(struct3D& s1, struct3D& s2){
    std::vector<int> common_points;
    set_intersection(s1.occupiedPositions.begin(),s1.occupiedPositions.end(),s2.occupiedPositions.begin(),s2.occupiedPositions.end(), std::back_inserter(common_points));
    return common_points;
}

// Does not say if yourself are touching, just say whether one neighbor touches.
// of course, if you are inside the structure, one of your neighbors will be inside the protein (except if it is an empty protein).
bool touch(struct3D& s1, int pos){
    // if(s1.sequence.size() == 0) return false;
    // if(s1.sequence.size() == 1) return (s1.startingPosition == pos);
    vector<int> nb = lattice::idNeighbors(pos);
    for(int i = 0; i < (int) nb.size(); ++i){
        if(contains(s1.occupiedPositions, nb[i])) return true;
    }
    return false;
}

int nbTouchPoints(struct3D& s1, struct3D & s2){
    if(s1.sequence.size() == 0) return 0;
    if(s2.sequence.size() == 0) return 0;
    int res = 0;
    std::set<int>::iterator it;
    for(it = s1.occupiedPositions.begin(); it != s1.occupiedPositions.end(); ++it){
        res += nbTouchPoints(s2, *it);
    }
    return res;
}

int nbTouchPoints(set<int>& occupPos, int pos){
    // if(s1.sequence.size() == 0) return false;
    // if(s1.sequence.size() == 1) return (s1.startingPosition == pos);
    vector<int> nb = lattice::idNeighbors(pos);
    int res = 0;
    for(int i = 0; i < (int) nb.size(); ++i){
        if(contains(occupPos, nb[i])) res++;
    }
    return res;
}

int nbTouchPoints(struct3D& s1, int pos){
    return nbTouchPoints(s1.occupiedPositions, pos);
}

//new way without needing 'points' :
set<int> neighborPositions(set<int> occupPos, set<int> toExclude){ //don't give to exclude as reference, will be modified inside
    // for each positions, add the neighbors in the list, and remove the positions of the protein itself.
    set<int> accuNeigh;
    for( set<int>::iterator sit = occupPos.begin(); sit != occupPos.end(); ++sit){
        vector<int> np = lattice::idNeighbors(*sit);
        // inserts it into the set
        std::copy( np.begin(), np.end(), std::inserter( accuNeigh, accuNeigh.end()));
    }
    // now removes the protein itself and excluded positions
    toExclude.insert(occupPos.begin(), occupPos.end());
    set<int> result;
    std::set_difference(accuNeigh.begin(), accuNeigh.end(), toExclude.begin(), toExclude.end(),
        std::inserter(result, result.end()));
    return result;
}
set<int> neighborPositions(struct3D& s1, set<int> toExclude){ //don't give to exclude as reference, will be modified inside
    return neighborPositions(s1.occupiedPositions, toExclude);
}

/* old way set<int> neighborPositions(struct3D& s1, set<int> toExclude){ //don't give to exclude as reference, will be modified inside
    // for each positions, add the neighbors in the list, but not those in the protein itself.
    int nbRes = s1.points.size();
    set<int> accuNeigh;
    for(int i = 0; i < nbRes; ++i){
        vector<int> np = lattice::idNeighbors(s1.points[i].IDposition);
        // inserts it into the set
        std::copy( np.begin(), np.end(), std::inserter( accuNeigh, accuNeigh.end()));
    }
    // now removes the protein itself and excluded positions
    toExclude.insert(s1.occupiedPositions.begin(), s1.occupiedPositions.end());
    set<int> result;
    std::set_difference(accuNeigh.begin(), accuNeigh.end(), toExclude.begin(), toExclude.end(),
        std::inserter(result, result.end()));
    return result;
} */



// takes whatever structure, and finds the group representative (by rotations) in the relative form
string normalizeAbsolute(string seq){
    // 1/ find which rotation to start by S
    if(seq.size() == 0) return seq;
    moveDirection axisToRotate = UnDefined;
    int nbTurns; // 0, 1 or 2
    switch(seq[0]){
        case 'S': { nbTurns = 0; break;}
        case 'R': { axisToRotate = Up; nbTurns = 1; break;}
        case 'L': { axisToRotate = Down; nbTurns = 1; break;}
        case 'U': { axisToRotate = Left; nbTurns = 1; break;}
        case 'D': { axisToRotate = Right; nbTurns = 1; break;}
        case 'B': { axisToRotate = Up; nbTurns = 2; break;} // any axis would work.
        default: { cerr << "ERR: normalizeAbsolute, unauthoried character " << seq[0] << endl; return string("");}
    }
    //cout << "Turning with axis" << intToMoveChar(axisToRotate) << endl;
    string seq2 = seq;
    if(nbTurns > 0) seq2 = easyRotate(seq, axisToRotate);
    //cout << "newSeq" << seq2 << endl;
    if(nbTurns > 1) seq2 = easyRotate(seq2, axisToRotate);
    //cout << "newSeq2" << seq2 << endl;

    // 2/ find which rotation to first move Up.
    int L = seq2.size();
    char firstTurn = '-';
    for(int i = 1; i < L; ++i){ // could have put a while;..
        if(seq2[i] != 'S') {
            firstTurn = seq2[i];
            break;
        }
    }
    //cout << "First turn is " << firstTurn << endl;

    // now just need to turn around the 'S' axis.
    switch(firstTurn){
        case '-': { return seq2;} // there is no turn, only SSSS....
        case 'S': { cerr << "Y'a une couille dans le potage !" << endl; return seq2;} // should not happen, but anyways,
        case 'U': { return seq2;}
        case 'R': { axisToRotate = Backwards; nbTurns = 1; break;}
        case 'L': { axisToRotate = Straight; nbTurns = 1; break;}
        case 'D': { axisToRotate = Straight; nbTurns = 2; break;}
        case 'B': { cerr << "ERR: normalizeAbsolute, Y'a une couille dans le potage !" << endl; return string("");} // any axis would work.
        default: { cerr << "ERR: normalizeAbsolute, unauthoried character " << seq[0] << endl; return string("");}
    }

    string seq3 = seq2;
    if(nbTurns > 0) seq3 = easyRotate(seq2, axisToRotate);
    if(nbTurns > 1) seq3 = easyRotate(seq3, axisToRotate);
    //cout << "After second turn" << seq3 << endl;

    return seq3;
}








//// Gives the sequence for the same 3D structure, but starting from an opposite Y vector. => not sure it works.
string opposite(string seq){
    int s = seq.size();
    string res = seq;
    for(int i = 1; i < s; ++i){
        switch(seq[i]){
            case 'U': res[i] = 'D'; break;
            case 'D': res[i] = 'U'; break;
            case 'R': res[i] = 'L'; break;
            case 'L': res[i] = 'R'; break;
            case 'S': res[i] = 'S'; break;
            default: cerr << "ERR: opposite(" << seq << "), B is forbidden except at position 0" << endl;
        }
    }
    return res;
}

char opp(char c){
    switch(c){
    case 'U': return 'D'; break;
    case 'D': return 'U'; break;
    case 'L': return 'R'; break;
    case 'R': return 'L'; break;
    case 'S': return 'B'; break;
    case 'B': return 'S'; break;
    default: cerr << "ERR: revert, unknown move " << c << endl;
    }
    return '#';
}

// Tool function that does the dirty job: takes independent steps in space, rotates them and reconstructs it as a list of moves
// Input : a list of INDEPENDENT absolute moves
// Output: the absolute sequence + list of absolute moves after rotation
std::pair<string,string> absoluteRotate(string seq, moveDirection absAxis, bool clockwise){
    // there are only 6 rotations possible = possible axis and clockwise.
    // unclockwise rotations are transformed into clockwise when the axis is opposite.
    if(absAxis == UnDefined){
        cerr << "ERR: absoluteRotate, Undefined is not an accepted rotation axis" << endl;
        return std::pair<string, string>(string(""), seq);
    }
    moveDirection newAxis = absAxis;
    if(!clockwise){
       switch(absAxis){ // note that these are absolute directions
           case Straight: {newAxis = Backwards; break;}
           case Right: {newAxis = Left; break;}
           case Left: {newAxis = Right; break;}
           case Up: {newAxis = Down; break;}
           case Down: {newAxis = Up; break;}
           case Backwards: {newAxis = Straight; break;}
           case UnDefined: {} // should not happen.
       }
    }

    // rotation straight
    bool loaded = false;
    static std::map<std::pair<moveDirection, char>, char> rotations;
    if(!loaded){
        rotations[std::pair<moveDirection, char> (Straight, 'S')] = 'S';
        rotations[std::pair<moveDirection, char> (Straight, 'R')] = 'D';
        rotations[std::pair<moveDirection, char> (Straight, 'L')] = 'U';
        rotations[std::pair<moveDirection, char> (Straight, 'U')] = 'R';
        rotations[std::pair<moveDirection, char> (Straight, 'D')] = 'L';
        rotations[std::pair<moveDirection, char> (Straight, 'B')] = 'B';
        rotations[std::pair<moveDirection, char> (Straight, '-')] = '-';
        rotations[std::pair<moveDirection, char> (Right, 'S')] = 'U';
        rotations[std::pair<moveDirection, char> (Right, 'R')] = 'R';
        rotations[std::pair<moveDirection, char> (Right, 'L')] = 'L';
        rotations[std::pair<moveDirection, char> (Right, 'U')] = 'B';
        rotations[std::pair<moveDirection, char> (Right, 'D')] = 'S';
        rotations[std::pair<moveDirection, char> (Right, 'B')] = 'D';
        rotations[std::pair<moveDirection, char> (Right, '-')] = '-';
        rotations[std::pair<moveDirection, char> (Up, 'S')] = 'L';
        rotations[std::pair<moveDirection, char> (Up, 'R')] = 'S';
        rotations[std::pair<moveDirection, char> (Up, 'L')] = 'B';
        rotations[std::pair<moveDirection, char> (Up, 'U')] = 'U';
        rotations[std::pair<moveDirection, char> (Up, 'D')] = 'D';
        rotations[std::pair<moveDirection, char> (Up, 'B')] = 'R';
        rotations[std::pair<moveDirection, char> (Up, '-')] = '-';
        rotations[std::pair<moveDirection, char> (Backwards, 'S')] = 'S';
        rotations[std::pair<moveDirection, char> (Backwards, 'R')] = 'U';
        rotations[std::pair<moveDirection, char> (Backwards, 'L')] = 'D';
        rotations[std::pair<moveDirection, char> (Backwards, 'U')] = 'L';
        rotations[std::pair<moveDirection, char> (Backwards, 'D')] = 'R';
        rotations[std::pair<moveDirection, char> (Backwards, 'B')] = 'B';
        rotations[std::pair<moveDirection, char> (Backwards, '-')] = '-';
        rotations[std::pair<moveDirection, char> (Left, 'S')] = 'D';
        rotations[std::pair<moveDirection, char> (Left, 'R')] = 'R';
        rotations[std::pair<moveDirection, char> (Left, 'L')] = 'L';
        rotations[std::pair<moveDirection, char> (Left, 'U')] = 'S';
        rotations[std::pair<moveDirection, char> (Left, 'D')] = 'B';
        rotations[std::pair<moveDirection, char> (Left, 'B')] = 'U';
        rotations[std::pair<moveDirection, char> (Left, '-')] = '-';
        rotations[std::pair<moveDirection, char> (Down, 'S')] = 'R';
        rotations[std::pair<moveDirection, char> (Down, 'R')] = 'B';
        rotations[std::pair<moveDirection, char> (Down, 'L')] = 'S';
        rotations[std::pair<moveDirection, char> (Down, 'U')] = 'U';
        rotations[std::pair<moveDirection, char> (Down, 'D')] = 'D';
        rotations[std::pair<moveDirection, char> (Down, 'B')] = 'L';
        rotations[std::pair<moveDirection, char> (Down, '-')] = '-';
        loaded = true;
    }

    string res = seq;
    int L = seq.size();
    for(int i = 0; i < L; ++i){
        res[i] = rotations[std::pair<moveDirection, char>(newAxis, seq[i])];
    }
    // this is the new sequence of absolute moves, but now needs to translate into a structure

    string rebuilt(res);
    moveDirection lastYAxis = initialYaxis((moveDirection) charMoveToInt(res[0]));
    moveDirection lastAbsDirection = (moveDirection) charMoveToInt(res[0]);

    int LS2 = res.size();
    for(int i = 1; i < LS2; ++i){
        moveDirection newRelativeMove = UnDefined;
        for(int j = 0; j < Nb_Moves_Relative; ++j){
            std::pair<moveDirection, moveDirection> possible = nextAbsoluteMove(lastAbsDirection,lastYAxis, (moveDirection) j);
            if(possible.first == (moveDirection) charMoveToInt(res[i])){
                newRelativeMove = (moveDirection) j;
                lastAbsDirection = (moveDirection) charMoveToInt(res[i]);
                lastYAxis = possible.second;
            }
        }
        if(newRelativeMove  == UnDefined){
            cerr << "ERR: absoluteRotate(seq=" << seq << ", absolute moves after rot. " << res << " -> being reconstructed " << rebuilt << endl;
            cerr << " at position " << i << "inside rot there is an opposite (absolute) direction than the last move. Forbidden to go backwards !" << endl;
            return std::pair<string, string> (string(), string());
        }
        rebuilt[i] = intToMoveChar(newRelativeMove);
    }
    return std::pair<string, string>(rebuilt, res);
}






// ACHTUNG!!! this function only rotates a sequences of independent, absolute directions (not moves !!)
// call it twice to rotate more than 90°. Will be as fast.
// don't need to give the starting position because it will not change, or need to translate according
// to the fix point later.

string easyRotate(string sequence, moveDirection absAxis, bool clockwise){
    // does the same as calling the struct3D constructor, but just decomposes into absolute separate moves
    string seqAbsMoves;
    if(sequence.size() == 0) return seqAbsMoves;
    moveDirection currentDir = (moveDirection) charMoveToInt(sequence[0]);
    moveDirection currentYaxis = initialYaxis((moveDirection) charMoveToInt(sequence[0]));
    seqAbsMoves.push_back(intToMoveChar(currentDir));
    for(int i = 1; i < (int) sequence.size(); ++i){
        std::pair<moveDirection, moveDirection> next = nextAbsoluteMove(currentDir, currentYaxis, (moveDirection) charMoveToInt(sequence[i]));
        currentDir = next.first;
        currentYaxis = next.second;
        seqAbsMoves.push_back(intToMoveChar(next.first));
    }
    return absoluteRotate(seqAbsMoves, absAxis, clockwise).first;
}

void sanityCheckEasyRotate(){
    string S = "SUSSDUSSDSSDLLUSSLSDDSLSSSDRLSUDDLURDSLDLUURRLUSR";

    // original way
    struct3D s3D = struct3D(S);
    moveDirection absAxis = Up;
    bool clockwise = true;
    string seqAbs = absoluteRotate(s3D.separatedSingleMoves, absAxis, clockwise).first;

    string seqAbsNew = easyRotate(S, absAxis, clockwise);
    if(seqAbs.compare(seqAbsNew) != 0) cerr << "ERR: sanityCheckEasyRotate failed" << endl;
    else cerr << "sanityCheckEasyRotate succeeded" << endl;
}


// can only rotate around the first move axis (can be relative or absolute)
string relativeRotate(string seq, bool clockwise){
    if(seq.size() == 0) return seq;
    return easyRotate(seq, (moveDirection) charMoveToInt(seq[0]), clockwise);






    /*
    struct3D s(seq);
    string target = absoluteRotate(s.seqAbsMoves, (moveDirection) charMoveToInt(seq[0]), clockwise);

    string rebuilt(seq);
    moveDirection lastYAxis = initialYaxis((moveDirection) charMoveToInt(target[0]));
    moveDirection lastAbsDirection = (moveDirection) charMoveToInt(target[0]);

    int LS2 = target.size();
    for(int i = 1; i < LS2; ++i){
        moveDirection newRelativeMove = UnDefined;
        for(int j = 0; j < Nb_Moves_Relative; ++j){
            std::pair<moveDirection, moveDirection> possible = nextAbsoluteMove(lastAbsDirection,lastYAxis, (moveDirection) j);
            if(possible.first == (moveDirection) charMoveToInt(target[i])){
                newRelativeMove = (moveDirection) j;
                lastAbsDirection = (moveDirection) charMoveToInt(target[i]);
                lastYAxis = possible.second;
            }
        }
        if(newRelativeMove  == UnDefined){
            cerr << "ERR: relativeRotate(" << seq << " -> absolutes " << s.seqAbsMoves;
            cerr << " at position " << i << "inside seq there is an opposite (absolute) direction than the last mobe. Forbidden to go backwards !" << endl;
            return string();
        }
        rebuilt[i] = intToMoveChar(newRelativeMove);
    }
    return rebuilt;*/
}





vector<string> allRotationsStruct(string base, moveDirection ForcingStart, moveDirection firstTurn);

// note : cuts all the structures.
// Does it actually work for absolute structures ??
// it seems that absolute and relative are the same: XSSSSSX[...], either X1=S and X2=U or just free.
// No, this is not true !!! §§§§§ change that
vector<AbsoluteStructure> allRotationsStruct(CompactStructure base, moveDirection ForcingStart, moveDirection firstTurn){
    if(!checkSyntaxSequence(intToRelative(base))) {
        cerr << "ERR: rotateStruct(" << base << ", ...) incorrect sequence\n";
        return vector<AbsoluteStructure>();
    }
    vector<string> res1 = allRotationsStruct(intToRelative(base), ForcingStart, firstTurn);
    vector<AbsoluteStructure> res2;
    int nr = res1.size();
    for(int i = 0; i < nr; ++i){
        res2.push_back(absoluteToInt(res1[i]));
    }
    return res2;
}

// I don't think this function works if forcingstart is different than the first move
// Seems this function is made for relative sequences. Does it work for absolute ?
// BACKWARDS IS INCLUDED. Output is absolute !
vector<string> allRotationsStruct(string base, moveDirection ForcingStart, moveDirection firstTurn){
    string S = cut(base);
    if(S.size() == 0) return vector<string>();

    //cout << "rotate, got " << S << endl;
    /// Enumerates all possible starting directions
    ///  => replaces the first move by an absolute direction (except if forced)
    vector<string> orientated;
    if(ForcingStart == Nb_Moves_Relative){
        for(int i = 0; i <= (int) Nb_Moves_Relative; ++i){
            string Scopy = S;
            Scopy[0] = intToMoveChar((i == (int) Nb_Moves_Relative) ? (int) Backwards : (moveDirection) i);
            //cout << "Push Init:" << Scopy << endl;
            orientated.push_back(Scopy);
        }
    } else {
        string Scopy = S;
        Scopy[0] = intToMoveChar(ForcingStart);
        //cout << "Push Init:" << Scopy << endl;
        orientated.push_back(Scopy);
    }

    /// From all oriented sequences, makes all possible rotations (the first non-starting S is replaced by anything. Note : it should be a U before replacing)
    vector<string> res;
    int nbOr = orientated.size();
    for(int i = 0; i < nbOr; ++i){
        string Scopy = orientated[i];
        // note:
        unsigned int k = 1; // starts at 1, to avoid the first letter
        while ((k < Scopy.size()) && (Scopy[k] == 'S')) {++k;} // finds the first non S
        if(k == Scopy.size()) {
            res.push_back(Scopy); // the sequence is SSSSSS, no change.
        } else {
            if(firstTurn == Nb_Moves_Relative){
                for(int j = 0; j < (int) Nb_Moves_Relative; ++j){
                    Scopy[k] = intToMoveChar(j); // Back not permitted (never permitted, moves are still relative to each-other).
                    //cout << "Push " << Scopy << " -> " << absoluteToInt(Scopy) << " -> " << intToAbsolute(absoluteToInt(Scopy)) << endl;
                    res.push_back(Scopy);
                }
            } else {
                Scopy[k] = firstTurn;
                res.push_back(Scopy);
            }
        }
    }
    return res;
}






void testRotate(){
    cout << "Testing rotation " << endl;
    vector<string> res;
    for(int i = 0; i < 10; ++i){
        if(i == 0) res = allRotationsStruct(string("SSURLDSUUSR"));
        if(i == 1) res = allRotationsStruct(string("SSURLDSUUSR"), Up);
        if(i == 2) res = allRotationsStruct(string("SSURLDSUUSR"), Down);
        if(i == 3) res = allRotationsStruct(string("SURLDSUUSR"), Down);
        if(i == 4) res = allRotationsStruct(string("SSURLDSUUSR"), Straight);
        if(i == 5) res = allRotationsStruct(string("SSSSSSS"));
        if(i == 6) res = allRotationsStruct(string("SSSSSSS"), Up);
        if(i == 7) res = allRotationsStruct(string("SSSSSSSSSSU"));
        for(unsigned int j = 0; j < res.size(); ++j){
            cout << res[j] << endl;
        }
        cout << " ===================================== " << endl;
    }
}












// functions from used from compact.h
// enum moveDirection{Straight, Right, Left, Up, Down, Backwards};
// unsigned long long int stringToInt(string prot);
// string intToString(unsigned long long protID);
// unsigned long long absoluteToInt(string prot);
// string intToAbsolute(unsigned long long protID);
// char intToMoveChar(int m);
// int charToInt(char c);

typedef long long unsigned int CompactStructure; // starts S(S*)U...
typedef long long unsigned int AbsoluteStructure; // whatever sequence (URLDS)

//struct SpaceStructure {
//    int gridStartPosition;
//    CompactStructure AbsoluteStructure;
//};

string printCompact(CompactStructure cs){
    return intToRelative(cs);
}
string printAbs(AbsoluteStructure as){
    return intToAbsolute(as);
}
string print(SpaceStructure &ss){
    stringstream res;
    vector<int> initPos = lattice::positionFromID(ss.gridStartPosition);
    res << "SpaceStruct, start " << ss.gridStartPosition << " (" << printVector(initPos) << ") " << printAbs(ss.AbsoluteStructure) << endl;
    return res.str();
}

CompactStructure compacte(AbsoluteStructure protID){
    string S = intToAbsolute(protID);
    //cout << "compacte got " << protID << " -> " << S << endl;
    //int originL = S.size();
    S = cut(S);
    int lg = S.size();
    S[0] = 'S';
    int k = 1;
    while((k < lg) && (S[k] == 'S')) {++k;}
    if(k < lg) S[k] = 'U';
    //S = string(originL - lg, '-') + S;
    //cout << fill(S) << endl;
    return relativeToInt(fill(S));
}



void testRotations(){
#ifdef ALLOW_GRAPHICS
    string S1 = "SLLURDULSD";
    string S2 = "BULSLLURDULSD";
    string S3 = "LSSUDLRL";
    string S4 = "SSSLRLRL";

    char *c[] = {(char*)"Hello",nullptr};
    glDisplay(0,c);

    for(int i = 1; i <= 4; ++i){
        string usedS;
        switch(i){
        case 1: usedS = S1; break;
        case 2: usedS = S2; break;
        case 3: usedS = S3; break;
        case 4: usedS = S4; break;
        }
        struct3D* emptySS = new struct3D("");
        addToDisplay(emptySS,false);
        struct3D* SS = new struct3D(usedS);
        addToDisplay(SS,false);
        /*for(int j = 0; j < (int) UnDefined; ++j){
            cout << "(absolute) Rotating sequence : " << usedS << " around axis " << intToMoveChar(j) << endl;
            addToDisplay(emptySS,false);
            addToDisplay(SS,false);
            string seqNew = absoluteRotate(SS->seqAbsMoves, (moveDirection) j, true).first;
            struct3D* SSN1 = new struct3D(seqNew); addToDisplay(SSN1, false);
            cout << "   -> axis " << intToMoveChar(j) << ",    clockwise, res=" << seqNew << endl;
            string seqNew2 = absoluteRotate(SSN1->seqAbsMoves, (moveDirection)j, true).first;
            struct3D* SSN2 = new struct3D(seqNew2); addToDisplay(SSN2, false);
            cout << "   -> axis " << intToMoveChar(j) << ", 2X clockwise, res=" << seqNew2 << endl;
            string seqNew3 = absoluteRotate(SSN2->seqAbsMoves, (moveDirection)j, true).first;
            struct3D* SSN3 = new struct3D(seqNew3); addToDisplay(SSN3, false);
            cout << "   -> axis " << intToMoveChar(j) << ", 3X clockwise, res=" << seqNew3 << endl;
            string seqNew4 = absoluteRotate(SS->seqAbsMoves, (moveDirection) j, false).first;
            struct3D* SSN4 = new struct3D(seqNew4); addToDisplay(SSN4, false);
            cout << "   -> axis " << intToMoveChar(j) << ", un-clockwise, res=" << seqNew4 << endl;
        }*/
        /*
        for(int j = 0; j < (int) UnDefined; ++j){
            cout << "(absolute) Rotating sequence : " << usedS << " around axis " << intToMoveChar(j) << endl;
            addToDisplay(emptySS,false);
            addToDisplay(SS,false);
            string seqNew = easyRotate(usedS, (moveDirection) j, true);
            struct3D* SSN1 = new struct3D(seqNew); addToDisplay(SSN1, false);
            cout << "   -> axis " << intToMoveChar(j) << ",    clockwise, res=" << seqNew << endl;
            string seqNew2 = easyRotate(seqNew, (moveDirection)j, true);
            struct3D* SSN2 = new struct3D(seqNew2); addToDisplay(SSN2, false);
            cout << "   -> axis " << intToMoveChar(j) << ", 2X clockwise, res=" << seqNew2 << endl;
            string seqNew3 = easyRotate(seqNew2, (moveDirection)j, true);
            struct3D* SSN3 = new struct3D(seqNew3); addToDisplay(SSN3, false);
            cout << "   -> axis " << intToMoveChar(j) << ", 3X clockwise, res=" << seqNew3 << endl;
            string seqNew4 = easyRotate(usedS, (moveDirection) j, false);
            struct3D* SSN4 = new struct3D(seqNew4); addToDisplay(SSN4, false);
            cout << "   -> axis " << intToMoveChar(j) << ", un-clockwise, res=" << seqNew4 << endl;
        }
        */
        /*
        cout << "(relative) Rotating sequence : " << usedS << endl;
        {
            addToDisplay(emptySS,false);
            addToDisplay(SS,false);
            string seqNew = relativeRotate(usedS, true);
            struct3D* SSN1 = new struct3D(seqNew); addToDisplay(SSN1, false);
            cout << "   ->    clockwise, res=" << seqNew << endl;
            string seqNew2 = relativeRotate(seqNew, true);
            struct3D* SSN2 = new struct3D(seqNew2); addToDisplay(SSN2, false);
            cout << "   -> 2X clockwise, res=" << seqNew2 << endl;
            string seqNew3 = relativeRotate(seqNew2, true);
            struct3D* SSN3 = new struct3D(seqNew3); addToDisplay(SSN3, false);
            cout << "   -> 3X clockwise, res=" << seqNew3 << endl;
            string seqNew4 = relativeRotate(usedS, false);
            struct3D* SSN4 = new struct3D(seqNew4); addToDisplay(SSN4, false);
            cout << "   -> un-clockwise, res=" << seqNew4 << endl;
        }*/

        cout << "Now normalize the sequence to a relative sequence starting with S(S*)U" << endl;
        {
            string seqNew = normalizeAbsolute(usedS);
            struct3D* SSN1 = new struct3D(seqNew); addToDisplay(SSN1, false);
            cout << "   ->    initial seq    =" << usedS << endl;
            cout << "      ...normalized, res=" << seqNew << endl;
        }
    }
    glutMainLoop();

    //   string normalizeAbsolute(string seq);


 //   string absoluteRotate(string seq, moveDirection absAxis, bool clockwise = true);
 //   string relativeRotate(string seq, bool clockwise = true);
 //   string normalizeAbsolute(string seq);
#endif
}


// describes the same structure from the other side, with the *exact same* configuration in space
// i.e. the absolute moves are the same but reverted and in opposite order.
// this function can not rotate nor change the absolute configuration.
// better to keep it simple because it is used heavily inside bigReceptors ...
string revert(string seq1){
    if(seq1.size() == 0) return seq1;
    struct3D s(seq1);
    string rev = s.separatedSingleMoves;
    std::reverse(rev.begin(), rev.end());
    string rebuilt(seq1);
    for(unsigned int i = 0; i < rev.size() ; ++i){
        rev[i] = opp(rev[i]);
    }
    rebuilt[0] = rev[0];
    moveDirection lastYAxis = initialYaxis((moveDirection) charMoveToInt(rev[0]));
    moveDirection lastAbsDirection = (moveDirection) charMoveToInt(rev[0]);

    int LS2 = seq1.size();
    for(int i = 1; i < LS2; ++i){
        moveDirection newRelativeMove = UnDefined;
        for(int j = 0; j < Nb_Moves_Relative; ++j){
            std::pair<moveDirection, moveDirection> possible = nextAbsoluteMove(lastAbsDirection,lastYAxis, (moveDirection) j);
            if(possible.first == (moveDirection) charMoveToInt(rev[i])){
                newRelativeMove = (moveDirection) j;
                lastAbsDirection = (moveDirection) charMoveToInt(rev[i]);
                lastYAxis = possible.second;
            }
        }
        if(newRelativeMove  == UnDefined){
            cerr << "ERR: revert(" << seq1 << " -> absolutes " << s.separatedSingleMoves;
            cerr << " at position " << i << "inside seq there is an opposite (absolute) direction than the last mobe. Forbidden to go backwards !" << endl;
            return string();
        }
        rebuilt[i] = intToMoveChar(newRelativeMove);
    }
    return rebuilt;
}

/* Big error with this function string revert(string seq1, moveDirection forceFirstAbs, moveDirection forceFirstyAxis){
    if(seq1.size() == 0) return seq1;
    struct3D s(seq1);
    string rev = s.seqAbsMoves;
    std::reverse(rev.begin(), rev.end());
    string rebuilt(seq1);
    for(unsigned int i = 0; i < rev.size() ; ++i){
        rev[i] = opp(rev[i]);
    }
    rebuilt[0] = rev[0];
    moveDirection lastYAxis = initialYaxis((moveDirection) charMoveToInt(rev[0]));
    if(lastYAxis != UnDefined) lastYAxis = forceFirstyAxis;
    moveDirection lastAbsDirection = (moveDirection) charMoveToInt(rev[0]);
    if(lastAbsDirection != UnDefined) lastAbsDirection = forceFirstAbs;

    int LS2 = seq1.size();
    for(int i = 1; i < LS2; ++i){
        moveDirection newRelativeMove = UnDefined;
        for(int j = 0; j < Nb_Moves_Relative; ++j){
            std::pair<moveDirection, moveDirection> possible = nextAbsoluteMove(lastAbsDirection,lastYAxis, (moveDirection) j);
            if(possible.first == (moveDirection) charMoveToInt(rev[i])){
                newRelativeMove = (moveDirection) j;
                lastAbsDirection = (moveDirection) charMoveToInt(rev[i]);
                lastYAxis = possible.second;
            }
        }
        if(newRelativeMove  == UnDefined){
            cerr << "ERR: revert(" << seq1 << " -> absolutes " << s.seqAbsMoves;
            cerr << " at position " << i << "inside seq there is an opposite (absolute) direction than the last mobe. Forbidden to go backwards !" << endl;
            return string();
        }
        rebuilt[i] = intToMoveChar(newRelativeMove);
    }
    return rebuilt;
}*/

string fuse(string seq1, string seq2){
    if(seq1.size() == 0) return seq2;
    if(seq2.size() == 0) return seq1;
    int LS2 = seq2.size();
    struct3D s(seq1);
    struct3D s2(seq2);
    string rebuilt(seq2);
    moveDirection lastYAxis = (moveDirection) charMoveToInt(s.listYAxis[s.sequence.size() - 1]);
    moveDirection lastAbsDirection = (moveDirection) charMoveToInt(s.separatedSingleMoves[s.sequence.size() - 1]);
    for(int i = 0; i < LS2; ++i){
        moveDirection newRelativeMove = UnDefined;
        for(int j = 0; j < Nb_Moves_Relative; ++j){
            std::pair<moveDirection, moveDirection> possible = nextAbsoluteMove(lastAbsDirection,lastYAxis, (moveDirection) j);
            if(possible.first == (moveDirection) charMoveToInt(s2.separatedSingleMoves[i])){
                newRelativeMove = (moveDirection) j;
                lastAbsDirection = (moveDirection) charMoveToInt(s2.separatedSingleMoves[i]);
                lastYAxis = possible.second;
            }
        }
        if(newRelativeMove  == UnDefined){
            cerr << "ERR: fuse(" << seq1 << " -> absolutes " << s.separatedSingleMoves << "," << seq2 << " -> absolutes " << s2.separatedSingleMoves;
            cerr << " at position " << i << "inside seq2 there is an opposite (absolute) direction than the last mobe. Forbidden to go backwards !" << endl;
            return string();
        }
        rebuilt[i] = intToMoveChar(newRelativeMove);
    }
    return seq1 + rebuilt;
/*
    moveDirection lastYAxis = (moveDirection) charMoveToInt(s.listYAxis[s.sequence.size() - 1]);
    moveDirection relativeMovementForFirstSeq2Move = Nb_Moves;
    moveDirection newYAxis = Nb_Moves;
    }

    int posTurn = 1;
    while((posTurn < seq2.size()) && (seq2[posTurn] == 'S')){
        posTurn++;
    }

    // Now, the structure should be inverted if necessary, in order to have the same end Yaxis of seq1 and the same starting axis for seq2.
    // do before modif!!!
    moveDirection firstAbsMoveSeq2 = (moveDirection) charMoveToInt(seq2[0]);
    struct3D s2(seq2);

    // now can modif
    seq2[0] = intToMoveChar(relativeMovementForFirstSeq2Move);
    if(posTurn >= seq2.size()){ // finishes in SSSS : no care about rotation
        return seq1 + opposite(seq2);
    }


    moveDirection yAxisAfterFirstMoveSeq2 =newYAxis ;
    moveDirection firstAbsTurnSeq2 = (moveDirection) charMoveToInt(s2.seqAbsMoves[posTurn]);
    moveDirection firstRelativeTurnSeq2 = Nb_Moves;
    cout << "The first absMove in seq2 is : Rel " << intToMoveChar(firstAbsMoveSeq2) << endl;
    cout << "The new Y axis after doing the first rel move " << intToMoveChar(relativeMovementForFirstSeq2Move) << " is |- = " << intToMoveChar(newYAxis) << endl;
    cout << "The first turn    in seq2 is : Rel " << seq2[posTurn] << endl;

    for(int i = 0; i < (int) Nb_Moves; ++i){
        if(nextAbsoluteMove(firstAbsMoveSeq2 , yAxisAfterFirstMoveSeq2, (moveDirection) i).first == firstAbsTurnSeq2){
            firstRelativeTurnSeq2 = (moveDirection) i;
        }
    }
    if(firstRelativeTurnSeq2 == Nb_Moves) cerr << "ERR: turning not found ??" << endl;
    seq2[posTurn] = intToMoveChar(firstRelativeTurnSeq2);
    return seq1 + seq2;
*/
}








string cut(string a){ // removes the ---. Note: should start with the ---, will not detect if SU---SSU
    int L = a.size();
    int i = 0;
    while((i < L) && (a[i] == '-')){
        ++i;
    }
    if(i == L) return string("");
    if(i == 0) return a;
    return a.substr(i, L-i);
}

string fill(string a){
    if(a.size() > MAXL) {
        cerr << "ERR: fill(" << a << "), exceeding the MAXL size " << MAXL << endl;
        return a;
    }
    return string(MAXL - a.size(), '-') + a;
}

// Function that generates all sequences (ULRDS)* (recursive)
// outputs are of size MAXL. Implemented by using char arrays: base[MAXL+1] and then converts to string;
vector<string> generateRandomSeqSize(int L){

    if(L > 9) cerr << "ERR: generateRandomSeqSize(" << L << "), too big number !" << endl;
    char base[MAXL+1];
    for(int i = 0; i < MAXL; ++i){
        base[i] = '-';
    }
    base[MAXL] = '\0';

    vector<string> l3m;
    if(L == 0){
        l3m.push_back(string(base));
        return l3m;
    }
    vector<string> smallerProts = generateRandomSeqSize(L-1);
    int smallS = smallerProts.size();
    for(int i = 0; i < smallS; ++i){
        char newWord[MAXL+1];
        strcpy(newWord,smallerProts[i].c_str()); // dest, source
        newWord[MAXL - L] = 'R';
        l3m.push_back(newWord);
        newWord[MAXL - L] = 'L';
        l3m.push_back(newWord);
        newWord[MAXL - L] = 'U';
        l3m.push_back(newWord);
        newWord[MAXL - L] = 'D';
        l3m.push_back(newWord);
        newWord[MAXL - L] = 'S';
        l3m.push_back(string(newWord));
    }
    return l3m;
}

// Function that generate all acceptable sequences for a relative structure -----[S]*U[ULRDS]*
// all outputs are strings of size MAXL
vector<string> generateRelativeSeqSizeL(int L){
    char base[MAXL+1];
    for(int i = 0; i < MAXL; ++i){
        base[i] = '-';
    }
    base[MAXL] = '\0';

    vector<string> l3m;
    if(L == 0){
        l3m.push_back(string(base));
        return l3m;
    }
    if(L == 1){
        base[MAXL-1] = 'S';
        l3m.push_back(string(base));
        return l3m;
    }
    // maybe L=2 is not necessary
    if(L == 2){
        base[MAXL-2] = 'S';
        base[MAXL-1] = 'S';
        l3m.push_back(string(base));
        base[MAXL-1] = 'U';
        l3m.push_back(string(base));
        return l3m;
    }

    // first choice : start with S (recursive), + prot
    vector<string> smallerProts = generateRelativeSeqSizeL(L-1);
    int smallS = smallerProts.size();
    for(int i = 0; i < smallS; ++i){
        char newWord[MAXL+1];
        strcpy(newWord,smallerProts[i].c_str()); // dest, source
        newWord[MAXL-L] = 'S';
        l3m.push_back(string(newWord));
    }

    // second choice : start with SU + random
    vector<string> smallerRandom = generateRandomSeqSize(L-2);
    int smallS2 = smallerRandom.size();
    for(int i = 0; i < smallS2; ++i){
        char newWord[MAXL+1];
        strcpy(newWord,smallerRandom[i].c_str()); // dest, source
        newWord[MAXL-L] = 'S';
        newWord[MAXL-L+1] = 'U';
        l3m.push_back(string(newWord));
    }
    // Number of acceptable sequences versus possible random sequences
    // 1 => 3  => 10  => 42
    // 5 => 25 => 125 => 625 => 3125 => 15625

    return l3m;
}

// Pools together all the possible sequences of different sizes up to L included.
vector<string> generateRelativeSeqSizeLAndLess(int L){
    if(L == 0) return generateRelativeSeqSizeL(L);
    vector<string> v1 = generateRelativeSeqSizeL(L);
    vector<string> v2 = generateRelativeSeqSizeLAndLess(L-1);
    v2.insert(v2.end(), v1.begin(), v1.end());
    return v2;
}

// Resizes the strings in a vector (of whatever size) to size L. Note: information might be lost if cut too much.
vector<string> resized(vector<string> v2, int L){
    for(unsigned int i = 0; i < v2.size() ; ++i){
        if((int) v2[i].size() >= L) v2[i] = v2[i].substr(v2[i].size()-L, L); // to check
    }
    return v2;
}


/* Example of output would be :
    l3m.push_back("-----"); // 0
    l3m.push_back("----S"); // 1
    l3m.push_back("---SS"); // 1 + 1
    l3m.push_back("---SU");
    l3m.push_back("--SSS"); // 1 + 1 + 5
    l3m.push_back("--SSU");
    l3m.push_back("--SUR");
    l3m.push_back("--SUL");
    l3m.push_back("--SUU");
    l3m.push_back("--SUD");
    l3m.push_back("--SUS");
    l3m.push_back("-SSSS"); // 1 + 1 + 5 + 25
    l3m.push_back("-SSSU");
    l3m.push_back("-SSUR");
    l3m.push_back("-SSUL");
    l3m.push_back("-SSUU");
    l3m.push_back("-SSUD");
    l3m.push_back("-SSUS");
    l3m.push_back("-SURR");
    l3m.push_back("-SURL");
    l3m.push_back("-SURU");
    l3m.push_back("-SURD");
    l3m.push_back("-SURS");
    l3m.push_back("-SULR");
    l3m.push_back("-SULL");
    l3m.push_back("-SULU");
    l3m.push_back("-SULD");
    l3m.push_back("-SULS");
    l3m.push_back("-SUUR");
    l3m.push_back("-SUUL");
    l3m.push_back("-SUUU");
    l3m.push_back("-SUUD");
    l3m.push_back("-SUUS");
    l3m.push_back("-SUDR");
    l3m.push_back("-SUDL");
    l3m.push_back("-SUDU");
    l3m.push_back("-SUDD");
    l3m.push_back("-SUDS");
    l3m.push_back("-SUSR");
    l3m.push_back("-SUSL");
    l3m.push_back("-SUSU");
    l3m.push_back("-SUSD");
    l3m.push_back("-SUSS");*/




// function that shows numbers of sequences and efficiency of encoding, depending on the number of bytes used.
void testEncoding(){
    int cpt = 0;
    for(int L = 0; L < 7; ++L){ // maximum possible is L < 13 here.

        //cout << "\n\n\n ================== L = " << L << " ======================== \n\n\n";
        vector<string> myTest = generateRelativeSeqSizeL(L);
        int NB = myTest.size();
        //cpt += NB;
        for(int i = 0; i < NB; ++i){
            cout << cpt++ << "\t" << myTest[i] << endl;
        }
        long int nb =  cpt + pow(5.0, L);

        //long int pM5 = (long int) cpt;
        int p2 = 1+ (long int) (log((double) cpt) / log(2.0));
//        cout << "L=" << L << ", 5^L =" << pM5 << ", binary " << p2 << " equiv " << (int) pow(2.0, p2) << "\t" << (double) pM5 / (double) (pow(2.0, p2)) << endl;

        //cout << "For coding in size" << L << ", Would make " << nb << "(" << cpt << " + " << pow(5.0, L) << ") vs " << log((double) nb) / log(2.0) << " instead of " << (long int) pow(6.0, L) << endl;
        cout << "For coding starting sequences size" << L << ", Would make " << nb << "(" << cpt << " + " << pow(5.0, L) << ") vs " << p2 << " => " << (long long int ) pow(2.0, (double) p2) << " instead of " << (long int) pow(6.0, L) << endl;

    }
    {
        //unsigned long long int Uim1 = 1;
        unsigned long long int  Ui = 1;
        unsigned long long int  Sum = 1; // for the case L = 1
        unsigned long long int pw2 = 1;
        int cpt2 = 0;

        for(int L = 1; L < 31; ++L){ // maximum possible here.
            unsigned long long int  Uip1 = (unsigned long long) pow(5.0, L-2) + Ui;
            //Uim1 = Ui;
            Ui = Uip1;
            Sum += Ui;
            while((cpt2 < 64) && (pw2 < Sum)){
                cpt2++;
                pw2 *= (unsigned long long int) 2;
            }
            cout << L << "\t" << Sum << "\t" << cpt2 << " => " << pw2-1 << "\tefficiency " << (double) Sum / (double) (pw2-1)  << endl;
        }
    }

    cout << "Efficiency of encoding the last residues in binary (which group to take ?)\n";
    for(int L = 0; L < 16; ++L){
        long int pM5 = (long int) pow(5.0, (double)L);
        int p2 = 1+ (long int) (log(pow(5.0, (double)L)) / log(2.0));
        cout << "L=" << L << ", 5^L =" << pM5 << ", binary " << p2 << " equiv " << (int) pow(2.0, p2) << "\t" << (double) pM5 / (double) (pow(2.0, p2)) << endl;
    }

    {
        int L = 4;
        cpt = 0;
        vector<string> myTest = generateRelativeSeqSizeLAndLess(L);
        int NB = myTest.size();
        for(int i = 0; i < NB; ++i){
            cout << cpt++ << "\t" << myTest[i] << endl;
        }
        cout << "For coding in size" << L << ", " << cpt << endl;
    }

    {
        unsigned long long int test5 = 1;
        unsigned long long int test2 = 1;

        int i2 = 0;
        for(int i = 0; i < 30; ++i){
            test5 = test5 * (unsigned long long int) 5;
            while(test2 < test5){
                i2++;
                test2 = test2 * (unsigned long long int) 2;
            }
            cout << i << "\t" << test5 << "\t" << i2 << "\t" << test2 << "\t" << (double) test5 / (double) test2 << endl;
        }
    }
}

//vector<string> fuse(vector<string> &v1, vector<string> &v2){

//}

#define ERRORSEQ 16000
unsigned long long int intToID(int part1, int part2, int part3, int part4, int part5){
    unsigned long long int res = part1;
    // note : ID=16000 is the only authorized ERROR structure.
    if((part1 < 0) || (part1 >= 200) ||
            (part2 < 0) || (part2 >= 15625+200) || // might be + 199 instead of +200
            (part3 < 0) || (part3 >= 15625+200) ||
            (part4 < 0) || (part4 >= 15625+200) ||
            (part5 < 0) || ((part5 >= 15625+200) && (part5 != 16000)) ) cerr << "ERR: intoToID(" << part1 << "," << part2 << "," << part3 << "," << part4 << "," << part5 << "), part1 should be in [0-255]; parts 2-5 should be in [0-16383]" << endl;
    res = res * (unsigned long long int) 16384 + (unsigned long long int) part2;
    res = res * (unsigned long long int) 16384 + (unsigned long long int) part3;
    res = res * (unsigned long long int) 16384 + (unsigned long long int) part4;
    res = res * (unsigned long long int) 16384 + (unsigned long long int) part5;
    return res;
}

string intToSring(int part1, int part2, int part3, int part4, int part5){
    static bool loaded = false;
    static vector<string>* starts5;
    static vector<string>* starts6; // with '-' at start
    if(!loaded){
        starts5 = new vector<string>(resized(generateRelativeSeqSizeLAndLess(5),5));
        starts6 = new vector<string>(resized(generateRelativeSeqSizeLAndLess(5),6)); // still size 5, but 1 more - before
        loaded = true;
    }

    // first part : starting sequence, on 8 bits forstarting  5-mers (200 possible ones),
    if(part1 >= (int) starts5->size()) {cerr << "Protein ID incorrect " << intToID(part1, part2, part3, part4, part5) << " from part 1=" << part1 << " that should lie within [0:" << starts5->size() << "]" <<endl; return string("ERROR");}
    string Spart1 = (*starts5)[part1];

    // all remaining parts : on 14 bits for 6-mers,
    //  -> either 5^6 mers of the inside of a protein
    //  -> either starting sequences of size 6 (already included in the previous line)
    //  -> either starting sequences of size 5 or less (-*****) => 200 possible cases
    for(int i = 0; i < 4; ++i){
        int partX;
        switch(i){
            case 0: {partX = part2; break;}
            case 1: {partX = part3; break;}
            case 2: {partX = part4; break;}
            case 3: {partX = part5;
                if(part5 == 16000) return string("ERROR");
            break;}
        }

        string SpartX;
        if(partX < 15625){
            int q1 = partX/3125;
            int r1 = partX - 3125*q1;
            int q2 = r1 / 625;
            int r2 = r1 - 625*q2;
            int q3 = r2 / 125;
            int r3 = r2 - 125*q3;
            int q4 = r3 / 25;
            int r4 = r3 - 25*q4;
            int q5 = r4 / 5;
            int r5 = r4 - 5*q5;
            int q6 = r5;
            if((q1 > 4) || (q2 > 4) || (q3 > 4) || (q4 > 4) || (q5 > 4) || (q6 > 4)){
                cerr << "Divisions error " << partX << " -> " << q1 << "," << q2 << "," << q3 << "," << q4 << "," << q5 << "," << q6 << endl;
            }
            char buf[7];
            buf[0] = intToMoveChar(q1);
            buf[1] = intToMoveChar(q2);
            buf[2] = intToMoveChar(q3);
            buf[3] = intToMoveChar(q4);
            buf[4] = intToMoveChar(q5);
            buf[5] = intToMoveChar(q6);
            buf[6] = '\0';
            SpartX = string(buf);
        } else {
            if((partX - 15625) >= (int) starts6->size()) {cerr << "Protein ID incorrect " << intToID(part1, part2, part3, part4, part5) << ", from part 1=" << partX-15625 << " that should lie within [0:" << starts5->size() << "]" <<endl; return string("ERROR");}
            SpartX = (*starts6)[partX - 15625];
            //cout << "got ID and text : " << partX << "\t" << SpartX << endl;
        }
        Spart1.append(SpartX);
    }
    return Spart1;
}


int String6ToInt(string S){
    static bool loaded = false;
    static map<string, int>* dicoStartSeqs;
    if(!loaded){
        dicoStartSeqs = new map<string, int>();
        vector<string> reverseCorresp = resized(generateRelativeSeqSizeLAndLess(5),6);
        int LS = reverseCorresp.size();
        for(int i = 0; i < LS; ++i){
            //cout << 15625 + i << "\t" << reverseCorresp[i] << endl;
            dicoStartSeqs->insert(std::pair<string, int> (reverseCorresp[i], i));
        }
        loaded = true;
    }

    std::map<string,int>::iterator it = (dicoStartSeqs->find(S));
    if (it == dicoStartSeqs->end()){
        int L = S.size();
        if(L != 6) {cerr << "ERR: String6ToInt only takes strings of size 6" << endl; return -1;}
        int res = 0;
        for(int i = 0; i < L; ++i){
            res += charMoveToInt(S[i]);
            if(i < L-1) res = res * 5;
        }
        return res;
    } else {
        if(it->second > (int) dicoStartSeqs->size()) cerr << "ERR: Why does this happen?" << endl;
        return 15625 + it->second;
    }
}

bool checkSyntaxSequence(string toCheck){
    int L= toCheck.size();
    if(L > 29) {cerr << "ERR: this program is only valid for proteins towards size 29, to fit inside 64 bits (prot=" << toCheck << ")" << endl; return false;}
    bool stillMinus = true;
    bool stillS = true;
    for(int i = 0; i < L; ++i){
        // first should only be '-' or 'S'
        if(stillMinus && (toCheck[i] != '-')){
            if(toCheck[i] == 'S') stillMinus = false;
            else {
                cerr << "ERR: incorrect protein structure : " << toCheck << endl;
                return false;
            }
        }
        // then should be either 'S' or 'U'
        if((!stillMinus) && (stillS) && (toCheck[i] != 'S')){
            if(toCheck[i] == 'U') stillS = false;
            else {
                cerr << "ERR: incorrect protein structure : " << toCheck << endl;
                return false;
            }
        }
        // if U then can do everything now
        if((!stillMinus) && (!stillS)){
            if(charMoveToInt(toCheck[i]) == Nb_Moves_Relative){ // in case it's a -
                cerr << "ERR: incorrect protein structure : " << toCheck << endl;
                return false;
            }    // just to raise an error if wrong characters

        }
    }
    return true;
}

// will be completed by '-'. so, can already contain '-' or not, no difference.
unsigned long long int relativeToInt(string prot){
    size_t L = prot.size();
    if(!checkSyntaxSequence(prot)) return ERRORSEQ; // 200 is the error sequence

    static bool loaded = false;
    static map<string, int>* dicoStartSeqs;
    if(!loaded){
        dicoStartSeqs = new map<string, int>();
        vector<string> reverseCorresp = resized(generateRelativeSeqSizeLAndLess(5),5);
        size_t LS = reverseCorresp.size();
        //cout << "note : " << LS << "starting sequences of size 5" << endl;
        for(size_t i = 0; i < LS; ++i){
            //cout << reverseCorresp[i] << "\t" << i << endl;
            dicoStartSeqs->insert(std::pair<string, int> (reverseCorresp[i], i));
        }
        loaded = true;
    }
    string enlarged = string(29-L,'-') + prot;
    string part1 = enlarged.substr(0,5);
    string part2 = enlarged.substr(5,6);
    string part3 = enlarged.substr(11,6);
    string part4 = enlarged.substr(17,6);
    string part5 = enlarged.substr(23,6);
    //cout << "Info : for prot " << prot << ", substrings and subIDs are:" << endl;
    std::map<string,int>::iterator it = (dicoStartSeqs->find(part1));
    if (it == dicoStartSeqs->end()){
        cerr << "ERR : Sequence (" << part1 << " not found in the list of starting sequences for size 5" << endl;
    }
    /*cout << part1 << " -> (start, lg5:)" << dicoStartSeqs->find(part1)->second << endl;
    cout << part2 << " -> (all, lg6) " << String6ToInt(part2) << endl;
    cout << part3 << " -> (all, lg6)" << String6ToInt(part3) << endl;
    cout << part4 << " -> (all, lg6)" << String6ToInt(part4) << endl;
    cout << part5 << " -> (all, lg6)" << String6ToInt(part5) << endl;*/



    int IDpart1 = (dicoStartSeqs->find(part1)->second);
    //cout << IDpart1 << "\t";
    unsigned long long int res = IDpart1;
    for(int i = 0; i < 4; ++i){
        //cout << res << " will get + " << String6ToInt(enlarged.substr(5+(i*6),6)) << endl;
        res = (res * (unsigned long long int) 16384) + String6ToInt(enlarged.substr(5+(i*6),6));
    }
    //cout << res << endl;
    return res;


}



string intToRelative(unsigned long long protID){
    unsigned long long part1 = protID; part1 >>= (4*14);
    unsigned long long part2 = protID; part2 >>= (3*14); part2 = part2 & (unsigned long long) 16383; // 32767 = 1...1 14 times
    unsigned long long part3 = protID; part3 >>= (2*14); part3 = part3 & (unsigned long long) 16383;
    unsigned long long part4 = protID; part4 >>= (1*14); part4 = part4 & (unsigned long long) 16383;
    unsigned long long part5 = protID;                   part5 = part5 & (unsigned long long) 16383;
    //cout << "Info : for ID " << protID << ", parts are " << part1 << "," << part2 << "," << part3 << "," << part4 << "," << part5 << endl;
    return intToSring(part1, part2, part3, part4, part5);
}

// when sequences have "-" inside ...
size_t seqLength(string sequence){
    size_t s = sequence.size();
    for(size_t i = 0; i < s; ++i){
       if(sequence[i] != '-') return s - i;
    }
    return 0;
}

string seqClean(string sequence){
    size_t L = seqLength(sequence);
    return sequence.substr(sequence.size() - L, L);
}




/// Note the intToAbsolute and vice-versa don't raise errors if the sequences are bullshit. Should be implemented.
string intToAbsolute(unsigned long long protID){
    string res = intToRelative(protID);
    //cout << res << " has length " << seqLength(res) << endl;
    if(!(seqClean(res).substr(0,2).compare(string("SU")))){ // start is different than Backwards
        size_t lg = seqLength(res);
        //return res.substr(2,(res.size()-2));
        return string(res.size() - lg + 2, '-') + res.substr(res.size() - lg + 2,(res.size()-2 - lg));
    } else { // starts with B. Syntax : SSU...
        if(res.size() < 3) cerr << "??? IntToAbsolute: too short sequence" << endl;
        size_t lg = seqLength(res);
        return string(res.size() - 1 - lg + 3, '-') + string("B") + res.substr(res.size() - lg + 3,(res.size()-3 - lg));
    }
}

unsigned long long absoluteToInt(string prot){
    if((prot.size() > 0) && (prot[0] == 'B')){
        prot[0] = 'U';
        return relativeToInt(string("SS")+prot);
    } else {
        return relativeToInt(string("SU")+prot);
    }
}

void testCoding(string totest){
    unsigned long long int t1 = relativeToInt(totest);
    string p2 = intToRelative(t1);

    if(totest.size() <= 29) {
        string enlarged = string(29-totest.size(),'-') + totest;
        if(p2.compare(enlarged)) cerr << "Test failed for sequence " << p2 << " != (init) " << enlarged << endl;
        cout << totest << " -> " << t1 << " -> " << p2 << endl;
        cout << endl << endl;
    }
}

void testCodingOnExampleSequences(){
    cout << "========= filling and cutting string sequences to MAXL size, with '-' " << endl;
    cout << fill(string("SULS")) << endl;
    cout << fill(string("")) << endl;
    cout << fill(string("-SSS")) << endl;
    cout << fill(string("SUSLSLLSLLSUSUUSSSSDDRDSDRRDSRDSRDS")) << endl;
    cout << intToRelative(15) << fill(intToRelative(15)) << endl;
    cout << cut("-----SUSLD") << endl;
    cout << cut("SUSLD") << endl;
    cout << cut("---------")<< endl;
    cout << cut("") << endl;
    cout << cut("-") << endl;
    // correct sequences
   testCoding(string("SSUDLRURLSURLDUSUSDDURDSSU"));
   testCoding(string("SSUDLRUSUUSUSDDURDSSU"));
   testCoding(string("SUUDRLDULDSDDURDLSSU"));
   testCoding(string("SSUDRLDULDSDDURDSSU"));
   testCoding(string("SSUDRLDLUSDDURDSSU"));
   testCoding(string("SSUDRLDUSDDURDSSU"));
   testCoding(string("SSUDLLDURDSSU"));
   testCoding(string("SSUDLURLU"));
   testCoding(string("---SSUDLURLU"));
   testCoding(string("SSSSSSSSSSSSSSSS"));
   testCoding(string("SU"));
   testCoding(string("S"));

   // Incorrect sequences
   testCoding(string("SSU-DLRURLSURLDUSUSDDURDSSU"));
   testCoding(string("SSADLRUSUUSUSDDURDSSU"));
   testCoding(string("UUDRLDULDSDDURDLSSU"));
   testCoding(string("SSDDRLDULDSDDURDSSU"));
   testCoding(string("SSSSSDUDRLDLUSDDURDSSU"));
   testCoding(string("SUDLRLRSSSLUDLSLLSSUDRLDUSDDURDSSU"));

}




string printVector(vector<string> v){
    //cerr << "Fin" << v.size() << endl;
    stringstream res;
    for(unsigned int i = 0; i < v.size(); ++i){
        res << i << "\t" << v[i] << endl;
    }
    //cerr << res.str() << endl;
    return res.str();
}

string printVector(vector<int> v){
    stringstream res;
    for(unsigned int i = 0; i < v.size(); ++i){
        if(i > 0) res << "\t";
        res << v[i];
    }
    return res.str();
}

string printVector(vector<double> v){
    stringstream res;
    for(unsigned int i = 0; i < v.size(); ++i){
        if(i > 0) res << "\t";
        res << v[i];
    }
    return res.str();
}

string printVector(vector<size_t> v){
    stringstream res;
    for(unsigned int i = 0; i < v.size(); ++i){
        if(i > 0) res << "\t";
        res << v[i];
    }
    return res.str();
}


set<int> stringToSet(string text, char sep){
    if(sep != '\t'){
        std::replace(text.begin(), text.end(), sep, '\t');
    }
    stringstream read(text);
    int position;
    set<int> res;
    while(read >> position){
        res.insert(position);
    }
    return res;
}

bool isIncluded(set<int>& thisSet, set<int>& intoThatOne){
    size_t S1 = thisSet.size();
    size_t S2 = intoThatOne.size();
    if(S1 > S2) return false;
    set<int> merged = union_sets(thisSet, intoThatOne);
    return(merged.size() == intoThatOne.size());
}

string setToString(set<int> mySet, char sep){
    stringstream res;
    for(set<int>::iterator it = mySet.begin(); it != mySet.end(); ++it){
        if(it != mySet.begin()) res << sep;
        res << *it;
    }
    return res.str();
}

set<int> vectorToSet(vector<int> v){
    return set<int>(v.begin(), v.end());
}

vector<int> setToVector(set<int> v){
    return vector<int>(v.begin(), v.end());
}


/*string convertMoveWithNewYAxis(string seq, moveDirection newYAxis){
    if(seq.size() <= 1) return seq; // first move is absolute => not affected by the y axis
    struct3D s(seq);
    int SL = seq.size();
    moveDirection currentDir = (moveDirection) charMoveToInt(seq[0]);

    for(int i = 0; i < SL; ++i){

    }
}*/





