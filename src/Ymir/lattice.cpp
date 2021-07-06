#include "lattice.h"
#include "proteins.h"

vector<int> lattice::positionFromID(int IDcoord){
    int z = IDcoord / (XWidth * YWidth);
    int y = (IDcoord - XWidth * YWidth * z)/YWidth;
    int x = (IDcoord - XWidth * YWidth * z - YWidth*y);
    //cout << x << "," << y << "," << z << endl;
    vector<int> res(3,0);
    res[0] = x;
    res[1] = y;
    res[2] = z;
    return res;
}

int lattice::idFromPosisition(vector<int> position){
    if(position.size() != 3){cerr << "ERR: getIdFromPosition(vector), the vector should have size 3;"; return -1;}
    if((position[0] < 0) || (position[1] < 0) || (position[2] < 0) || (position[0] >= XWidth) || (position[1] >= YWidth) || (position[2] >= ZWidth))
        cerr << "ERR: lattice::idFromPosition, Incorrect position : " << position[0] << "," << position[1] << "," << position[2] << endl;
    return position[0] + YWidth*position[1] + YWidth*ZWidth*position[2];
}

int lattice::idFromPosisition(int x, int y, int z){
    if((x < 0) || (x >= XWidth) || (y < 0) || (y >= YWidth) ||(z < 0) || (z >= ZWidth)){
        cerr << "ERR: lattice::idFromPosisition(" << x << "," << y << "," << z << "), out of bounds " << endl;
    }
    return x + YWidth* y + YWidth*ZWidth* z;
}

int lattice::idFromPosisition(double x, double y, double z){
    if((x < 0) || (x >= XWidth) || (y < 0) || (y >= YWidth) ||(z < 0) || (z >= ZWidth)){
        cerr << "ERR: lattice::idFromPosisition(" << x << "," << y << "," << z << "), out of bounds " << endl;
    }
    return static_cast<int>(x) + YWidth* static_cast<int>(y) + YWidth*ZWidth* static_cast<int>(z);
}

vector<int> lattice::idNeighbors(int IDcoord){
    //int z = IDcoord / 1000000;
    //int y = (IDcoord - 1000000*z)/1000;
    //int x = (IDcoord - 1000000*z - 1000*y);
    vector<int> res;
    res.push_back(IDcoord+1);
    res.push_back(IDcoord-1);
    res.push_back(IDcoord+YWidth);
    res.push_back(IDcoord-YWidth);
    res.push_back(IDcoord+YWidth*ZWidth);
    res.push_back(IDcoord-YWidth*ZWidth);
    return res;
}

vector<int> lattice::idFreeNeighbors(int IDcoord, set<int>& listPositions){ // according to the already occupied space.
    vector<int> toTest = idNeighbors(IDcoord);
    int nb = toTest.size();
    vector<int> res;
    for(int i = 0; i < nb; ++i){
         if((listPositions.find(toTest[i])) == listPositions.end()){ //i.e. not found
             res.push_back(toTest[i]);
         }
    }
    return res;
}

vector<int> lattice::idLargeNeighbors(int IDcoord){
    //int z = IDcoord / 1000000;
    //int y = (IDcoord - 1000000*z)/1000;
    //int x = (IDcoord - 1000000*z - 1000*y);
    vector<int> res;
    res.push_back(IDcoord+1);
    res.push_back(IDcoord-1);
    res.push_back(IDcoord+YWidth);
    res.push_back(IDcoord-YWidth);
    res.push_back(IDcoord+YWidth*ZWidth);
    res.push_back(IDcoord-YWidth*ZWidth);

    res.push_back(IDcoord+1 + YWidth);
    res.push_back(IDcoord+1 - YWidth);
    res.push_back(IDcoord-1 + YWidth);
    res.push_back(IDcoord-1 - YWidth);
    res.push_back(IDcoord+1 +YWidth*ZWidth);
    res.push_back(IDcoord+1 -YWidth*ZWidth);
    res.push_back(IDcoord-1 +YWidth*ZWidth);
    res.push_back(IDcoord-1 -YWidth*ZWidth);
    res.push_back(IDcoord+YWidth*ZWidth + YWidth);
    res.push_back(IDcoord+YWidth*ZWidth - YWidth);
    res.push_back(IDcoord-YWidth*ZWidth + YWidth);
    res.push_back(IDcoord-YWidth*ZWidth - YWidth);
    return res;
}

vector<int> lattice::idFreeLargeNeighbors(int IDcoord, set<int>& listPositions){ // according to the already occupied space.
    vector<int> toTest = idLargeNeighbors(IDcoord);
    int nb = toTest.size();
    vector<int> res;
    for(int i = 0; i < nb; ++i){
         if((listPositions.find(toTest[i])) == listPositions.end()){ //i.e. not found
             res.push_back(toTest[i]);
         }
    }
    return res;
}

int lattice::getIdFromAbsoluteMove(int IDcurrent, moveDirection MoveToNext){
    if(MoveToNext == UnDefined) cerr << "ERR: lattice::getIdFromAsoluteMove(), Undefined is not an acceptable direction" << endl;
    vector<int> pos = lattice::positionFromID(IDcurrent);
    vector<int> actualMove = moveVector(MoveToNext);
    pos[0] += actualMove[0];
    pos[1] += actualMove[1];
    pos[2] += actualMove[2];
    return lattice::idFromPosisition(pos);
}

bool lattice::testPos(int x, int y, int z){
    if((x < 0) || (x >= XWidth)) return false;
    if((y < 0) || (y >= YWidth)) return false;
    if((z < 0) || (z >= ZWidth)) return false;
    return true;
}

bool lattice::testPos(vector<int> v){
    if(v.size() != 3) return false;
    return testPos(v[0], v[1], v[2]);
}

int lattice::centralPosition(){ // ususally 31, 31, 31
    return idFromPosisition(XWidth / 2, YWidth / 2, ZWidth / 2);
}

bool lattice::areNeighbors(int pos1, int pos2){
    if(abs(pos1 - pos2) == 1) return true;
    if(abs(pos1 - pos2) == XWidth) return true;
    if(abs(pos1 - pos2) == XWidth*YWidth) return true;
    return false;
}

//bool lattice::isNeighbor(int pos1, int pos2){
//    return ((abs(pos1 - pos2) == 1) || (abs(pos1 - pos2) == XWidth) || (abs(pos1 - pos2) == XWidth*YWidth));
//}

void testLattice(){
    cout << "Testing the lattice class" << endl;
    cout << "predefined sizes : " << XWidth << "," << YWidth << "," << ZWidth << endl;
    vector<int> toTest = {0, 15, 200, 1000, 15000};
    for(unsigned int i = 0; i < toTest.size(); ++i){
        cout << "positionFromID\t" << toTest[i] << "->" << printVector(lattice::positionFromID(toTest[i])) << " -> " << lattice::idFromPosisition(lattice::positionFromID(toTest[i])) << endl;
    }
    cout << "Neighbors " << endl;
    // Note: don't use neighbors from border values, will become negative ...
    toTest = {15000, 15250, 18500};
    for(unsigned int i = 0; i < toTest.size(); ++i){
        cout << "Neighbors for " << toTest[i] << ": " << printVector(lattice::positionFromID(toTest[i])) << endl;
        vector<int> idNeigh = lattice::idNeighbors(toTest[i]);
        for(unsigned int j = 0; j < idNeigh.size(); ++j){
            cout << "   ... " << idNeigh[j] << ": " << printVector(lattice::positionFromID(idNeigh[j])) << endl;
        }

    }
 /*   cout << "Free Neighbors for " << toTest[i] << ": " << printVector(lattice::positionFromID(toTest[i])) << endl;
    vector<int> idNeigh = lattice::idNeighbors(toTest[i]);
    for(int j = 0; j < idNeigh.size(); ++j){
        cout << "   ... " << idNeigh[j] << ": " << printVector(lattice::positionFromID(idNeigh[j])) << endl;
    }  */
}
