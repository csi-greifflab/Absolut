#include "discretize.h"
#include <fstream>
#include "../Ymir/plot3d.h"
#include "../Ymir/lattice.h"
#include "../latFit/biu/Point.hh"
#include "../latFit/latFit.h"
#include "../Ymir/lattice.h"
#include "topology.h"
#include <sstream>
#include <map>

bool isVector3D(vector3D &v){return v.size() == 3;}
bool isMat3D(mat3D & M){
    if(M.size() != 3) return false;
    for(size_t i = 0; i < 3; ++i){
        if(M[i].size() != 3) return false;
    }
    return true;
}

string print(vector3D v){
    stringstream res;
    res << "[ " << v[0] << ", " << v[1] << ", " << v[2] << "]";
    return res.str();
}

vector3D operator *(double d, vector3D v){
    if(!isVector3D(v)) return vector3D();
    vector3D res = {d*v[0], d*v[1], d*v[2]};
    return res;
}

vector3D operator *(mat3D M, vector3D v){
    vector3D res = vector<double>(3, NAN);
    if(!isMat3D(M)) return res;
    if(!isVector3D(v)) return res;
    for(size_t i = 0; i < 3; ++i){
        res[i] = 0;
        for(size_t j = 0; j < 3; ++j){
            res[i] +=  M[i][j] * v[j];
        }
    }
    return res;
}

mat3D emptyMat(){
    mat3D m = vector< vector<double>>();
    m.resize(3);
    m[0].resize(3,NAN);
    m[1].resize(3,NAN);
    m[2].resize(3,NAN);
    return m;
}

mat3D Rotx(double a){
    mat3D m = emptyMat();
    m[0][0] = 1        ;     m[0][1] = 0        ;    m[0][2] = 0        ;
    m[1][0] = 0        ;     m[1][1] = cos(a)   ;    m[1][2] = -sin(a)  ;
    m[2][0] = 0        ;     m[2][1] = sin(a)   ;    m[2][2] = cos(a)   ;
    return m;
}

mat3D Roty(double a){
    mat3D m = emptyMat();
    m[0][0] = cos(a)   ;     m[0][1] = 0        ;    m[0][2] = sin(a)   ;
    m[1][0] = 0        ;     m[1][1] = 1        ;    m[1][2] = 0        ;
    m[2][0] = -sin(a)  ;     m[2][1] = 0        ;    m[2][2] = cos(a)   ;
    return m;
}

mat3D Rotz(double a){
    mat3D m = emptyMat();
    m[0][0] = cos(a)   ;     m[0][1] = -sin(a)  ;    m[0][2] = 0        ;
    m[1][0] = sin(a)   ;     m[1][1] = cos(a)   ;    m[1][2] = 0        ;
    m[2][0] = 0        ;     m[2][1] = 0        ;    m[2][2] = 1        ;
    return m;
}
vector3D baseXAxis(){vector3D v = {1., 0., 0.}; return v;}
vector3D baseYAxis(){vector3D v = {0., 1., 0.}; return v;}
vector3D baseZAxis(){vector3D v = {0., 0., 1.}; return v;}

// Starting from the z vector, phi is around the y axis (first), -90 to +90, and then theta is around the z axis (0-360),
vector3D rotate (vector3D v, double theta, double phi){
    mat3D m1 = Roty(phi);
    mat3D m2 = Rotz(theta);
    return (m2 * (m1 * v));
}

double PI(){return 3.14159265358979323846;}

// equally spaced directions
vector<polarCoord> getDirections(int resulution){
    vector<polarCoord> res = vector<polarCoord>();
    if(resulution < 1) return res;
    //if((n < 0) || (n >= resulution)) {cerr << "ERR: getDirection(" << n << "," << totalN << "), incorrect n)"; return res;}
    // choice of latitude
    for(int n = 1; n <= resulution-1; ++n){
        //n = 2;
        double thetaJ = PI() * n / resulution - PI() / 2.0; // j = 0 to n-1 latitude
        // choice of points
        int nj = floor(0.5 + sqrt(3) * resulution * cos(thetaJ));
        //cout << nj << " positions " << endl;
        double shift = 0;
        if((n % 2) != 0) shift += 0.5;
        for(int j = 0; j < nj; ++j){
            // their theta and phi are inverted, and their phi is from equator, our phi is from north pole (z axis)
            res.push_back(polarCoord(-PI() + 2.0 * PI() *(((double) j + shift) / ((double) nj)), (PI() * 0.5 - thetaJ)));
            //cout << "theta = " << -PI() + 2.0 * PI() *(((double) j + shift) / ((double) nj)) << "\tPhi = " <<  (PI() * 0.5 - thetaJ) << endl;
        }
        //n = 1000;
    }
    return res;
}
// algo 1
// https://projecteuclid.org/download/pdf_1/euclid.em/1067634731
// algo 2
// https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf






void testRotate3D(){

#ifdef ALLOW_GRAPHICS
    vector<vector3D>* equallyDistributedPositions = new vector<vector3D>();
    for(unsigned int i = 0; i < 10; ++i){
        vector3D v = 10*rotate(baseZAxis(), 0, i*PI()/10.);
        cout << print(v) << endl;
        equallyDistributedPositions->push_back(v);
    }
    for(unsigned int i = 0; i < 10; ++i){
        vector3D v = 10*rotate(baseZAxis(), i*PI()/10., PI()/4.);
        cout << print(v) << endl;
        equallyDistributedPositions->push_back(v);
    }
    glDisplay();
    addToDisplay(equallyDistributedPositions);
    glutMainLoop();
#endif
}

void testEquallyDistributedPositions(){

#ifdef ALLOW_GRAPHICS

    vector<polarCoord> dirs = getDirections(20);
    vector<vector3D>* equallyDistributedPositions = new vector<vector3D>();
    for(unsigned int i = 0; i < dirs.size(); ++i){
        vector3D v = 10*rotate(baseZAxis(), dirs[i].first, dirs[i].second);
        //cout << print(v) << endl;
        equallyDistributedPositions->push_back(v);
    }
    glDisplay();
    addToDisplay(equallyDistributedPositions);
    glutMainLoop();
#endif
}

vector3D direction(double x1, double y1, double z1, double x2, double y2, double z2){
    vector3D res = {x2 - x1, y2 - y1, z2 - z1};
    return res;
}

vector3D direction(vector3D pos1, vector3D pos2){
    if((!isVector3D(pos1)) || (!isVector3D(pos2))) cerr << "ERR: wrong length of vectors in operator ^" << endl;
    return direction(pos1[0], pos1[1], pos1[2], pos2[0], pos2[1], pos2[2]);
}

vector3D operator ^ (vector3D u, vector3D v){
    if((!isVector3D(u)) || (!isVector3D(v))) cerr << "ERR: wrong length of vectors in operator ^" << endl;
    vector3D res = {u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0] * v[2], u[0]*v[1] - u[1] * v[0]};
    return res;
}

vector3D operator + (vector3D u, vector3D v){
    if((!isVector3D(u)) || (!isVector3D(v))) cerr << "ERR: wrong length of vectors in operator +" << endl;
    vector3D res = {u[0] + v[0], u[1] + v[1], u[2] + v[2]};
    return res;
}

double operator * (vector3D u, vector3D v){
    if((!isVector3D(u)) || (!isVector3D(v))) cerr << "ERR: wrong length of vectors in operator +" << endl;
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

double norm(vector3D v){
    if(!isVector3D(v)) cerr << "ERR: wrong length of vectors in operator ^" << endl;
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}



std::pair<vector3D, vector3D> moveVector(vector3D Xaxis, vector3D Yaxis, char moveLetter){
    if(!isVector3D(Xaxis)) cerr << "ERR: MoveVector, input vector with wrong size, should be 3" << endl;
    vector3D resX = {NAN, NAN, NAN};
    vector3D resY = {NAN, NAN, NAN};
    double distance = norm(Xaxis);
    if(distance < 1e-16) cerr << "ERR: moveVector, Xaxis is infinitesimal/null" << endl;
    if(fabs(distance - norm(Yaxis)) > 0.1*distance){
        cout << "ERR: moveVector: Please provide vectors of same length for Xaxis and Yaxis, such that the return has same size as well." << endl;
        cout << "X=" << print(Xaxis) << ", Y=" << print(Yaxis) << endl;
    }
    switch(moveLetter){
        case 'R': {resX = - 1.0 * Yaxis;            resY = Xaxis; break;}
        case 'L': {resX = Yaxis;                    resY = -1.0 * Xaxis; break;}
        case 'U': {resX = Xaxis ^ Yaxis;            resY = Yaxis; break;}
        case 'D': {resX = -1.0 * (Xaxis ^ Yaxis);   resY = Yaxis; break;}
        case 'S': {resX = Xaxis;                    resY = Yaxis; break;}
        case 'B': {resX = -1.0 * Xaxis;             resY = -1.0 * Yaxis; break;}
        default:  {
        /*cerr << "ERR: moveVector(UnDefined) is not an allowed direction" << endl This error always happens during init; */
        break;}
    }
    resX = (distance / norm(resX)) * resX;
    resY = (distance / norm(resY)) * resY;
    if(fabs(norm(resX) - distance) > 0.0001 * distance) cerr << "ERR: moveVector, failed to produce a vector X with same length. Why?" << endl;
    if(fabs(norm(resY) - distance) > 0.0001 * distance) cerr << "ERR: moveVector, failed to produce a vector Y with same length. Why?" << endl;
    return std::pair<vector3D, vector3D>(resX, resY);
}

bool pointsAreEquall(vector3D v1, vector3D v2){
    return (norm(direction(v1, v2)) < 0.1);
}

void latFitToLattice::transform(){

    // initialize variables
    string structure;
    int position = lattice::centralPosition();
    structure.push_back('S');
    initXAxis = {NAN, NAN, NAN};
    initYAxis = {NAN, NAN, NAN};
    listIntPositions.clear();
    listRotatedOriginPositions.clear();
    IDeachStartingResidue.clear();

    if(IDresidues.size() != listPositions.size()){
        cerr << "ID residues has wrong size " << IDresidues.size() << endl;
        return;
    }
    if(listPositions.size() < 3) {
        return;
    }
    // Now we have at least 3 points.
    // 1- find the first turn.
    int IDfirstNonAlignedPoint = 2; // it's at least the 3rd point (position 2)
    bool firstPointFound = false;
    vector3D lastDirection = direction(listPositions[0], listPositions[1]);
    IDeachStartingResidue.push_back(IDresidues[0]);
    vector3D lastYaxis = {NAN, NAN, NAN}; // only starts at first turn
    //cout << listPositions.size() << " positions to discretize \n";
    for(unsigned int i = 2; i < listPositions.size(); ++i){
        vector3D newDirection = direction(listPositions[i-1], listPositions[i]);
        //cout << "oldX: " << print(lastDirection) << ", oldY " << print(lastYaxis) << ", newDir " << print(newDirection) << ", norms "
        //     << norm(lastDirection) << "," << norm(lastYaxis) << "," << norm(newDirection) << endl;

        if(firstPointFound){ // i.e. started turning
            vector<char> tryMoves = {'S', 'U', 'D', 'L', 'R', 'B'};
#define nMovesRel 5
#define nMovesAbs 6
            bool found = false;
            //if(IDresidues[i] != IDresidues[i-1] + 1) cerr << "Detected jump!" << endl;
            //cout << "Order is " << IDresidues[i] << "\t" << IDresidues[i-1] << endl;
            for(unsigned int j = 0; j < nMovesRel; ++j){
                if(pointsAreEquall(listPositions[i], listPositions[i-1] + moveVector(lastDirection, lastYaxis, tryMoves[j]).first)){
                    if(!found) {/*cout << tryMoves[j]; cout << "Dont get why should COUT here" << endl;*/}
                    else cout << "Euh... two moves give same position??" << endl;
                    if(!found) structure.push_back(tryMoves[j]);
                    found = true;
                    lastYaxis = moveVector(lastDirection, lastYaxis, tryMoves[j]).second;
                }
            }
            if(!found) {



                //cout << "? at position " << i << "\t, arriving at: " << print(listPositions[i]) << endl;
                //cout << "oldXaxis: " << print(lastDirection) << ", oldYaxis " << print(lastYaxis) << ", newDir " << print(newDirection) << ", norms "
                //     << norm(lastDirection) << "," << norm(lastYaxis) << "," << norm(newDirection) << endl;
                //cout << "Starting a new absolute structure from this position" << endl;

                if(true){
                    // Need to initiate a new structure, but with absolute coordinates now.
                    structures.push_back(structure);
                    positions.push_back(position);
                    structure = string("");

                    // project the new position:
                    double distance = norm(initXAxis);
                    vector3D xAxis = initXAxis;
                    vector3D yAxis = initYAxis;
                    if(fabs(norm(xAxis) - norm(yAxis)) > 0.001 * norm(xAxis)) cerr << "ERR: x and y axis with different sizes" << endl;
                    vector3D zAxis = xAxis ^ yAxis;
                    zAxis = (1. / norm(xAxis)) * zAxis;

                    //cout << "Init Axis: X" << print(xAxis) << ", Y" << print(yAxis) << ", Z" << print(zAxis) << endl;
                    double x = (1. / norm(xAxis)) * (direction(listPositions[0], listPositions[i]) * xAxis) / distance;
                    double y = (1. / norm(yAxis)) * (direction(listPositions[0], listPositions[i]) * yAxis) / distance;
                    double z = (1. / norm(zAxis)) * (direction(listPositions[0], listPositions[i]) * zAxis) / distance;
                    // In theory those positions should be integers

                    int xi = static_cast<int>(floor(x + 0.499) );
                    int yi = static_cast<int>(floor(y + 0.499) );
                    int zi = static_cast<int>(floor(z + 0.499) );

                    if((fabs(x - (static_cast<double>(xi))) > 0.01) ||
                       (fabs(y - (static_cast<double>(yi))) > 0.01) ||
                       (fabs(z - (static_cast<double>(zi))) > 0.01)) cout << "ERR: coordinates of jumping point are not integers: x=" << x << "/" << xi << ", y=" << y << "/" << yi << ", z=" << z << "/" << zi << "\n It means the compartments of the proteins are not on the same lattice?" << endl;

                    int x0 = lattice::positionFromID(lattice::centralPosition())[0];
                    int y0 = lattice::positionFromID(lattice::centralPosition())[1];
                    int z0 = lattice::positionFromID(lattice::centralPosition())[2];
                    position = lattice::idFromPosisition(x0 + xi, y0 + yi, z0 + zi);

                    if(i == listPositions.size()-1){
                        //cerr << "WRN: single AA at the end in the structure." << endl;
                        structure.push_back('?');                    }
                    else {
                        i = i + 1;
                        // lastDirection and lastYAxis are the same, just changed the reference start, and this time can start with B.
                        for(unsigned int j = 0; j < nMovesAbs; ++j){
                            if(pointsAreEquall(listPositions[i], listPositions[i-1] + moveVector(initXAxis, initYAxis, tryMoves[j]).first)){
                                //if(!found) cout << "Yeah! Start move is " << tryMoves[j] << endl;
                                if(found) cout << "Euh... two moves give same position??" << endl;
                                if(!found) structure.push_back(tryMoves[j]);
                                found = true;
                                //lastDirection = moveVector(lastDirection, lastYaxis, tryMoves[j]).first;
                                //lastYaxis = moveVector(lastDirection, lastYaxis, tryMoves[j]).second;
                                newDirection = moveVector(initXAxis, initYAxis, tryMoves[j]).first;
                                lastYaxis = moveVector(initXAxis, initYAxis, tryMoves[j]).second;
                                IDeachStartingResidue.push_back(IDresidues[i]);
                            }
                        }
                        if(!found) {
                            cout << "ERR: A jump is followed by points that are not at integer positions according to the same grid" << endl;
                            cout << "   ... from " << print(listPositions[i-1]) << " to " << print(listPositions[i]) << endl;
                        }
                    }
                }
            }
        }

        if(!firstPointFound){
            if(norm(newDirection ^ lastDirection) < 0.1* norm(lastDirection)){ // shit loads of imprecisions in latfit output
                structure.push_back('S');
            } else {
                IDfirstNonAlignedPoint = i;
                firstPointFound = true;
                lastYaxis = (+1.0 / norm(lastDirection)) * (newDirection ^ lastDirection); // first move will be SU, the Y axis is conserved by Up, and it happens to be new ^ last
                structure.push_back('U');
                // saves the coordinates of the protein, for jumps
                initXAxis = lastDirection;
                initYAxis = lastYaxis;
                //cout << "Turn strength = " << norm(newDirection ^ lastDirection) << ", Found starting Y" << print(lastYaxis) << endl;


                // check the full structure to be on this grid
                //cout << "Check the position of all points on the grid made by this X,Y starting axis" << endl;
                double distance = norm(initXAxis);
                vector3D xAxis = initXAxis;
                vector3D yAxis = initYAxis;
                if(fabs(norm(xAxis) - norm(yAxis)) > 0.001 * norm(xAxis)) cerr << "ERR: x and y axis with different sizes" << endl;
                vector3D zAxis = xAxis ^ yAxis;
                zAxis = (1. / norm(xAxis)) * zAxis;

                //if(listOriginPositions.size() != listPositions.size()) cerr << "ERR: Original positions of size " << listOriginPositions.size() << " instead of " << listPositions.size() << endl;
                for(unsigned int j = 0; j < listPositions.size(); ++j){
                    double x = (1. / norm(xAxis)) * (direction(listPositions[0], listPositions[j]) * xAxis) / distance;
                    double y = (1. / norm(yAxis)) * (direction(listPositions[0], listPositions[j]) * yAxis) / distance;
                    double z = (1. / norm(zAxis)) * (direction(listPositions[0], listPositions[j]) * zAxis) / distance;

                    int xi = static_cast<int>(floor(x + 0.499) );
                    int yi = static_cast<int>(floor(y + 0.499) );
                    int zi = static_cast<int>(floor(z + 0.499) );

                    if((fabs(x - (static_cast<double>(xi))) > 0.01) ||
                       (fabs(y - (static_cast<double>(yi))) > 0.01) ||
                       (fabs(z - (static_cast<double>(zi))) > 0.01)) cout << "ERR: coordinates of jumping point are not integers: x=" << x << "/" << xi << ", y=" << y << "/" << yi << ", z=" << z << "/" << zi << "\n It means the compartments of the proteins are not on the same lattice?" << endl;

                    //vector3D intPos = {(double)xi, (double)yi, (double)zi};
                    vector3D intPos = {(double)x, (double)y, (double)z};
                    listIntPositions.push_back(intPos);

                    if(j < listOriginPositions.size()){
                        x = (1. / norm(xAxis)) * (direction(listPositions[0], listOriginPositions[j]) * xAxis) / distance;
                        y = (1. / norm(yAxis)) * (direction(listPositions[0], listOriginPositions[j]) * yAxis) / distance;
                        z = (1. / norm(zAxis)) * (direction(listPositions[0], listOriginPositions[j]) * zAxis) / distance;

                        vector3D doublePos = {x, y, z};
                        //cout << "xyz origin " << x << "," << y << "," << z << endl;
                        listRotatedOriginPositions.push_back(doublePos);
                    }
                }
            }
        }

        lastDirection = newDirection;
    }

    structures.push_back(structure);
    positions.push_back(position);

    //cout << "first turning point is " << IDfirstNonAlignedPoint << endl;
}

superProtein* latFitToLattice::asSuperProtein(){

    superProtein* inLattice = nullptr;                                           //new superProtein("",-1);
    for(unsigned int i = 0; i < structures.size(); ++i){

        struct3D toAdd = struct3D(structures[i], UnDefined,  positions[i]);
        if(i == 0){ inLattice = new superProtein(toAdd);
           //cout << "Starting protein: " << endl;
           //cout << print(*inLattice) << endl;
        }
        else {inLattice = new superProtein(insert(*inLattice, toAdd, IDeachStartingResidue[i]));
            //cout << "Adding the structure: " << endl;
            //cout << print(toAdd) << endl;
            //cout << "Resulting protein: " << endl;
            //cout << print(*inLattice) << endl;
        }
        //cout << "\t" << a.positions[i] << "\t" << a.structures[i] << endl;
    }

    if(!inLattice){
        cerr << "ERR: latFitToLattice::asSuperProtein(), Could not transform discretized PDB structure into the lattice " << endl;
        return inLattice;
    }

    inLattice->setAAs(sequence);

    return inLattice;
}



size_t findAAFrolListPos(vector< vector3D >& list, vector3D positionAA){
    if(positionAA.size() != 3) {cerr << "ERR: findAAFrolListPos, input position has wrong size " << positionAA.size() << endl; return -1;}
    vector<int> res;
    size_t N = list.size();
    for(size_t i = 1; i < N; ++i){
        if(list[i].size() != 3) cerr << "ERR: findAAFrolListPos, the list of positions has a position of size " << list[i].size() << " at index " << i << " instead of 3" << endl;
        if(fabs(list[i][0] - positionAA[0]) + fabs(list[i][1] - positionAA[1]) + fabs(list[i][2] - positionAA[2]) < 0.01){
            cout << "Found AA at ID " << i << endl;
            res.push_back(i);
        }
    }
    if(res.size() > 1) cerr << "ERR: found the same position on multiple AAs (?) at indices " << printVector(res) << endl;
    if(res.size() == 0) return -1;
    return res[0];
}



bool compForDist(std::pair<int, double> a, std::pair<int, double> b){
    return (a.second < b.second);
}

vector<std::pair<int, double>> sortPositionsByDistance(set<int> possiblePositions, vector<double> actualFloatingPosition){
    vector<std::pair<int, double> > toSort; // position - distance
    if(actualFloatingPosition.size() != 3){
        cerr << "ERR: sortPositionsByDistance, the actualFloatingPosition should have size 3" << endl;
        return toSort;
    }
    for(set<int>::iterator it = possiblePositions.begin(); it != possiblePositions.end(); ++it){
        vector3D doublePos = {(double) lattice::positionFromID(*it)[0],(double) lattice::positionFromID(*it)[1], (double) lattice::positionFromID(*it)[2]};
        toSort.push_back(std::pair<int, double>(*it, sqrt(norm(direction(doublePos, actualFloatingPosition)))));
    }
    sort(toSort.begin(), toSort.end(), compForDist);
//    for(size_t i = 0; i < toSort.size(); ++i){
//        res.push_back(toSort[i].first);
//    }
    return toSort;
}



vector3D toolPDBtoLattice(vector3D realPos, vector3D initXAxis, vector3D initYAxis, vector3D OriginPos){
    if(realPos.size() != 3) {cerr << "PDBtoLattice received vector of size " << realPos.size() << endl; return vector3D();}
    double distance = norm(initXAxis);
    vector3D zAxis = initXAxis ^ initYAxis;
    zAxis = (1. / norm(initXAxis)) * zAxis;

    double x = (1. / norm(initXAxis)) * (direction(OriginPos, realPos) * initXAxis) / distance;
    double y = (1. / norm(initYAxis)) * (direction(OriginPos, realPos) * initYAxis) / distance;
    double z = (1. / norm(zAxis)) *     (direction(OriginPos, realPos) * zAxis) / distance;

    vector3D res = {x,y,z};
    return res;
}

vector< vector3D> pooledPDBtoLattice(vector<vector3D> &toTransform, vector3D initXAxis, vector3D initYAxis, vector3D OriginPos){
    size_t N = toTransform.size();
    vector<vector3D> res;
    for(size_t i = 0; i < toTransform.size(); ++i){
        vector3D tf = toolPDBtoLattice(toTransform[i], initXAxis, initYAxis, OriginPos);
        res.push_back(tf);
    }
    return res;
}

vector3D latFitToLattice::PDBtoLattice(vector3D realPos){
    return toolPDBtoLattice(realPos, initXAxis, initYAxis, listPositions[0]);
}

vector3D latFitToLattice::LatticetoPDB(vector3D intPos){
    if(intPos.size() != 3) {cerr << "LatticetoPDB received vector of size " << intPos.size() << endl; return vector3D();}
    double distance = norm(initXAxis);
    vector3D zAxis = initXAxis ^ initYAxis;
    zAxis = (1. / norm(initXAxis)) * zAxis;
    vector3D realPos = ((intPos[0] * initXAxis) + (intPos[1] * initYAxis) + (intPos[2] * zAxis)); // no need to multiply by distance
    return realPos;
}

void latFitToLattice::testConversions(){
    vector3D zAxis = initXAxis ^ initYAxis;
    zAxis = (1. / norm(initXAxis)) * zAxis;
    cout << "Ref: XAxis = " << print(initXAxis) << ", YAxis = " << print(initYAxis) << endl;
    cout << "The center is at position " << print(listPositions[0]) << endl;
    cout << "R=>L: XAxis becomes " << print(PDBtoLattice(listPositions[0] + initXAxis)) << endl;
    cout << "R=>L: YAxis becomes " << print(PDBtoLattice(listPositions[0] + initYAxis)) << endl;
    cout << "R=>L: ZAxis becomes " << print(PDBtoLattice(listPositions[0] + zAxis)) << endl;
    cout << "R=>L: 2*XAxis + 1*Zaxis becomes " << print(PDBtoLattice(listPositions[0] + (2*initXAxis) + zAxis)) << endl;

    cout << "L=>R: (1,0,0) becomes " << print(listPositions[0] + LatticetoPDB({1.,0,0})) << endl;
    cout << "L=>R: (1,0,0) becomes " << print(listPositions[0] + LatticetoPDB({0,1.,0})) << endl;
    cout << "L=>R: (1,0,0) becomes " << print(listPositions[0] + LatticetoPDB({0,0,1.})) << endl;
    cout << "L=>R: (0.9,1.2,-0.3) becomes " << print(listPositions[0] + LatticetoPDB({0.9,1.2,-0.3})) << endl;

    // note: we stay double so never rounded to integer, so can go back
    cout << "R=>L=>R: (0.9,1.2,-0.3) becomes " << print(LatticetoPDB(PDBtoLattice(listPositions[0] + vector3D({0.9,1.2,-0.3})))) << endl;
    cout << "L=>R=>L: (0.9,1.2,-0.3) becomes " << print(PDBtoLattice(listPositions[0] + LatticetoPDB({0.9,1.2,-0.3}))) << endl;
}


// possible types of Discretization: FuC, CoM, CA, or None
vector< vector3D > getPDBChainCoarseGrainedPositions(string PDB, string chains, string typeDiscretization){

    bool addNANsBetweenChains = true;

    vector< vector3D > res;
    vector< vector3D > cAlphas;
    vector<string> AAs; // this is just for debug

    // We will only include the full backbone in the case FuC, or None where we take all atoms
    bool includeBackboneInMassCenter = (!typeDiscretization.compare("FuC")) || (!typeDiscretization.compare("None"));
    bool onlyCA = !typeDiscretization.compare("CA");
    string prevResidueID = "?"; // We will include the insertion
    char prevChain = '?';

    ifstream f;
    string lineType;
    cout << "Transforming chains " << chains << " as coarse grained AAs in " << PDB << endl;
    f.open(PDB.c_str());
    if(!f) {
        cerr << "ERR: latFitToLattice::parseLatFitPDB, file not found:\n" << PDB << endl; return res;
    }

    biu::DPointVec sideChain;
    biu::DblPoint cAlpha(0.0,0.0,0.0);
    bool cAlphaFound = false;

    while(f.good()){
        char buf[10001];
        f.getline(buf, 10000);
        string full(buf);
        string lineType = full.substr(0, 5);



        if(!lineType.compare(string("ATOM "))){
            // http://www.wwpdb.org/documentation/file-format-content/format33/sect6.html#LINK

            string Atom = full.substr(13,3); // sometimes we can use 4 cahracters (12,4), don't know why. Here do like larfit with 3 characters
            string Residue = full.substr(17,3);
            char Chain = full.substr(21,1)[0];

            //cout << "Line " << lineType << " Atom |" << Atom << "| Chain |" << Chain << "| Residue |" << Residue << "|, Chains " << chains << endl;

            if(chains.find(Chain) != std::string::npos){

                string ResidueID = full.substr(22,5); // this includes the insertion tag
                //cout << "Line " << lineType << " Atom |" << Atom << "| Chain |" << Chain << "| Residue |" << Residue << "|, ID" << ResidueID << ", Chains " << chains << endl;

                string insertion = full.substr(26,1);
                stringstream ResDescr1;
                ResDescr1 << Chain << ResidueID; // << Residue;

                if ((ResidueID .compare(prevResidueID)) || (prevChain != Chain)){
                    if(cAlphaFound){
                        cAlphaFound = false;
                        biu::DblPoint pt = getCentroid(sideChain);
                        sideChain.clear();
                        res.push_back(vector3D({pt.getX(), pt.getY(), pt.getZ()}));
                        //cout << "Side Chain Finished" << printVector(vector3D({pt.getX(), pt.getY(), pt.getZ()})) << endl;

                        if((prevChain != Chain) && (addNANsBetweenChains))  res.push_back(vector3D({NAN, NAN, NAN}));
                    } else {
                        cout << "A Communist was here" << endl;
                        //it should be the first line I guess...
                    }
                }
                prevResidueID = ResidueID;
                prevChain = Chain;

                // get cAlpha atom
                if (!Atom.compare("CA ")) {
                    if(cAlphaFound) cerr << "ERR: found twice the C alpha, mayve chain changed??" << endl;
                    double x = atof(full.substr(30,8).c_str());
                    double y = atof(full.substr(38,8).c_str());
                    double z = atof(full.substr(46,8).c_str());
                    cAlpha = biu::DblPoint(x,y,z);
                    cAlphaFound = true;
                    cAlphas.push_back(vector3D({x,y,z}));
                    cout << "CA " << printVector(vector<double>({x,y,z})) << endl;
                }
                // only add atoms not belonging to the backbone to the side chain
                if ((!onlyCA) && (includeBackboneInMassCenter ||
                                  ( Atom.compare("CA ") != 0
                                    &&	Atom.compare("N  ") != 0
                                    &&	Atom.compare("C  ") != 0
                                    &&	Atom.compare("O  ") != 0)
                                  ))
                {
                    // coordinates
                    double x = atof(full.substr(30,8).c_str());
                    double y = atof(full.substr(38,8).c_str());
                    double z = atof(full.substr(46,8).c_str());
                    //cerr << "Residue " <<  aaName << " seqNum " << seqNum << " got atom " << line.substr(13,3) << endl;
                    sideChain.push_back(biu::DblPoint(x,y,z));
                    //cout << "Push Side Chain" << endl;
                }
            }
        }
    }

    if(sideChain.size() > 0){
        biu::DblPoint pt = getCentroid(sideChain);
        res.push_back(vector3D({pt.getX(), pt.getY(), pt.getZ()}));
        //cout << "Final side chain" << endl;
    }

    f.close();

    if(onlyCA) return cAlphas;
    else return res;
}
string AA_PDBName(char AA_name){
    switch(AA_name){
        case 'C': return "CYS";
        case 'M': return "MET";
        case 'F': return "PHE";
        case 'I': return "ILE";
        case 'L': return "LEU";
        case 'V': return "VAL";
        case 'W': return "TRP";
        case 'Y': return "TYR";
        case 'A': return "ALA";
        case 'G': return "GLY";
        case 'T': return "THR";
        case 'S': return "SER";
        case 'N': return "ASN";
        case 'Q': return "GLN";
        case 'D': return "ASP";
        case 'E': return "GLU";
        case 'H': return "HIS";
        case 'R': return "ARG";
        case 'K': return "LYS";
        case 'P': return "PRO";
        case '?': return "XXX";
        default: {return "XXX";}
    }
}

// only for positive numbers
string fixedDigitsInteger(size_t val, size_t nbChar){
    stringstream res;
    if((val < 10) && (nbChar > 0))  res << string(nbChar-1, ' ') << val;
    else if((val < 100) && (nbChar > 1)) res << string(nbChar-2, ' ') << val;
    else if((val < 1000) && (nbChar > 2)) res << string(nbChar-3, ' ') << val;
    else if (nbChar > 3) res << string(nbChar-4, ' ') << val;
    else {cerr << "ERR: fixedDigitsInteger " << val  << "," << nbChar << " doen't fit in the given number of characters" << endl; return string(nbChar, '?');}
    return res.str();
}

string fixedDigitsFloat(double val, size_t nbCharBeforeInclMinus, size_t decimals){
    stringstream toChop; toChop << val;
    if(toChop.str().find(".") == string::npos){
        toChop << "." << string('0', decimals);
    } else {
        toChop << string('0', decimals);
    }
    string longerString = string(' ', nbCharBeforeInclMinus) + toChop.str();
    size_t posComa = longerString.find(".");
    string res = longerString.substr(posComa - nbCharBeforeInclMinus, posComa + decimals);
    return res;
}

void testInOutLatticePDB(){

    // test 1:
    superProtein* P1 = new superProtein("BRUDLUSSS");
    P1 ->setAAs("KRTDLASSRP");
    superProtein* P2 = new superProtein("URUDLUSSS");
    P2 ->setAAs("PRIGLASSFW");

    string file1 = generatePDBfromLattice({3.0, 0.0, 0.0}, {0.0, 3.0, 0.0}, {50, 50, 50}, P1, P2);
    ofstream f1("ManualGenerate.txt"); f1 << file1; f1.close();

    // test 2:
    // read a PDB from latfit, then read the same PDB but with another chain, include it in visualization
    latFitToLattice a = latFitToLattice();
    a.parseLatFitPDB("C:/Users/pprobert/Desktop/Main/B-CurrentZapotect/Zapotec/build-Absolut-Desktop_Qt_5_12_5_MinGW_64_bit-Release/1ADQ_AdiscretizedFuC5.25.pdb");
    a.transform();
    superProtein* P3 = a.asSuperProtein();

    vector<vector<double> > positionsAB = getPDBChainCoarseGrainedPositions("C:/Users/pprobert/Desktop/Main/B-CurrentZapotect/Zapotec/build-Absolut-Desktop_Qt_5_12_5_MinGW_64_bit-Release/1ADQ.pdb", "LH", "FuC");
    vector<vector<double> >* transformed = new vector<vector<double> >(pooledPDBtoLattice(positionsAB, a.initXAxis, a.initYAxis, a.listPositions[0]));


#ifdef ALLOW_GRAPHICS
    glDisplay();
    addToDisplay(P3, false);
    addToDisplay(transformed);
    glutMainLoop();
#endif

}

string generatePDBfromLattice(vector3D xAxis, vector3D yAxis, vector3D originPos, superProtein* p1, superProtein* p2, superProtein* p3){
    int nChains = 1;
    if(p1 == nullptr) {cerr << "generatePDBfromLattice, null firts protein" << endl; return string("ERROR");}
    if(p2 != nullptr) nChains = 2;
    if(p3 != nullptr) nChains = 3;

    stringstream PDBout;
    PDBout << "HEADER    LATTICE PROTEIN STRUCTURE               01-JAN-D0   9ABC              "
              "TITLE     FIT OF THE PDB STRUCTURE 1ADQ ONTO A LATTICE                          "
              "COMPND    MOL_ID: 1;                                                            "
              "COMPND   2 MOLECULE: LATTICE PROTEIN;                                           ";
    if(nChains == 1) PDBout << "COMPND   3 CHAIN: A;                                                            ";
    if(nChains == 2) PDBout << "COMPND   3 CHAIN: A, B;                                                         ";
    if(nChains == 3) PDBout << "COMPND   3 CHAIN: A, B, C;                                                      ";
    PDBout << "COMPND   4 ENGINEERED: YES                                                      "
              "SOURCE    MOL_ID: 1                                                             "
              "KEYWDS    LATTICE EXPORT                                                        "
              "EXPDTA    YMIR THEORETICAL MODEL                                                "
              "AUTHOR    ABSOLUT SOFTWARE                                                      ";



    for(int nC = 0; nC < nChains; ++nC){
        string AAsChain;
        char chainID = '?';
        if(nC == 0){
            AAsChain = p1->getAAseq();
            chainID = 'A';
        }
        if(nC == 1){
            AAsChain = p2->getAAseq();
            chainID = 'B';
        }
        if(nC == 2){
            AAsChain = p3->getAAseq();
            chainID = 'B';
        }

        // formatting in exactly 4 characters


        size_t j = 0;
        for(size_t i = 0; i < (AAsChain.size()/13) + 1; ++i){
            // formatting the line ID from 1 to X with exactly the number of characters
            PDBout << "SEQRES" << fixedDigitsInteger(i+1, 4) << " " << chainID << " " << fixedDigitsInteger(AAsChain.size(),4) << "  ";

            for(int k = 0; k < 13; ++k){
                if(j < AAsChain.size()){
                    PDBout << AA_PDBName(AAsChain[j]) << " ";
                } else {
                    PDBout << "    ";
                }
                j++;
            }
            PDBout << "\n";
        }
    }

    // here, we need to transform back
    for(int nC = 0; nC < nChains; ++nC){

        string AAsChain;
        char chainID = '?';
        superProtein* currentProt = nullptr;

        if(nC == 0){
            AAsChain = p1->getAAseq();
            chainID = 'A';
            currentProt = p1;
        }
        if(nC == 1){
            AAsChain = p2->getAAseq();
            chainID = 'B';
            currentProt = p2;
        }
        if(nC == 2){
            AAsChain = p3->getAAseq();
            chainID = 'B';
            currentProt = p3;
        }


        if(p1->size() != AAsChain.size()) cerr << "ERR: inconsistent AAseq size and AAs size inside the superProtein" << endl;
        for(size_t i = 0; i < static_cast<size_t>(p1->size()); ++i){
            int posLattice = p1->points[i].IDposition;
            vector<int> intPos = lattice::positionFromID(posLattice);
            vector3D doublePos = {static_cast<double>(intPos[0]), static_cast<double>(intPos[1]), static_cast<double>(intPos[2])};
            // Heum, this is probably the other function, latticeToPDB...
            vector3D realPos = toolPDBtoLattice(doublePos, xAxis, yAxis, originPos);
            PDBout << "HETATM" << fixedDigitsInteger(i+1,5) << "  CA  " << AA_PDBName(AAsChain[i]) << " " << chainID << " " << fixedDigitsInteger(i+1,5) << "    " << fixedDigitsFloat(realPos[0], 4, 3) <<  fixedDigitsFloat(realPos[1], 4, 3) << fixedDigitsFloat(realPos[2], 4, 3) << "  1.00  0.00           C  " << endl;
        }
    }
    //HETATM    1  CA  PRO L 238       4.069 -25.860   5.700  1.00  0.00           C
    //HETATM    2  CA  SER L 239       7.694 -22.367   7.191  1.00  0.00           C
    //HETATM    3  CA  VAL L 240       4.945 -18.532   4.890  1.00  0.00           C
    return PDBout.str();
}

// Now I need:
// - new data structure for vector<vector3D> (more than one?) - option to hide
// - test: take PDB, discretize and go back to PDB
// - test: take same PDB as latfit, and check if we have the same center of masses.


set<int> latFitToLattice::addGlycans(readGlycans& glycansFromOriginalPDB, set<int>& occupiedPositions){

    set<int> atmosphere = getAtmosphere(occupiedPositions);

    //cout << "Got atmosphere: " << print(atmosphere) << endl;
    vector< set<int> > possiblePositionsEachGlycan;
    vector<int> shiftAsInt = lattice::positionFromID(lattice::centralPosition());
    vector<double> shift = vector<double>(shiftAsInt.begin(), shiftAsInt.end());
    // So, the original 'glycans' contains the


//        readGlycans(string filename, string chain);
//        size_t nGlycans;
//        vector< vector3D > listPositionGlycans;
//        vector< std::pair<char, int> > AAboundByGlycans;
//        vector< vector<double> > PositionsAAboundByGlycans;
//        string print();

    size_t N = glycansFromOriginalPDB.nGlycans;
    if(glycansFromOriginalPDB.listPositionGlycans.size() != N) cerr << "ERR: latFitToLattice::addGlycans, inconsistent size of listPositionGlycans" << endl;

    double distance = norm(initXAxis);
    if(fabs(norm(initXAxis) - norm(initYAxis)) > 0.001 * norm(initYAxis)) cerr << "ERR: x and y axis with different sizes" << endl;
    vector3D zAxis = initXAxis ^ initYAxis;
    zAxis = (1. / norm(initXAxis)) * zAxis;


    //glycansFromOriginalPDB
    vector< vector< std::pair<int, double > > > rankedPositionsForEachGlycan;

    for(size_t i = 0; i < N; ++i){

        double x = (1. / norm(initXAxis)) * (direction(listPositions[0], glycansFromOriginalPDB.listPositionGlycans[i]) * initXAxis) / distance;
        double y = (1. / norm(initYAxis)) * (direction(listPositions[0], glycansFromOriginalPDB.listPositionGlycans[i]) * initYAxis) / distance;
        double z = (1. / norm(zAxis)) *     (direction(listPositions[0], glycansFromOriginalPDB.listPositionGlycans[i]) * zAxis) / distance;

        cout << endl;
        cout << "Glycan" << i+1 << "/" << N << ": real " << glycansFromOriginalPDB.listPositionGlycans[i][0] << "," << glycansFromOriginalPDB.listPositionGlycans[i][1] << "," << glycansFromOriginalPDB.listPositionGlycans[i][2] << " => rescaled at " << x + shift[0] << " " << y + shift[1] << " " << z  + shift[2] << endl;

        size_t pos_aa = findAAFrolListPos(listOriginPositions, glycansFromOriginalPDB.PositionsAAboundByGlycans[i]);
        if(pos_aa > 0){
            double xAA = (1. / norm(initXAxis)) * (direction(listPositions[0], listPositions[pos_aa]) * initXAxis) / distance;
            double yAA = (1. / norm(initYAxis)) * (direction(listPositions[0], listPositions[pos_aa]) * initYAxis) / distance;
            double zAA = (1. / norm(zAxis)) *     (direction(listPositions[0], listPositions[pos_aa]) * zAxis) / distance;

            // DANGER: Make sure that (int) XX is always getting a positive number, if not it becomes the rounding opposite => Adds shift before doing (int). Shound use floor or so...
            vector<int> posAAinLattice = {(int) (xAA+0.1 + shift[0]), (int) (yAA+0.1 + shift[1]), (int) (zAA+0.1 + shift[2])};
            cout << "   Found matching AA, original pos " << printVector(listOriginPositions[pos_aa]) << ", afterr latfit " << printVector(listPositions[pos_aa]) << " and rescaled " << xAA + shift[0] << "," << yAA + shift[1] << "," << zAA + shift[2] << " => " << printVector(posAAinLattice) << endl;
            int pos = lattice::idFromPosisition(posAAinLattice);
            vector<int> possible = lattice::idNeighbors(pos);
            set<int> remaining;
            cout << "   ID possible neighbors are: " << endl;

            // Here we check the contraints where a glycan should be:
            // 1: Not inside an occupied position of the protein
            // 2: - AND - should belong to the atmosphere

            for(vector<int>::iterator it = possible.begin(); it != possible.end(); ++it){
                cout << "  -> " << *it << ", " << printVector(lattice::positionFromID(*it)) << endl;
                if(occupiedPositions.find(*it) == occupiedPositions.end()){
                    if(atmosphere.find(*it) != atmosphere.end()){
                        remaining.insert(*it);
                    }
//                    else {
//                        set<int> neighPos = vectorToSet(lattice::idNeighbors(*it));
//                        if(intersection_sets(atmosphere,neighPos).size() > 0){
//                            remaining.insert(*it);
//                            cout << "Exception for position " << *it << endl;
//                        }
//                    }
                }
            }
            if(remaining.size() == 0){
                cout << "WRN: The AA that is bound by this glycan is NOT A SURFACE AA => Taking the neighbors of the glycan instead!!" << endl;
                vector<int> posGlyinLattice = {(int) (x+0.1 + shift[0]), (int) (y+0.1 + shift[1]), (int) (z+0.1 + shift[2])};
                vector<int> possible2 = lattice::idNeighbors(lattice::idFromPosisition(posGlyinLattice));
                ////// NOOOOOOOO!! Here we need to also include the glycan => DOne now, need to redo the AGS with glycan
                possible2.push_back(lattice::idFromPosisition(posGlyinLattice));

                // Here we check the contraints where a glycan should be:
                // 1: Not inside an occupied position of the protein
                // 2: - AND - should belong to the atmosphere
                // 3: - OR - could be a forbidden point BUT touches the atmosphere => distance one from atmosphere.
                for(vector<int>::iterator it = possible2.begin(); it != possible2.end(); ++it){
                    cout << "  -> " << *it << ", " << printVector(lattice::positionFromID(*it)) << endl;
                    if(occupiedPositions.find(*it) == occupiedPositions.end()){
                        if(atmosphere.find(*it) != atmosphere.end()){
                            remaining.insert(*it);
                        }
                    }
                }
            }
            if(remaining.size() == 0){
                cout << "WRN: direct neighbors of the glycan are not free => Taking nonnearest neighbors!!" << endl;
                vector<int> posGlyinLattice = {(int) (x+0.1 + shift[0]), (int) (y+0.1 + shift[1]), (int) (z+0.1 + shift[2])};
                vector<int> possible3 = lattice::idLargeNeighbors(lattice::idFromPosisition(posGlyinLattice));

                for(vector<int>::iterator it = possible3.begin(); it != possible3.end(); ++it){
                    cout << "  -> " << *it << ", " << printVector(lattice::positionFromID(*it)) << endl;
                    if(occupiedPositions.find(*it) == occupiedPositions.end()){
                        if(atmosphere.find(*it) != atmosphere.end()){
                            remaining.insert(*it);
                        }
                    }
                }
            }
            cout << "   => ID remaining free positions for this glycan are: " << print(remaining) << endl;
            possiblePositionsEachGlycan.push_back(remaining);
            rankedPositionsForEachGlycan.push_back(sortPositionsByDistance(remaining,{x+shift[0],y+shift[0],z+shift[0]}));
        } else {
            possiblePositionsEachGlycan.push_back(set<int>());
            rankedPositionsForEachGlycan.push_back(sortPositionsByDistance(set<int>(),{x+shift[0],y+shift[0],z+shift[0]}));
        }
    }

    cout << "List of priorities for each glycan" << endl;
    for(size_t i = 0; i < rankedPositionsForEachGlycan.size(); ++i){
        vector<double> posReal = glycansFromOriginalPDB.listPositionGlycans[i];
        cout << "Glycan " << i << "\t";
        for(size_t j = 0; j < rankedPositionsForEachGlycan[i].size(); ++j){
               cout << "\t" << rankedPositionsForEachGlycan[i][j].first << "->" << printVector(lattice::positionFromID(rankedPositionsForEachGlycan[i][j].first)) <<",d=" << rankedPositionsForEachGlycan[i][j].second;
        }
        cout << endl;
    }

    size_t NG = rankedPositionsForEachGlycan.size();




    vector<int> posPerGlycan(NG, -1);
    size_t posFound = 0;
    int cpt = 0;
    while((posFound < NG) && (cpt < 10000)){
        for(size_t i = 0; i < NG; ++i){

            // completely suboptimal to put this part inside the loop...
            // possiblePositionsEachGlycan
            // makes the list of first positions
            std::map<int, bool> positionsRequestedTwice;
            // 1: to be in the list = at least one glycan => true, to be twice => False
            for(size_t i = 0; i < NG; ++i){
                if(rankedPositionsForEachGlycan[i].size() > 0){ // should not be empty, just avoid sgfault
                    // just for output
//                    cout << "GLycan " << i;
//                    for(size_t j = 0; j < rankedPositionsForEachGlycan[i].size(); ++j){
//                           cout << "\t" << rankedPositionsForEachGlycan[i][j].first << ",d=" << rankedPositionsForEachGlycan[i][j].second;
//                    }
//                    cout << endl;

                    int bestPos = rankedPositionsForEachGlycan[i][0].first;
                    if(positionsRequestedTwice.find(bestPos) == positionsRequestedTwice.end()){
                        positionsRequestedTwice[bestPos] = false;
                    } else {
                        positionsRequestedTwice[bestPos] = true;
                    }
                } else {
                    // Once a glycan has been attributed, it will habe no possible positions
                    //cerr << "ERR: Glycan " << i << " Has no possible positions!!" << endl;
                }
            }


            if(rankedPositionsForEachGlycan[i].size() > 0){ // empty means the glycan has been solved
                int bestPos = rankedPositionsForEachGlycan[i][0].first;
                std::map<int,bool>::iterator isFree = positionsRequestedTwice.find(bestPos);
                if(isFree == positionsRequestedTwice.end()){cerr << "ERR: BestPositions are not all in the 'positionsRequestedTwice" << endl;}
                else{
                    // If the position is only claimed by this glycan, gives it
                    if(isFree->second){
                        rankedPositionsForEachGlycan[i].clear();
                        posPerGlycan[i] = bestPos;
                        posFound++;
                        cout << "Glycan " << i << " attributed to position " << bestPos << endl;
                    } else {
                        // Now there are other conflicting glycans wanting the same position.
                        // A glycan with only this possibility wins. If all glycans have two choices at least,
                        // the glycan with minimum distance wins. Then, this position is removed as possible for all glycans
                        size_t bestGlycan = i;
                        double bestDist = rankedPositionsForEachGlycan[i][0].second;
                        // Now solving the problem for all glycans requesting this position. So j > i.
                        for(size_t j = i+1; j < NG; ++j){
                            if(rankedPositionsForEachGlycan[j][0].first == bestPos){
                                if(rankedPositionsForEachGlycan[j].size() == 0){
                                    bestDist = -1; // then he wins,
                                    bestGlycan = j;
                                }
                                if(rankedPositionsForEachGlycan[j][0].second < bestDist){
                                    bestDist = rankedPositionsForEachGlycan[j][0].second;
                                    bestGlycan = j;
                                }
                            }
                        }
                        // The glycan bestGlycan wins,
                        rankedPositionsForEachGlycan[bestGlycan].clear();
                        posPerGlycan[bestGlycan] = bestPos;
                        posFound++;
                        cout << "Glycan " << bestGlycan << " attributed to position " << bestPos << endl;

                        // Now removing this position for all glycans.
                        for(size_t j = i+1; j < NG; ++j){
                            // use iterator to remove elements easily
                            for(vector<std::pair<int,double>>::iterator itk = rankedPositionsForEachGlycan[j].begin(); itk != rankedPositionsForEachGlycan[j].end(); ++itk){
                                if(itk->first == bestPos){
                                    itk = rankedPositionsForEachGlycan[j].erase(itk);
                                    break;
                                }
                                // Note, in theory we should do itk-- because the new elements is at the same (old position)
                                // but anyways the same number should not happen twice so we don't care;
                            }
                        }

                        // and updating the list of free positions is done each time


                        // Might be conflicts, like if two glycans want only the same position. In this case one less position will
                        // be blocked... hopefully it doesn't happen too much, or we could extend to the second order neighbors, just a bit
                        // more complex...

                    }
                }
            }
        }
    }
    cout << "Found " << posPerGlycan.size() << " positions for " << NG << " Glycans " << endl;
    return vectorToSet(posPerGlycan);

}

char identifyAA(string typeResidue){
    static std::map<string, char> AAchar;
    static bool loaded = false;
    if(!loaded){
        AAchar["CYS"] = 'C';
        AAchar["MET"] = 'M';
        AAchar["PHE"] = 'F';
        AAchar["ILE"] = 'I';
        AAchar["LEU"] = 'L';
        AAchar["VAL"] = 'V';
        AAchar["TRP"] = 'W';
        AAchar["TYR"] = 'Y';
        AAchar["ALA"] = 'A';
        AAchar["GLY"] = 'G';
        AAchar["THR"] = 'T';
        AAchar["SER"] = 'S';
        AAchar["ASN"] = 'N';
        AAchar["GLN"] = 'Q';
        AAchar["ASP"] = 'D';
        AAchar["GLU"] = 'E';
        AAchar["HIS"] = 'H';
        AAchar["ARG"] = 'R';
        AAchar["LYS"] = 'K';
        AAchar["PRO"] = 'P';
        loaded = true;
    }
    std::map<string,char>::iterator it = AAchar.find(typeResidue);
    if(it == AAchar.end()) {return '?';}
    else return it->second;
}




readGlycans::readGlycans(string filename, string chainFilter){
    // Parse lines LINK for glycans will store: Chain+IDAA+AA in a string, and match the found glycan ID
    // example: LINK         ND2 ASN A  88                 C1  NAG A1016     1555   1555  1.44
    //          will become  "A88ASN" => "A1016"
    std::map<string, string> contactsNAG_NAN;

    // Parse HETATM and stores the position of each of them by name:
    // HETATM17723  C1  NAG A1016     173.605 169.074 172.021  1.00125.52           C
    // becomes  "A1016" => {173.605, 169.074, 172.021};
    std::map<string, vector<double>> posHETATM;

    //std::map<string, vector<double> > positionsAllCAs;

    ifstream f;
    string lineType;
    cout << "Looking for glycans in " << filename << endl;
    f.open(filename.c_str());
    if(!f) {
        cerr << "ERR: latFitToLattice::parseLatFitPDB, file not found:\n" << filename << endl; return;
    }
    while(f.good()){
        char buf[10001];
        f.getline(buf, 10000);
        string full(buf);
        string lineType = full.substr(0, 5);

        //cout << "Line " << lineType << endl;
        if(!lineType.compare(string("LINK "))){
            // http://www.wwpdb.org/documentation/file-format-content/format33/sect6.html#LINK

            string Atom = full.substr(12,3);
            string Residue = full.substr(17,3);
            string Chain = full.substr(21,1);
            string ResidueID = full.substr(22,4);
                stringstream transform1(ResidueID);
                int resIDint;
                transform1 >> resIDint;
            string insertion = full.substr(26,1);
                stringstream ResDescr1;
                ResDescr1 << Chain << resIDint; // << Residue;

            string Atom2 = full.substr(42,3);
            string Residue2 = full.substr(47,3);
            string Chain2 = full.substr(51,1);
            string ResidueID2 = full.substr(52,4);
                stringstream transform2(ResidueID2);
                int resIDint2;
                transform2 >> resIDint2;
            string insertion2 = full.substr(56,1);
            string linkDist = full.substr(73,5);
                stringstream ResDescr2;
                ResDescr2 << Chain2 << resIDint2; // << Residue2;


            // will be encoded as: AA+Residue->GlyID, as:  "A88ASN" => "A1016"
            // Case 1: first residue is an AA
            char getAA1 = identifyAA(Residue);
            if((getAA1 != '?') && (chainFilter.find(Chain[0]) != std::string::npos) && ((!Residue2.compare(string("NAN"))) || (!Residue2.compare(string("NAG"))))){
                //cout << "Identified AA linked to glycan at: " << Chain << resIDint << Residue << " and towards " << Chain2 << resIDint2 << Residue2 <<  endl;
                if(contactsNAG_NAN.find(ResDescr1.str()) == contactsNAG_NAN.end()){
                    contactsNAG_NAN.insert(std::pair<string, string>(ResDescr1.str() + Residue, ResDescr2.str()));
                } else {
                    cerr << "ERR: The PDB file has more than one glycan for an AA, see line LINK for ID " << ResDescr1.str() << endl;
                }
            }

            // Case 2: second residue is an AA
            char getAA2 = identifyAA(Residue2);
            if((getAA2 != '?') && (chainFilter.find(Chain2[0]) != std::string::npos) && ((!Residue.compare(string("NAN"))) || (!Residue.compare(string("NAG"))))){
                //cout << "Identified AA linked to glycan at: " << Chain2 << resIDint2 << Residue2 << " and towards " << Chain << resIDint << Residue <<   endl;
                if(contactsNAG_NAN.find(ResDescr2.str()) == contactsNAG_NAN.end()){
                    contactsNAG_NAN.insert(std::pair<string, string>(ResDescr2.str() + Residue2, ResDescr1.str()));
                } else {
                    cerr << "ERR: The PDB file has more than one glycan for an AA, see line LINK for ID " << ResDescr2.str() << endl;
                }
            }
        }

        // In case we use pdb-tools to remove insertions, the ID of residues will change, so we need to also save the XYZ position of residues
        // of interest, to find later their ID in the modified PDB
        // will store: string(Chain-IDres-AA) to position(X,Y,Z), ex: A88ASN=> {10.2, 15.3, 8.8}
        if(!lineType.compare(string("ATOM "))){

            string serial = full.substr(6, 5);
            string Atom = full.substr(12, 4);
            string Residue = full.substr(17,3);
            string Chain = full.substr(21,1);
            string ResidueID = full.substr(22,4);
                stringstream transform1(ResidueID);
                int resIDint;
                transform1 >> resIDint;
            string insertion = full.substr(26,1);

            string posX = full.substr(30,8);
            string posY = full.substr(38,8);
            string posZ = full.substr(46,8);
                stringstream PX(posX);
                stringstream PY(posY);
                stringstream PZ(posZ);
                double posXdouble = 0, posYdouble = 0, posZdouble = 0;
                PX >> posXdouble;
                PY >> posYdouble;
                PZ >> posZdouble;

            //cout << Atom << endl;
            //if(!Atom.compare(" CA ")){
                if(!insertion.compare(" ")){
                    stringstream thisAA; thisAA << Chain << resIDint << Residue;
                    //cout << thisAA.str() << " => " << posXdouble << " " << posYdouble << " " << posZdouble << endl;
                    std::map<string, vector<double> >::iterator it = positionsAllCAs.find(thisAA.str());
                    if(it == positionsAllCAs.end()){
                           positionsAllCAs[thisAA.str()] = {posXdouble, posYdouble, posZdouble};
                    } else {
                           it->second.push_back(posXdouble);
                           it->second.push_back(posYdouble);
                           it->second.push_back(posZdouble);
                        //cerr << "ERR: readGlycans(), In this PDB (" << filename << "), an AA is defined twice?? (Residue " << thisAA.str() << endl;
                    }
                }
                // we do not accept glycans binding to insertions, this is too tedious for now!!
            //}
        }



        // Stupid PDB uses different nb of chars for the line type... so HETATM => HETAT with the 5 first chars.
        if(!lineType.compare(string("HETAT"))){
            // http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#HETATM

            string serial = full.substr(6, 5);
            string Atom = full.substr(12, 3);
            string Residue = full.substr(17,3);
            string Chain = full.substr(21,1);
            string ResidueID = full.substr(22,4);
                stringstream transform1(ResidueID);
                int resIDint;
                transform1 >> resIDint;
            string insertion = full.substr(26,1);

            string posX = full.substr(30,8);
            string posY = full.substr(38,8);
            string posZ = full.substr(46,8);
                stringstream PX(posX);
                stringstream PY(posY);
                stringstream PZ(posZ);
                double posXdouble = 0, posYdouble = 0, posZdouble = 0;
                PX >> posXdouble;
                PY >> posYdouble;
                PZ >> posZdouble;

           if((!Residue.compare("NAG")) || (!Residue.compare("NAN"))){
               stringstream ResDescr1;
               ResDescr1 << Chain << resIDint; // << Residue;

               std::map<string, vector<double>>::iterator it = posHETATM.find(ResDescr1.str());
               if(it == posHETATM.end()){
                   vector<double> posVec = {posXdouble, posYdouble, posZdouble};
                   posHETATM.insert(std::pair<string, vector<double>>(ResDescr1.str(), posVec)) ;
               } else {
                   it->second.push_back(posXdouble);
                   it->second.push_back(posYdouble);
                   it->second.push_back(posZdouble);
                   //cerr << "ERR: The PDB file has two times the line HETATM for ID " << ResDescr1.str() << endl;
               }
               //cout << ResDescr1.str() << ", Got position: " << posXdouble << "," << posYdouble << "," << posZdouble << endl;
           }
        }
    }
    //cout << "\n\n" << listPositions.size() << endl;
    //cout << "End of file: have read " << listPositions.size() << " Discretized positions, and " << listOriginPositions.size() << " original positions" << endl;
    f.close();


    // transform the list of AA positions into Fused Center

    std::map<string, vector<double> >::iterator it4;
    for(it4 = positionsAllCAs.begin(); it4 != positionsAllCAs.end(); ++it4){

        vector<double> positions = it4->second;
        if((positions.size() % 3) != 0){
            cerr << "ERR: atom positions for " << it4->first << " is not a multiple of 3" << endl;
        }
        if(positions.size() == 0) cerr << "ERR: Empty Residue positions ??" << endl;

        biu::DPointVec centralPosition;
        for(size_t i = 0; i < positions.size()/3; ++i){
            biu::DblPoint oneAtom(positions[3*i], positions[3*i+1], positions[3*i+2]);
            centralPosition.push_back(oneAtom);
        }
        biu::DblPoint center = getCentroid(centralPosition);

        it4->second = {center.getX(), center.getY(), center.getZ()};
    }


    // Now retrieving the glycans one by one
    // std::map<string, string> contactsNAG_NAN;
    nGlycans = contactsNAG_NAN.size();
    //listPositionGlycans.resize(nGlycans);
    for(std::map<string, string>::iterator it = contactsNAG_NAN.begin(); it != contactsNAG_NAN.end(); ++it){
        cout << "AA" << it->first << " interacts with " << it->second << "\t";

        // I am using the latfit vector structures to use the exact same function to get the center of masses (getCentroid). I could also have done it manually...
        biu::DPointVec glycan;

        std::map<string, vector<double> >::iterator it2 = posHETATM.find(it->second);
        if(it2 == posHETATM.end()){
            cerr << "ERR: couldn't find the positions of the associated glycan" << it->second << endl;
            listPositionGlycans.push_back({0,0,0});
        } else {
            vector<double> positions = it2->second;
            if((positions.size() % 3) != 0){
                cerr << "ERR: glycan vector position for " << it->second << " is not a multiple of 3" << endl;
             }
            cout << "found " << positions.size() << " positions\t";
            if(positions.size() == 0) cerr << "Empty glycan positions ??" << endl;

            for(size_t i = 0; i < positions.size()/3; ++i){

                 biu::DblPoint oneAtom(positions[3*i], positions[3*i+1], positions[3*i+2]);
                 glycan.push_back(oneAtom);
            }
            biu::DblPoint center = getCentroid(glycan);
            cout << "Glycan " << it->second << " Has center at " << center.getX() << "," << center.getY() << "," << center.getZ() << "\t";
            listPositionGlycans.push_back({center.getX(), center.getY(), center.getZ()});
        }

        if(it->first.size() < 2) cerr << "ERR: parse glycans, Empty/wrong identifier for AA" << it->first << ", should be chain + IDresidue" << endl;
        char chainAA = it->first[0];
        int IDres;
        stringstream read(it->first.substr(1, it->first.size()-1));
        read >> IDres;

        AAboundByGlycans.push_back(std::pair<char,int>(chainAA, IDres));
        std::map<string, vector<double> >::iterator it3 = positionsAllCAs.find(it->first);
        if(it3 == positionsAllCAs.end()){
            cerr << "ERR: readGlycans(), could not retrieve the position of " << it->first << endl;
            PositionsAAboundByGlycans.push_back({NAN, NAN, NAN});
            cout << "- didn't find AA position" << endl;
        } else {
            PositionsAAboundByGlycans.push_back(it3->second);
            cout << "AAposition is " << printVector(it3->second) << endl;
        }


//        cout << endl;
    }
    cout << "Found " << nGlycans << "Glycans" << endl;

    if(listPositionGlycans.size() != nGlycans) cerr << "ERR: Glycans, inconsistency size between number found " << nGlycans<< " and size of listPositionGlycans " << listPositionGlycans.size() << endl;
    if(AAboundByGlycans.size() != nGlycans) cerr << "ERR: Glycans, inconsistency size between number found " << nGlycans<< " and size of AAboundByGlycans " << AAboundByGlycans.size() << endl;

    // now, need to fill these infos in the mother class:
//    int nGlycans;
//    vector< vector3D > listPositionGlycans;
//    vector< std::pair<char, int> > AAboundByGlycans;






//        chain.rbegin()->idBegin = seqNum;
//        chain.rbegin()->atom = "CoM";
//        cAlpha = biu::DblPoint(0.0,0.0,0.0);
//        chain.rbegin()->idBegin = seqNum;
//        chain.rbegin()->atom = "CoM";

//        chain.rbegin()->points.push_back(cAlpha);
//        or
//        chain.rbegin()->points.push_back(getCentroid(sideChain));
//        chain.rbegin()->seq.push_back(aaName);


//    listPositionGlycans

}


// Careful, this function parses the output of LATFIT, where discretized atoms appear as HETATM, and glycans are removed...
void latFitToLattice::parseLatFitPDB(string filename, string typeAtoms){



    std::map<string, char> AAchar;
    AAchar["CYS"] = 'C';
    AAchar["MET"] = 'M';
    AAchar["PHE"] = 'F';
    AAchar["ILE"] = 'I';
    AAchar["LEU"] = 'L';
    AAchar["VAL"] = 'V';
    AAchar["TRP"] = 'W';
    AAchar["TYR"] = 'Y';
    AAchar["ALA"] = 'A';
    AAchar["GLY"] = 'G';
    AAchar["THR"] = 'T';
    AAchar["SER"] = 'S';
    AAchar["ASN"] = 'N';
    AAchar["GLN"] = 'Q';
    AAchar["ASP"] = 'D';
    AAchar["GLU"] = 'E';
    AAchar["HIS"] = 'H';
    AAchar["ARG"] = 'R';
    AAchar["LYS"] = 'K';
    AAchar["PRO"] = 'P';

    listPositions.clear();
    listOriginPositions.clear();
    sequence.clear();
    IDresidues.clear();
    ifstream f;
    string lineType;
    f.open(filename.c_str());
    if(!f) {
        cerr << "ERR: latFitToLattice::parseLatFitPDB, file not found:\n" << filename << endl; return;
    }
    while(f.good()){
        f >> lineType;
        if(!lineType.compare(typeAtoms)){ // HETATM by default
            int ID = 0;
            string typeAtom;
            string typeResidue;
            //string typeChain;
            int IDinProt = 0;
            double x = NAN;
            double y = NAN;
            double z = NAN;
            f >> ID >> typeAtom >> typeResidue;
            char typeChain;
            f.get(typeChain); // stupid PDB format where chain and ID can touch each other
            while(typeChain == ' '){
                f.get(typeChain);
            }
            //cout << "got" << typeChain << endl;
            f >> IDinProt >> x >> y >> z;
            //cout << ID << " " << typeAtom << " " << typeResidue << " " << typeChain << " " << IDinProt << " " << x << " " << y << " " << z;
            if(typeChain == 'L'){
            //if(!typeChain.compare(string("L"))){
                vector3D res = {x,y,z};
                //cout << typeResidue << "\t" << IDinProt << "\t" << print(res) << endl;
                listPositions.push_back(res);

                char typeAA = '?';
                std::map<string,char>::iterator it = AAchar.find(typeResidue);
                if(it == AAchar.end()) {cerr << "WRN: Amino Acid " << typeResidue << " not known, encoded as '?'" << endl;}
                else typeAA = it->second;
                sequence.push_back(typeAA);
                IDresidues.push_back(IDinProt);
            }
            else { // original PDB coordinates
                vector3D res = {x,y,z};
                //cout << typeResidue << "\t" << IDinProt << "\t" << print(res) << endl;
                listOriginPositions.push_back(res);
            }
            char buf[100001];
            f.getline(buf, 100000);
        } else {
            char buf[100001];
            f.getline(buf, 100000);
        }
    }
    //cout << "\n\n" << listPositions.size() << endl;
    //cout << "End of file: have read " << listPositions.size() << " Discretized positions, and " << listOriginPositions.size() << " original positions" << endl;
    f.close();
}


void TestPDBtoLat() {
    // latTest -pdbFile="C:/Qt/PDBsMat/latpack-1.9.1/build-testLat-Desktop_Qt_5_5_1_MinGW_32bit-Release/6nqd_remove_ins.pdb" -pdbAtom=CA -pdbChain=A -pdbChainGaps -lat=CUB -outMode=PDB -outFile=out.pdb
    latFitToLattice a = latFitToLattice();
    //a.parseLatFitPDB("C:/Users/Philippe/Desktop/Discrete/6nqs_out_remIns.pdb");
    a.parseLatFitPDB("C:/Qt/1NCA_NdiscretizedFuC5.25.pdb");//Zapotec/PDB/3ECA/outEC3.pdb");
    a.transform();
#ifdef ALLOW_GRAPHICS

    a.testConversions();
    glDisplay();

    // This piece of code is now redundant with asSuperProtein
    if(a.structures.size() != a.positions.size()) cerr << "ERR: incompatible positions and structure numbers" << endl;
    for(unsigned int i = 0; i < a.structures.size(); ++i){
        cout << a.positions[i] << "\t" << a.structures[i] << endl;
        struct3D* test = new struct3D(a.structures[i], UnDefined, a.positions[i]);
        addToDisplay(test, false);
        //addDissected(test, UnDefined, a.positions[i]);
    }
    addToDisplay(&(a.listIntPositions));
    addToDisplay(&(a.listRotatedOriginPositions));

    glutMainLoop();
#endif
}





