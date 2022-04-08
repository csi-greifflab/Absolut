#ifndef DLAB_H
#define DLAB_H


struct superProtein;

class dlab
{
public:
    dlab();
};

#include <vector>
#include <string>
#include <set>

struct voxelGrid {
    std::vector< std::vector< std::vector<int>>> content;
    size_t sizeX, sizeY, sizeZ;
    voxelGrid(size_t x, size_t y, size_t z);
    voxelGrid(const voxelGrid& toCopy);
    voxelGrid(superProtein* P);
    voxelGrid(std::string structure, int position, std::string AAs);
    voxelGrid operator + (const voxelGrid& toAdd);

    std::vector<double> getCentre();
    double access(int x, int y, int z);
    std::string asText(bool dim2 = true);
    std::vector<voxelGrid*> PossibleCenterings(std::vector<double> c, size_t newSizeX = 0, size_t newSizeY = 0, size_t newSizeZ = 0);
    voxelGrid reshapeAroundCenter(size_t newSizeX, size_t newSizeY, size_t newSizeZ);

    voxelGrid interfaceWith(voxelGrid& other);
};

std::vector < std::pair<std::string, std::string> > latticesComplexes(superProtein* antigen, int startPos, std::string structure, std::string AAseq, size_t latticeSize = 6);

void testVoxelGrid();
void testCentering();


std::set<std::string> split5(std::string s);
int fnat(std::string code1, std::string code2);

#endif // DLAB_H
