#ifndef PLOT3D_H
#define PLOT3D_H

/// \file
/// \brief Home-made functions to visualize multiple structures in space
/// \defgroup OpenGL 3D visualization of structures (plot3D.h/cpp)
/// Usage:
///
/// char *c[] = {(char*)"Hello",NULL};
/// glDisplay(1, c);
/// addToDisplay(...);
///  ...
/// addToDisplay(...);
/// addToDisplay(...);
/// glutMainLoop();
///
/// Make sure that the pointers given to addToDisplay are not deleted before displaying!
/// Control: via mouse rolling, and via keyboard 'q' or 'd': next or previous structure.
///
/// Use #define ALLOW_GRAPHICS here in plot3D.h in order to be able to use OpenGL functions. Do not define it in order to compile without openGL (no other library required then)
/// \date 9th October 2019 \author Philippe A. Robert

// It is easier to use DEFINES += "ALLOW_GRAPHICS" into the .pro or .pri files
//#define ALLOW_GRAPHICS




#ifdef ALLOW_GRAPHICS

#ifdef _WIN32
#include <GL/glut.h>
#else
#include <GL/glut.h>
#endif

#include "proteins.h"
#include <vector>
using namespace std;



/// \brief Initializes OpenGL. Append=true allows to add to previous display if existed \ingroup OpenGL
void glDisplay(int argc = 0, char** argv = nullptr, bool append = false);

/// \brief Write text at a position in space. \ingroup OpenGL
void drawBitmapText(string s,float x,float y,float z);

/// \brief Add a new protein to the display. if always_visible is true, it will be present each time 'q' or 'd' are presses during visualization.
/// The pointer is stored in a list, but no copy is made. Make sure the given pointer will not be deleted. \ingroup OpenGL
void addToDisplay(superProtein* new_s, bool always_visible = true);

/// \brief Add a new structure to the display. \ingroup OpenGL
void addToDisplay(struct3D* new_s, bool always_visible = true);

/// \brief Add a new set of positions (green bowls) to the display. \ingroup OpenGL
void addToDisplay(set<int>* new_s);

/// \brief Add a new set of positions (green bowls) to the display. \ingroup OpenGL
void addToDisplay(vector<vector<double>>* listpoints);

void addToDisplay(vector< std::pair<int, double> > * newHeatmap);

void setHeatmapColorZero(vector<double> v);

void setHeatmapColorOne(vector<double> v);

/// \brief Remove all objects from the display. \ingroup OpenGL
void clearDisplay();

/// \brief Decompose a structures into each step of its development, from the first move to the complete structure. \ingroup OpenGL
void addDissected(struct3D* new_s,  moveDirection startDir = UnDefined, int startPos = -1);

void savePicture(string fileName);

void setPictureFileNamePrefix(string newPrefix);

#endif
#endif // PLOT3D_H
