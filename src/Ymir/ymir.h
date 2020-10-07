#ifndef YMIR_H
#define YMIR_H

#include "lattice.h"
#include "compact.h"
#include "proteins.h"
#include "receptorligand.h"
#include "fastaffinity.h"

// Think about activating/deactivating the OpenGL graphical display, by #define ALLOW_GRAPHICS or not inside plot3d.h

// Windows: If using the provided qt .pro file, or the provided makefile, the libraries are already linked to the Ymir subfolders.
// => nothing to do for compiling
// => Note : for running OpeenGL displays on windows, put the 3 dll files inside the run folder of the exe file 

// Linux/MAC: If ALLOW_GRAPHICS is defined, you need to have glut and SOIL libraries installed/link.
// In linux, apt-get install freeglut (or so). The provided Makefile should work.
// If you wish to compile manually, think of including the library folders inside Ymir and use the flags for compilation, and to link by -lglut -lGLU -lGL -lSOIL

#endif

