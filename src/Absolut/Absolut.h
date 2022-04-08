#include "../Ymir/ymir.h"

#include "../Tools/md5.h"
#include "../Tools/zaprandom.h"
#include "../Tools/dirent.h"
#include "../Tools/stopwatch.h"

#include "motifFeatures.h"
#include "quality.h"
#include "selfEvo.h"
#include "antigenLib.h"
#include "importrepertoire.h"
#include "poolstructs.h"
#include "epitope.h"
#include "html.h"
#include "fileformats.h"
#include "topology.h"

// for absolut this is 10 and 11
#define DefaultReceptorSizeBonds 10
#define DefaultContactPoints 11

// for sleep
#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

#include <cstdlib>

#ifndef NO_LIBS
#include "discretize.h"
#endif

#ifndef NOQT
#include <QApplication>
#include "pdb.h"
#endif

#ifdef USE_MPI
#include "mpi.h"
#endif
