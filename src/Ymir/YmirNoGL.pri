# this file allows to do include Ymir library from another project 
# (without needing to include/compile all files one by one in your project)

# IMPORTANT: two folders are required: Ymir/ and Tools/, side by side in the same folder, because Ymir will look for ../Tools/blabla.h

# use: 		write 'include("blablavla/Ymir/Ymir.pri")' in a .pro file, this is enough. 
#			and add include "blabla/Ymir/Ymir.h" into your c++ files that need Ymir.

INCLUDEPATH += $$PWD
DEPENDPATH += $$PWD
message("Info: Ymir Library included from folder :")
message($$PWD)

QMAKE_CXXFLAGS +=  -Wreorder
QMAKE_CFLAGS += -Wno-unused-local-typedefs -Wno-unused-parameter -Wno-oldstyle-case -Wno-shorten-64-to-32 -Wno-sign-conversion -Wno-ignored-qualifiers -Wno-type-limits -Wno-misleading-indentation -Wno-sign-compare -Wno-implicit-float-conversion

#  necessary for the syntax vector<int> bla = {1, 2, 3};
QMAKE_CXXFLAGS += -std=c++11

HEADERS += \
	$$PWD/compact.h \
    $$PWD/lattice.h \
    $$PWD/proteins.h \
    $$PWD/receptorligand.h \
    $$PWD/fastaffinity.h \
    $$PWD/plot3d.h \
    $$PWD/../Tools/stopwatch.h \
    $$PWD/../Tools/md5.h \
    $$PWD/ymir.h \
    $$PWD/zaptrackball.h \
    $$PWD/../Tools/zaprandom.h \
    $$PWD/../Tools/nucleotides.h \
    $$PWD/../Tools/distribution.h

SOURCES += \
	$$PWD/compact.cpp \
    $$PWD/lattice.cpp \
    $$PWD/proteins.cpp \
    $$PWD/receptorligand.cpp \
    $$PWD/fastaffinity.cpp \
    $$PWD/plot3d.cpp \
    $$PWD/../Tools/stopwatch.cpp \
    $$PWD/../Tools/md5.cpp \
    $$PWD/zaptrackball.cpp \
    $$PWD/../Tools/zaprandom.cpp \
    $$PWD/../Tools/nucleotides.cpp \
    $$PWD/../Tools/distribution.cpp

QMAKE_CXXFLAGS += -std=c++11 -Wno-extra

# Note: the files using openGL (plot 3D) are still included into compiling, but because 'ALLOW_GRAPHICS' is not defined,
# the GL code is removed in those files by #ifdefs

