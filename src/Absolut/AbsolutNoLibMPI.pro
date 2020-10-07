include("../Ymir/YmirNoGL.pri")

DEFINES += "USE_MPI"
DEFINES += "NOQT"
DEFINES += "NO_LIBS"
CONFIG -= QT

HEADERS += \
    antigenLib.h \
	epitope.h \
	fileformats.h \
	html.h \
    importrepertoire.h \
    motifFeatures.h \
    selfEvo.h \
	poolStructs.h \
    quality.h \
    topology.h

SOURCES += \
    antigenLib.cpp \
	epitope.cpp \
    fileformats.cpp \
    html.cpp \
    importrepertoire.cpp \
    motifFeatures.cpp \
    selfEvo.cpp \
    quality.cpp \
	poolStructs.cpp \
    delimain.cpp \
    topology.cpp


win32: TARGET = AbominationMPI
unix: TARGET = AbominationMPI

QMAKE_CXXFLAGS += "-Wno-old-style-cast" "-Wno-shorten-64-to-32" "-Wno-sign-conversion" "-Wno-old-style-cast" "-Wno-implicit-int-float-conversion"

LIBS += -pthread

win32:{
    LIBS += -LC:/MyPrograms/MSMPI_SDK/Lib/x64/ -lmsmpi
    INCLUDEPATH += C:/MyPrograms/MSMPI_SDK/Include/
    DEPENDPATH += C:/MyPrograms/MSMPI_SDK/Lib/x64/
}

unix: {
	QMAKE_CXX = mpicxx
	QMAKE_CXX_RELEASE = $$QMAKE_CXX
	QMAKE_CXX_DEBUG = $$QMAKE_CXX
	QMAKE_LINK = $$QMAKE_CXX
	QMAKE_CC = mpicc

	QMAKE_CFLAGS += $$system(mpicc --showme:compile)
	QMAKE_LFLAGS += $$system(mpicxx --showme:link)
	QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
	QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
}

# This also sets the correct compile and link flags in addition to changing the linker to mpicxx as well. The

#https://gcc.gnu.org/onlinedocs/gcc-4.9.2/gcc/Warning-Options.html
CONFIG += warn_off
#https://gcc.gnu.org/onlinedocs/gcc-4.9.2/gcc/Warning-Options.html these are different
QMAKE_CXXFLAGS += -std=c++11 -Wno-extra
# -Wno-oldstyle-case -Wno-unused-parameter -Wno-shorten-64-to-32 -Wno-sign-conversion -Wno-ignored-qualifiers -Wno-type-limits -Wno-misleading-indentation -Wno-sign-compare -Wno-implicit-float-conversion
# QMAKE_CXXFLAGS += -Wall -Wextra -Wunsafe-loop-optimizations -pedantic -Wfloat-equal -Wundef -Wpointer-arith -Wcast-align -Wunreachable-code

