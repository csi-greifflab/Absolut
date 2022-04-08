include("../Ymir/YmirNoGL.pri")
DEFINES += "NOQT"

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
	dlab.h \
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
	dlab.cpp \
    topology.cpp



win32: TARGET = AbsolutNoLib
unix: TARGET = AbsolutNoLib

QMAKE_CXXFLAGS += "-Wno-old-style-cast" "-Wno-shorten-64-to-32" "-Wno-sign-conversion" "-Wno-old-style-cast" "-Wno-implicit-int-float-conversion"

#https://gcc.gnu.org/onlinedocs/gcc-4.9.2/gcc/Warning-Options.html
CONFIG += warn_off
#https://gcc.gnu.org/onlinedocs/gcc-4.9.2/gcc/Warning-Options.html these are different
QMAKE_CXXFLAGS += -std=c++11 -Wno-extra
# -Wno-oldstyle-case -Wno-unused-parameter -Wno-shorten-64-to-32 -Wno-sign-conversion -Wno-ignored-qualifiers -Wno-type-limits -Wno-misleading-indentation -Wno-sign-compare -Wno-implicit-float-conversion
# QMAKE_CXXFLAGS += -Wall -Wextra -Wunsafe-loop-optimizations -pedantic -Wfloat-equal -Wundef -Wpointer-arith -Wcast-align -Wunreachable-code

CONFIG -= QT