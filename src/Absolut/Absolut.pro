include("../Ymir/Ymir.pri")
include("../latFit/latFit.pri")

HEADERS += \
    antigenLib.h \
    dlab.h \
    epitope.h \
    fileformats.h \
    generatemutants.h \
    html.h \
    importrepertoire.h \
    motifFeatures.h \
    poolstructs.h \
    runListTasks.h \
    selfEvo.h \
    discretize.h \
    pdb.h \
    ../Tools/dirent.h \
    quality.h \
    topology.h

SOURCES += \
    antigenLib.cpp \
    dlab.cpp \
    epitope.cpp \
    fileformats.cpp \
    generatemutants.cpp \
    html.cpp \
    importrepertoire.cpp \
    motifFeatures.cpp \
    poolstructs.cpp \
    runListTasks.cpp \
    selfEvo.cpp \
    discretize.cpp \
    pdb.cpp \
    quality.cpp \
    delimain.cpp \
    topology.cpp

FORMS += \
    pdb.ui

win32: TARGET = Absolut
unix: TARGET = Absolut

QMAKE_CXXFLAGS += "-Wno-old-style-cast" "-Wno-shorten-64-to-32" "-Wno-sign-conversion" "-Wno-old-style-cast" "-Wno-implicit-int-float-conversion"

QT += core gui widgets

#https://gcc.gnu.org/onlinedocs/gcc-4.9.2/gcc/Warning-Options.html
CONFIG += warn_off
#https://gcc.gnu.org/onlinedocs/gcc-4.9.2/gcc/Warning-Options.html these are different
QMAKE_CXXFLAGS += -std=c++11 -Wno-extra
# -Wno-oldstyle-case -Wno-unused-parameter -Wno-shorten-64-to-32 -Wno-sign-conversion -Wno-ignored-qualifiers -Wno-type-limits -Wno-misleading-indentation -Wno-sign-compare -Wno-implicit-float-conversion
# QMAKE_CXXFLAGS += -Wall -Wextra -Wunsafe-loop-optimizations -pedantic -Wfloat-equal -Wundef -Wpointer-arith -Wcast-align -Wunreachable-code

