INCLUDEPATH += $$PWD
DEPENDPATH += $$PWD
message("Info: Absolut library included from folder :")
message($$PWD)


include("../Ymir/Ymir.pri")
include("../latFit/latFit.pri")

HEADERS += \
    $$PWD/antigenLib.h \
    $$PWD/epitope.h \
    $$PWD/fileformats.h \
    $$PWD/generatemutants.h \
    $$PWD/html.h \
    $$PWD/importrepertoire.h \
    $$PWD/motifFeatures.h \
    $$PWD/poolstructs.h \
    $$PWD/selfEvo.h \
    $$PWD/discretize.h \
    $$PWD/pdb.h \
    $$PWD/../Tools/dirent.h \
    $$PWD/topology.h

SOURCES += \
    $$PWD/antigenLib.cpp \
    $$PWD/epitope.cpp \
    $$PWD/fileformats.cpp \
    $$PWD/generatemutants.cpp \
    $$PWD/html.cpp \
    $$PWD/importrepertoire.cpp \
    $$PWD/motifFeatures.cpp \
    $$PWD/oldScripts.cpp \
    $$PWD/poolstructs.cpp \
    $$PWD/selfEvo.cpp \
    $$PWD/discretize.cpp \
    $$PWD/pdb.cpp \
    $$PWD/topology.cpp

FORMS += \
    $$PWD/pdb.ui

QMAKE_CXXFLAGS += "-Wno-old-style-cast" "-Wno-shorten-64-to-32" "-Wno-sign-conversion" "-Wno-old-style-cast" "-Wno-implicit-int-float-conversion"

QT += core gui widgets

#https://gcc.gnu.org/onlinedocs/gcc-4.9.2/gcc/Warning-Options.html
CONFIG += warn_off
#https://gcc.gnu.org/onlinedocs/gcc-4.9.2/gcc/Warning-Options.html these are different
QMAKE_CXXFLAGS += -std=c++11 -Wno-extra
# -Wno-oldstyle-case -Wno-unused-parameter -Wno-shorten-64-to-32 -Wno-sign-conversion -Wno-ignored-qualifiers -Wno-type-limits -Wno-misleading-indentation -Wno-sign-compare -Wno-implicit-float-conversion
# QMAKE_CXXFLAGS += -Wall -Wextra -Wunsafe-loop-optimizations -pedantic -Wfloat-equal -Wundef -Wpointer-arith -Wcast-align -Wunreachable-code

