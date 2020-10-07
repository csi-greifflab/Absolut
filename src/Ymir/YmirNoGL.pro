SOURCES += \
    compact.cpp \
    lattice.cpp \
    proteins.cpp \
    receptorligand.cpp \
    fastaffinity.cpp \
    plot3d.cpp \
    zaptrackball.cpp \
    ../Tools/stopwatch.cpp \
    ../Tools/md5.cpp \
    ../Tools/zaprandom.cpp \
    ../Tools/distribution.cpp \
    ../Tools/nucleotides.cpp \
    YmirMain.cpp

HEADERS += \
    compact.h \
    lattice.h \
    proteins.h \
    receptorligand.h \
    fastaffinity.h \
    plot3d.h \
    zaptrackball.h \
    ymir.h \
    ../Tools/md5.h \
    ../Tools/stopwatch.h \
    ../Tools/zaprandom.h \
    ../Tools/nucleotides.h \
    ../Tools/distribution.h
	


QMAKE_CXXFLAGS += -std=c++11 -Wno-unused-parameter

win32: TARGET = YmirNoGL.exe
unix: TARGET = YmirNoGL

