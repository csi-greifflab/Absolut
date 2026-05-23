DEFINES += "ALLOW_GRAPHICS"

include("soil/SOIL.pri")

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
    ../Tools/distribution.h

QMAKE_CXXFLAGS += -std=c++11
#-Wno-unused-parameter

unix: LIBS += -lglut -lGLU -lGL -lSOIL

win32: {
    LIBS += -LC:/Programs/freeglut/bin/ -lfreeglut
    #INCLUDEPATH += C:/Programs/freeglut/include/
    #DEPENDPATH += C:/Programs/freeglut/lib/
    LIBS +=  -lopengl32 -lglu32 -lpthread -lglut32
}

win32: TARGET = Ymir.exe
unix: TARGET = Ymir

