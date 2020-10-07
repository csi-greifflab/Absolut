# In order to compile yout C++ files using Moonfit, two ways:
# Option 1: design a .pro file containing your files and add the following command inside : include(Moonfit/moonfit.pri).
# Option 2: add all C++ and h files from moonfit folder and subfolders into your .pro (together with your own files including moonfit).
#			in that case, please make sure to use the following commands inside your .pro 

# A - Options that you can put inside your .pro file (just uncomment there)

#name of the executable file generated
#TARGET = LatFit.exe

#put += console to run in a separate terminal
#CONFIG += console

#bundles might be required for MAC OS 
#CONFIG -= app_bundle

#TEMPLATE = app

# this function allows to do include(Moonfit.pri) from another project
INCLUDEPATH += $$PWD
DEPENDPATH += $$PWD
message("Info: LatFit Library included from folder :")
message($$PWD)

# B - options for including Moonfit.

# note: the options widgets and core gui are necessary for the production of the ui_....h automatically from interface files (.ui).
# QT += core gui
# greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport


# ======================== BOOST PART ==================================
# boost is not required anymore for the code.
# win32: INCLUDEPATH += C:\Qt\Boost\boost_1_62_0\boost_1_62_0
#in windows, boosts raises the following error : "boost unable to find numeric literal operator operator Q" => need to add the following option
# win32: QMAKE_CXXFLAGS += -fext-numeric-literals
# to avoid lots of useless warnings from boost library files (wtf !)
QMAKE_CXXFLAGS +=  -Wreorder
QMAKE_CFLAGS += -Wno-unused-local-typedefs -Wno-unused-parameter -Wno-oldstyle-case -Wno-unused-parameter -Wno-shorten-64-to-32 -Wno-sign-conversion -Wno-ignored-qualifiers -Wno-type-limits -Wno-misleading-indentation -Wno-sign-compare -Wno-implicit-float-conversion


# Necessary only for boost - would also be necessary for the syntax vector<int> bla = {1, 2, 3};
QMAKE_CXXFLAGS += -std=c++11

HEADERS += \
	$$PWD/latFit.h  \
    $$PWD/biu/Rotator3D.hh \
    $$PWD/biu/SuperPos_Kabsch.hh \
    $$PWD/biu/LatticeProteinUtil.hh \
    $$PWD/biu/LatticeDescriptorCKW.hh \
    $$PWD/biu/LatticeDescriptorCUB.hh \
    $$PWD/biu/LatticeDescriptorFCC.hh \
    $$PWD/biu/LatticeDescriptorSQR.hh \
    $$PWD/biu/LatticeModel.hh \
    $$PWD/biu/Matrix.hh \
    $$PWD/biu/OptionParser.hh \
    $$PWD/biu/Point.hh \
    $$PWD/biu/LatticeDescriptor.hh



SOURCES += \
    $$PWD/latFit.cc \
    $$PWD/biu/Rotator3D.cc \
    $$PWD/biu/SuperPos_Kabsch.cc \
    $$PWD/biu/LatticeProteinUtil.cc \
    $$PWD/biu/LatticeDescriptorCKW.cc \
    $$PWD/biu/LatticeDescriptorCUB.cc \
    $$PWD/biu/LatticeDescriptorFCC.cc \
    $$PWD/biu/LatticeDescriptorSQR.cc \
    $$PWD/biu/LatticeModel.cc \
    $$PWD/biu/OptionParser.cc \
    $$PWD/biu/LatticeDescriptor.cc

DISTFILES += \
    $$PWD/biu/LatticeModel.icc \
    $$PWD/biu/Matrix.icc \
    $$PWD/biu/LatticeDescriptor.icc

win32 {
    INCLUDEPATH += C:/MyPrograms/gsl-2.6/
  # INCLUDEPATH += $$PWD/gsl/gsl-2.2.1-shared/include/
    INCLUDEPATH += C:/MyPrograms/gsl-2.6/.libs/
  # INCLUDEPATH += $$PWD/gsl/gsl-2.2.1-shared/lib
    LIBS += -LC:/MyPrograms/gsl-2.6/.libs/  -llibgsl-25
    LIBS += -LC:/MyPrograms/gsl-2.6/cblas/.libs -llibgslcblas-0
  # LIBS += -L$$PWD/gsl/gsl-2.2.1-shared/bin  -llibgsl-19 -llibgslcblas-0
}

unix {
    LIBS += -lgsl -lgslcblas
}


