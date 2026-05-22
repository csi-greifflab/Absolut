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


QMAKE_CXXFLAGS += -std=c++11

HEADERS += \
    $$PWD/kabsch_nogsl.hh \
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
    $$PWD/kabsch_nogsl.cc \
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
