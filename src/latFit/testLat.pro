HEADERS += \
	kabsch_nogsl.hh \
    biu/Rotator3D.hh \
    biu/SuperPos_Kabsch.hh \
    biu/LatticeProteinUtil.hh \
    biu/LatticeDescriptorCKW.hh \
    biu/LatticeDescriptorCUB.hh \
    biu/LatticeDescriptorFCC.hh \
    biu/LatticeDescriptorSQR.hh \
    biu/LatticeModel.hh \
    biu/Matrix.hh \
    biu/OptionParser.hh \
    biu/Point.hh \
    biu/LatticeDescriptor.hh



SOURCES += \
	kabsch_nogsl.cc \
    latFit.cc \
    biu/Rotator3D.cc \
    biu/SuperPos_Kabsch.cc \
    biu/LatticeProteinUtil.cc \
    biu/LatticeDescriptorCKW.cc \
    biu/LatticeDescriptorCUB.cc \
    biu/LatticeDescriptorFCC.cc \
    biu/LatticeDescriptorSQR.cc \
    biu/LatticeModel.cc \
    biu/OptionParser.cc \
    biu/LatticeDescriptor.cc

DISTFILES += \
    biu/LatticeModel.icc \
    biu/Matrix.icc \
    biu/LatticeDescriptor.icc
