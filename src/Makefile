MAKEFILE      = Makefile

EQ            = =

####### Compiler, tools and options

CC            = gcc
CXX           = g++

MPICC         = mpicc
MPICXX        = mpicxx

DEFINES       = -DNOQT -DNO_LIBS
DEFINES_MPI   = -DUSE_MPI -DNOQT -DNO_LIBS

CFLAGS        = -pipe -Wno-unused-local-typedefs -Wno-unused-parameter -Wno-oldstyle-case -Wno-shorten-64-to-32 -Wno-sign-conversion -Wno-ignored-qualifiers -Wno-type-limits -Wno-misleading-indentation -Wno-sign-compare -Wno-implicit-float-conversion -O2 -w -D_REENTRANT -fPIC $(DEFINES)
CXXFLAGS      = -pipe -Wreorder -std=c++11 -Wno-extra -Wno-old-style-cast -Wno-shorten-64-to-32 -Wno-sign-conversion -Wno-implicit-int-float-conversion -O2 -w -D_REENTRANT -fPIC $(DEFINES)

CXXFLAGSMPI   = -pipe -Wreorder -std=c++11 -Wno-extra -Wno-old-style-cast -Wno-shorten-64-to-32 -Wno-sign-conversion -Wno-implicit-int-float-conversion -O2 -w -D_REENTRANT -fPIC -DMPICH_IGNORE_CXX_SEEK $(DEFINES_MPI)

INCPATH       = -IDelicab -I. -IYmir -I. -isystem /usr/include/libdrm 
DEL_FILE      = rm -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p
COPY          = cp -f
COPY_FILE     = cp -f
COPY_DIR      = cp -f -R
INSTALL_FILE  = install -m 644 -p
INSTALL_PROGRAM = install -m 755 -p
INSTALL_DIR   = cp -f -R
DEL_FILE      = rm -f
SYMLINK       = ln -f -s
DEL_DIR       = rmdir
MOVE          = mv -f
TAR           = tar -cf
COMPRESS      = gzip -9f
DISTNAME      = DelicabNoLib1.0.0
DISTDIR = /home/pprobert/Desktop/TempSpace/testsDelicab/Bindings/.tmp/DelicabNoLib1.0.0

LINK          = g++
LINKMPICXX    = mpicxx

LFLAGS        = -Wl,-O1
LIBS          = $(SUBLIBS) -lpthread   
AR            = ar cqs
RANLIB        = 
SED           = sed
STRIP         = strip

####### Output directory

OBJECTS_DIR   = ./

####### Files

SOURCES       = Ymir/compact.cpp \
		Ymir/lattice.cpp \
		Ymir/proteins.cpp \
		Ymir/receptorligand.cpp \
		Ymir/fastaffinity.cpp \
		Ymir/plot3d.cpp \
		Tools/stopwatch.cpp \
		Tools/md5.cpp \
		Ymir/zaptrackball.cpp \
		Tools/zaprandom.cpp \
		Tools/nucleotides.cpp \
		Tools/distribution.cpp \
		Delicab/antigenLib.cpp \
		Delicab/epitope.cpp \
		Delicab/fileformats.cpp \
		Delicab/html.cpp \
		Delicab/importrepertoire.cpp \
		Delicab/motifFeatures.cpp \
		Delicab/selfEvo.cpp \
		Delicab/quality.cpp \
		Delicab/poolstructs.cpp \
		Delicab/delimain.cpp \
		Delicab/topology.cpp
OBJECTS       = compact.o \
		lattice.o \
		proteins.o \
		receptorligand.o \
		fastaffinity.o \
		plot3d.o \
		stopwatch.o \
		md5.o \
		zaptrackball.o \
		zaprandom.o \
		nucleotides.o \
		distribution.o \
		antigenLib.o \
		epitope.o \
		fileformats.o \
		html.o \
		importrepertoire.o \
		motifFeatures.o \
		selfEvo.o \
		quality.o \
		poolstructs.o \
		delimain.o \
		topology.o
DIST          = Delicab/.qmake.stash \
		Delicab/Ymir/YmirNoGL.pri \
		Delicab/Delicab/DelicabNoLib.pro Ymir/compact.h \
		Ymir/lattice.h \
		Ymir/proteins.h \
		Ymir/receptorligand.h \
		Ymir/fastaffinity.h \
		Ymir/plot3d.h \
		Tools/stopwatch.h \
		Tools/md5.h \
		Ymir/ymir.h \
		Ymir/zaptrackball.h \
		Tools/zaprandom.h \
		Tools/nucleotides.h \
		Tools/distribution.h \
		Delicab/antigenLib.h \
		Delicab/epitope.h \
		Delicab/fileformats.h \
		Delicab/html.h \
		Delicab/importrepertoire.h \
		Delicab/motifFeatures.h \
		Delicab/selfEvo.h \
		Delicab/poolstructs.h \
		Delicab/quality.h Ymir/compact.cpp \
		Ymir/lattice.cpp \
		Ymir/proteins.cpp \
		Ymir/receptorligand.cpp \
		Ymir/fastaffinity.cpp \
		Ymir/plot3d.cpp \
		Tools/stopwatch.cpp \
		Tools/md5.cpp \
		Ymir/zaptrackball.cpp \
		Tools/zaprandom.cpp \
		Tools/nucleotides.cpp \
		Tools/distribution.cpp \
		Delicab/antigenLib.cpp \
		Delicab/epitope.cpp \
		Delicab/fileformats.cpp \
		Delicab/html.cpp \
		Delicab/importrepertoire.cpp \
		Delicab/motifFeatures.cpp \
		Delicab/selfEvo.cpp \
		Delicab/quality.cpp \
		Delicab/poolstructs.cpp \
		Delicab/delimain.cpp
QMAKE_TARGET  = DelicabNoLib
DESTDIR       = 
TARGET        = AbsolutNoLib
TARGETMPI     = AbominationMPI



first: all
####### Build rules

DelicabNoLib:  $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)

DelicabNoLibMPI:  $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)


Makefile: Delicab/DelicabNoLib.pro .qmake.stash \
		Ymir/YmirNoGL.pri \
		Delicab/DelicabNoLib.pro
	$(QMAKE) -o Makefile Delicab/DelicabNoLib.pro

	
.qmake.stash:
Ymir/YmirNoGL.pri:
Delicab/DelicabNoLib.pro:

all: clean Makefile DelicabNoLib 

MPIcxx: clean
MPIcxx: CC = $(MPICC)
MPIcxx: CXX = $(MPICXX)
MPIcxx: LINK = $(LINKMPICXX)
MPIcxx: CXXFLAGS = $(CXXFLAGSMPI)
MPIcxx: TARGET = $(TARGETMPI)
MPIcxx: DelicabNoLibMPI

MPIc++: clean
MPIc++: CC = $(MPICC)
MPIc++: CXX = mpic++
MPIc++: LINK = mpic++
MPIc++: CXXFLAGS = $(CXXFLAGSMPI)
MPIc++: TARGET = $(TARGETMPI)
MPIc++: DelicabNoLibMPI

MPIgxx: clean
MPIgxx: CC = mpigcc
MPIgxx: CXX = mpigxx
MPIgxx: LINK = mpigxx
MPIgxx: CXXFLAGS = $(CXXFLAGSMPI)
MPIgxx: TARGET = $(TARGETMPI)
MPIgxx: DelicabNoLibMPI

dist: distdir FORCE
	(cd `dirname $(DISTDIR)` && $(TAR) $(DISTNAME).tar $(DISTNAME) && $(COMPRESS) $(DISTNAME).tar) && $(MOVE) `dirname $(DISTDIR)`/$(DISTNAME).tar.gz . && $(DEL_FILE) -r $(DISTDIR)

distdir: FORCE
	@test -d $(DISTDIR) || mkdir -p $(DISTDIR)
	$(COPY_FILE) --parents $(DIST) $(DISTDIR)/
	$(COPY_FILE) --parents Ymir/compact.h Ymir/lattice.h Ymir/proteins.h Ymir/receptorligand.h Ymir/fastaffinity.h Ymir/plot3d.h Tools/stopwatch.h Tools/md5.h Ymir/ymir.h Ymir/zaptrackball.h Tools/zaprandom.h Tools/nucleotides.h Tools/distribution.h Delicab/antigenLib.h Delicab/epitope.h Delicab/fileformats.h Delicab/html.h Delicab/importrepertoire.h Delicab/motifFeatures.h Delicab/selfEvo.h Delicab/poolstructs.h Delicab/quality.h $(DISTDIR)/
	$(COPY_FILE) --parents Ymir/compact.cpp Ymir/lattice.cpp Ymir/proteins.cpp Ymir/receptorligand.cpp Ymir/fastaffinity.cpp Ymir/plot3d.cpp Tools/stopwatch.cpp Tools/md5.cpp Ymir/zaptrackball.cpp Tools/zaprandom.cpp Tools/nucleotides.cpp Tools/distribution.cpp Delicab/antigenLib.cpp Delicab/epitope.cpp Delicab/fileformats.cpp Delicab/html.cpp Delicab/importrepertoire.cpp Delicab/motifFeatures.cpp Delicab/selfEvo.cpp Delicab/quality.cpp Delicab/poolstructs.cpp Delicab/delimain.cpp $(DISTDIR)/


clean: compiler_clean 
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~ core *.core


distclean: clean 
	-$(DEL_FILE) $(TARGET) 
	-$(DEL_FILE) .qmake.stash
	-$(DEL_FILE) Makefile


####### Sub-libraries

mocclean: compiler_moc_header_clean compiler_moc_objc_header_clean compiler_moc_source_clean

mocables: compiler_moc_header_make_all compiler_moc_objc_header_make_all compiler_moc_source_make_all

check: first

benchmark: first

compiler_rcc_make_all:
compiler_rcc_clean:
compiler_moc_predefs_make_all: moc_predefs.h
compiler_moc_predefs_clean:
	-$(DEL_FILE) moc_predefs.h

compiler_moc_header_make_all:
compiler_moc_header_clean:
compiler_moc_objc_header_make_all:
compiler_moc_objc_header_clean:
compiler_moc_source_make_all:
compiler_moc_source_clean:
compiler_yacc_decl_make_all:
compiler_yacc_decl_clean:
compiler_yacc_impl_make_all:
compiler_yacc_impl_clean:
compiler_lex_make_all:
compiler_lex_clean:
compiler_clean: compiler_moc_predefs_clean 

####### Compile

compact.o: Ymir/compact.cpp Ymir/compact.h \
		Ymir/lattice.h \
		Tools/zaprandom.h \
		Ymir/plot3d.h \
		Ymir/proteins.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o compact.o Ymir/compact.cpp

lattice.o: Ymir/lattice.cpp Ymir/lattice.h \
		Ymir/compact.h \
		Ymir/proteins.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o lattice.o Ymir/lattice.cpp

proteins.o: Ymir/proteins.cpp Tools/zaprandom.h \
		Ymir/proteins.h \
		Ymir/compact.h \
		Ymir/lattice.h \
		Ymir/plot3d.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o proteins.o Ymir/proteins.cpp

receptorligand.o: Ymir/receptorligand.cpp Ymir/receptorligand.h \
		Ymir/compact.h \
		Ymir/lattice.h \
		Ymir/proteins.h \
		Ymir/plot3d.h \
		Tools/md5.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o receptorligand.o Ymir/receptorligand.cpp

fastaffinity.o: Ymir/fastaffinity.cpp Ymir/fastaffinity.h \
		Ymir/receptorligand.h \
		Ymir/compact.h \
		Ymir/lattice.h \
		Ymir/proteins.h \
		Tools/distribution.h \
		Tools/zaprandom.h \
		Tools/md5.h \
		Ymir/plot3d.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o fastaffinity.o Ymir/fastaffinity.cpp

plot3d.o: Ymir/plot3d.cpp Ymir/plot3d.h \
		Ymir/proteins.h \
		Ymir/compact.h \
		Ymir/lattice.h \
		Ymir/zaptrackball.h \
		Tools/zaprandom.h \
		Ymir/soil/SOIL.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o plot3d.o Ymir/plot3d.cpp

stopwatch.o: Tools/stopwatch.cpp 
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o stopwatch.o Tools/stopwatch.cpp

md5.o: Tools/md5.cpp Tools/md5.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o md5.o Tools/md5.cpp

zaptrackball.o: Ymir/zaptrackball.cpp Ymir/zaptrackball.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o zaptrackball.o Ymir/zaptrackball.cpp

zaprandom.o: Tools/zaprandom.cpp Tools/zaprandom.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o zaprandom.o Tools/zaprandom.cpp

nucleotides.o: Tools/nucleotides.cpp Tools/nucleotides.h \
		Ymir/proteins.h \
		Ymir/compact.h \
		Tools/zaprandom.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o nucleotides.o Tools/nucleotides.cpp

distribution.o: Tools/distribution.cpp Tools/zaprandom.h \
		Tools/distribution.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o distribution.o Tools/distribution.cpp

antigenLib.o: Delicab/antigenLib.cpp Delicab/antigenLib.h \
		Ymir/ymir.h \
		Ymir/lattice.h \
		Ymir/compact.h \
		Ymir/proteins.h \
		Ymir/receptorligand.h \
		Ymir/fastaffinity.h \
		Tools/distribution.h \
		Tools/zaprandom.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o antigenLib.o Delicab/antigenLib.cpp

epitope.o: Delicab/epitope.cpp Delicab/epitope.h \
		Delicab/fileformats.h \
		Ymir/compact.h \
		Ymir/plot3d.h \
		Ymir/proteins.h \
		Ymir/receptorligand.h \
		Ymir/lattice.h \
		Delicab/motifFeatures.h \
		Ymir/ymir.h \
		Ymir/fastaffinity.h \
		Delicab/importrepertoire.h \
		Delicab/poolstructs.h \
		Delicab/antigenLib.h \
		Tools/distribution.h \
		Tools/zaprandom.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o epitope.o Delicab/epitope.cpp

fileformats.o: Delicab/fileformats.cpp Delicab/fileformats.h \
		Delicab/antigenLib.h \
		Ymir/ymir.h \
		Ymir/lattice.h \
		Ymir/compact.h \
		Ymir/proteins.h \
		Ymir/receptorligand.h \
		Ymir/fastaffinity.h \
		Tools/distribution.h \
		Tools/zaprandom.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o fileformats.o Delicab/fileformats.cpp

html.o: Delicab/html.cpp Delicab/html.h \
		Ymir/ymir.h \
		Ymir/lattice.h \
		Ymir/compact.h \
		Ymir/proteins.h \
		Ymir/receptorligand.h \
		Ymir/fastaffinity.h \
		Tools/md5.h \
		Delicab/discretize.h \
		Delicab/motifFeatures.h \
		Delicab/importrepertoire.h \
		Delicab/quality.h \
		Delicab/selfEvo.h \
		Tools/zaprandom.h \
		Tools/dirent.h \
		Delicab/antigenLib.h \
		Tools/distribution.h \
		Delicab/poolstructs.h \
		Delicab/pdb.h \
		Ymir/plot3d.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o html.o Delicab/html.cpp

importrepertoire.o: Delicab/importrepertoire.cpp Delicab/importrepertoire.h \
		Delicab/antigenLib.h \
		Ymir/ymir.h \
		Ymir/lattice.h \
		Ymir/compact.h \
		Ymir/proteins.h \
		Ymir/receptorligand.h \
		Ymir/fastaffinity.h \
		Tools/distribution.h \
		Tools/zaprandom.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o importrepertoire.o Delicab/importrepertoire.cpp

motifFeatures.o: Delicab/motifFeatures.cpp Delicab/motifFeatures.h \
		Ymir/ymir.h \
		Ymir/lattice.h \
		Ymir/compact.h \
		Ymir/proteins.h \
		Ymir/receptorligand.h \
		Ymir/fastaffinity.h \
		Delicab/importrepertoire.h \
		Ymir/plot3d.h \
		Tools/dirent.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o motifFeatures.o Delicab/motifFeatures.cpp

selfEvo.o: Delicab/selfEvo.cpp Delicab/selfEvo.h \
		Ymir/ymir.h \
		Ymir/lattice.h \
		Ymir/compact.h \
		Ymir/proteins.h \
		Ymir/receptorligand.h \
		Ymir/fastaffinity.h \
		Tools/nucleotides.h \
		Tools/zaprandom.h \
		Delicab/motifFeatures.h \
		Delicab/importrepertoire.h \
		Delicab/antigenLib.h \
		Tools/distribution.h \
		Ymir/plot3d.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o selfEvo.o Delicab/selfEvo.cpp

quality.o: Delicab/quality.cpp Delicab/quality.h \
		Delicab/antigenLib.h \
		Ymir/ymir.h \
		Ymir/lattice.h \
		Ymir/compact.h \
		Ymir/proteins.h \
		Ymir/receptorligand.h \
		Ymir/fastaffinity.h \
		Tools/distribution.h \
		Tools/zaprandom.h \
		Tools/stopwatch.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o quality.o Delicab/quality.cpp

poolstructs.o: Delicab/poolstructs.cpp Delicab/poolstructs.h \
		Ymir/ymir.h \
		Ymir/lattice.h \
		Ymir/compact.h \
		Ymir/proteins.h \
		Ymir/receptorligand.h \
		Ymir/fastaffinity.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o poolstructs.o Delicab/poolstructs.cpp

delimain.o: Delicab/delimain.cpp Ymir/ymir.h \
		Ymir/lattice.h \
		Ymir/compact.h \
		Ymir/proteins.h \
		Ymir/receptorligand.h \
		Ymir/fastaffinity.h \
		Tools/md5.h \
		Delicab/motifFeatures.h \
		Delicab/importrepertoire.h \
		Delicab/quality.h \
		Delicab/selfEvo.h \
		Tools/zaprandom.h \
		Tools/dirent.h \
		Delicab/antigenLib.h \
		Tools/distribution.h \
		Delicab/poolstructs.h \
		Delicab/epitope.h \
		Delicab/fileformats.h \
		Delicab/html.h \
		Delicab/discretize.h \
		Delicab/pdb.h \
		Ymir/plot3d.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o delimain.o Delicab/delimain.cpp

topology.o: Delicab/topology.cpp Delicab/topology.h \
		Ymir/ymir.h \
		Ymir/lattice.h \
		Ymir/compact.h \
		Ymir/proteins.h \
		Ymir/receptorligand.h \
		Ymir/fastaffinity.h \
		Ymir/plot3d.h \
		Delicab/fileformats.h \
		Delicab/motifFeatures.h \
		Delicab/importrepertoire.h \
		Delicab/antigenLib.h \
		Tools/distribution.h \
		Tools/zaprandom.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o topology.o Delicab/topology.cpp

####### Install

install:  FORCE

uninstall:  FORCE

FORCE:

