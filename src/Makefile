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

INCPATH       = -IAbsolute -I. -IYmir -I. -isystem /usr/include/libdrm 
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
DISTNAME      = AbsolutNoLib1.0.0

LINK          = g++
LINKMPICXX    = mpicxx

LFLAGS        = -Wl
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
		Absolut/antigenLib.cpp \
		Absolut/epitope.cpp \
		Absolut/fileformats.cpp \
		Absolut/html.cpp \
		Absolut/importrepertoire.cpp \
		Absolut/motifFeatures.cpp \
		Absolut/selfEvo.cpp \
		Absolut/quality.cpp \
		Absolut/poolstructs.cpp \
		Absolut/delimain.cpp \
		Absolut/topology.cpp \
		Absolut/dlab.cpp
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
		topology.o \
		dlab.o
DIST          = Absolut/.qmake.stash \
		Absolut/Ymir/YmirNoGL.pri \
		Ymir/compact.h \
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
		Absolut/antigenLib.h \
		Absolut/epitope.h \
		Absolut/fileformats.h \
		Absolut/html.h \
		Absolut/importrepertoire.h \
		Absolut/motifFeatures.h \
		Absolut/selfEvo.h \
		Absolut/poolstructs.h \
		Absolut/quality.h \
		Absolut/topology.h \
		Absolut/dlab.h \
		Ymir/compact.cpp \
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
		Absolut/antigenLib.cpp \
		Absolut/epitope.cpp \
		Absolut/fileformats.cpp \
		Absolut/html.cpp \
		Absolut/importrepertoire.cpp \
		Absolut/motifFeatures.cpp \
		Absolut/selfEvo.cpp \
		Absolut/quality.cpp \
		Absolut/poolstructs.cpp \
		Absolut/topology.cpp \
		Absolut/dlab.cpp \
		Absolut/delimain.cpp
QMAKE_TARGET  = AbsolutNoLib
DESTDIR       = 
TARGET        = AbsolutNoLib
TARGETMPI     = AbominationMPI
TARGETFULL    = Absolut



first: all
####### Build rules

AbsolutNoLib:  $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)

AbsolutNoLibMPI:  $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)


Makefile: Absolut/AbsolutNoLib.pro .qmake.stash \
		Ymir/YmirNoGL.pri \
		Absolut/AbsolutNoLib.pro
	$(QMAKE) -o Makefile Absolut/AbsolutNoLib.pro
Makefile: TARGET = $(TARGETFULL)

	
.qmake.stash:
Ymir/YmirNoGL.pri:
Absolut/AbsolutNoLib.pro:

all: clean Makefile AbsolutNoLib 

full: clean
full: CC = $(CC)
full: CXX = $(CXX)
full: LINK = $(LINKCXX)
full: CXXFLAGS = $(CXXFLAGS)
full: TARGET = $(TARGETFULL)
full: Absolut

MPIcxx: clean
MPIcxx: CC = $(MPICC)
MPIcxx: CXX = $(MPICXX)
MPIcxx: LINK = $(LINKMPICXX)
MPIcxx: CXXFLAGS = $(CXXFLAGSMPI)
MPIcxx: TARGET = $(TARGETMPI)
MPIcxx: AbsolutNoLibMPI

MPIc++: clean
MPIc++: CC = $(MPICC)
MPIc++: CXX = mpic++
MPIc++: LINK = mpic++
MPIc++: CXXFLAGS = $(CXXFLAGSMPI)
MPIc++: TARGET = $(TARGETMPI)
MPIc++: AbsolutNoLibMPI

MPIgxx: clean
MPIgxx: CC = mpigcc
MPIgxx: CXX = mpigxx
MPIgxx: LINK = mpigxx
MPIgxx: CXXFLAGS = $(CXXFLAGSMPI)
MPIgxx: TARGET = $(TARGETMPI)
MPIgxx: AbsolutNoLibMPI

dist: distdir FORCE
	(cd `dirname $(DISTDIR)` && $(TAR) $(DISTNAME).tar $(DISTNAME) && $(COMPRESS) $(DISTNAME).tar) && $(MOVE) `dirname $(DISTDIR)`/$(DISTNAME).tar.gz . && $(DEL_FILE) -r $(DISTDIR)

distdir: FORCE
	@test -d $(DISTDIR) || mkdir -p $(DISTDIR)
	$(COPY_FILE) --parents $(DIST) $(DISTDIR)/
	$(COPY_FILE) --parents Ymir/compact.h Ymir/lattice.h Ymir/proteins.h Ymir/receptorligand.h Ymir/fastaffinity.h Ymir/plot3d.h Tools/stopwatch.h Tools/md5.h Ymir/ymir.h Ymir/zaptrackball.h Tools/zaprandom.h Tools/nucleotides.h Tools/distribution.h Absolut/antigenLib.h Absolut/epitope.h Absolut/fileformats.h Absolut/html.h Absolut/importrepertoire.h Absolut/motifFeatures.h Absolut/selfEvo.h Absolut/poolstructs.h Absolut/quality.h $(DISTDIR)/
	$(COPY_FILE) --parents Ymir/compact.cpp Ymir/lattice.cpp Ymir/proteins.cpp Ymir/receptorligand.cpp Ymir/fastaffinity.cpp Ymir/plot3d.cpp Tools/stopwatch.cpp Tools/md5.cpp Ymir/zaptrackball.cpp Tools/zaprandom.cpp Tools/nucleotides.cpp Tools/distribution.cpp Absolut/antigenLib.cpp Absolut/epitope.cpp Absolut/fileformats.cpp Absolut/html.cpp Absolut/importrepertoire.cpp Absolut/motifFeatures.cpp Absolut/selfEvo.cpp Absolut/quality.cpp Absolut/poolstructs.cpp Absolut/delimain.cpp $(DISTDIR)/


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

antigenLib.o: Absolut/antigenLib.cpp Absolut/antigenLib.h \
		Ymir/ymir.h \
		Ymir/lattice.h \
		Ymir/compact.h \
		Ymir/proteins.h \
		Ymir/receptorligand.h \
		Ymir/fastaffinity.h \
		Tools/distribution.h \
		Tools/zaprandom.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o antigenLib.o Absolut/antigenLib.cpp

epitope.o: Absolut/epitope.cpp Absolut/epitope.h \
		Absolut/fileformats.h \
		Ymir/compact.h \
		Ymir/plot3d.h \
		Ymir/proteins.h \
		Ymir/receptorligand.h \
		Ymir/lattice.h \
		Absolut/motifFeatures.h \
		Ymir/ymir.h \
		Ymir/fastaffinity.h \
		Absolut/importrepertoire.h \
		Absolut/poolstructs.h \
		Absolut/antigenLib.h \
		Tools/distribution.h \
		Tools/zaprandom.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o epitope.o Absolut/epitope.cpp

fileformats.o: Absolut/fileformats.cpp Absolut/fileformats.h \
		Absolut/antigenLib.h \
		Ymir/ymir.h \
		Ymir/lattice.h \
		Ymir/compact.h \
		Ymir/proteins.h \
		Ymir/receptorligand.h \
		Ymir/fastaffinity.h \
		Tools/distribution.h \
		Tools/zaprandom.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o fileformats.o Absolut/fileformats.cpp

html.o: Absolut/html.cpp Absolut/html.h \
		Ymir/ymir.h \
		Ymir/lattice.h \
		Ymir/compact.h \
		Ymir/proteins.h \
		Ymir/receptorligand.h \
		Ymir/fastaffinity.h \
		Tools/md5.h \
		Absolut/discretize.h \
		Absolut/motifFeatures.h \
		Absolut/importrepertoire.h \
		Absolut/quality.h \
		Absolut/selfEvo.h \
		Tools/zaprandom.h \
		Tools/dirent.h \
		Absolut/antigenLib.h \
		Tools/distribution.h \
		Absolut/poolstructs.h \
		Absolut/pdb.h \
		Ymir/plot3d.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o html.o Absolut/html.cpp

importrepertoire.o: Absolut/importrepertoire.cpp Absolut/importrepertoire.h \
		Absolut/antigenLib.h \
		Ymir/ymir.h \
		Ymir/lattice.h \
		Ymir/compact.h \
		Ymir/proteins.h \
		Ymir/receptorligand.h \
		Ymir/fastaffinity.h \
		Tools/distribution.h \
		Tools/zaprandom.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o importrepertoire.o Absolut/importrepertoire.cpp

motifFeatures.o: Absolut/motifFeatures.cpp Absolut/motifFeatures.h \
		Ymir/ymir.h \
		Ymir/lattice.h \
		Ymir/compact.h \
		Ymir/proteins.h \
		Ymir/receptorligand.h \
		Ymir/fastaffinity.h \
		Absolut/importrepertoire.h \
		Ymir/plot3d.h \
		Tools/dirent.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o motifFeatures.o Absolut/motifFeatures.cpp

selfEvo.o: Absolut/selfEvo.cpp Absolut/selfEvo.h \
		Ymir/ymir.h \
		Ymir/lattice.h \
		Ymir/compact.h \
		Ymir/proteins.h \
		Ymir/receptorligand.h \
		Ymir/fastaffinity.h \
		Tools/nucleotides.h \
		Tools/zaprandom.h \
		Absolut/motifFeatures.h \
		Absolut/importrepertoire.h \
		Absolut/antigenLib.h \
		Tools/distribution.h \
		Ymir/plot3d.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o selfEvo.o Absolut/selfEvo.cpp

quality.o: Absolut/quality.cpp Absolut/quality.h \
		Absolut/antigenLib.h \
		Ymir/ymir.h \
		Ymir/lattice.h \
		Ymir/compact.h \
		Ymir/proteins.h \
		Ymir/receptorligand.h \
		Ymir/fastaffinity.h \
		Tools/distribution.h \
		Tools/zaprandom.h \
		Tools/stopwatch.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o quality.o Absolut/quality.cpp

poolstructs.o: Absolut/poolstructs.cpp Absolut/poolstructs.h \
		Ymir/ymir.h \
		Ymir/lattice.h \
		Ymir/compact.h \
		Ymir/proteins.h \
		Ymir/receptorligand.h \
		Ymir/fastaffinity.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o poolstructs.o Absolut/poolstructs.cpp

delimain.o: Absolut/delimain.cpp Ymir/ymir.h \
		Ymir/lattice.h \
		Ymir/compact.h \
		Ymir/proteins.h \
		Ymir/receptorligand.h \
		Ymir/fastaffinity.h \
		Tools/md5.h \
		Absolut/motifFeatures.h \
		Absolut/importrepertoire.h \
		Absolut/quality.h \
		Absolut/selfEvo.h \
		Tools/zaprandom.h \
		Tools/dirent.h \
		Absolut/antigenLib.h \
		Tools/distribution.h \
		Absolut/poolstructs.h \
		Absolut/epitope.h \
		Absolut/fileformats.h \
		Absolut/html.h \
		Absolut/discretize.h \
		Absolut/pdb.h \
		Ymir/plot3d.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o delimain.o Absolut/delimain.cpp

dlab.o: Absolut/dlab.cpp Absolut/dlab.h \
		Ymir/ymir.h \
		Ymir/lattice.h \
		Ymir/compact.h \
		Ymir/proteins.h \
		Ymir/receptorligand.h \
		Ymir/fastaffinity.h \
		Tools/md5.h \
		Absolut/motifFeatures.h \
		Absolut/importrepertoire.h \
		Absolut/quality.h \
		Absolut/selfEvo.h \
		Tools/zaprandom.h \
		Tools/dirent.h \
		Absolut/antigenLib.h \
		Tools/distribution.h \
		Absolut/poolstructs.h \
		Absolut/epitope.h \
		Absolut/fileformats.h \
		Absolut/html.h \
		Absolut/topology.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o dlab.o Absolut/dlab.cpp

topology.o: Absolut/topology.cpp Absolut/topology.h \
		Ymir/ymir.h \
		Ymir/lattice.h \
		Ymir/compact.h \
		Ymir/proteins.h \
		Ymir/receptorligand.h \
		Ymir/fastaffinity.h \
		Ymir/plot3d.h \
		Absolut/fileformats.h \
		Absolut/motifFeatures.h \
		Absolut/importrepertoire.h \
		Absolut/antigenLib.h \
		Tools/distribution.h \
		Tools/zaprandom.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o topology.o Absolut/topology.cpp

####### Install

install:  FORCE

uninstall:  FORCE

FORCE:

