
# TeX sources of the documentation
DOC_SOURCE =	\
				manual-LatFit.tex \
				manual-LatFold.tex \
				manual-LatFoldVec.tex \
				manual-LatVec.tex 

# TeX sources source files
DOC_SOURCE_SRC =	\
					fit_both.png \
					fit_lat_sc.png \
					fit_pdb.png \
					algorithmicx.sty \
					algpseudocode.sty \
					square-stucked.png \
					version-info.tex.in


################################################################################
# IF DOCUMENTATION HAS TO BE COMPILED
if build_docs
################################################################################

# documentation target files 
DOC_TARGETS = $(DOC_SOURCE:.tex=.pdf)

# temporary files created during documentation compilation
DOC_TARGET_TMP_FILES =	$(DOC_SOURCE:.tex=.aux) \
						$(DOC_SOURCE:.tex=.bbl) \
						$(DOC_SOURCE:.tex=.blg) \
						$(DOC_SOURCE:.tex=.log) 

# files to delete with clean
MOSTLYCLEANFILES =	$(DOC_TARGETS) $(DOC_TARGET_TMP_FILES)

# compilation of TeX file via pdflatex and removal of temporary files
.tex.pdf: $(DOC_SOURCE) $(DOC_SOURCE_SRC)
	pdflatex $<
	pdflatex $<
	$(RM) $(DOC_TARGET_TMP_FILES)
	
all : $(DOC_SOURCE) $(DOC_SOURCE_SRC) $(DOC_TARGETS)

EXTRA_DIST = $(DOC_SOURCE) $(DOC_SOURCE_SRC) $(DOC_TARGETS) 

################################################################################
else
################################################################################

# DO NOTHING

EXTRA_DIST = $(DOC_SOURCE) $(DOC_SOURCE_SRC) 

################################################################################
endif
################################################################################

