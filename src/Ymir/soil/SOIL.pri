# SOIL is a library to export or load pictures from openGL. This can be used to make pictures or automated videos of receptors.
# this function allows to do include(SOIL.pri) from another project
INCLUDEPATH += $$PWD
DEPENDPATH += $$PWD
message("Info: SOIL Library included from folder :")
message($$PWD)

# Necessary only for boost - would also be necessary for the syntax vector<int> bla = {1, 2, 3};
QMAKE_CXXFLAGS += -std=c++11

HEADERS += \
	$$PWD/image_DXT.h  \
	$$PWD/image_helper.h\
	$$PWD/SOIL.h\
	$$PWD/stb_image_aug.h\
	$$PWD/stbi_DDS_aug.h\
	$$PWD/stbi_DDS_aug_c.h
	
SOURCES += \
    $$PWD/image_DXT.c \
	$$PWD/image_helper.c\
	$$PWD/SOIL.c\
	$$PWD/stb_image_aug.c\

unix: LIBS += -lGL

win32: {
    LIBS += -LC:/MyPrograms/freeglut-3.2.1/bin/ -lfreeglut
    INCLUDEPATH += C:/MyPrograms/freeglut-3.2.1/include/
    DEPENDPATH += C:/MyPrograms/freeglut-3.2.1/lib/
    LIBS +=  -lopengl32
}
