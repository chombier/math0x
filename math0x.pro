# place customizations in .qmake.cache

TEMPLATE = app
TARGET = bin/test

# depends
CONFIG -= qt
CONFIG += link_pkgconfig
PKGCONFIG += eigen3

# compiler flags
DEFINES = 
QMAKE_CXXFLAGS += -std=c++11 
QMAKE_CXXFLAGS_RELEASE += -DNDEBUG

# sources, headers
INCLUDEPATH += ..
DEPENDPATH += . .. 

SOURCES = \
    main.cpp \
#		test.cpp 

# install targets
headers.path = $$PREFIX/include/math0x
headers.files = *.h tuple func
header.depends = clean

QMAKE_CLEAN += tuple/*~ func/*~

INSTALLS += headers
