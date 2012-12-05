
TEMPLATE = app
TARGET = bin/test

INCLUDEPATH += .. /usr/include/eigen3
DEPENDPATH += . .. 

headers.path = $$PREFIX/include/math0x
headers.files = *.h tuple func
header.depends = clean

QMAKE_CLEAN += tuple/*~ func/*~

INSTALLS += headers

# QMAKE_CXX = /usr/lib/gcc-snapshot/bin/g++
# QMAKE_CXX = clang++

QMAKE_CXXFLAGS += -std=c++11

QMAKE_CXXFLAGS_RELEASE += -DNDEBUG

SOURCES = \
    main.cpp \
#    test.cpp \
    
# CONFIG += debug
