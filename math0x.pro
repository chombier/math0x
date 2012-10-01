
TEMPLATE = app
TARGET = bin/test

INCLUDEPATH += .. /usr/include/eigen3
DEPENDPATH += . .. 

# QMAKE_CXX = /usr/lib/gcc-snapshot/bin/g++
# QMAKE_CXX = clang++

QMAKE_CXXFLAGS += -std=c++11

SOURCES = \
    main.cpp \
#    test.cpp \
    
# CONFIG += debug
