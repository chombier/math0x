
TEMPLATE = app
TARGET = bin/test

INCLUDEPATH += .. /usr/include/eigen3
DEPENDPATH += . .. 

QMAKE_CXXFLAGS += -std=c++11

SOURCES = main.cpp

