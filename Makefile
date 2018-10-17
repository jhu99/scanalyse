OS = $(shell uname -s)
CXX = g++
DEBUG = yes

# Default mode is "Release"
DEFAULT_MODE  = Release
MODE         ?= $(DEFAULT_MODE)
# If mode is something other than "Debug" or "Release",throw a fit
ifneq ($(MODE),Debug)
ifneq ($(MODE),Release)
$(error MODE must be one of {Debug,Release})
endif
endif
# Flags for compiling the code
ifeq ($(MODE),Debug)
CXXFLAGS = -Wall -g3 -DDEBUG -std=c++0x -DVERBOSE -Ilib/
else
CXXFLAGS = -Wall -O3 -ffast-math -Ilib/ -std=c++0x -DNDEBUG
endif
# all: cellTest funTest ArgParserTest
cellTest: tests/cellTest.cpp ./lib/cell/cell.o
	${CXX} ${CXXFLAGS} -o $@ $^
funTest: tests/funTest.cpp lib/cell/cell.o lib/fun/fun.o
	${CXX} ${CXXFLAGS} -o $@ $^
ArgParserTest: test/testargparser.cpp lib/argparser.o
	${CXX} ${CXXFLAGS} -o $@ $^
clean:
	rm lib/*.o
