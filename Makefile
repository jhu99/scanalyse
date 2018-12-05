OS = $(shell uname -s)
CXX = g++
H5CXX = h5c++
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
CXXFLAGS = -Wall -g3 -DDEBUG -std=c++0x -DVERBOSE -Ilib/ -I/usr/local/include/ -I/usr/include/hdf5/serial/
CXXGSLFLAGS = -Wall -g3 -DDEBUG -std=c++0x -DVERBOSE -Ilib/ -I/usr/local/include/ -lgsl -lgslcblas -lpthread
else
CXXFLAGS = -Wall -O3 -ffast-math -Ilib/ -I/usr/include/hdf5/serial/ -std=c++0x -DNDEBUG
H5CXXFLAGS = -std=c++0x -Ilib/
CXXGSLFLAGS = -Wall -O3 -ffast-math -Ilib/ -I/usr/local/include/ -std=c++0x -DNDEBUG -lgsl -lgslcblas -lpthread
endif
all: cellTest funTest linearRegressionTest argParserTest linearRegressionParameterTest qqNormTest HDF5ReaderTest geneTopsTest
#ArgParserTest
cellTest: tests/cellTest.cpp lib/cell/cell.o
	${CXX} $^ ${CXXFLAGS} -o $@ 
funTest: tests/funTest.cpp lib/cell/cell.o lib/fun/fun.o
	${CXX} $^ ${CXXFLAGS} -o $@ 
linearRegressionTest: tests/linearRegressionTest.cpp lib/cell/cell.o lib/linearRegression/linearRegressionParameter.o lib/linearRegression/linearRegression.o
	${CXX} $^ ${CXXGSLFLAGS} -o $@ 
linearRegressionParameterTest: tests/linearRegressionParameterTest.cpp lib/linearRegression/linearRegressionParameter.o
	${CXX} $^ ${CXXGSLFLAGS} -o $@
argParserTest: tests/argparsertest.cpp lib/argparser/argparser.o
	${CXX} $^ ${CXXFLAGS} -o $@ 
qqNormTest:tests/qqNormTest.cpp lib/qqNorm/qqNorm.o lib/qqNorm/caculateInterface.o
	${CXX} $^ ${CXXGSLFLAGS} -o $@ 
HDF5Reader.o:lib/HDF5Reader/HDF5Reader.cpp 
	${H5CXX} $^ ${H5CXXFLAGS} -c -o $@
HDF5ReaderTest:tests/testHDF5Reader.cpp HDF5Reader.o
	${H5CXX} $^ ${CXXFLAGS} -o $@  
geneTopsTest:tests/geneTopsTest.cpp lib/geneTop/geneExpressionTop.o lib/geneTop/geneCount.o HDF5Reader.o
	${H5CXX} $^ ${CXXFLAGS} -o $@ 
rankTest:tests/rankTest.cpp lib/rank/rankNormalize.o lib/rank/geneInfo.o HDF5Reader.o
	${H5CXX} $^ ${CXXFLAGS} -o $@ 
clean:
	rm lib/cell/*.o lib/fun/*.o lib/linearRegression/*.o lib/argparser/*.o lib/qqNorm/*.o lib/HDF5Reader/*.o lib/geneTop/*.o lib/rank/*.o *.o cellTest funTest linearRegressionTest argParserTest qqNormTest geneTopsTest HDF5ReaderTest rankTest
