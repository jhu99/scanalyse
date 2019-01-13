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
all: cellTest funTest linearRegressionTest argParserTest linearRegressionParameterTest qqNormTest SparseMatrixTest geneTopsTest fetch_batchTest FiltrationTest
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
SparseMatrix.o:lib/SparseMatrix/SparseMatrix.cpp 
	${H5CXX} $^ ${H5CXXFLAGS} -c -o $@
SparseMatrixTest:tests/SparseMatrixTest.cpp SparseMatrix.o lib/rank/geneInfo.o lib/rank/rankNormalize.o lib/qqNorm/qqNorm.o lib/qqNorm/caculateInterface.o
	${H5CXX} $^ ${CXXGSLFLAGS} -o $@  
geneTopsTest:tests/geneTopsTest.cpp lib/geneTop/geneExpressionTop.o lib/geneTop/geneCount.o SparseMatrix.o lib/qqNorm/qqNorm.o
	${H5CXX} $^ ${CXXGSLFLAGS} -o $@ 
rankTest:tests/rankTest.cpp lib/rank/rankNormalize.o lib/rank/geneInfo.o SparseMatrix.o lib/qqNorm/qqNorm.o
	${H5CXX} $^ ${CXXGSLFLAGS} -o $@ 
fetch_batchTest:tests/fetch_batchTest.cpp SparseMatrix.o lib/qqNorm/qqNorm.o
	${H5CXX} $^ ${CXXGSLFLAGS} -o $@ 
geneVariationTest:tests/geneVariationTest.cpp lib/geneVariationTop/geneVariation.o lib/geneVariationTop/geneVariationTop.o lib/geneTop/geneExpressionTop.o lib/geneTop/geneCount.o SparseMatrix.o lib/qqNorm/qqNorm.o
	${H5CXX} $^ ${CXXGSLFLAGS} -o $@ 
FiltrationTest:tests/FiltrationTest.cpp lib/Filtration/Filtration.o lib/geneVariationTop/geneVariation.o lib/geneVariationTop/geneVariationTop.o lib/geneTop/geneCount.o SparseMatrix.o lib/qqNorm/qqNorm.o
	${H5CXX} $^ ${CXXGSLFLAGS} -o $@
clean:
	rm lib/cell/*.o lib/fun/*.o lib/linearRegression/*.o lib/argparser/*.o lib/qqNorm/*.o lib/SparseMatrix/*.o lib/geneTop/*.o lib/rank/*.o *.o cellTest funTest linearRegressionTest argParserTest qqNormTest geneTopsTest HDF5ReaderTest rankTest fetch_batchTest FiltrationTest
