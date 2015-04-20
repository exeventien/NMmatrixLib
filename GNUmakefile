all:	solver generator

LIB_SOURCE = NMmatrix.cpp Nvector.cpp MatrixMethods.cpp MatrixReader.cpp Generator.cpp
SOLVER_SOURCE = MainApp.cpp
GENERATOR_SOURCE  = GeneratorApp.cpp

SHARED_FLAGS = -fPIC -shared 
LINKER_FLAGS = -L. -lnmmatlib -Wl,-rpath,. 
COMMON_FLAGS = -g -I. 


lib:
	g++ $(COMMON_FLAGS) $(SHARED_FLAGS) $(LIB_SOURCE) -o libnmmatlib.so

solver:	lib
	g++ $(COMMON_FLAGS) $(SOLVER_SOURCE) $(LINKER_FLAGS) -o matrixSolver

generator: lib
	g++ $(COMMON_FLAGS) $(GENERATOR_SOURCE) $(LINKER_FLAGS) -o matrixGenerator



