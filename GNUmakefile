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
	g++ $(COMMON_FLAGS) $(LINKER_FLAGS) $(SOLVER_SOURCE) -o matrixSolver

generator: lib
	g++ $(COMMON_FLAGS) $(LINKER_FLAGS) $(GENERATOR_SOURCE) -o matrixGenerator



