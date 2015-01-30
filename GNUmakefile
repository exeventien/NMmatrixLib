all:	solver generator

solver_source = MainApp.cpp NMmatrix.cpp Nvector.cpp MatrixMethods.cpp MatrixReader.cpp

solver:
	g++ -g -I. $(solver_source) -o matrixSolver

generator_source = Generator.cpp GeneratorApp.cpp Nvector.cpp

generator:
	g++ -I. $(generator_source) -o matrixGenerator



