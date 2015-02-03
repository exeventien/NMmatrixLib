#include <cstdio>
#include <Nvector.h>
#include <NMmatrix.h>

#ifndef MATRIX_READER
#define MATRIX_READER
class MatrixReader{
	private:
	NMmatrix* imat;
	Nvector* ivec;
	FILE* file;
	int n, m;	
	
	public:
	MatrixReader(const char*);
	~MatrixReader();	

	int getN();
	int getM();

	void readMatrix(NMmatrix&);
	void readVector(Nvector&);
};
#endif
