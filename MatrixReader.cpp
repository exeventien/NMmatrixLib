#include <cstring>
#include <cstdlib>
#include <MatrixReader.h>

MatrixReader::MatrixReader(const char* name){
	char strN[5], strM[5];
	file = fopen(name, "r");
	fscanf(file, "%s %s", strN, strM);
	n = atoi(strN);
	m = atoi(strM);
	imat = new NMmatrix(n,m);
	ivec = new Nvector(n);
}

MatrixReader::~MatrixReader(){
	delete imat;
	delete ivec;
	fclose(file);
}

int MatrixReader::getN(){
	return n;
}

int MatrixReader::getM(){
	return m;
}

void MatrixReader::readMatrix(NMmatrix& mat){
	static bool matrixRead = 0;
	char value[20];
	if(matrixRead == 0){
		for(int i = 0; i < m; i++)
			for(int j = 0; j < n; j++){
				fscanf(file, "%s", value);
				imat->set(j, i, atof(value));
			}
		matrixRead = 1;
	}
	mat.copy(*imat);
}

void MatrixReader::readVector(Nvector& vec){
	static bool vectorRead = 0;
	char value[20];
	if(vectorRead == 0){
		for(int i = 0; i < n; i++){
			fscanf(file, "%s", value);
			ivec->set(i, atof(value));
		}
		vectorRead = 1;
	}
	vec.copy(*ivec);
}
