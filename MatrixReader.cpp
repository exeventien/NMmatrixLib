#include <cstring>
#include <cstdlib>
#include <MatrixReader.h>

MatrixReader::MatrixReader(const char* name){
	char strN[5], strM[5];
	file = fopen(name, "r");
	fscanf(file, "%s %s", strN, strM);
	n = atoi(strN);
	m = atoi(strM);
}

MatrixReader::~MatrixReader(){
	fclose(file);
}

int MatrixReader::getN(){
	return n;
}

int MatrixReader::getM(){
	return m;
}

void MatrixReader::readMatrix(NMmatrix& mat){
	char value[20];
	for(int i = 0; i < m; i++)
		for(int j = 0; j < n; j++){
			fscanf(file, "%s", value);
			mat.set(j, i, atof(value));
		}
}

void MatrixReader::readVector(Nvector& vec){
	char value[20];
	for(int i = 0; i < n; i++){
		fscanf(file, "%s", value);
		vec.set(i, atof(value));
	}
}
