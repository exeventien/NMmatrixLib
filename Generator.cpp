#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <Generator.h>

static bool RNG_SEEDED = 0;

Generator::Generator(int nSize, int mSize){
	if(!RNG_SEEDED){
		srand(time(0));
		RNG_SEEDED = 1;
	}
	n = nSize;
	m = mSize;
	rows = new Nvector*[m+1];
}

Generator::~Generator(){
	for(int i = 0; i < m+1; i++)	
		delete rows[i];
	delete[] rows;
}

Nvector* Generator::generateRow(){
	Nvector* temp = new Nvector(n);
	for(int i = 0; i < n; i++)
		temp->set(i, (double)((rand()%10000)/100.00 * (rand()%2 ? 1.0 : -1.0))); 	
	return temp;
}

void Generator::createDiagonalElement(int row, Nvector* vec){
	double rowSum = 0;
	for(int i = 0; i < n; i++)
		rowSum += fabs(vec->get(i));
	vec->set(row, rowSum);
}

void Generator::createMatrix(){
	for(int i = 0; i < m; i++){
		Nvector* entry = generateRow();	
		createDiagonalElement(i, entry);
		rows[i] = entry;
	}
}

void Generator::printFile(const char* path){
	FILE* file = fopen(path, "w+");
	fprintf(file, "%i\n%i\n", n, m);
	rows[m] = generateRow();  // vector b appended to the end of the file
	for(int i = 0; i < m+1; i++){
		for(int j = 0; j < n; j++){
			fprintf(file, "%f ", rows[j]->get(i));
		}
		fprintf(file, "\n");
	}
}
