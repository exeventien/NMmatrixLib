#include <cstdio>
#include <cstring>
#include <Nvector.h>
#include <NMmatrix.h>

NMmatrix::NMmatrix(int nSize, int mSize){
	n = nSize;
	m = mSize;
	mat = new double[n*m];		
}

NMmatrix::NMmatrix(const NMmatrix& copy){
	n = copy.n;
	m = copy.m;
	mat = new double[n*m];
	memcpy(mat, copy.mat, sizeof(double)*n*m);
}

NMmatrix::~NMmatrix(){
	delete[] mat;
}

void NMmatrix::vectorTransform(Nvector& vec){
	Nvector temp(vec);
	for(int i = 0; i < m; i++){
		double rowM = 0;
		for(int j = 0; j < n; j++){
			rowM += mat[(n*i)+j] * temp.get(j);
		}
		vec.set(i, rowM);
	}
}

void NMmatrix::diagonalInverse(){
	if(n == m){
		for(int i = 0; i < m; i++){
			for(int j = 0; j < n; j++){
				if(i != j)
					mat[(n*i)+j] = 0.0;
				else if(mat[(n*i)+j] != 0.0)
					mat[(n*i)+j] = 1.0 / mat[(n*i)+j];
			}
		}
	}
}


void NMmatrix::diagonalRemainder(){
	if(n == m){
		for(int i = 0; i < m; i++){
			for(int j = 0; j < n; j++){
				if(i == j)
					mat[(n*i)+j] = 0.0;
			}
		}
	}
}

void NMmatrix::inverseDiagonalBelow(){
	for(int i = 0; i < m; i++)
		for(int j = 0; j < n; j++){
			if(j > i)
				mat[(n*i)+j] = 0.0;
			else if(mat[(n*i)+j] != 0.0)
				mat[(n*i)+j] = 1.0/mat[(n*i)+j];
		}
}

void NMmatrix::aboveDiagonal(){
        for(int i = 0; i < m; i++)
                for(int j = 0; j < n; j++)
                        if(j <= i)
                                mat[(n*i)+j] = 0.0;
}

void NMmatrix::getRow(int row, Nvector& vec){
	for(int i = 0; i < n; i++)
		vec.set(i, mat[(row*n)+i]);
}

void NMmatrix::set(int nElement, int mElement, double value){
	mat[(mElement*n)+nElement] = value;
}

void NMmatrix::copy(const NMmatrix& copy){
        n = copy.n;
        m = copy.m;
	delete[] mat;
        mat = new double[n*m];
        memcpy(mat, copy.mat, sizeof(double)*n*m);
}

void NMmatrix::print() const{
	for(int i = 0; i < m; i++){
		printf("\n[");
		for(int j = 0; j < n; j++){
			printf(" %2.4f", mat[(n*i)+j]);
		}
		printf(" ]");
	}
	printf("\n");
}
