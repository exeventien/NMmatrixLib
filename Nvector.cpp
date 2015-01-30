#include <cstdio>
#include <cstring>
#include <Nvector.h>

Nvector::Nvector(int nSize){
	n = nSize;
	vec = new double[n];
}

Nvector::Nvector(const Nvector& copy){
	n = copy.n;
	vec = new double[n];
	memcpy(vec, copy.vec, sizeof(double)*n);
}

Nvector::~Nvector(){
	delete[] vec;
}

void Nvector::scalarMultiply(double scalar){
	for(int i = 0; i < n; i++)
		vec[i] *= scalar;
}

void Nvector::addVector(const Nvector& bvec){
	for(int i = 0; i < n; i++)
		vec[i] += bvec.vec[i];	
}

double Nvector::dotProduct(const Nvector& b){
	double sum = 0.0;
	for(int i = 0; i < n; i++)
		sum += vec[i] * b.vec[i];
	return sum;
}

int Nvector::getSize() const{
	return n;
}

void Nvector::set(int nElement, double value){
	vec[nElement] = value;
}

double Nvector::get(int nElement) const{
	return vec[nElement];
}

void Nvector::copy(const Nvector& copy){
	delete[] vec;
	n = copy.n;
	vec = new double[n];
	memcpy(vec, copy.vec, sizeof(double)*n);
}

void Nvector::print(){
	printf("\n[");
	for(int i = 0; i < n; i++)
		printf(" %0.4f", vec[i]);
	printf(" ]\n");
}
