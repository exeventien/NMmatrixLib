#include <Nvector.h>

#ifndef GENERATOR
#define GENERATOR

class Generator{
private:
	int n, m;	
	Nvector** rows;

        Nvector* generateRow();
        void createDiagonalElement(int, Nvector*);
public:
	Generator(int, int);
	~Generator();

	void printFile(const char*);
	void createMatrix();
};

#endif

