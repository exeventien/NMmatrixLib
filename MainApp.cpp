#include <cstdlib>
#include <cstdio>
#include <Nvector.h>
#include <NMmatrix.h>
#include <MatrixMethods.h>
#include <MatrixReader.h>

static const char JACOBI[] = "1.) Jacobi Method";
static const char GAUSS_SEIDEL[] = "2.) Gauss-Seidel Method";


void options(int, MatrixReader&);

int main(int argv, const char** argc){
	if(argv >= 2){
		MatrixReader reader(argc[1]);
		int option;
		while(1){
			printf("\n\nChoose a method to solve the matrix.\n%s\n%s\n\n>", JACOBI, GAUSS_SEIDEL);
			scanf("%i", &option);		
			options(option, reader);
		}
	}
	return 0;
}

void options(int method, MatrixReader& reader){
        NMmatrix mat(reader.getN(), reader.getM());
        reader.readMatrix(mat);

        Nvector vec(reader.getN());
        reader.readVector(vec);

       	Nvector x(reader.getN());
	
	switch(method){
		case 1:
			MatrixMethods::jacobiMethod(mat, vec, x);
			x.print();
			break;
		case 2:
			MatrixMethods::gaussSeidelMethod(mat, vec, x);
			x.print();
			break;
		case 3:
		
			break;
		case 4:

			break;			
		case 5:
		
			break;
		default:
			exit(0);
	}

}