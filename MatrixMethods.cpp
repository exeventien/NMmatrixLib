#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <Nvector.h>
#include <NMmatrix.h>
#include <MatrixMethods.h>

void MatrixMethods::jacobiMethod(const NMmatrix& mat, const Nvector& b, Nvector& x1){
	NMmatrix d(mat), r(mat);
	mat.print();
	Nvector x0(b.getSize());
	d.diagonalInverse();
	d.print();
	r.diagonalRemainder();
	r.print();
	randomizeVector(x1);
	x1.print();
	while(!vectorsConverge(x0, x1)){
		x0.copy(x1);
		Nvector temp(x0);
		r.vectorTransform(temp);
		temp.scalarMultiply(-1.0);
		temp.addVector(b);
		d.vectorTransform(temp);
		x1.copy(temp);
		x1.print();
	}		
}


void MatrixMethods::gaussSeidelMethod(const NMmatrix& mat, const Nvector& b, Nvector& x1){
	Nvector x0(b.getSize());
	randomizeVector(x1);
	do{
		x0.copy(x1);
		for(int i = 0; i < b.getSize(); i++){
			double dot = 0;
			for(int j = 0; j < b.getSize(); j++)
				if(j != i)
					dot += mat.get(j, i)*x0.get(j);					
			x1.set(i, (b.get(i) - dot)/mat.get(i, i));
		}
	}while(!vectorsConverge(x0, x1));
}

void MatrixMethods::SORmethod(const NMmatrix& mat, const Nvector& b, Nvector& x1){
	Nvector x0(b.getSize());
	randomizeVector(x1);
	do{
		x0.copy(x1);
		for(int i = 0; i < b.getSize(); i++){
			double dot = 0;
			for(int j = 0; j < b.getSize(); j++)
				if(j != i)
					dot += mat.get(j, i)*x0.get(j);	
			x1.set(i, (1.0 - relaxation)*x0.get(i) + ((relaxation/mat.get(i, i))*(b.get(i) - dot)));	
		}			
	}while(!vectorsConverge(x0, x1));
}

void MatrixMethods::gaussianEliminationMethod(const NMmatrix& mat, const Nvector& b, Nvector& x1){
	x1.copy(b);
	NMmatrix cmat(mat);
	for(int diagonal = 0; diagonal < b.getSize(); diagonal++){
		for(int colIndex = diagonal+1; colIndex < b.getSize(); colIndex++){
			if(cmat.get(diagonal, colIndex) != 0.0){
				Nvector opRow(b.getSize());
				double scalar = cmat.get(diagonal, colIndex)/cmat.get(diagonal, diagonal);
				cmat.getRow(diagonal, opRow);
				opRow.scalarMultiply(scalar);
				for(int j = 0; j < opRow.getSize(); j++){
					cmat.set(j, colIndex, cmat.get(j, colIndex)-opRow.get(j));
				}
				cmat.set(diagonal, colIndex, 0.0);
				x1.set(colIndex, x1.get(colIndex)-(x1.get(diagonal)*scalar));
			}
			cmat.print();
		}
	}
	for(int diagonal = b.getSize()-1; diagonal >= 0; diagonal--){
		for(int rowIndex = b.getSize()-1; rowIndex >= diagonal; rowIndex--){
			if(rowIndex != diagonal){
				x1.set(diagonal, x1.get(diagonal) - (cmat.get(rowIndex, diagonal) * x1.get(rowIndex)));
			}
			else{
				x1.set(diagonal, x1.get(diagonal)/cmat.get(diagonal, diagonal));
				cmat.set(diagonal, diagonal, 1.0);
			}
		}
	}
}

void MatrixMethods::decompositionLUMethod(const NMmatrix& mat, const Nvector& b, Nvector& x1){
	NMmatrix lMat(mat), uMat(mat);
	initLUMatrices(lMat, uMat);

	printf("\nL Matrix:\n");
	lMat.print();
	printf("\nU Matrix:\n");
	uMat.print();
	printf("\n");	

	for(int diagonal = 0; diagonal < b.getSize(); diagonal++){
                if(diagonal != 0){
                        printf("L Column Operations\n");
                        for(int lRow = diagonal; lRow < b.getSize(); lRow++){
                             double element = mat.get(diagonal, lRow);
                             printf("%0.4f ", element);
                             for(int i = 0; i < b.getSize(); i++)
                                if(i != diagonal){
                                    element -= (uMat.get(diagonal, i)*lMat.get(i, lRow));
                                    printf("- (%0.4f * %0.4f) ", uMat.get(diagonal, i), lMat.get(i, lRow));
                                }
                             lMat.set(diagonal, lRow, element/uMat.get(diagonal, diagonal));
                             printf("/%0.4f = %0.4f\n", uMat.get(diagonal, diagonal), element/uMat.get(diagonal, diagonal));
                        }
                        lMat.print();
                        printf("\n");
                }

		if(diagonal != b.getSize()-1){
                        printf("U Row Operations\n");
                        for(int uColumn = diagonal+1; uColumn < b.getSize(); uColumn++){
                                double element = mat.get(uColumn, diagonal);
                                printf("%0.4f ", element);
                                for(int i = 0; i < b.getSize(); i++)
                                    if(i != diagonal){
                                        element -= (lMat.get(i, diagonal)*uMat.get(uColumn, i));
                                        printf("- (%0.4f * %0.4f) ", lMat.get(i, diagonal), uMat.get(uColumn, i));      
                                    }
                                uMat.set(uColumn, diagonal, element/lMat.get(diagonal, diagonal));
                                printf("/%0.4f = %0.4f\n", lMat.get(diagonal, diagonal), element/lMat.get(diagonal, diagonal));
                        }
                        uMat.print();
                        printf("\n");
                }
	}

	Nvector y(b);
	for(int diagonal = 0; diagonal < b.getSize(); diagonal++){
		for(int rowIndex = 0; rowIndex <= diagonal; rowIndex++){
			if(rowIndex != diagonal)
				y.set(diagonal, y.get(diagonal) - (lMat.get(rowIndex, diagonal)*y.get(rowIndex)));
			else
	 			y.set(diagonal, y.get(diagonal)/lMat.get(diagonal, diagonal));	
		}
	}

	x1.copy(y);
        for(int diagonal = b.getSize()-1; diagonal >= 0; diagonal--){
                for(int rowIndex = b.getSize()-1; rowIndex >= diagonal; rowIndex--){
                        if(rowIndex != diagonal)
                                x1.set(diagonal, x1.get(diagonal) - (uMat.get(rowIndex, diagonal) * x1.get(rowIndex)));
                        else
                                x1.set(diagonal, x1.get(diagonal)/uMat.get(diagonal, diagonal));
                }
        }
}

void MatrixMethods::randomizeVector(Nvector& vec){
	srand(time(0));
	for(int i = 0; i < vec.getSize(); i++)
		vec.set(i, (double)(rand()%20));
}

bool MatrixMethods::vectorsConverge(const Nvector& vec0, const Nvector& vec1){
	if(vec0.getSize() != vec1.getSize())
		return 0;
	for(int i = 0; i < vec0.getSize(); i++)
		if(fabs(vec0.get(i)-vec1.get(i)) > MARGIN_OF_ERROR)
			return 0;
	return 1;
}

void MatrixMethods::initLUMatrices(NMmatrix& lMat, NMmatrix& uMat){
	if(uMat.getN() != lMat.getN() || uMat.getM() != lMat.getM())
		return;

	for(int i = 0; i < uMat.getN(); i++)
		for(int j = 0; j < uMat.getM(); j++){
			if(j == 0 && i != 0)	
				lMat.set(i, j, 0.0);
				uMat.set(i, j, 0.0);
			if(i == 0)
				uMat.set(i, j, 0.0);
			if(j != 0 && i != 0){	
				uMat.set(i, j, 0.0);
				lMat.set(i, j, 0.0);
			}
			if(j == i)
				uMat.set(i, j, 1.0);
		}	
}
