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
	randomizeVector(x0);
	NMmatrix uMat(mat), lMat(mat);
	uMat.aboveDiagonal();
	lMat.diagonalBelow();
	int row = 0;
	while(!vectorsConverge(x1, x0)){
		Nvector currentRow(b.getSize()), nextRow(b.getSize());
		uMat.getRow(row, currentRow);
		lMat.getRow((row==b.getSize()-1) ? 0 : row+1, nextRow);
		double k0 = x1.dotProduct(currentRow);
		double k1 = x0.dotProduct(nextRow); 
		x1.set(row, (b.get(row) - k1 - k0)/mat.get(row, row));
		x1.print();
		if(row == b.getSize())
			row = 0;
		else
			row++;
	}

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
	cmat.print();
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

void MatrixMethods::randomizeVector(Nvector& vec){
	srand(time(0));
	for(int i = 0; i < vec.getSize(); i++)
		vec.set(i, (double)(rand()%50));
}

bool MatrixMethods::vectorsConverge(const Nvector& vec0, const Nvector& vec1){
	if(vec0.getSize() != vec1.getSize())
		return 0;
	for(int i = 0; i < vec0.getSize(); i++)
		if(fabs(vec0.get(i)-vec1.get(i)) > MARGIN_OF_ERROR)
			return 0;
	return 1;
}
