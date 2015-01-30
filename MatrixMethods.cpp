#include <cstdlib>
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
	NMmatrix l(mat), u(mat);
	l.inverseDiagonalBelow();
	l.print();
	u.aboveDiagonal();
	u.print();
	Nvector x0(b.getSize());
	randomizeVector(x1);
	int row = 0;
	while(!vectorsConverge(x1, x0)){
		Nvector uRow(b.getSize()), lRow(b.getSize());
		u.getRow(row, uRow);
		l.getRow(row, lRow);
		double ux0 = x0.dotProduct(uRow);
		x0.set(row, (b.get(row) - ux0));
		double lbux0 = x0.dotProduct(lRow);
		x1.set(row, lbux0);
		x0.set(row, lbux0);
		if(row == b.getSize())
			row = 0;
		else
			row++;
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
