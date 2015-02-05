#ifndef NM_MATRIX
#define NM_MATRIX
class NMmatrix{
private:
	int n, m;
	double* mat;

public:
	NMmatrix(int, int);
	NMmatrix(const NMmatrix&);
	~NMmatrix();

	void vectorTransform(Nvector&);
	void diagonalInverse();
	void diagonalRemainder();

	void diagonalBelow();
	void aboveDiagonal();
	void getRow(int, Nvector&);
	double get(int, int) const;

	void set(int, int, double);
	void copy(const NMmatrix&);
	void print() const;
};
#endif
