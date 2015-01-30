#ifndef MATRIX_METHODS
#define MATRIX_METHODS

class MatrixMethods{
	private:
	static const double MARGIN_OF_ERROR = 0.00001;
	static void randomizeVector(Nvector&);
        static bool vectorsConverge(const Nvector&, const Nvector&);


	public:
	static void jacobiMethod(const NMmatrix&, const Nvector&, Nvector&);
	
	static void gaussSeidelMethod(const NMmatrix&, const Nvector&, Nvector&);
};

#endif
