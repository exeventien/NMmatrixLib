#ifndef MATRIX_METHODS
#define MATRIX_METHODS

class MatrixMethods{
	private:
	static const double MARGIN_OF_ERROR = 0.00001;
	
	static const double relaxation = 1.4;

	static void randomizeVector(Nvector&);
    static bool vectorsConverge(const Nvector&, const Nvector&);
	static void initLUMatrices(NMmatrix&, NMmatrix&);

	public:
	static void jacobiMethod(const NMmatrix&, const Nvector&, Nvector&);
	
	static void gaussSeidelMethod(const NMmatrix&, const Nvector&, Nvector&);

	static void SORmethod(const NMmatrix&, const Nvector&, Nvector&);

	static void gaussianEliminationMethod(const NMmatrix&, const Nvector&, Nvector&);
	
	static void decompositionLUMethod(const NMmatrix&, const Nvector&, Nvector&);
};

#endif
