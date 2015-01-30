#ifndef N_VECTOR
#define N_VECTOR
class Nvector{
private:
	int n;
	double* vec;

public:
	Nvector(int);
	Nvector(const Nvector&);
	~Nvector();

	void scalarMultiply(double);
	void addVector(const Nvector&);
	double dotProduct(const Nvector&);
	int getSize() const;
	void set(int, double);
	double get(int) const;
	void copy(const Nvector&);
	void print();
};
#endif
