#include <cstdlib>
#include <Generator.h>

int main(int argc, char** argv){
	if(argc >= 4){
		Generator generator(atoi(argv[1]), atoi(argv[2]));
		generator.createMatrix();
		generator.printFile(argv[3]);
	}
	return 0;
}
