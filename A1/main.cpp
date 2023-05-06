#include "library.hpp"
#include "multiply.hpp"
  
int main(int argc,char* argv[]){
    assert(argc==3);
    matrix_square(string(argv[1]),string(argv[2]));
    return 0;
}