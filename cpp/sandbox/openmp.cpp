#include <omp.h>
#include <iostream>

int main(void){
    std::cout << "Num threads: " << omp_get_num_procs() << std::endl;
    #pragma omp parallel
    {
        std::cout << "Hello, OpenMP!\n";
    }
    return 0;
}