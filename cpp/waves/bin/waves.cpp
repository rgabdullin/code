#include <iostream>

typedef struct SimulationParams{
    int Nx;
    int Ny;
    int Nz;
    int Nt;
    double Lx;
    double Ly;
    double Lz;
    double Lt;
    double hx;
    double hy;
    double hz;
    double ht;
} Params;

class U{
    double data;
    SimulationParams p;

};

int main(void){

}