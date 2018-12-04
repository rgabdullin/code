#include "GridFunc.h"

// GridParams
GridParams::GridParams() {}
GridParams::GridParams(double L, int N){
    Lx = Ly = Lz = L;
    Nx = Ny = Nz = N-1;

    hx = Lx/(Nx+1);
    hy = Ly/(Ny+1);
    hz = Lz/(Nz+1);
}

// GridFunc
double*     GridFunc::GetData() { return &data[0]; }
int         GridFunc::GetSize() { return size; }
GridParams& GridFunc::GetGrid() { return grid; }

GridFunc::GridFunc(){}
GridFunc::GridFunc(const GridParams& grid){
    this->grid = grid;

    size = (grid.Nx+1) * (grid.Ny+1) * (grid.Nz+1);
    
    if(size > 0){
        data.resize(size);
        std::fill(data.begin(),data.end(),'\0');
    }
}

double& GridFunc::operator()(int i, int j,int k) { return data[k + (grid.Nz+1)*(j + (grid.Ny+1) * i)]; }


