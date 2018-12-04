#pragma once
#include <vector>

struct GridParams{
    double Lx, Ly, Lz;
    double hx, hy, hz;
    int Nx, Ny, Nz;

    GridParams();
    GridParams(double L, int N);
};

class GridFunc{
    GridParams grid;
    std::vector<double> data;
    int size;
public:
    double*         GetData();
    int             GetSize();
    GridParams&     GetGrid();

    double& operator()(int i, int j,int k);

    GridFunc();
    GridFunc(const GridParams& grid);
};