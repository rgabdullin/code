#include <iostream>
#include <fstream>

#include <cmath>
#include <queue>
#include <vector>

#include <omp.h>

#include "GridFunc.h"

struct SimulationParams{
    GridParams grid;
    int Nt;
    int fps;
    double Lt;
    double ht;

    bool isPeriodicX, isPeriodicY, isPeriodicZ;

    SimulationParams(){}
    SimulationParams(double L, int N, int T, int fps, bool isPeriodicX = false, bool isPeriodicY = false, bool isPeriodicZ = false):
        grid(L,N)
    {
        Nt = T * fps;
        Lt = T;
        ht = Lt/Nt;

        this->fps = fps;

        this->isPeriodicX = isPeriodicX;
        this->isPeriodicY = isPeriodicY;
        this->isPeriodicZ = isPeriodicZ;
    }
} params;

inline double phi(int i, int j, int k){
    return sin(
        fmin(
            4*sqrt(
                pow(params.grid.hx * i - params.grid.Lx/2, 2)
                + pow(params.grid.hy * j - params.grid.Ly/2, 2)
                + pow(params.grid.hz * k - params.grid.Lz/2, 2)
            ),
            M_PI
        )
    );
}

int main(int argc, char* argv[]){
    int batch_size = 512;

    params = SimulationParams(M_PI,96,5,120,true,true,true);

    std::cout << "Hello, Wave Equation!" << std::endl;
    std::cout << "sizeof(GridParams) = " << sizeof(GridParams) << std::endl;
    std::cout << "sizeof(SimulationParams) = " << sizeof(SimulationParams) << std::endl;
    std::cout << "sizeof(bool) = " << sizeof(bool) << std::endl;


    std::ofstream fout("/home/ruslixag/Jupyter/ParallelComputing/cache/result.bin");
    fout.write((const char*)&params,sizeof(params));

    std::vector<GridFunc> u;
    for(int t = 0; t < params.Nt; ++t){
        if((t == 0) || ((t+1)%8 == 0))
            std::cout << "\rComputing " << t+1 << "/" << params.Nt << std::flush;

        u.push_back(GridFunc(params.grid));
        int n = u.size()-1;
        int Nx = params.grid.Nx;
        int Ny = params.grid.Ny;
        int Nz = params.grid.Nz;
        
        if(n == 0){
            int size = (Nx+1)*(Ny+1)*(Nz+1);
            #pragma omp parallel for
            for(int idx = 0; idx < size; ++idx){
                int k = idx % (Nz+1),
                    j = (idx / (Nz+1)) % (Ny+1),
                    i = (idx / ((Nz+1)*(Ny+1)));
                u[n](i,j,k) = phi(i,j,k);
            }
        }
        else if(n == 1){
            int size = (Nx+1)*(Ny+1)*(Nz+1);
            #pragma omp parallel for
            for(int idx = 0; idx < size; ++idx){
                int k = idx % (Nz+1),
                    j = (idx / (Nz+1)) % (Ny+1),
                    i = (idx / ((Nz+1)*(Ny+1)));
                if((i>0)&(i<Nx)&(j>0)&(j<Ny)&(k>0)&(k<Nz))
                    u[n](i,j,k) = u[n-1](i,j,k) + pow(params.ht,2)/2 * (
                        (phi(i+1,j,k) - 2*phi(i,j,k) + phi(i-1,j,k))/pow(params.grid.hx,2)+
                        (phi(i,j+1,k) - 2*phi(i,j,k) + phi(i,j-1,k))/pow(params.grid.hy,2)+
                        (phi(i,j,k+1) - 2*phi(i,j,k) + phi(i,j,k-1))/pow(params.grid.hz,2)
                    );
            }
        }
        else {
            int size = (Nx+1)*(Ny+1)*(Nz+1);
            #pragma omp parallel for
            for(int idx = 0; idx < size; ++idx){
                int k = idx % (Nz+1),
                    j = (idx / (Nz+1)) % (Ny+1),
                    i = (idx / ((Nz+1)*(Ny+1)));
                if((i>0)&(i<Nx)&(j>0)&(j<Ny)&(k>0)&(k<Nz))
                    u[n](i,j,k) = 2*u[n-1](i,j,k) - u[n-2](i,j,k) + pow(params.ht,2) * (
                        (u[n-1](i+1,j,k) - 2*u[n-1](i,j,k) + u[n-1](i-1,j,k))/pow(params.grid.hx,2)+
                        (u[n-1](i,j+1,k) - 2*u[n-1](i,j,k) + u[n-1](i,j-1,k))/pow(params.grid.hy,2)+
                        (u[n-1](i,j,k+1) - 2*u[n-1](i,j,k) + u[n-1](i,j,k-1))/pow(params.grid.hz,2)
                    );
            }
        }
        //сохранение
        
        if(u.size() == batch_size){
            for(int n = 0; n < batch_size-2; ++n)
                fout.write((const char*)u[n].GetData(),u[n].GetSize() * sizeof(double));
            fout.flush();
            
            u.erase(u.begin(),u.begin() + (batch_size-2));
        }    
    }
    for(int n = 0; n < u.size(); ++n)
        fout.write((const char*)u[n].GetData(),u[n].GetSize() * sizeof(double));

    return 0;
}