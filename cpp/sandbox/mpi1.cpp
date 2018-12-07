#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>

class MPIGridFunc{
    int Nx, Ny, Nz;
    int Dx, Dy, Dz, Ox, Oy, Oz;
    int Bx, By, Bz, Gx, Gy, Gz;
    
    std::vector<double> data;
    std::vector<double> extern_data[6];

    void PrepareExternData(){
        for(int j = 0; j < Dy; ++j)
            for(int k = 0; k < Dz; ++k){
                int i = 0, target = (Gx > 1? 0: 1);
                extern_data[target][k + Dz*j] = data[k + Dz*(j + Dy*i)];
                i = Dx-1; target = (Gx > 1? 1: 0);
                extern_data[target][k + Dz*j] = data[k + Dz*(j + Dy*i)];
            }
        for(int i = 0; i < Dx; ++i)
            for(int k = 0; k < Dz; ++k){
                int j = 0, target = (Gy > 1? 2: 3);
                extern_data[target][k + Dz*i] = data[k + Dz*(j + Dy*i)];
                j = Dy-1; target = (Gy > 1? 3: 2);
                extern_data[target][k + Dz*i] = data[k + Dz*(j + Dy*i)];
            }
        for(int i = 0; i < Dx; ++i)
            for(int j = 0; j < Dy; ++j){
                int k = 0, target = (Gz > 1? 4: 5);
                extern_data[target][j + Dy*i] = data[k + Dz*(j + Dy*i)];
                k = Dz-1; target = (Gz > 1? 5: 4);
                extern_data[target][j + Dy*i] = data[k + Dz*(j + Dy*i)];
            }
    }

public:
    MPIGridFunc(int Nx, int Ny, int Nz, int Gx, int Gy, int Gz, int Bx, int By, int Bz){
        this->Nx = Nx;
        this->Ny = Ny;
        this->Nz = Nz;

        this->Gx = Gx;
        this->Gy = Gy;
        this->Gz = Gz;

        this->Bx = Bx;
        this->By = By;
        this->Bz = Bz;

        if((Nx % Gx)||(Ny % Gy)||(Nz % Gz))
            std::cerr << "Плохой размер сетки" << std::endl;

        Dx = Nx / Gx;
        Dy = Ny / Gy;
        Dz = Nz / Gz;

        Ox = Dx * Bx;
        Oy = Dy * By;
        Oz = Dz * Bz;

        data.resize(Dx * Dy * Dz);

        std::fill(data.begin(), data.end(), nan(""));

        extern_data[0].resize(Dy*Dz);
        extern_data[1].resize(Dy*Dz);

        extern_data[2].resize(Dx*Dz);
        extern_data[3].resize(Dx*Dz);
        
        extern_data[4].resize(Dx*Dy);
        extern_data[5].resize(Dx*Dy);

        /*std::cout << "[Ox,Oy,Oz;Dx,Dy,Dz] = " << "[" 
            << Ox << "," << Oy << "," << Oz << ";"
            << Dx << "," << Dy << "," << Dz << "]" << std::endl;/**/
    }

    void SyncMPI(MPI_Comm comm){
        this->PrepareExternData();

        int crd[3];

        int my_rank;

        MPI_Comm_rank(comm, &my_rank);

        crd[0] = Bx;crd[1] = By;crd[2] = Bz;

        int target[6];
        int delta[6][3] = {
            {-1,0,0},{1,0,0},
            {0,-1,0},{0,1,0},
            {0,0,-1},{0,0,1}
        };

        for(int i = 0; i < 6; i++){
            crd[0] = Bx + delta[i][0];
            crd[1] = By + delta[i][1];
            crd[2] = Bz + delta[i][2];
            
            MPI_Cart_rank(comm,crd,&target[i]);            
        }

        MPI_Status status;

        for(int axis = 0; axis < 3; axis++){
            int tp = (axis == 0? Bx : (axis == 1? By : Bz)) % 2;
            for(int tmp = 0; tmp < 2; tmp++){
                tp = 1 - tp;

                int target_idx = 2 * axis + (1 - tp);

                int send_tag = my_rank * 100 + axis * 10 + tp;
                int recv_tag = target[target_idx] * 100 + axis * 10 + (1-tp);

                if(my_rank != target[target_idx]){
                    MPI_Sendrecv_replace(&extern_data[target_idx][0],extern_data[target_idx].size(),
                        MPI_DOUBLE,target[target_idx],send_tag,target[target_idx],recv_tag,
                        comm,&status);
                }
            }
        }
    }

    double Get(int i, int j, int k){
        i -= Ox; j -= Oy; k -= Oz;
        if((j >= 0)&&(j<Dy)&&(k>=0)&&(k<Dz)){
            if(i == -1)
                return extern_data[0][k + Dz*j];
            if((i >= 0)&&(i < Dx))
                return data[k + Dz*(j + Dy*i)];
            if(i == Dx)
                return extern_data[1][k + Dz*j];
        }
        if((i >= 0)&&(i<Dx)&&(k>=0)&&(k<Dz)){
            if(j == -1)
                return extern_data[2][k + Dz*i];
            if(j == Dy)
                return extern_data[3][k + Dz*i];
        }
        if((i >= 0)&&(i<Dx)&&(j>=0)&&(j<Dy)){
            if(k == -1)
                return extern_data[4][j + Dy*i];
            if(k == Dz)
                return extern_data[5][j + Dy*i];
        }
        return nan("");
    }
    bool Set(int i, int j, int k, double v){
        i -= Ox; j -= Oy; k -= Oz;
        if((i < 0)||(i >= Dx)||(j < 0)||(j >= Dy)||(k < 0)||(k >= Dz))
            return false;

        data[k + Dz*(j + Dy*i)] = v;
        return true;
    }
};

int main(int argc, char* argv[]){
    int rank, size;
    MPI_Comm comm;

    int dim[3], period[3], reorder;
    int coord[3];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // int Gsize = int(pow(size,1.0/3));
    // for(;Gsize*Gsize*Gsize < size; ++Gsize);
    // dim[0]=Gsize; dim[1]=Gsize; dim[2]=Gsize;
    
    dim[0]=1; dim[1]=2; dim[2]=3;

    int req_size = dim[0] * dim[1] * dim[2];

    if(req_size != size){
        if(rank == 0)
            std::cout << "Please run with cubic thread number (" << size << "," << req_size << ")" << std::endl;

        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    period[0]=1; period[1]=1; period[2]=1;
    reorder=1;

    MPI_Cart_create(MPI_COMM_WORLD, 3, dim, period, reorder, &comm);

    MPI_Cart_coords(comm, rank, 3, coord);

    int Nx, Ny, Nz;
    if(rank == 0){
        Nx = 2 * dim[0];
        Ny = 3 * dim[1];
        Nz = 3 * dim[2];
    }

    MPI_Bcast(&Nx,1,MPI_INT,0,comm);
    MPI_Bcast(&Ny,1,MPI_INT,0,comm);
    MPI_Bcast(&Nz,1,MPI_INT,0,comm);

    int Bx = coord[0], By = coord[1], Bz = coord[2];
    int Gx = dim[0], Gy = dim[1], Gz = dim[2];

    MPIGridFunc u(Nx,Ny,Nz,Gx,Gy,Gz,Bx,By,Bz);

    int count = 0;
    for(int i = -1; i <= Nx; i++)
        for(int j = -1; j <= Ny; j++)
            for(int k = -1; k <= Nz; k++)
                count += u.Set(i,j,k, k + Nz*(j + Ny * i));
    u.SyncMPI(comm);

    for(int i = -1; i <= Nx; i++){
        if(rank == 0)
            std::cout << "[" << std::endl;
        for(int j = -1; j <= Ny; j++){
            std::cout << "\t";
            for(int k = -1; k <= Nz; k++){
                double loc_value = u.Get(i,j,k);
                loc_value = (isnan(loc_value)? -1e300: loc_value);

                double glob_value = loc_value;
                // MPI_Reduce(&loc_value,&glob_value,1,MPI_DOUBLE,MPI_MAX,0,comm);

                if(rank == 0){
                    if(glob_value == -1e300)
                        std::cout << "?" << "\t";
                    else
                        std::cout << glob_value << "\t";
                }
            }
            if(rank == 0)
                std::cout << std::endl;
        }
        if(rank == 0)
            std::cout << "]" << std::endl;
    }

    if(rank == 0){
        std::cout << "Result:" << std::endl;
    }

    for(int i = 0; i < Nx; i++){
        if(rank == 0)
            std::cout << "[" << std::endl;
        for(int j = 0; j < Ny; j++){
            std::cout << "\t";
            for(int k = 0; k < Nz; k++){
                double loc_value = u.Get(i,j,k);
                loc_value = (isnan(loc_value)? -1e300: loc_value);

                double glob_value = loc_value;
                MPI_Reduce(&loc_value,&glob_value,1,MPI_DOUBLE,MPI_MAX,0,comm);

                if(rank == 0){
                    if(glob_value == -1e300)
                        std::cout << "?" << "\t";
                    else
                        std::cout << glob_value << "\t";
                }
            }
            if(rank == 0)
                std::cout << std::endl;
        }
        if(rank == 0)
            std::cout << "]" << std::endl;
    }

    MPI_Finalize();
    return 0;
}