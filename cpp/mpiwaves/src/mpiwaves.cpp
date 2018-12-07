#include <mpi.h>
#include <omp.h>

#include <iomanip>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>

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
    void Init(int Nx, int Ny, int Nz, int Gx, int Gy, int Gz, int Bx, int By, int Bz){
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
            std::cerr << "Bad grid size" << std::endl;

        Dx = Nx / Gx;
        Dy = Ny / Gy;
        Dz = Nz / Gz;

        Ox = Dx * Bx;
        Oy = Dy * By;
        Oz = Dz * Bz;

        data.resize(Dx * Dy * Dz);

        extern_data[0].resize(Dy*Dz);
        extern_data[1].resize(Dy*Dz);

        extern_data[2].resize(Dx*Dz);
        extern_data[3].resize(Dx*Dz);
        
        extern_data[4].resize(Dx*Dy);
        extern_data[5].resize(Dx*Dy);
    }

    MPIGridFunc() {}
    MPIGridFunc(int Nx, int Ny, int Nz, int Gx, int Gy, int Gz, int Bx, int By, int Bz){
        this->Init(Nx,Ny,Nz, Gx,Gy,Gz, Bx,By,Bz);
    }
    MPIGridFunc(int N, int G, int Bx, int By, int Bz){
        this->Init(N,N,N, G,G,G, Bx,By,Bz);
    }

    void SyncMPI(MPI_Comm comm){
        this->PrepareExternData();

        int crd[3];
        int my_rank;
        MPI_Status status;

        MPI_Comm_rank(comm, &my_rank);

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

    double GetLocalIndex(int i, int j, int k){
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

    bool SetLocalIndex(int i, int j, int k, double v){
        if((i < 0)||(i >= Dx)||(j < 0)||(j >= Dy)||(k < 0)||(k >= Dz))
            return false;

        data[k + Dz*(j + Dy*i)] = v;
        return true;
    }


    int GetN(int i) {return (i == 0? Nx : (i == 1? Ny : Nz));}
    int GetB(int i) {return (i == 0? Bx : (i == 1? By : Bz));}
    int GetG(int i) {return (i == 0? Gx : (i == 1? Gy : Gz));}
    int GetD(int i) {return (i == 0? Dx : (i == 1? Dy : Dz));}
    int GetO(int i) {return (i == 0? Ox : (i == 1? Oy : Oz));}

    double Get(int i, int j, int k) {return GetLocalIndex(i - Ox, j - Oy, k - Oz);}
    bool Set(int i, int j, int k, double v) {return SetLocalIndex(i - Ox, j - Oy, k - Oz, v);}    
};

void PrintMPIGridFunc(MPI_Comm comm, MPIGridFunc& u, int print_rank = -1, bool extended = false){
    bool full_print = (print_rank == -1);
    
    print_rank = (full_print? 0 : print_rank);
    int ext = int(extended);

    int rank;
    MPI_Comm_rank(comm, &rank);

    for(int i = 0 - ext; i < u.GetN(0) + ext; i++){
        if(print_rank == rank)
            std::cout << "[" << std::endl;
        for(int j = 0 - ext; j < u.GetN(1) + ext; j++){
            if(print_rank == rank)
                std::cout << "\t";
            for(int k = 0 - ext; k < u.GetN(2) + ext; k++){
                double loc_value = u.Get(i,j,k);
                loc_value = (isnan(loc_value)? -1e300: loc_value);

                double glob_value = loc_value;
                if(full_print)
                    MPI_Reduce(&loc_value,&glob_value,1,MPI_DOUBLE,MPI_MAX,0,comm);

                if(print_rank == rank){
                    if(glob_value == -1e300)
                        std::cout << "?" << "\t";
                    else
                        std::cout << std::setprecision(3) << glob_value << "\t";
                }
            }
            if(print_rank == rank)
                std::cout << std::endl;
        }
        if(print_rank == rank)
            std::cout << "]" << std::endl;
    }
}

double UAnalytics(double x, double y, double z, double t){
    return sin(x) * sin (y) * cos(z) * cos(sqrt(3) * t);
}

double Phi(double x, double y, double z){
    return sin(x) * sin (y) * cos(z);
}

int main(int argc, char* argv[]){
    int rank, size;
    MPI_Comm comm;

    int dim[3], period[3], reorder;
    int coord[3];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int Gsize = int(pow(size,1.0/3));
    for(;Gsize*Gsize*Gsize < size; ++Gsize);

    dim[0]=Gsize; dim[1]=Gsize; dim[2]=Gsize;
    period[0]=1; period[1]=1; period[2]=1;

    int req_size = dim[0] * dim[1] * dim[2];
    if(req_size != size){
        if(rank == 0)
            std::cout << "Please run with cubic thread number (got " << req_size << " instead of " << size << ")" << std::endl;

        MPI_Abort(MPI_COMM_WORLD, 1);
    }


    MPI_Cart_create(MPI_COMM_WORLD, 3, dim, period, 1, &comm);
    MPI_Cart_coords(comm, rank, 3, coord);

    if(argc <= 1){
        if(rank == 0)
            std::cout << "Need at least 1 argument" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int N;
    if(rank == 0){
        N = atoi(argv[1]);
        std::cout << "N = " << N << std::endl;

        if(N % Gsize){
            std::cout << N << " %% " << Gsize << " != 0" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    MPI_Bcast(&N,1,MPI_INT,0,comm);

    int Nx = N, Ny = N, Nz = N;
    double Lx = 2 * M_PI, Ly = 2 * M_PI, Lz = 2 * M_PI;
    double hx = Lx / (Nx-1), hy = Ly / (Ny-1), hz = Lz / (Nz-1);

    bool is_periodic_x = false, is_periodic_y = false, is_periodic_z = true;

    int T = 1, fps = 120;
    // int frames = T * fps;

    int frames = 120;
    double ht = 1.0 / fps;

    std::vector<MPIGridFunc> u;

    bool compute_metrics = true;

    if(rank == 0){
        std::cout << "Running" << std::endl;
        std::cout << "Num procs: " << omp_get_num_procs() << std::endl;
        std::cout << "Num threads: " << omp_get_num_threads() << std::endl;
    }
    for(int n = 0; n < frames; n++){
        if(rank == 0){
            std::cout << "Frame " << n << std::endl;
        }

        u.push_back(MPIGridFunc(Nx,Ny,Nz,dim[0],dim[1],dim[2],coord[0],coord[1],coord[2]));

        int Dx = u[n].GetD(0);
        int Dy = u[n].GetD(1);
        int Dz = u[n].GetD(2);

        int Ox = u[n].GetO(0);
        int Oy = u[n].GetO(1);
        int Oz = u[n].GetO(2);

        double linf_metrics = 0;

        #pragma omp parallel for
        for(int idx = 0; idx < Dx * Dy * Dz; ++idx){
            int k = idx % Dz,
                j = (idx / Dz) % Dy,
                i = idx / (Dz * Dy);
            double x = (Ox + i) * hx, 
                y = (Oy + j) * hy,
                z = (Oz + k) * hz;
            double value;
            if(n == 0)
                value = Phi(x,y,z);
            if(n == 1)
                value = u[n-1].GetLocalIndex(i,j,k) + pow(ht,2)/2 * (
                    (Phi(x+hx,y,z) - 2*Phi(x,y,z) + Phi(x-hx,y,z))/pow(hx,2)+
                    (Phi(x,y+hy,z) - 2*Phi(x,y,z) + Phi(x,y-hy,z))/pow(hy,2)+
                    (Phi(x,y,z+hz) - 2*Phi(x,y,z) + Phi(x,y,z-hz))/pow(hz,2)
                );
            if(n >= 2)
                value = 2*u[n-1].GetLocalIndex(i,j,k) - u[n-2].GetLocalIndex(i,j,k) + pow(ht,2) * (
                    (u[n-1].GetLocalIndex(i+1,j,k) - 2*u[n-1].GetLocalIndex(i,j,k) + u[n-1].GetLocalIndex(i-1,j,k))/pow(hx,2)+
                    (u[n-1].GetLocalIndex(i,j+1,k) - 2*u[n-1].GetLocalIndex(i,j,k) + u[n-1].GetLocalIndex(i,j-1,k))/pow(hy,2)+
                    (u[n-1].GetLocalIndex(i,j,k+1) - 2*u[n-1].GetLocalIndex(i,j,k) + u[n-1].GetLocalIndex(i,j,k-1))/pow(hz,2)
                );

            if(!is_periodic_x && ((x == 0)||(x == Lx)))
                value = 0;
            if(!is_periodic_y && ((y == 0)||(y == Ly)))
                value = 0;
            if(!is_periodic_z && ((z == 0)||(z == Lz)))
                value = 0;
            u[n].SetLocalIndex(i,j,k,value);

            if(compute_metrics)
                linf_metrics = fmax(linf_metrics, fabs(u[n].GetLocalIndex(i,j,k) - UAnalytics(x,y,z,n*ht)));
        }

        u[n].SyncMPI(comm);
        
        if(compute_metrics){
            double result = 0;
            MPI_Reduce(&linf_metrics,&result,1,MPI_DOUBLE,MPI_MAX,0,comm);
            if(rank == 0)
                std::cout << "L_inf = " << result << std::endl;
        }      
        // PrintMPIGridFunc(comm,u[n],-1,false);
    }

    MPI_Finalize();
    return 0;
}