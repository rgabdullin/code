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