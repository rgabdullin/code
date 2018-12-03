for(int i = 1; i < params.grid.Nx; ++i){
    std::cout << "[" << std::endl;
    for(int j = 1; j < params.grid.Ny; ++j){
        std::cout << "\t";
        for(int k = 1; k < params.grid.Nz; ++k){
            double value = phi(i,j,k);
            u[n](i,j,k) = value;
            std::cout << u[n](i,j,k) << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "]" << std::endl;
}