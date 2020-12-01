#include <iostream>
#include <fstream>
#include "bifiltration_cell.h"

// Contributed by Michael Kerber
template<typename FT>
void bifiltration_to_firep(std::vector<BifiltrationCell<FT>> bf, int dim, std::ofstream& ofile) {
    ofile << "firep" << std::endl;
    ofile << "Critical value (r)" << std::endl;
    ofile << "Order (k)" << std::endl;
    long no_dim_minus_1=0, no_dim=0, no_dim_plus_1=0;
    std::map<int,int> new_id_d_minus_1;
    std::map<int,int> new_id_d;
    for(const auto& bc : bf) {
        if(bc.d==dim-1) {
          new_id_d_minus_1[bc.id]=no_dim_minus_1;
          no_dim_minus_1++;
        }
        if(bc.d==dim) {
          new_id_d[bc.id]=no_dim;
          no_dim++;
        }
        if(bc.d==dim+1) {
          no_dim_plus_1++;
        }
    }
    ofile << no_dim_plus_1 << " " << no_dim << " " << no_dim_minus_1 << std::endl;
    std::cout << no_dim_plus_1 << " " << no_dim << " " << no_dim_minus_1 << std::endl;

    for(const auto& bc : bf) {
        if(bc.d==dim+1) {
            ofile << std::fixed << std::setprecision(12) << CGAL::to_double(bc.r) << " " << -bc.k << " ; ";
            std::vector <int> bd;
            for(int idx : bc.boundary) {
                bd.push_back(new_id_d[idx]);
            }
            std::sort(bd.begin(),bd.end());
            for(int idx : bd) {
                ofile << idx << " ";
            }
            ofile << std::endl;
        }
    }
    
    for(const auto& bc : bf) {
        if(bc.d==dim) {
            ofile << std::fixed << std::setprecision(12) << CGAL::to_double(bc.r) << " " << -bc.k << " ; ";
            std::vector <int> bd;
            for(int idx : bc.boundary) {
                bd.push_back(new_id_d_minus_1[idx]);
            }
            std::sort(bd.begin(),bd.end());
            for(int idx : bd) {
                ofile << idx << " ";
            }
            ofile << std::endl;
        }
    }
}
