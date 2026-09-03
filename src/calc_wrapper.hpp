#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <iterator>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>   // std::pair
#include <array>
#include <cstdio>
#include <filesystem>

#include <gsl/gsl_errno.h>

#include <chrono>   //for execution time.

#include "calculation_core.hpp"
#include "helper.hpp"
//#include "morse_fitting.hpp"

inline int calc_wrapper(std::vector<double>& V, std::vector<double>& r, std::vector<double>& V_v, std::vector<double>& r_v, std::vector<double> &v, AllParams& p, CalcParam &cp) {

    
   //create the output folder for all outputs.
    try {
      //  std::filesystem::create_directories("./output");
    } catch (const std::exception& e) {
        std::cerr << "Failed to create output folder: " << e.what() << "\n";
    }

    //in-process variables for calculation.
    SplineNaK::Spline s;

    RKRContext rkr = CreateRKRContext(p, cp);     //rkr procedure.
    //std::cout<< p << std::endl;


    std::vector<double> r_raw;
    std::vector<double> V_raw;
    std::vector<double> r_spl;
    std::vector<double> V_spl;

    std::vector<double> r_l;
    std::vector<double> V_l;

    std::vector<double> r_i;
    std::vector<double> V_i;

    std::vector<double> r_o;
    std::vector<double> V_o;

    //std::vector<double> r;
    //std::vector<double> V;

    std::vector<double> d1V_l;
    std::vector<double> d2V_l;

    std::cout<< "Kaiser correction for vmin: " << rkr.ve.kaiser_correction_vmin() << std::endl;
    
    calc_rawRKR(V_raw, r_raw,nullptr,nullptr, nullptr, p, cp, rkr);

    setEquilibriumPoint(V_raw, r_raw, cp, p);

    if(!cp.detectInwardCurv){
        cp.detectOutwardCurv= detectOutwardCurvature(V_raw, r_raw,  cp);
        //std::cout<< cp.detectOutwardCurv<<std::endl;
    }
    s=SplineFit_to_RKR(V_raw, r_raw,V_spl, r_spl, cp);
    estimateInnerWall_params(s, r_raw, V_raw, cp, rkr);
    determineExtraPIndices(V_raw, r_raw, cp);
    construct_ExtraPCurve(rkr, cp, V_i, r_i, V_o, r_o);
   // std::cout <<"v_ex :" <<cp.v_ex <<" v_ex_level: "<< cp.v_ex_level << '\n';

    combine_inner_outer_rkr(V_i, r_i, V_o, r_o, V_raw, r_raw, V, r, cp);


    //need to add the discrete energy levels to the output later.
    std::cout<<"v_min = "<<cp.vmin<<", v_max = "<<cp.vmax<<", v_step = "<<cp.v_step<<std::endl;
    std::cout<<"v_ex_level = "<<cp.v_ex_level<<std::endl;
    std::cout<<cp<<std::endl;
    return 0;

}

inline int calc_morseV(std::vector<double>& V, std::vector<double>& r, AllParams& p) {
    // Calculate Morse potential values based on the provided parameters
    double De = p.ex.De; // Dissociation energy
    double ke = p.ex.ke; // Force constant
    p.morse.De = De;
    p.morse.ke = ke;
    double a = p.morse.a(); // Morse parameter
    double re = p.morse.Re; // Equilibrium bond length

    std::cout << "Calculating Morse potential with parameters: De = " << De << ", a = " << a << ", re = " << re << std::endl;
    if (r.empty()) {
        std::cout << "Error: r vector is empty." << std::endl;
        return 1;
    }

    for (size_t i = 0; i < r.size(); ++i) {
        double r_i = r[i];
        std::cout << "Calculating Morse potential for r[" << i << "] = " << r_i << std::endl;
        V.push_back(De * (1 - std::exp(-a * (r_i - re))) * (1 - std::exp(-a * (r_i - re))));
    }
    std::cout << "Morse potential calculation completed." << std::endl;

    return 0;
}

