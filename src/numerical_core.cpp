#include "numerical_core.hpp"
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <iterator>

#include <vector>
#include "rkr_procedure.h"
#include <iostream>

#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

#include <gsl/gsl_spline.h> // for cubic spline fitting.
#include <gsl/gsl_errno.h>
#include <utility>   // std::pair
#include <limits>
#include <array>
#include <cstdio>


//Adding the not-a-knot cspline c++ impl using eigen library.

#include "SplineNaK2.h"


class RGWriter {
private:
    std::ofstream fout;
    bool twoColumn;    // true = r & E, false = single-column

public:
    // --- Constructor for 2-column mode ---
    RGWriter(const std::string& filename)
        : twoColumn(true)
    {
        fout.open(filename, std::ios::out | std::ios::trunc);
        if (!fout.is_open()) {
            throw std::runtime_error("RGWriter: cannot open file '" + filename + "'");
        }

        fout << std::scientific << std::setprecision(8);
    }

    // --- Constructor for single-column output ---
    RGWriter(const std::string& filename, bool singleColumnMode)
        : twoColumn(!singleColumnMode)
    {
        fout.open(filename, std::ios::out | std::ios::trunc);
        if (!fout.is_open()) {
            throw std::runtime_error("RGWriter: cannot open file '" + filename + "'");
        }

        fout << std::scientific << std::setprecision(8);
    }

    // --- Write pair OR single value ---
    void writePoint(double r, double E) {
        if (!fout.is_open()) return;

        if (twoColumn) {
            fout << std::setw(15) << r
                 << "  "
                 << std::setw(15) << E
                 << '\n';
        } else {
            // Single column → only write "r"
            fout << std::setw(15) << r << '\n';
        }
    }

    // Overload for single-column writer:
    void writeValue(double v) {
        if (!fout.is_open()) return;
        fout << std::setw(15) << v << '\n';
    }

    // Destructor closes file
    ~RGWriter() {
        if (fout.is_open())
            fout.close();
    }
};


void plotXY(const std::vector<double>& x, const std::vector<double>& y, const std::string& title = "Plot") {
    if (x.size() != y.size() || x.empty()) {
        std::cerr << "Error: x and y must have same nonzero length\n";
        return;
    }

    FILE* gp = _popen("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\" -persist", "w");
    if (!gp) {
        std::cerr << "Error: Could not open gnuplot\n";
        return;
    }

    fprintf(gp, "set terminal wxt size 1000,700\n");
    fprintf(gp, "set title 'RKR Potential Well'\n");
    fprintf(gp, "set xlabel 'Internuclear distance r (Å)'\n");
    fprintf(gp, "set ylabel 'Energy (cm^{-1})'\n");
    fprintf(gp, "set grid\n");
   
    fprintf(gp, "plot '-'  with points pt 7 ps 0.5 lc rgb 'blue' title 'Potential Curve'\n");

    for (size_t i = 0; i < x.size(); ++i) {
        if (std::isfinite(x[i]) && std::isfinite(y[i])) {
            fprintf(gp, "%f %f\n", x[i], y[i]);
        }
    }
    fprintf(gp, "e\n");
    fflush(gp);

    _pclose(gp);
}
 

void calc(std::vector<double> &E, std::vector<double>& r,
          std::vector<double> &E_v, std::vector<double>& r_v,
          const AllParams& p)
{
    //don't abort the program if error is received.
    gsl_set_error_handler_off();


    // Constants (your originals)
    const double Cu = 16.857629206;      // hbar^2/2 in amu Å^2 cm^-1
    const double me = 0.000548579909;    // electron mass (amu)
    const double netcharge = 0.0;        

    const double space = p.ex.space;

    const double m1 = p.ex.m1;
    const double m2 = p.ex.m2;
    const double mu = (m1 * m2) / (m1 + m2 - netcharge * me);

    rkr_procedure ve(
        p.gv.we, p.gv.xwe, p.gv.ywe, p.gv.zwe, p.gv.awe,
        p.gv.bwe, p.gv.cwe, p.gv.dwe, p.gv.ewe, p.gv.fwe,
        p.bv.Be, p.bv.ae, p.bv.ye, p.bv._1e, p.bv._2e,
        p.ex.Te, p.ex.De, p.ex.ke, p.ex.Re,
        p.ex.Vmax, p.ex.space, p.ex.ladderspace, p.ex.UseKaiser
    );

    const double vmax_trunc = ve.checkAndTruncateVmax();
    const double vmax = std::min(p.ex.Vmax, vmax_trunc);
    const double vmin = ve.kaiser_correction_vmin();

    const int vmax_int = std::min(p.ex.Vmax + 1e-12, vmax_trunc);

    r_v.reserve(vmax_int + 1);
    E_v.reserve(vmax_int + 1);

    // Always include v=0 (zero-point)
    {
        int v0 = 0;
        double G0 = ve.G(v0);         // ve.G includes Y00_0 already in your class
        double E0 = p.ex.Te + G0;
        r_v.push_back(v0);
        E_v.push_back(E0);
    }

    // Integer ladder levels
    for (int v = 1; v <= vmax_int; ++v) {
        double Gv = ve.G(v);
        if (!std::isfinite(Gv)) break;

        r_v.push_back(v);
        ////printf("Te %f\n", p.ex.Te);
        E_v.push_back(p.ex.Te + Gv);
    }

    // --- Local helpers (unchanged) ---
    auto r1_from_fg = [](double fF, double gF) {
        const double s = std::sqrt(fF/gF + fF*fF);
        return s - fF;
    };
    auto r2_from_fg = [](double fF, double gF) {
        const double s = std::sqrt(fF/gF + fF*fF);
        return s + fF;
    };

    //external functions
    auto calculate_r1_n = [&](double vib){

        // Fleming integrals
        const double fF_2 = ve.fFleming(vib, Cu, mu, vmin, 1e-12);
        const double gF_2 = ve.gFleming(vib, Cu, mu, vmin, 1e-12);

        if (!std::isfinite(fF_2) || !std::isfinite(gF_2) || gF_2 == 0.0){//printf("The integrals are not finite!!! Try another vibrational level.\n");
        }

        return r1_from_fg(fF_2, gF_2);

    };

    auto calculate_G_n = [&](double vib){
        return p.ex.Te + ve.G(vib);
    };


    auto fitAB = [](double n, double d1, double d2, double E1, double E2,
                    double& A, double& B, double v_detect)
    {
        double denom = (1.0/std::pow(d1,n) - 1.0/std::pow(d2,n));
       // if (std::abs(denom) < 1e-24) { //printf("The faulty return\n");A = 0.0; B = E1; return; }
        A = (E1 - E2) / denom;
        B = E2 - A / std::pow(d2,n);
        //printf("The estimated A and B are %f and %f", A, B);
        // Example MATLAB vs C++:
        // 669782.9284 and -388.461  (MATLAB)
        // 675222      and -446.016  (old C++)
    };

    
    // --- PASS 1: build raw series and detect inward curvature ---
    std::vector<double> v_s, G_s, r1_s, r2_s, f_s;


    //std vector for energy levels and internuclear distances for discrete vibrational levels.
    std::vector<double> r1_n, r2_n; 


    //v_s.reserve(1000000); G_s.reserve(1000000);
    //r1_s.reserve(1000000); r2_s.reserve(1000000); f_s.reserve(1000000);

    bool inwardDetected = false;
    double v_detect = vmin;
    //printf("\nThe v-minimum and v-maximum are: %f and %f\n", vmin, vmax);
    double r_detect = 0.0;
    double last_r1 = std::numeric_limits<double>::infinity();

    for (double v = vmin+1e-6; v <= vmax; v += space) {
        const double Gv = ve.G(v) + p.ex.Te;

        if (ve.dG(v) <= 0.0){//printf("This is where it breaks!!!"); 
            break; // extra safety
        }   

        // Fleming integrals
        const double fF = ve.fFleming(v, Cu, mu, vmin, 1e-12);
        const double gF = ve.gFleming(v, Cu, mu, vmin, 1e-12);

        if (!std::isfinite(fF) || !std::isfinite(gF) || gF == 0.0){//printf("This is where it breaks 2!!!"); 
                                                                    break;}

        const double r1v = r1_from_fg(fF, gF);

        const double r2v = r2_from_fg(fF, gF);

        ////printf("The r1v and r2v values %f %f\n", r1v, r2v);


        if (!std::isfinite(r1v) || !std::isfinite(r2v) || !std::isfinite(Gv)){
            //printf("This is the place %f\n", v);
            break;
            }
        // inward curvature detection (inner wall should keep decreasing as v increases)
        if (r1v > last_r1) {
            inwardDetected = true;
            v_detect = std::round(0.65 * v);
            //r_detect = last_r1;

            // Fleming integrals
            const double fF2 = ve.fFleming(v_detect, Cu, mu, vmin, 1e-12);
            const double gF2 = ve.gFleming(v_detect, Cu, mu, vmin, 1e-12);

            if (!std::isfinite(fF2) || !std::isfinite(gF2) || gF2 == 0.0) break;
            r_detect = r1_from_fg(fF2, gF2);


            //printf("v_detect and r_detect %f %f\n", v_detect, r_detect);
            break; // stop collecting raw tail
        }

        v_s.push_back(v);
        G_s.push_back(Gv);
        r1_s.push_back(r1v);
        r2_s.push_back(r2v);
        f_s.push_back(fF);

        last_r1 = r1v;
       
    }

    //export important parameters of this calculation.
    WriteColumnsToFile("test_data_cpp.dat", v_s,G_s,r1_s, r2_s,f_s);
    //this checks out.. the core calculations are perfectly matching with Matlab's.

    //printf("\nTHIS IS A LOCATION OF BREAK 1 %f\n", r1_s.size());

    // If no inward curvature: return exactly like you intended (inner then outer; duplicated energies)
    if (!inwardDetected) {
         //
        // Build final axes exactly like code expects: inner then outer, energies duplicated
        r1_s.erase(r1_s.begin()); //erases sigularity induced problematic points at the beginning.
        r2_s.erase(r2_s.begin());
        G_s.erase(G_s.begin());
        //erase the end values too
        r1_s.erase(r1_s.end()); //erases sigularity induced problematic points at the end
        r2_s.erase(r2_s.end());
        G_s.erase(G_s.end());


        std::vector<double> r_axis = r1_s;
        r_axis.insert(r_axis.end(), r2_s.begin(), r2_s.end());
        std::vector<double> E_axis = G_s;
        E_axis.insert(E_axis.end(), E_axis.begin(), E_axis.end()); // duplicate to match two branches

        E.insert(E.end(), E_axis.begin(), E_axis.end());
        r.insert(r.end(), r_axis.begin(), r_axis.end());



        //save in file for further analysis
        printf("The file is being generated....\n");
        WriteRGToFile("Evsr_complete.dat",r, E );
        printf("The file generation is complete, Evsr_complete.dat\n");

        return;
    }

    // 1) Build turningPoints-like arrays exactly
    std::vector<double> r_axis;
    std::vector<double> E_axis;
    //r_axis.reserve(2 * r1_s.size());
    //E_axis.reserve(2 * r1_s.size());

    for (std::size_t i = 0; i < r1_s.size(); ++i) {
        r_axis.push_back(r1_s[i]);   // r1(v)
        E_axis.push_back(G_s[i]);    // Te + G(v)

        r_axis.push_back(r2_s[i]);   // r2(v)
        E_axis.push_back(G_s[i]);    // Te + G(v)
    }

    // 2) Sort by r (like [turningPointsOrdered, I] = sort(turningPoints))
    SortPaired(E_axis, r_axis);
    gsl_interp_accel* acc_d   = gsl_interp_accel_alloc();
    gsl_spline* spl_d         = gsl_spline_alloc(gsl_interp_steffen, r_axis.size());
    gsl_spline_init(spl_d, r_axis.data(), E_axis.data(), r_axis.size());


    std::vector<double> V_vals, dV_vals, d2V_vals;
    //spline_eval_cpp(r_spline, E_spline, distanceDomain, V_vals, dV_vals, d2V_vals);


    //test splining######################################################
    std::vector<double> r_branch = r_axis;
    std::vector<double> E_branch = E_axis;

    std::vector<double> r_spl_branch;
    std::vector<double> E_spl_branch;
    
    //std::reverse(r1_branch.begin(), r1_branch.end());
    //std::reverse(E1_branch.begin(), E1_branch.end());
    
    //keeps r stictly monotonic.
    //SortPaired(E1_branch, r1_branch);



    SplineNaK::Spline s;
    s.setPoints(r_branch, E_branch);
    
    int NN=500000;
    //domain
    //double xx = *std::min_element(std::begin(r1_branch), std::end(r1_branch));
    //xx = xx;

    
    // sanity check
    if (E_branch.size() != r_branch.size()) {
        throw std::runtime_error("E1_branch and r1_branch must have same size");
    }

    // iterator to max element in E1_branch
    auto xx_rmin = *std::min_element(
        std::begin(r_branch),
        std::end(r_branch)
    );

     // iterator to max element in E1_branch
    auto xx_rmax = *std::max_element(
        std::begin(r_branch),
        std::end(r_branch)
    );

    //printf("The minimum and max of r1_branch: %f %f\n", xx_rmin, xx_rmax);

    double dx = 1e-5;      // like distanceDomain step in MATLAB
    

    //find different point for xx_rmax;
    //double xx_rmax2 =  (xx_rmax + 15.0*xx_rmin)/16.0;
    double xx_rmax2 =  (xx_rmax + xx_rmin)/2.0;

    //printf("The xx_rmin is %f and the new xx_rmax position is %f and r_detect is %f\n", xx_rmin, xx_rmax, r_detect);



    for(double xx=  xx_rmax-2.0*dx; xx>=xx_rmin-dx ; xx=xx-dx){
    
        r_spl_branch.push_back(xx);
        E_spl_branch.push_back(s(xx));
    }


    WriteRGToFile("inner_wall_before_spline.dat", r_branch, E_branch);
    WriteRGToFile("inner_wall_after_spline.dat", r_spl_branch, E_spl_branch);


    //now we have a segmented spline.

    //lets estimate derivatives, first and second at few points back from the tail.
    if (r_spl_branch.size() < 4)
        throw std::runtime_error("r1_branch must have at least 4 elements");

    //the last element of the array
    std::size_t idx_target = r_spl_branch.size() - 32;
    int dummy_idx= 1;
    double r_eval = r_spl_branch[idx_target];

    double dV   = s.deriv(r_eval);
    double d2V = s.deriv2(r_eval);

    // Now power estimate.
    // Given: rEval, Vp (= dV/dr), Vpp (= d2V/dr2).

    double n_estimate_rounded = 12.0;  // default fallback

    if (std::isfinite(r_eval) && std::isfinite(dV) && std::isfinite(d2V) && std::abs(dV) >= 1e-14) {
        double n_estimate = -(d2V / dV) * r_eval - 1.0;

        if (std::isfinite(n_estimate)) {
            n_estimate_rounded = std::round(n_estimate);

            //printf("The n estimate rounded power for inner wall is: %f . dV = %f. d2V= %f\n. r_eval = %f at %d\n", n_estimate_rounded, dV, d2V, r_eval, dummy_idx);
        }
    }


    double n_sum = 0.0;
    int    n_cnt = 0;


    //canonical stl alrogrith
    auto it = std::min_element(r_spl_branch.begin(), r_spl_branch.end(),
    [&](double a, double b) {
        return std::abs(a - r_detect) < std::abs(b - r_detect);
    });

    std::size_t idx = std::distance(r_spl_branch.begin(), it);
    double closest_value = *it;

    //printf("The closest value in the r_branch is %f , indexed at %d\n", closest_value, idx);



    int pivot_idx = idx;
    int start_idx = pivot_idx - 100;
    int end_idx   = pivot_idx + 100;
    
    int size_frac = int(r_spl_branch.size()/200.0);
    //printf("The 0.02 of the array is : %d", size_frac);

    for (int k = start_idx; k < end_idx; ++k) {
        r_eval = r_spl_branch[k];
    
        dV   = s.deriv(r_eval);
        d2V = s.deriv2(r_eval);

        if (!std::isfinite(r_eval) ||
            !std::isfinite(dV) ||
            !std::isfinite(d2V) ||
            std::abs(dV) < 1e-14)
            continue;

        double n_i = -(d2V / dV) * r_eval - 1.0;

        if (!std::isfinite(n_i))
            continue;

        n_sum += n_i;
        ++n_cnt;

        ////printf("idx = %d, n_i = %f  (r=%f, dV=%f, d2V=%f)\n", k, n_i, r_eval, dV, d2V);
    }

    // fallback if nothing valid
    n_estimate_rounded = 12.0;

    if (n_cnt > 0) {
        double n_avg = n_sum / static_cast<double>(n_cnt);
        n_estimate_rounded = std::round(n_avg);
    }

    //printf("Averaged inner-wall power n = %f (from %d points)\n", n_estimate_rounded, n_cnt);


    //std:://printf("power n = %f (r=%g, V'=%g, V''=%g)\n", n_rounded_e, rEval, Vp, Vpp);

    //###################################################################

    double dr = 0.0001;      // like distanceDomain step in MATLAB
    
    std::vector<double> r_spline, E_spline;

    CubicSplineFit(r_axis, E_axis, dr, r_spline, E_spline);
    //NotAKnotSplineFit(r_axis, E_axis, dr, r_spline, E_spline);


    WriteRGToFile("Evsr_data_cpp_spline.dat", r_spline, E_spline);

    // --- PASS 2: rebuild tail with inverse-power inner wall ---

    const size_t N = r1_s.size();

    // Convert negative concavity v-level (v_detect) into a vibrational index
    int vibIndex_nc = static_cast<int>(std::round((v_detect - vmin) / space));
    if (vibIndex_nc < 0) vibIndex_nc = 0;
    if (vibIndex_nc >= static_cast<int>(r1_s.size()))
        vibIndex_nc = static_cast<int>(r1_s.size()) - 1;

    //printf("The vibrational quantum number detected: %f\n", vibIndex_nc);

    
    RGWriter writer("first_derivatives.txt");
    RGWriter writer2("first_derivatives_feed-insE.dat");
    RGWriter writer3("first_derivatives_feed-insr.dat");
    // New estimatePowerN: use MATLAB spline derivatives at r1_s via dV_vals/d2V_vals
    auto estimatePowerN = [&](const std::vector<double>& r1_levels,
                          int vibIndex_nc_local) -> double
{
    if (r1_levels.empty()) {
        std::cout << "This1\n";
        return 12.0;
    }

    // --- make range match MATLAB ---
    // MATLAB uses 1-based negativeConcavityVibLevel.
    // Our vibIndex_nc_local is 0-based → add 1 to mimic MATLAB index.
    int nc_mat = vibIndex_nc_local + 1;   // 1-based equivalent

    
    int iStart_mat = static_cast<int>(std::floor(0.5 * v_detect));
    int iEnd_mat   = static_cast<int>(std::floor(0.9 * v_detect));

    //printf("The start and end vib index are %d and %d at %f\n", iStart_mat, iEnd_mat, v_detect);

    // MATLAB indices are at least 1
    //if (iStart_mat < 1) iStart_mat = 1;

    // convert back to 0-based for C++
    int iStart = iStart_mat ;
    int iEnd   = iEnd_mat;

   // iStart = std::max(iStart, 0);
//    iEnd   = std::min(iEnd, static_cast<int>(r1_levels.size()) - 1);

    if (iStart > iEnd)
        std::swap(iStart, iEnd);
    // --- end range logic ---

    double acc = 0.0;
    int cnt = 0;

    //printf("The start and end vib index for c++ are %d and %d\n", iStart, iEnd, v_detect);

    for (int i = iStart; i <= iEnd; ++i) {

              // Fleming integrals
        const double fF_1 = ve.fFleming(i, Cu, mu, vmin, 1e-12);
        const double gF_1 = ve.gFleming(i, Cu, mu, vmin, 1e-12);

        if (!std::isfinite(fF_1) || !std::isfinite(gF_1) || gF_1 == 0.0){//printf("The integrals are not finite anymore!!!"); 
                                                                          break;}

        const double r1v_n = r1_from_fg(fF_1, gF_1);
        //const double r2v = r2_from_fg(fF, gF);


        double rEval = r1v_n;//r1_levels[static_cast<std::size_t>(i)];

        //old steffen-like cspline from gsl
        //double Vp  = gsl_spline_eval_deriv (spl_d, rEval, acc_d);
        //double Vpp = gsl_spline_eval_deriv2(spl_d, rEval, acc_d);

        //new Not-A-Knot cspline from github.
        double Vp   = s.deriv(rEval);
        double Vpp = s.deriv2(rEval);

        //printf("Vp and Vpp at %d at %f are %f and %f\n", i,rEval, Vp, Vpp);

        writer.writeValue(Vp);

        if (!std::isfinite(Vp) || !std::isfinite(Vpp) || std::abs(Vp) < 1e-14){
            //printf("This is where it breaks!!!!\n");
            continue;

        }
        double n_est = -(Vpp / Vp) * rEval - 1.0;
        if (std::isfinite(n_est)) {
            acc += n_est;
            ++cnt;
        }
    }

    if (cnt == 0) {
        std::cout << "This2\n";
        return 12.0;
    }

    double n_avg     = acc / static_cast<double>(cnt);
    double n_rounded = std::round(n_avg);

    //if (n_rounded >= 25.0)
        //n_rounded = 25.0;
    
    //if (n_rounded<0){
       // n_rounded = abs(n_rounded);
    //}


    //std::printf("MATLAB-like power n = %f\n", std::round(n_avg));
    return n_rounded;
};


    // Estimate power n using spline-smoothed potential and r1_s
    const double n = estimatePowerN(r1_s, vibIndex_nc);

    vibIndex_nc = std::clamp(vibIndex_nc, 0, (int)r1_s.size() - 1);

    // Compute smoothing indices from MATLAB rule: floor(0.85 * v_n)



    int i2_vib = static_cast<int>(std::floor(0.85 * vibIndex_nc));
    int i1_vib = i2_vib - 1;

    i1_vib = std::clamp(i1_vib, 0, (int)r1_s.size() - 1);
    i2_vib = std::clamp(i2_vib, 0, (int)r1_s.size() - 1);

    // Extract the two turning-point radii and energies
    double r1_i1 = r1_s[i1_vib];
    double r1_i2 = r1_s[i2_vib];

    double E1_i1 = G_s[i1_vib];  // already includes Te
    double E1_i2 = G_s[i2_vib];

    // Fit A and B exactly like MATLAB
    double A = 0.0, B = 1.0;

    const double d1 = calculate_r1_n(std::floor(v_detect*0.85)-1.0);
    const double d2 = calculate_r1_n(std::floor(v_detect*0.85));

    const double E1 =  calculate_G_n(std::floor(v_detect*0.85)-1.0);
    const double E2 =  calculate_G_n(std::floor(v_detect*0.85));

    fitAB(n, d1, d2, E1, E2, A, B, v_detect);
    //printf("\n %e and %e\n", A, B);

    auto inner_r_from_G = [&](double Gval) -> double {
        double denom = (Gval - B);
        //if (denom <= 0.0){//printf("Am I here?\n"); return std::numeric_limits<double>::quiet_NaN();}
        return std::pow(A / denom, 1.0 / n);
    };

    // Start rebuilding from ~v_detect (mimic MATLAB)
    const double v_tail_start = std::floor(v_detect / space) * space;
    //printf("v_tail_start %f\n", v_tail_start);

    // Emit “head”: all raw samples up to v_tail_start
    std::vector<double> r1_all, r2_all, G_all;
    for (size_t i = 0; i < v_s.size(); ++i) {
        if (v_s[i] > v_tail_start) break;
        r1_all.push_back(r1_s[i]);
        r2_all.push_back(r2_s[i]);
        G_all.push_back(G_s[i]);
    }

    // Rebuild tail using analytic inner wall and Fleming distance 2*fF(v)
    double last_inner = r1_all.empty()
                        ? std::numeric_limits<double>::infinity()
                        : r1_all.back();

    //printf("The vtail and v max are %f  and %f\n",  vmin, vmax);
    for (double v = v_tail_start + space; v <= vmax; v += space) {
    
            const double Gv = ve.G(v) + p.ex.Te;
            ////printf("%f : %f\n", v, Gv);
            if (!std::isfinite(Gv)) {
                ////printf("Skip v=%g: Gv nonfinite\n", v);
                continue;
            }

            const double fac = A/(Gv - B);
            if (!(fac > 0.0) || !std::isfinite(fac)) {
                ////printf("Skip v=%g: denom=%g\n", v, denom);
                continue;
            }

            double r1v = inner_r_from_G(Gv);
            if (!std::isfinite(r1v)) {
                ////printf("Skip v=%g: r1v nonfinite\n", v);
                continue;
            }

            if (r1v >= last_inner) r1v = std::nextafter(last_inner, 0.0);
            last_inner = r1v;

            const double fF = ve.fFleming(v, Cu, mu, vmin, 1e-8);
            if (!std::isfinite(fF)) {
                ////printf("Skip v=%g: fF nonfinite\n", v);
                continue;
            }

            const double dr = std::max(2.0 * fF, 1e-9);  // now safe
            double r2v = r1v + dr;
            if (!std::isfinite(r2v)) {
                ////printf("Skip v=%g: r2v nonfinite\n", v);
                continue;
            }
           // //printf("######### %f : %f\n", v, Gv);

            r1_all.push_back(r1v);
            r2_all.push_back(r2v);
            G_all.push_back(Gv);
        
    }

    //printf("The nans are skipped\n");

    // Build final axes exactly like code expects: inner then outer, energies duplicated
    r1_all.erase(r1_all.begin()); //erases sigularity induced problematic points at the beginning.
    r2_all.erase(r2_all.begin());
    G_all.erase(G_all.begin());
    //erase the end values too
    r1_all.erase(r1_all.end()); //erases sigularity induced problematic points at the end
    r2_all.erase(r2_all.end());
    G_all.erase(G_all.end());

    std::vector<double> r_axis_all = r1_all;
    std::vector<double> E_axis_all = G_all;

    r_axis_all.insert(r_axis_all.end(), r2_all.begin(), r2_all.end());
    E_axis_all.insert(E_axis_all.end(), E_axis_all.begin(), E_axis_all.end());

    E.insert(E.end(), E_axis_all.begin(), E_axis_all.end());
    r.insert(r.end(), r_axis_all.begin(), r_axis_all.end());


    //save in file for further analysis
    printf("The file is being generated....\n");
    WriteRGToFile("Evsr_complete.dat",r, E );
    printf("The file generation is complete, Evsr_complete.dat\n");

    gsl_spline_free(spl_d);
    gsl_interp_accel_free(acc_d);


}


//generate first derivative from fitted Not-A-Knot c-spline.
void generate_dVspline(std::vector<double> & dVspline, std::vector<double> r_range){
    



}




//calculate the De

void estimate_morse_params(double & a, double &De, double &re, double &ke, const AllParams& p){
    const double ke_i =  p.ex.ke;
    const double De_i = p.ex.De;

    //the three parameters for morse potential
    //from inpput file.
    De =  De_i;
    ke =  ke_i;
    //calculated.
    a = std::pow(ke_i/(2.0*De_i),0.5);

    //use RKR-estinate.

}


// Sorts r ascending and reorders E accordingly (in-place),
// then removes points that create unphysical large jumps.
// Criteria:
//
//  (A) If dr is tiny but dE is huge  -> drop point (vertical spike).
//  (B) If slope |dE/dr| is huge     -> drop point (outlier jump).
//
// Tune thresholds below if needed.
void SortPaired(std::vector<double>& E, std::vector<double>& r)
{
    if (r.size() != E.size()) {
        throw std::runtime_error("SortPaired: r and G must have the same size.");
    }
    if (r.empty()) return;

    // ---- 1) Sort by r ----
    std::vector<size_t> idx(r.size());
    for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;

    std::sort(idx.begin(), idx.end(),
              [&](size_t i1, size_t i2) { return r[i1] < r[i2]; });

    std::vector<double> r_sorted, E_sorted;
    //r_sorted.reserve(r.size());
    //E_sorted.reserve(E.size());

    for (size_t i = 0; i < idx.size(); ++i) {
        r_sorted.push_back(r[idx[i]]);
        E_sorted.push_back(E[idx[i]]);
    }

    // ---- 2) Jump-based cleanup ----
    // You can play with these:
    const double drTiny   = 1e-8;   // "same r" scale (Å)
    const double dELarge  = 1e-2;   // big energy jump (cm^-1); your bad case was ~4e-2
    const double slopeMax = 1e6;   // max |dE/dr| allowed (cm^-1/Å)

    std::vector<double> r_clean, E_clean;
    //r_clean.reserve(r_sorted.size());
    //E_clean.reserve(E_sorted.size());

    // keep first point
    r_clean.push_back(r_sorted[0]);
    E_clean.push_back(E_sorted[0]);

    for (size_t i = 1; i < r_sorted.size(); ++i) {

        double dr = r_sorted[i] - r_clean.back();
        double dE = E_sorted[i] - E_clean.back();

        // guard against NaNs
        if (!std::isfinite(dr) || !std::isfinite(dE) ||
            !std::isfinite(r_sorted[i]) || !std::isfinite(E_sorted[i])) {
            continue; // drop non-finite points
        }

        // (A) vertical spike: almost same r but large E jump
        if (std::abs(dr) < drTiny && std::abs(dE) > dELarge) {
            continue; // drop this point
        }

        // (B) slope spike: huge derivative jump
        if (std::abs(dr) > 0.0) {
            double slope = std::abs(dE / dr);
            if (slope > slopeMax) {
                continue; // drop this point
            }
        }

        // keep point
        r_clean.push_back(r_sorted[i]);
        E_clean.push_back(E_sorted[i]);
    }

    // ---- 3) Replace original vectors ----
    r.swap(r_clean);
    E.swap(E_clean);
}



//################################################################//
//################################################################//


// Build a cubic spline E(r) and evaluate on a fine grid.
void CubicSplineFit(const std::vector<double>& r_in,
                    const std::vector<double>& E_in,
                    double dr,
                    std::vector<double>& r_out,
                    std::vector<double>& E_out)
{
    if (r_in.size() != E_in.size()) {
        throw std::runtime_error("CubicSplineFit: r_in and E_in must have the same size.");
    }
    if (r_in.size() < 2) {
        throw std::runtime_error("CubicSplineFit: need at least 2 points for spline.");
    }
    if (dr <= 0.0) {
        throw std::runtime_error("CubicSplineFit: dr must be positive.");
    }

    for (std::size_t i = 1; i < r_in.size(); ++i) {
        if (r_in[i] <= r_in[i-1]) {
            throw std::runtime_error("CubicSplineFit: r_in must be strictly increasing.");
        }
    }

    std::size_t n = r_in.size();

    gsl_interp_accel* acc   = gsl_interp_accel_alloc();
    //gsl_spline* spline      = gsl_spline_alloc(gsl_interp_cspline, n);
    
    const gsl_interp_type* T = (n < gsl_interp_steffen->min_size)
                         ? gsl_interp_linear
                         : gsl_interp_steffen;

    gsl_spline* spline = gsl_spline_alloc(T, n);


    gsl_spline_init(spline, r_in.data(), E_in.data(), n);

    double rmin = r_in.front();
    double rmax = r_in.back();

    std::size_t Nout = static_cast<std::size_t>((rmax - rmin)/dr) + 1;

    r_out.clear();
    E_out.clear();
    r_out.reserve(Nout);
    E_out.reserve(Nout);

    for (std::size_t i = 0; i < Nout; ++i) {
        double r = rmin + i * dr;
        if (r > rmax) r = rmax;

        double E = gsl_spline_eval(spline, r, acc);

        r_out.push_back(r);
        E_out.push_back(E);
    }

    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
}


void LinearResample(const std::vector<double>& r_in,
                    const std::vector<double>& E_in,
                    double dr,
                    std::vector<double>& r_out,
                    std::vector<double>& E_out)
{
    if (r_in.size() != E_in.size())
        throw std::runtime_error("LinearResample: size mismatch.");
    if (r_in.size() < 2)
        throw std::runtime_error("LinearResample: need at least 2 points.");
    if (dr <= 0.0)
        throw std::runtime_error("LinearResample: dr must be > 0.");

    for (std::size_t i = 1; i < r_in.size(); ++i) {
        if (r_in[i] <= r_in[i-1])
            throw std::runtime_error("LinearResample: r_in must be strictly increasing.");
    }

    double rmin = r_in.front();
    double rmax = r_in.back();

    std::size_t Nout = static_cast<std::size_t>((rmax - rmin)/dr) + 1;
    r_out.clear();
    E_out.clear();
    r_out.reserve(Nout);
    E_out.reserve(Nout);

    auto interp = [&](double r) {
        auto it = std::upper_bound(r_in.begin(), r_in.end(), r);
        if (it == r_in.begin()) return E_in.front();
        if (it == r_in.end())   return E_in.back();

        std::size_t i1 = std::size_t(it - r_in.begin());
        std::size_t i0 = i1 - 1;

        double x0 = r_in[i0], x1 = r_in[i1];
        double y0 = E_in[i0], y1 = E_in[i1];

        double t = (r - x0) / (x1 - x0);
        return y0 + t*(y1 - y0);
    };

    for (std::size_t i = 0; i < Nout; ++i) {
        double r = rmin + i*dr;
        if (r > rmax) r = rmax;

        r_out.push_back(r);
        E_out.push_back(interp(r));
    }
}


void clean_up_data(const std::vector<double>& r_raw,
                   const std::vector<double>& E_raw,
                   std::vector<double>& r_clean,
                   std::vector<double>& E_clean){
    for (size_t i = 0; i < r_raw.size(); ++i) {
        if (!std::isfinite(r_raw[i]) || !std::isfinite(E_raw[i])) continue;

        if (r_clean.empty() || r_raw[i] > r_clean.back() + 1e-8) {
            r_clean.push_back(r_raw[i]);
            E_clean.push_back(E_raw[i]);
        }
    }
}


//for morse potential
double MorseParams::a() const {
    if (De > 0.0 && ke > 0.0)
        return std::sqrt(ke / (2.0 * De));
    return 0.0;
}

double morseV(double r, const MorseParams& m) {
    const double a = m.a();
    const double e = std::exp(-a * (r - m.Re));
    return m.Te + m.De * (1.0 - e) * (1.0 - e);
}

void makeMorseCurve(const MorseParams& m,
                    double rmin, double rmax, double dr,
                    std::vector<double>& r_out,
                    std::vector<double>& V_out)
{
    r_out.clear();
    V_out.clear();
    if (dr <= 0.0 || rmax <= rmin) return;

    std::size_t N = static_cast<std::size_t>((rmax - rmin) / dr) + 1;
    r_out.reserve(N);
    V_out.reserve(N);

    double r = rmin;
    for (std::size_t i = 0; i < N; ++i, r += dr) {
        r_out.push_back(r);
        V_out.push_back(morseV(r, m));
    }
}

void makeMorseOnGrid(const MorseParams& m,
                     const std::vector<double>& r_in,
                     std::vector<double>& V_out)
{
    V_out.resize(r_in.size());
    for (std::size_t i = 0; i < r_in.size(); ++i)
        V_out[i] = morseV(r_in[i], m);
}


// Writes paired arrays (r, G) to a text file.
void WriteRGToFile(const std::string& filename,
                   const std::vector<double>& r,
                   const std::vector<double>& G)
{
    if (r.size() != G.size()) {
        throw std::runtime_error("WriteRGToFile: r and G must have the same size.");
    }

    std::ofstream fout(filename);
    if (!fout.is_open()) {
        throw std::runtime_error("WriteRGToFile: cannot open file '" + filename + "'");
    }

    fout << std::scientific << std::setprecision(8);
    for (size_t i = 0; i < r.size(); ++i) {
        fout << std::setw(15) << r[i]
             << "  "
             << std::setw(15) << G[i]
             << '\n';
    }

    fout.close();
}
