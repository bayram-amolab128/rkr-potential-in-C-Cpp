#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include "helper.hpp"
#include "SplineNaK2.h"
#include "rkr_procedure.h"


extern bool inwardDetected;
extern bool outwardDetected;

//holds all the in claculation parameters.
struct CalcParam{

    //inner wall suspect 
    size_t i_s{0};       // suspect Index in the original rkr vector
    double r_s{0.0};       // suspect Turning point distance (Å)
    double V_s{0.0};       // suspect Potential energy (cm^-1)
    double v_s{0.0};       // suspected continuous vibrational level where inner wall curves inward or outward.
    size_t v_s_level{0}; // suspected nearest integer v.


    //rkr calculation parameters.
    double vmin{0.0};
    double vmax{0.0};
    double v_step{0.01};

    //well minimum
    double V_re{0.0};
    double re{0.0};
    size_t i_re{0};

    //inverse-power wall extrapolation parameters.
    double A{0.0};
    double B{0.0};
    double n{1.0};
    double v_ex{1.0};
    double V_ex{0.0};
    size_t v_ex_level{0};
    double r_ex{1.0};
    size_t i_ex_inner{0};
    size_t i_ex_outer{0};   //index of outer wall extrapolation point

    bool detectInwardCurv{false};
    bool detectOutwardCurv{false};

    //
    double abserrtol{0.0};    //fixed abs. error tolerance.
    double relerrtol{1e-9};    //fixed rel. error tolerance.

    //
    double rmin{0.0};
    double rmax{0.0};

};


struct RKRContext
{
    rkr_procedure ve;

    double Cu{};
    double me{};
    double netcharge{};
    double mu{};
    double Te{};

    double r1_from_fg(double fF, double gF) const
    {
        const double s = std::sqrt(fF / gF + fF * fF);
        return s - fF;
    }

    double r2_from_fg(double fF, double gF) const
    {
        const double s = std::sqrt(fF / gF + fF * fF);
        return s + fF;
    }
};



AllParams instantiate_params();

void print_execution_time(auto t1, auto t2);


RKRContext CreateRKRContext(const AllParams& p, CalcParam &cp);
void calc_rawRKR(std::vector<double> &V, std::vector<double>& r, const AllParams& p, CalcParam &cp, RKRContext &rkr);
void determineExtraPIndices(std::vector<double> &V, std::vector<double> &r, CalcParam &cp);
SplineNaK::Spline SplineFit_to_RKR(std::vector<double>& V, std::vector<double>& r, std::vector<double>& V_spl, std::vector<double>& r_spl, const CalcParam &cp);

void estimateInnerWall_params( SplineNaK::Spline &s, const std::vector<double>& raw_V,const std::vector<double>& raw_r,CalcParam &cp, RKRContext &rkr);

bool detectOutwardCurvature(const std::vector<double>&V, const std::vector<double>&r,  CalcParam &cp);
void setEquilibriumPoint( std::vector<double>& V,  std::vector<double>& r, CalcParam &cp);
void construct_ExtraPCurve(RKRContext &rkr, const CalcParam &cp,
                         std::vector<double>& V_inner, std::vector<double>& r_inner,
                         std::vector<double>& V_outer, std::vector<double>& r_outer);
void combine_inner_outer_rkr(std::vector<double> &V_inner, std::vector<double>& r_inner,
                             std::vector<double> &V_outer, std::vector<double>& r_outer, 
                             std::vector<double> &raw_V, std::vector<double>& raw_r, 
                             std::vector<double> &V, std::vector<double>& r, 
                             const CalcParam &cp);

//added morse curve generators.
double morseV(double r, const MorseParams& m);
void makeMorseCurve(const MorseParams& m,
                    double rmin, double rmax, double dr,
                    std::vector<double>& r_out,
                    std::vector<double>& V_out);


std::ostream& operator<<(std::ostream& os, const CalcParam& p);
std::ostream& operator<<(std::ostream& os, const AllParams& p);