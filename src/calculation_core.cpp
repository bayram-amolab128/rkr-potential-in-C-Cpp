
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
#include <limits>    //for infs.
#include <array>
#include <cstdio>
#include <filesystem>
#include <cstddef>
#include <stdexcept>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_integration.h>

//external libraries.
#include "SplineNaK2.h"

//user defined headers.
//#include "rkr_procedure.h"


//header 
#include "calculation_core.hpp"


bool inwardDetected = false;
bool outwardDetected = true;
bool useFleming = true; //true: use Fleming integrals, false: use Klein integrals.

RKRContext CreateRKRContext(const AllParams& p, CalcParam &cp)
{
    RKRContext rkr{
        rkr_procedure(
            p.gv.we, p.gv.xwe, p.gv.ywe, p.gv.zwe, p.gv.awe,
            p.gv.bwe, p.gv.cwe, p.gv.dwe, p.gv.ewe, p.gv.fwe,p.gv.gwe, p.gv.hwe, p.gv.iwe, p.gv.jwe, p.gv.kwe, p.gv.lwe,
            p.bv.Be, p.bv.ae, p.bv.ye, p.bv._1e, p.bv._2e, p.bv._3e,p.bv._4e, p.bv._5e, p.bv._6e, p.bv._7e, p.bv._8e, p.bv._9e, p.bv._10e, p.bv._11e, p.bv._12e,p.bv._13e,
            p.ex.Te, p.ex.De, p.ex.ke, p.ex.Re,
            p.ex.Vmax, p.ex.space, p.ex.ladderspace, p.ex.UseKaiser
        )
    };

    rkr.Cu = 16.857629206;
    rkr.me = 0.000548579909;
    rkr.netcharge = 0.0;

    const double m1 = p.ex.m1;
    const double m2 = p.ex.m2;

    rkr.mu = (m1 * m2) / (m1 + m2 - rkr.netcharge * rkr.me);
    rkr.Te = p.ex.Te;

    //const double vmax_trunc = rkr.ve.checkAndTruncateVmax();
    const double vmax = p.ex.Vmax;//std::min(p.ex.Vmax, vmax_trunc);
    const double vmin = rkr.ve.kaiser_correction_vmin();

    cp.v_step = p.ex.space;
    cp.vmin   = vmin;
    cp.vmax   = vmax;

    return rkr;
}


////calculates raw RKR data (V,r) without worrying about unphysical behaviour in the inner wall.
void calc_rawRKR(std::vector<double> &V, std::vector<double>& r,
          const AllParams& p, CalcParam &cp, RKRContext &rkr)
{
    //don't abort the program if error is received.
    //gsl_set_error_handler_off();
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);

  
    //control variables.
    double abserrtol      = cp.abserrtol;    //fixed abs. error tolerance.
    double relerrtol      = cp.relerrtol;    //fixed rel. error tolerance.
    //errors.
    double errfF         = 0.0;      //absolute error from gsl_integration_qags.
    double errgF         = 0.0;

    double errfF2         = 0.0;      //absolute error from gsl_integration_qags.
    double errgF2         = 0.0;


    //in-loop variables.
    double fF  =0.0;
    double gF  =0.0;                                                  
    double r1v =0.0;                                               
    double r2v =0.0;
    double last_r1 = std::numeric_limits<double>::infinity();
    double last_Gv =  std::numeric_limits<double>::infinity();
    double last_v =  std::numeric_limits<double>::infinity();
    double Gv  =0.0;
    double v_s  =0.0, r_s;

    //stl containers for calculated values. 
    std::vector<double> v_;
    std::vector<double> fF_;
    std::vector<double> gF_; 
    std::vector<double> r1;
    std::vector<double> r2;
    std::vector<double> G;
    //for errors
    std::vector<double> err_fF;
    std::vector<double> err_gF;

    std::vector<double> err_fF2;
    std::vector<double> err_gF2;


    const double vmax     =  cp.vmax;
    const double vmin     =  cp.vmin-0.01; //Small numerical offset to smoothen bottom of the potential well.
    const double v_step   =  cp.v_step;
    cp.V_re  =  rkr.Te;
    //calculates the turning points, Gv, i.e., raw RKR.
    for (double v = vmin; v <= vmax; v += v_step) {//-rkr_step-abserrtol
        Gv = rkr.ve.G(v) + rkr.Te;
        if (rkr.ve.dG(v) <= 0.0){printf("dG(v) is less than 0.0!\n"); 
            break; // extra safety
        }   

        // Fleming integrals called from rkr_procedure.h
        fF = rkr.ve.f(useFleming,v, rkr.Cu, rkr.mu, vmin, abserrtol,relerrtol,&errfF,w);
        gF = rkr.ve.g(useFleming,v, rkr.Cu, rkr.mu, vmin, abserrtol,relerrtol,&errgF,w);

        //evaluation errors.
        if (!std::isfinite(fF) || !std::isfinite(gF) || gF == 0.0){
            //printf("Undefined or nan for fF=%g, gF=%g, or gF=%g is 0.0\n", fF, gF, gF); 
            //break;
            continue; // extra safety, skip this v and continue.
        }
 

        v_.push_back(v);
        fF_.push_back(fF);
        gF_.push_back(gF);
        err_fF.push_back(errfF);
        err_gF.push_back(errgF);


        //calculating turning points.                                                           
        r1v = rkr.r1_from_fg(fF, gF);
        r2v = rkr.r2_from_fg(fF, gF);

        if (!std::isfinite(r1v) || !std::isfinite(r2v) || !std::isfinite(Gv)){
            //printf("nan r1v, r2v, or Gv at v=%f\n",v);
            break;
        }


        // inward curvature detection (inner wall should keep decreasing as v increases)
        if (r1v > last_r1) {
            cp.detectInwardCurv = true;
            cp.detectOutwardCurv = false;   // because if inward curvature, outward curvature wouldn't occur.
           
            //collect inward curvature suspect point.
            cp.v_s =  last_v;
            cp.v_s_level =  std::round(last_v);
            cp.r_s =  last_r1;
            cp.V_s =  last_Gv;
            cp.i_s =  0; //it's at the beginning of sorted data.

            const double v_ex = 0.65 * v;        //go down to 65% of suspect v to avoid curvature.
            const double v_ex_level = std::round(v_ex);

            fF = rkr.ve.f(useFleming,v_ex_level, rkr.Cu, rkr.mu, vmin, abserrtol,relerrtol,&errfF,w);
            gF = rkr.ve.g(useFleming,v_ex_level, rkr.Cu, rkr.mu, vmin, abserrtol,relerrtol,&errgF,w);

            cp.v_ex =  v_ex;
            cp.v_ex_level =  v_ex_level;

            //collect 65% extrapolation point.//not necessary.
            cp.r_ex =  rkr.r1_from_fg(fF, gF);       //set at v_ex_level*0.85 later.
            cp.V_ex =  rkr.ve.G(v_ex_level) + rkr.Te;//set at v_ex_level*0.85 later.

            //printf("Inward curvature detected at v = %f, select v_ex_level at 0.65 of v for extrapolation\n", v);
            break; // stop collecting raw tail
        }
        else {
            cp.detectInwardCurv = false;
            cp.detectOutwardCurv = true;
        }

        last_r1 = r1v;
        last_Gv =  Gv;
        last_v = v;

        //output.
        G.push_back(Gv);
        r1.push_back(r1v);
        r2.push_back(r2v);
    
    }

    //write errors as the function of v.
    //WriteRGToFile("./output/err_fF_vs_v.dat", v_, err_fF);
    //WriteRGToFile("./output/err_gF_vs_v.dat", v_, err_gF);

    //WriteRGToFile("./output/fFleming_vs_v.dat", v_, fF_);
    //WriteRGToFile("./output/gFleming_vs_v.dat", v_, gF_);


    //printf("Combining inner and outer branches together for raw RKR.\n");
    for (std::size_t i = 0; i < 2*r1.size(); ++i) { //r1.size() = r2.size() =G.size()
        if(i<r1.size()){
            r.push_back(r1[i]);      // r1(v)
            V.push_back(G[i]);    // Te + G(v)
        }
        else{
            r.push_back(r2[i-r1.size()]);      // r2(v)
            V.push_back(G[i-r1.size()]);    // Te + G(v)
        }
    }

    //printf("Sorting r in increasing order, and rearranging V(r) accordingly for raw RKR.\n");
    SortPaired(V, r);
    //also find the v_ex, r_ex index.


    //WriteRGToFile("./output/rawRKR.dat", r, V);

    //printf("The program exits raw rkr calculation loop.\n");
    gsl_integration_workspace_free(w);
}//end of calc.

void setEquilibriumPoint( std::vector<double>& V,  std::vector<double>& r, CalcParam &cp)
{
    if (V.size() != r.size()) {
        throw std::runtime_error("determineEquilibriumParameters: V and r must have the same size.");
    }

    std::vector<size_t> i_re = findMinimumIndices(V, r);
    cp.i_re = average(i_re);

    double re_sum = 0.0;
    for (size_t idx : i_re) {
        re_sum += r[idx];
    }
    cp.re = re_sum / static_cast<double>(i_re.size());
    //cp.V_re previously stored, Te. 
    V.push_back(cp.V_re);
    r.push_back(cp.re);

    SortPaired(V, r); // Ensure V and r are sorted after adding the equilibrium point.

    //WriteRGToFile("./output/rawRKR.dat", r, V);
    //printf("Equilibrium parameters determined: re = %f, V_re = %f\n", cp.re, cp.V_re);
}

//calculates the calc partameters such as v_ex, re, ...
void determineExtraPIndices(std::vector<double> &V, std::vector<double> &r, CalcParam &cp){
    
    const size_t i_ex_inner = findClosestIndex(r, cp.r_ex);
    cp.i_ex_inner =   i_ex_inner; //still need to calculate outer.
    //
    //std::cout<<"i_ex_inner "<<cp.i_ex_inner<<std::endl;


    size_t i_ex_outer = splitandfindClosestIndex(V, cp.i_re, cp.V_ex) + cp.i_re;
    cp.i_ex_outer =  i_ex_outer;

}


//sets the spline, also pre-calculates the spline fitted r_spl and V_spl
SplineNaK::Spline SplineFit_to_RKR(std::vector<double>& V, std::vector<double>& r,
                      std::vector<double>& V_spl, std::vector<double>& r_spl, const CalcParam &cp){
    
    SplineNaK::Spline s;
    //spline containers
    s.setPoints(r, V);
    
    V_spl.clear();
    r_spl.clear();
    // sanity check, if V and r have same size.
    if (V.size() != r.size()) {
        throw std::runtime_error("V and r must have same size\n");
    }

    // iterator to min element in r
    //auto rmin = *std::min_element(
    //    std::begin(r),
    //    std::end(r)
    //);

     // iterator to max element in r
   //auto rmax = *std::max_element(
   //    std::begin(r),
   //    std::end(r)
   //);
   double rmin = cp.rmin;
   double rmax = cp.rmax;


    //printf("The minimum and maximum of r: %f %f\n", rmin, rmax);

    double dx    = 1e-3;               // spline step size.
    double rmax2 =  (rmax + rmin)/2.0; //midpoint

    //generate the spline curve.
    for(double x = rmin-dx; x <= rmax-2.0*dx; x += dx){
        r_spl.push_back(x);
        V_spl.push_back(s(x));
    }
    //WriteRGToFile("./output/spline_rawRKR.dat", r_spl, V_spl);
    //printf("Not-A-Knot cubic spline is fitted through the raw RKR data\n");

    return s; //returns the evaluated spline.
}


bool detectOutwardCurvature(const std::vector<double>& V,
                       const std::vector<double>& r,
                       CalcParam &cp)
{
    std::vector<size_t> i_res =  findMinimumIndices(V,r);
    cp.i_re = average(i_res);
    //std::cout<<"The index for re = "<<cp.i_re<<std::endl;
    double lastCurvature = 1.0;

    const size_t i_re = cp.i_re;
    if (V.size() != r.size())
        throw std::runtime_error("detectOutwardCurvature: V and r size mismatch.");

    if (V.size() < 3){
        //std::cout<<"Not enough data points for curvature detection, returning no curvature."<<std::endl;
        return false;
    }
    if (i_re >= V.size())
        throw std::runtime_error("detectOutwardCurvature: i_re out of range.");

    // inner branch is from i_re going backward toward smaller r
    for (size_t i = i_re; i >= 2; --i)
    {
        double r0 = r[i];
        double r1 = r[i - 1]; //center.
        double r2 = r[i - 2];
        //std::cout<<i<<std::endl;
        double V0 = V[i];
        double V1 = V[i - 1];
        double V2 = V[i - 2];

        double dr1 = r0 - r1; //r1 needs to be smaller than r0, so dr1 should be positive.
        double dr2 = r1 - r2;

        double slope1 = (V0 - V1) / dr1;
        double slope2 = (V1 - V2) / dr2;

        double curvature = slope1 - slope2;
        //std::cout<<"curvature= "<<curvature<<std::endl;
        
        if(curvature<0.0 && lastCurvature>0.0){
            //std::cout<<"Warning: outward curvature detected at r = "<<r1<<std::endl;
            //calculate v-level
            
            cp.v_s=cp.vmax- i*cp.v_step;
            cp.v_s_level=std::floor(cp.v_s);
            cp.r_s=r1;
            cp.V_s=V1;
            cp.i_s=i;  
            cp.v_ex = cp.v_s;//if outward curvature detected, v_ex is set to v_s.
            cp.v_ex_level = std::floor(0.65*cp.v_ex);

            ///not necessary to calculate r_ex and V_ex here, as they will be calculated later.
            cp.V_ex = V1; //set at v_ex_level*0.85 later.
            cp.r_ex = r1; //set at v_ex_level*0.85 later.
            return true;  // indicates outward curvature detected, which is unphysical.
        }

        lastCurvature = curvature;
    }
    //printf("No outward curvature detected in the inner branch.\n");
    return false;           // no outward curvature detected
}

void estimateInnerWall_params( SplineNaK::Spline &s, 
                              const std::vector<double>& raw_V,
                              const std::vector<double>& raw_r,
                              CalcParam &cp, RKRContext &rkr)
{

    gsl_integration_workspace* w =gsl_integration_workspace_alloc(limit);
    //control variables.
    double abserrtol      = cp.abserrtol;    //fixed abs. error tolerance.
    double relerrtol      = cp.relerrtol;    //fixed rel. error tolerance.
    //errors.
    double errfF         = 0.0;      //absolute error from gsl_integration_qags.
    double errgF         = 0.0;

    double errfF2         = 0.0;      //absolute error from gsl_integration_qags.
    double errgF2         = 0.0;
    //in-loop variables.
    double fF  =0.0;
    double gF  =0.0;   
    double vmin =  cp.vmin;

    //estimate the n
    double n_i  =  0.0;
    //double n   =  double(r_slice.size());
    int n_avg   =  1;

    double v_step = cp.v_step;

    double vmax   = cp.vmax;

    //       
    //std::ofstream out("./output/n_est_data.dat");
    //if (!out) {
     //   throw std::runtime_error("Could not open output.dat");
    //}
    //out << std::fixed << std::setprecision(14);

    double d2V=0.0;
    double d1V=0.0;
    double r_eval=0.0;
    double n=0.0;


    double v_ex = cp.v_ex; 
    double v_ex_level = static_cast<double> (cp.v_ex_level); 

    double v_init =  std::floor(0.5*v_ex_level);
    double v_final = std::floor(0.9*v_ex_level);

    for(double v = v_init; v<=v_final; ++v){
        n=n+1.0;
        
        // Fleming integrals called from rkr_procedure.h


        fF = rkr.ve.f(useFleming,v, rkr.Cu, rkr.mu, vmin, abserrtol,relerrtol,&errfF,w);
        gF = rkr.ve.g(useFleming,v, rkr.Cu, rkr.mu, vmin, abserrtol,relerrtol,&errgF,w);
   
        //evaluation errors.
        if (!std::isfinite(fF) || !std::isfinite(gF) || gF == 0.0){
           // printf("Undefined or nan for fF=%g, gF=%g, or gF=%g is 0.0\n", fF, gF, gF); 
            break;
        }
 
        r_eval= rkr.r1_from_fg(fF, gF);
        d2V  = s.deriv2(r_eval);
        d1V  = s.deriv(r_eval);

         //n_i = n_i - d2V_slice[i]/d1V_slice[i] * r_slice[i] - 1.0;
        n_i = n_i - d2V/d1V * r_eval - 1.0;
        //out <<v<< '\t' <<r_eval << '\t' << d1V<< '\t' << d2V<<'\n';
    }


    //out.close();
    //average the n_is.
    n_avg = std::round(n_i/n);
    //std::cout<<"n_i "<< n_i << " n "<<n<<std::endl;


    double v_1 = std::floor(v_ex_level*0.85)-1.0;
    // Fleming integrals called from rkr_procedure.h

    fF = rkr.ve.f(useFleming,v_1, rkr.Cu, rkr.mu, vmin, abserrtol,relerrtol,&errfF,w);
    gF = rkr.ve.g(useFleming,v_1, rkr.Cu, rkr.mu, vmin, abserrtol,relerrtol,&errgF,w);

    double r1= rkr.r1_from_fg(fF, gF);
    double V1= rkr.ve.G(v_1) + rkr.Te;


    double v_2 = std::floor(v_ex_level*0.85);
    // Fleming integrals called from rkr_procedure.h

    fF = rkr.ve.f(useFleming,v_2, rkr.Cu, rkr.mu, vmin, abserrtol,relerrtol,&errfF,w);
    gF = rkr.ve.g(useFleming,v_2, rkr.Cu, rkr.mu, vmin, abserrtol,relerrtol,&errgF,w);

    double r2=rkr.r1_from_fg(fF, gF);
    double V2= rkr.ve.G(v_2) + rkr.Te;


    double A = (V1 - V2)/((1.0/std::pow(r1,n_avg)) - (1.0/std::pow(r2, n_avg))) ;
    //std::cout<<"A "<<A<<" "<< V1<<" "<<V2<<std::endl;
    double B = V2 - A/(std::pow(r2, n_avg));

    //construct_InnerWall( A,B, n_avg);
    //(Asmooth/r^(integerPower)) + Bsmooth;
    //printf("r1=%f r2=%f and V1=%f  V2=%f.\n Estimated n-power law value, A and B, %d, %g, and %g\n",r1, r2, V1, V2, n_avg, A, B);

    cp.A = A;
    cp.B = B;
    cp.n = double(n_avg);

    //also set extraploation r and V
    cp.r_ex = r2;
    cp.V_ex = V2;
    //std::cout<<"r_ex "<<cp.r_ex<<" V_ex "<<cp.V_ex<<std::endl;

    gsl_integration_workspace_free(w);

}


void construct_ExtraPCurve(RKRContext &rkr, const CalcParam &cp,
                         std::vector<double>& V_inner, std::vector<double>& r_inner,
                         std::vector<double>& V_outer, std::vector<double>& r_outer){
    gsl_integration_workspace* w =gsl_integration_workspace_alloc(limit);
    //(Asmooth/r^(integerPower)) + Bsmooth;
    double Gv=0.0;
    double rv_=0.0;
    double r1=1.0;
    double r2_=1.0;
    double r2=1.0;
    double fF=0.0;
    double gF=0.0;

    double vmin = cp.vmin;
    double vmax = cp.vmax;
    double v_step = cp.v_step;
    double abserrtol      = cp.abserrtol;    //fixed abs. error tolerance.
    double relerrtol      = cp.relerrtol;    //fixed rel. error
    double v_start = 0.85*cp.v_ex_level;     //start extrapolation from v_ex_level.

    double n = double(cp.n);
    double A = cp.A;
    double B = cp.B;

    //std::cout<<A<<" "<<B<<" "<<n<<std::endl;

    for (double v = v_start; v <= vmax; v +=v_step) {
        Gv = rkr.ve.G(v) + rkr.Te;
        
        rv_ = cp.A/(Gv-cp.B);

        r1 = std::pow(rv_, 1.0/n);

   
        fF = rkr.ve.f(useFleming,v, rkr.Cu, rkr.mu, vmin, abserrtol,relerrtol, nullptr, w);


        //this is where Near Dissociation extrapolation need to be applied.
        r2_ = 2.0*fF;
        r2 = r2_+r1;
       

        //fill the inner and outer wall vectors
        r_inner.push_back(r1);
        V_inner.push_back(Gv);
        r_outer.push_back(r2);
        V_outer.push_back(Gv);
    }

    SortPaired(V_inner, r_inner);
    SortPaired(V_outer, r_outer);
    //WriteRGToFile("./output/extraP_Inner.dat", r_inner, V_inner);
    //WriteRGToFile("./output/extraP_Outer.dat", r_outer, V_outer);
    gsl_integration_workspace_free(w);
}

void combine_inner_outer_rkr(std::vector<double>& V_inner, std::vector<double>& r_inner,
                             std::vector<double>& V_outer, std::vector<double>& r_outer,
                             std::vector<double>& raw_V,   std::vector<double>& raw_r,
                             std::vector<double>& V,       std::vector<double>& r,
                             const CalcParam& cp)
{
    V.clear();
    r.clear();

    // Inner wall
    for (size_t i = 0; i < r_inner.size(); ++i) {
        r.push_back(r_inner[i]);
        V.push_back(V_inner[i]);
    }

    // Raw RKR, starting from patch index inP.i_s
    size_t start = static_cast<size_t>(cp.i_ex_inner);
    size_t end = static_cast<size_t>(cp.i_ex_outer);

    if (start > raw_r.size())
        start = raw_r.size();

    for (size_t i = start; i < end; ++i) {
        r.push_back(raw_r[i]);
        V.push_back(raw_V[i]);
    }

    // Outer wall
    for (size_t i = 0; i < r_outer.size(); ++i) {
        r.push_back(r_outer[i]);
        V.push_back(V_outer[i]);
    }

    //WriteRGToFile("./output/totalPEC.dat", r, V);
}


//===morse curve for overlay.. not yet complete.
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



std::ostream& operator<<(std::ostream& os, const CalcParam& p)
{
    os << std::fixed << std::setprecision(16);

    os << "========== CalcParam ==========\n";

    os << "\n--- Inner-wall suspect ---\n";
    os << "i_s       : " << p.i_s << '\n';
    os << "r_s       : " << p.r_s << '\n';
    os << "V_s       : " << p.V_s << '\n';
    os << "v_s       : " << p.v_s << '\n';
    os << "v_s_level : " << p.v_s_level << '\n';

    os << "\n--- RKR parameters ---\n";
    os << "vmin      : " << p.vmin << '\n';
    os << "vmax      : " << p.vmax << '\n';
    os << "v_step    : " << p.v_step << '\n';

    os << "\n--- Well minimum ---\n";
    os << "V_re       : " << p.V_re << '\n';
    os << "re         : " << p.re << '\n';
    os << "i_re       : " << p.i_re << '\n';

    os << "\n--- Inverse-power extrapolation ---\n";
    os << "A          : " << std::scientific << std::setprecision(20)<<p.A << '\n';
    os << "B          : " << std::scientific << std::setprecision(20)<< p.B << '\n';
    os << "n          : " << std::fixed << std::setprecision(16)<< p.n << '\n';
    os << "v_ex       : " << p.v_ex << '\n';
    os << "V_ex       : " << p.V_ex << '\n';
    os << "v_ex_level : " << p.v_ex_level << '\n';
    os << "r_ex       : " << p.r_ex << '\n';
    os << "i_ex_inner : " << p.i_ex_inner << '\n';
    os << "i_ex_outer : " << p.i_ex_outer << '\n';

    os << "===============================\n";

    return os;
}


std::ostream& operator<<(std::ostream& os, const AllParams& p)
{
    os << std::scientific << std::setprecision(16);

    os << "========== AllParams ==========\n";

    os << "\n--- General ---\n";
    os << "state_name   : " << p.state_name << '\n';

    os << "\n--- G(v) coefficients ---\n";
    os << "we           : " << p.gv.we  << '\n';
    os << "xwe          : " << p.gv.xwe << '\n';
    os << "ywe          : " << p.gv.ywe << '\n';
    os << "zwe          : " << p.gv.zwe << '\n';
    os << "awe          : " << p.gv.awe << '\n';
    os << "bwe          : " << p.gv.bwe << '\n';
    os << "cwe          : " << p.gv.cwe << '\n';
    os << "dwe          : " << p.gv.dwe << '\n';
    os << "ewe          : " << p.gv.ewe << '\n';
    os << "fwe          : " << p.gv.fwe << '\n';
    os << "gwe          : " << p.gv.gwe << '\n';
    os << "hwe          : " << p.gv.hwe << '\n';
    os << "iwe          : " << p.gv.iwe << '\n';
    os << "jwe          : " << p.gv.jwe << '\n';
    os << "kwe          : " << p.gv.kwe << '\n';
    os << "lwe          : " << p.gv.lwe << '\n';

    os << "\n--- B(v) coefficients ---\n";
    os << "Be           : " << p.bv.Be  << '\n';
    os << "ae           : " << p.bv.ae  << '\n';
    os << "ye           : " << p.bv.ye  << '\n';
    os << "_1e          : " << p.bv._1e << '\n';
    os << "_2e          : " << p.bv._2e << '\n';
    os << "_3e          : " << p.bv._3e << '\n';
    os << "_4e          : " << p.bv._4e << '\n';
    os << "_5e          : " << p.bv._5e << '\n';
    os << "_6e          : " << p.bv._6e << '\n';
    os << "_7e          : " << p.bv._7e << '\n';
    os << "_8e          : " << p.bv._8e << '\n';
    os << "_9e          : " << p.bv._9e << '\n';
    os << "_10e         : " << p.bv._10e << '\n';
    os << "_11e         : " << p.bv._11e << '\n';
    os<< "_12e         : " << p.bv._12e << '\n';
    os<< "_13e         : " << p.bv._13e << '\n';

    os << "\n--- Extra parameters ---\n";
    os << "m1           : " << p.ex.m1 << '\n';
    os << "m2           : " << p.ex.m2 << '\n';
    //os << "netcharge    : " << p.ex.netcharge << '\n';
    os << "space        : " << p.ex.space << '\n';
    os << "ladderspace  : " << p.ex.ladderspace << '\n';
    os << "Vmax         : " << p.ex.Vmax << '\n';
    //os << "errortol     : " << p.ex.errortol << '\n';
    os << "UseKaiser    : " << p.ex.UseKaiser << '\n';
    os << "Te           : " << p.ex.Te << '\n';
    os << "De           : " << p.ex.De << '\n';
    os << "Ke           : " << p.ex.ke << '\n';
    os << "re           : " << p.ex.Re << '\n';

    os << "===============================\n";

    return os;
}

