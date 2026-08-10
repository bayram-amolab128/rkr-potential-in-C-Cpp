#ifndef RKR_PROCEDURE_H
#define RKR_PROCEDURE_H


#include <cmath>
#include <utility>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <algorithm>
#include <gsl/gsl_integration.h>


const size_t limit = 100000;
//double resabs = 0.0;
//double resasc = 0.0;

class rkr_procedure {
private:
    // coefficients
    double we, xwe, ywe, zwe, awe, bwe, cwe, dwe, ewe, fwe, gwe, hwe, iwe, jwe, kwe, lwe;
    double be, ae, ye, _1e, _2e, _3e, _4e, _5e, _6e,_7e, _8e, _9e, _10e, _11e, _12e, _13we;
    double Te, De, ke, re;
    double vmax, space, ladderspace, UseKaiser;

public:
    rkr_procedure(double we_,  double xwe_,  double ywe_, double zwe_, double awe_,
                      double bwe_, double cwe_,  double dwe_, double ewe_, double fwe_, double gwe_, double hwe_, double iwe_, double jwe_, double kwe_, double lwe_,
                      double be_,  double ae_,   double ye_,  double _1e_,  double _2e_, double _3e_, double _4e_, double _5e_, double _6e_, double _7e_, double _8e_, double _9e_, double _10e_, double _11e_, double _12e_, double _13we_,
                      double Te_,  double De_,   double ke_,  double re_,
                      double vmax_, double space_, double ladderspace_, double UseKaiser_);

    double G(double vib) const;
    double dG(double vib) const;
    double d2G(double vib) const;

    double kaiser_correction_Y00() const;
    double kaiser_correction_vmin() const;
    double checkAndTruncateVmax(double step = 0.0001, double abstol = 1e-6) const;

    double fFleming(double vib, double Cu, double mu, double vmin, double abstol = 1e-14, double reltol = 1e-10, double *err = nullptr, gsl_integration_workspace* w = nullptr) const;
    double gFleming(double vib, double Cu, double mu, double vmin, double abstol = 1e-14, double reltol = 1e-10, double *err =nullptr, gsl_integration_workspace* w = nullptr) const;

    double fKlein(double vib,double Cu,double mu,double vmin,double abstol = 1e-14,double reltol = 1e-10,double* err = nullptr, gsl_integration_workspace* w = nullptr) const;
    double gKlein(double vib,double Cu,double mu,double vmin,double abstol = 1e-14,double reltol = 1e-10,double* err = nullptr, gsl_integration_workspace* w = nullptr) const;

    double f(bool useFleming, double vib, double Cu, double mu, double vmin, double abstol = 1e-14, double reltol = 1e-10, double* err = nullptr, gsl_integration_workspace* w = nullptr) const;
    double g(bool useFleming, double vib, double Cu, double mu, double vmin, double abstol = 1e-14, double reltol = 1e-10, double* err = nullptr, gsl_integration_workspace* w = nullptr) const;

    double B(double vib) const;
    double dB(double vib) const;

    double Y00_0;
    double vmax_trunc;
};


// ---------------- Constructor ----------------
inline rkr_procedure::rkr_procedure(double we_,  double xwe_,  double ywe_, double zwe_, double awe_,
                                     double bwe_, double cwe_,  double dwe_, double ewe_, double fwe_, double gwe_, double hwe_, double iwe_, double jwe_, double kwe_, double lwe_,
                                     double be_,  double ae_,   double ye_,  double _1e_,  double _2e_, double _3e_, double _4e_, double _5e_, double _6e_, double _7e_, double _8e_, double _9e_, double _10e_, double _11e_, double _12e_, double _13we_,
                                     double Te_,  double De_,   double ke_,  double re_,
                                     double vmax_, double space_, double ladderspace_, double UseKaiser_)
    : we(we_), xwe(xwe_), ywe(ywe_), zwe(zwe_), awe(awe_), bwe(bwe_),
      cwe(cwe_), dwe(dwe_), ewe(ewe_), fwe(fwe_), gwe(gwe_), hwe(hwe_), iwe(iwe_), jwe(jwe_), kwe(kwe_), lwe(lwe_),
      be(be_), ae(ae_), ye(ye_), _1e(_1e_), _2e(_2e_), _3e(_3e_), _4e(_4e_), _5e(_5e_), _6e(_6e_), _7e(_7e_), _8e(_8e_), _9e(_9e_), _10e(_10e_), _11e(_11e_), _12e(_12e_), _13we(_13we_),
      Te(Te_), De(De_), ke(ke_), re(re_),
      vmax(vmax_), space(space_), ladderspace(ladderspace_), UseKaiser(UseKaiser_)
{
    Y00_0      = kaiser_correction_Y00();
    vmax_trunc = checkAndTruncateVmax();
    // no other side effects
}


// ---------------- Kaiser pieces ----------------
inline double rkr_procedure::kaiser_correction_Y00() const {
    double Y10 = we;
    double Y01 = be;
    double Y20 = -1.0 * xwe;
    double Y11 = -1.0 * ae;

    double y00 = 0.250 * (Y01 + Y20)
                 - ((Y11 * Y10) / (12.0 * Y01))
                 + (1.0 / Y01) * std::pow(((Y11 * Y10) / (12.0 * Y01)), 2);

    return y00;
}

inline double rkr_procedure::kaiser_correction_vmin() const {
    if (UseKaiser == 1.0) {
        double Y10 = we;
        double y00 = kaiser_correction_Y00();
        //std::cout << "Kaiser correction for Y00: " << y00 << std::endl;
        return -0.5 - y00 / Y10;
    }
    return -0.5;
}


// ---------------- G(v), dG, d2G ----------------
//inline double rkr_procedure::G(double vib) const {
//    double x = vib + 0.5;
//    return we * std::pow(x, 1) + xwe * std::pow(x, 2) + ywe * std::pow(x, 3)
//         + zwe * std::pow(x, 4) + awe * std::pow(x, 5) + bwe * std::pow(x, 6)
//         + cwe * std::pow(x, 7) + dwe * std::pow(x, 8) + ewe * std::pow(x, 9)
//         + fwe * std::pow(x,10) + gwe * std::pow(x,11) + hwe * std::pow(x,12)
//         + iwe * std::pow(x,13) + jwe * std::pow(x,14) + kwe * std::pow(x,15)
//         + lwe * std::pow(x,16) + Y00_0;
//}
//horner style for G(v) and dG(v)
inline double rkr_procedure::G(double vib) const
{
    const double x = vib + 0.5;

    double y = lwe;
    y = y*x + kwe;
    y = y*x + jwe;
    y = y*x + iwe;
    y = y*x + hwe;
    y = y*x + gwe;
    y = y*x + fwe;
    y = y*x + ewe;
    y = y*x + dwe;
    y = y*x + cwe;
    y = y*x + bwe;
    y = y*x + awe;
    y = y*x + zwe;
    y = y*x + ywe;
    y = y*x + xwe;
    y = y*x + we;

    return y*x + Y00_0;
}


//inline double rkr_procedure::dG(double vib) const {
//    double x = vib + 0.5;
//    return we
//         + 2.0 * xwe * std::pow(x, 1)
//         + 3.0 * ywe * std::pow(x, 2)
//         + 4.0 * zwe * std::pow(x, 3)
//         + 5.0 * awe * std::pow(x, 4)
//         + 6.0 * bwe * std::pow(x, 5)
//         + 7.0 * cwe * std::pow(x, 6)
//         + 8.0 * dwe * std::pow(x, 7)
//         + 9.0 * ewe * std::pow(x, 8)
//         +10.0 * fwe * std::pow(x, 9)
//         +11.0 * gwe * std::pow(x,10)
//         +12.0 * hwe * std::pow(x,11)
//         +13.0 * iwe * std::pow(x,12)
//         +14.0 * jwe * std::pow(x,13)
//         +15.0 * kwe * std::pow(x,14)
//         +16.0 * lwe * std::pow(x,15);
//}

inline double rkr_procedure::dG(double vib) const
{
    const double x = vib + 0.5;

    double y = 16.0*lwe;
    y = y*x + 15.0*kwe;
    y = y*x + 14.0*jwe;
    y = y*x + 13.0*iwe;
    y = y*x + 12.0*hwe;
    y = y*x + 11.0*gwe;
    y = y*x + 10.0*fwe;
    y = y*x +  9.0*ewe;
    y = y*x +  8.0*dwe;
    y = y*x +  7.0*cwe;
    y = y*x +  6.0*bwe;
    y = y*x +  5.0*awe;
    y = y*x +  4.0*zwe;
    y = y*x +  3.0*ywe;
    y = y*x +  2.0*xwe;

    return y*x + we;
}

inline double rkr_procedure::d2G(double vib) const
{
    const double x = vib + 0.5;

    double y = 240.0 * lwe;
    y = y*x + 210.0 * kwe;
    y = y*x + 182.0 * jwe;
    y = y*x + 156.0 * iwe;
    y = y*x + 132.0 * hwe;
    y = y*x + 110.0 * gwe;
    y = y*x +  90.0 * fwe;
    y = y*x +  72.0 * ewe;
    y = y*x +  56.0 * dwe;
    y = y*x +  42.0 * cwe;
    y = y*x +  30.0 * bwe;
    y = y*x +  20.0 * awe;
    y = y*x +  12.0 * zwe;
    y = y*x +   6.0 * ywe;

    return y*x + 2.0*xwe;
}


// ---------------- B(v), dB ----------------
//inline double rkr_procedure::B(double vib) const {
//    double x = vib + 0.5;
//    return be + ae * x + ye * std::pow(x,2) + _1e * std::pow(x,3) + _2e * std::pow(x,4) + _3e * std::pow(x,5) + _4e * std::pow(x,6) + _5e * std::pow(x,7) + _6e * std::pow(x,8)+ _7e * std::pow(x,9) + _8e * std::pow(x,10) + _9e * std::pow(x,11) + _10e * std::pow(x,12) + _11e * std::pow(x,13) + _12e * std::pow(x,14) +_13we * std::pow(x,15);
//}
inline double rkr_procedure::B(double vib) const
{
    const double x = vib + 0.5;

    double y = _13we;
    y = y*x + _12e;
    y = y*x + _11e;
    y = y*x + _10e;
    y = y*x + _9e;
    y = y*x + _8e;
    y = y*x + _7e;
    y = y*x + _6e;
    y = y*x + _5e;
    y = y*x + _4e;
    y = y*x + _3e;
    y = y*x + _2e;
    y = y*x + _1e;
    y = y*x + ye;
    y = y*x + ae;

    return y*x + be;
}

//inline double rkr_procedure::dB(double vib) const {
//    double x = vib + 0.5;
//    return +ae
//           + 2.0 * ye  * std::pow(x,1)
//           + 3.0 * _1e * std::pow(x,2)
//           + 4.0 * _2e * std::pow(x,3)
//           + 5.0 * _3e * std::pow(x,4)
//           + 6.0 * _4e * std::pow(x,5)
//           + 7.0 * _5e * std::pow(x,6)
//           + 8.0 * _6e * std::pow(x,7)
//           + 9.0 * _7e * std::pow(x,8)
//           + 10.0 * _8e * std::pow(x,9)
//           + 11.0 * _9e * std::pow(x,10)
//           + 12.0 * _10e * std::pow(x,11)
//           + 13.0 * _11e * std::pow(x,12)
//           + 14.0 * _12e * std::pow(x,13)
//           + 15.0 * _13we * std::pow(x,14);
//}

inline double rkr_procedure::dB(double vib) const
{
    const double x = vib + 0.5;

    double y = 15.0 * _13we;
    y = y*x + 14.0 * _12e;
    y = y*x + 13.0 * _11e;
    y = y*x + 12.0 * _10e;
    y = y*x + 11.0 * _9e;
    y = y*x + 10.0 * _8e;
    y = y*x +  9.0 * _7e;
    y = y*x +  8.0 * _6e;
    y = y*x +  7.0 * _5e;
    y = y*x +  6.0 * _4e;
    y = y*x +  5.0 * _3e;
    y = y*x +  4.0 * _2e;
    y = y*x +  3.0 * _1e;
    y = y*x +  2.0 * ye;

    return y*x + ae;
}


// ---------------- Vmax truncation ----------------
inline double rkr_procedure::checkAndTruncateVmax(double step, double abstol) const {
    std::vector<double> roots;

    for (double v = 0.0; v <= vmax - step; v += step) {
        double f1 = dG(v);
        double f2 = dG(v + step);

        if (f1 * f2 < 0.0) {
            double a = v, b = v + step;

            for (int it = 0; it < 100; ++it) {
                double mid = 0.5 * (a + b);
                double fm  = dG(mid);

                if (std::fabs(fm) < abstol) break;

                if (dG(a) * fm < 0.0) b = mid;
                else a = mid;
            }

            double root = 0.5 * (a + b);
            if (root > 0.0) roots.push_back(root);
        }
    }

    if (!roots.empty()) {
        std::sort(roots.begin(), roots.end());
        double Gmax = roots.front();
        return (vmax > Gmax) ? Gmax : vmax;
    }

    //std::cerr << "Warning: no extremum (dG=0) found up to v=" << vmax << "\n";
    return vmax;
}


// ---------------- Fleming integrals ----------------
struct KleinParams {
    const rkr_procedure* ve;
    double Gv;
    double dGv;
    double Bv;
};


static double integrand_f_klein(double x, void* params)
{
    KleinParams* p =static_cast<KleinParams*>(params);

    const double diffE = p->Gv - p->ve->G(x);

    // Guard against roundoff.
    // QAGS does not normally evaluate exactly at the endpoint.
    if (diffE <= 0.0)return 0.0;

    return 1.0 / std::sqrt(diffE);
}

static double integrand_g_klein(double x, void* params)
{
    KleinParams* p =static_cast<KleinParams*>(params);

    const double diffE =p->Gv - p->ve->G(x);

    if (diffE <= 0.0)return 0.0;

    return p->ve->B(x) /std::sqrt(diffE);
}

// integrand for first fleming integral
static double integrand_f_fleming(double x, void* params) {
    KleinParams* p = static_cast<KleinParams*>(params);

    double diffE = p->Gv - p->ve->G(x);
    if (diffE <= 0.0) diffE = 0.0;               // keep real like MATLAB expects
    double denom = std::sqrt(diffE);
    if (denom == 0.0) return 0.0;

    double diff = p->dGv - p->ve->dG(x);
    return diff / denom;
}

// integrand for second fleming integral
static double integrand_g_fleming(double x, void* params) {
    KleinParams* p = static_cast<KleinParams*>(params);

    double diffE = p->Gv - p->ve->G(x);
    if (diffE <= 0.0) diffE = 0.0;
    double denom = std::sqrt(diffE);
    if (denom == 0.0) return 0.0;

    double diff = p->ve->B(x) * p->dGv - p->Bv * p->ve->dG(x);
    return diff / denom;
}



inline double rkr_procedure::fFleming(double vib, double Cu, double mu, double vmin, double abstol, double reltol, double *err, gsl_integration_workspace* w) const {

    KleinParams params{this, G(vib), dG(vib), 0.0};

    double factor = std::sqrt(Cu / mu) * (1.0 / params.dGv);

    // MATLAB term: 2*sqrt(G(v)+Y00), where MATLAB's G already has Y00.
    double argTerm = params.Gv + Y00_0;
    if (argTerm < 0.0) argTerm = 0.0;
    double term = 2.0 * std::sqrt(argTerm);

    gsl_function F;
    F.function = &integrand_f_fleming;
    F.params   = &params;

    double result = 0.0, local_err=0.0;

    //                           AbsTol=abstol, RelTol=1e-6      
    //gsl_integration_qag(&F, vmin, vib, abstol, reltol, limit,GSL_INTEG_GAUSS15, w, &result, &local_err);
    //gsl_integration_qag(&F, vmin, vib, abstol, reltol, limit,GSL_INTEG_GAUSS31, w, &result, &local_err);
    //gsl_integration_qag(&F, vmin, vib, abstol, reltol, limit,GSL_INTEG_GAUSS61, w, &result, &local_err);
    //printf("vmin=%f, vib=%f, abstol=%e, reltol=%e\n", vmin, vib, abstol, reltol);
    gsl_integration_qags(&F, vmin, vib, abstol, reltol, limit, w, &result, &local_err);
    double resabs = 0.0;
    double resasc = 0.0;   
    //gsl_integration_qk61(&F,vmin,vib,&result,&local_err,&resabs,&resasc);
    if(err) *err =  local_err;
    //if(local_err>=1e-10) std::cerr << "Warning: the absolute error approached by the fFleming integral is " << local_err << "\n";
  
    
    double final_result =  (term + result);
   
    //because local_err is the absolute error of the integral.
    double rel_err =local_err / std::abs(final_result);
    //if(rel_err>=1e-10) std::cerr << "Warning: the relative error approached by the gFleming integral is " << rel_err << "\n";

    
    return factor *final_result;
}

inline double rkr_procedure::gFleming(double vib, double Cu, double mu, double vmin, double abstol, double reltol, double *err,  gsl_integration_workspace* w) const {

    KleinParams params{this, G(vib), dG(vib), B(vib)};

    double factor = std::sqrt(mu / Cu) * (1.0 / params.dGv);

    // MATLAB term: 2*B(v)*sqrt(G(v)+Y00)
    double argTerm = params.Gv + Y00_0;
    if (argTerm < 0.0) argTerm = 0.0;
    double term = 2.0 * params.Bv * std::sqrt(argTerm);

    gsl_function F;
    F.function = &integrand_g_fleming;
    F.params   = &params;

    double result = 0.0, local_err=0.0;
                                   //abserr, relerr
    // gsl_integration_qag(&F, vmin, vib, abstol, reltol, limit,GSL_INTEG_GAUSS15, w, &result, &local_err);
    //gsl_integration_qag(&F, vmin, vib, abstol, reltol, limit,GSL_INTEG_GAUSS31, w, &result, &local_err);
    //gsl_integration_qag(&F, vmin, vib, abstol, reltol, limit,GSL_INTEG_GAUSS61, w, &result, &local_err);
    gsl_integration_qags(&F, vmin, vib, abstol, reltol, limit, w, &result, &local_err);
    double resabs = 0.0;
    double resasc = 0.0;
    //gsl_integration_qk61(&F,vmin,vib,&result,&local_err,&resabs,&resasc);
    //if(local_err>=1e-10) std::cerr << "Warning: the absolute error approached by the gFleming integral is " << local_err << "\n";
    
    double final_result =  (term + result);

    double rel_err =local_err / std::abs(final_result);
    //if(rel_err>=1e-10) std::cerr << "Warning: the relative error approached by the gFleming integral is " << rel_err << "\n";

    
    return factor *final_result;
}


inline double rkr_procedure::fKlein(double vib,double Cu,double mu,double vmin,double abstol,double reltol,double* err, gsl_integration_workspace* w) const
{
    KleinParams params{this,G(vib),0.0,0.0};

    gsl_function F;
    F.function = &integrand_f_klein;
    F.params   = &params;

    double result = 0.0;
    double local_err = 0.0;

    /*
        QAGS is used here because the ordinary Klein integral
        retains the integrable square-root singularity at x=vib.
    */
    gsl_integration_qags(&F,vmin,vib,abstol,reltol,limit,w,&result,&local_err);

    if (err) *err = local_err;

    const double factor =std::sqrt(Cu / mu);

    return factor * result;
}


inline double rkr_procedure::gKlein(double vib,double Cu,double mu,double vmin,double abstol,double reltol,double* err, gsl_integration_workspace* w) const
{
    KleinParams params{this,G(vib),0.0,0.0};

    gsl_function F;
    F.function = &integrand_g_klein;
    F.params   = &params;

    double result = 0.0;
    double local_err = 0.0;

    /*
        QAGS is used here because the ordinary Klein integral
        retains the integrable square-root singularity at x=vib.
    */
    gsl_integration_qags(&F,vmin,vib,abstol,reltol,limit,w,&result,&local_err);

    if (err) *err = local_err;

    const double factor =std::sqrt(mu / Cu);

    return factor * result;
}

inline double rkr_procedure::f(bool useFleming, double vib, double Cu, double mu, double vmin, double abstol, double reltol, double* err, gsl_integration_workspace* w) const
{
    if (useFleming) {
        return fFleming(vib, Cu, mu, vmin, abstol, reltol, err, w);
    } else {
        return fKlein(vib, Cu, mu, vmin, abstol, reltol, err, w);
    }
}

inline double rkr_procedure::g(bool useFleming, double vib, double Cu, double mu, double vmin, double abstol, double reltol, double* err, gsl_integration_workspace* w) const
{
    if (useFleming) {
        return gFleming(vib, Cu, mu, vmin, abstol, reltol, err, w);
    } else {
        return gKlein(vib, Cu, mu, vmin, abstol, reltol, err, w);
    }
}

#endif
