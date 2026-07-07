#ifndef RKR_PROCEDURE_H
#define RKR_PROCEDURE_H


#include <cmath>
#include <utility>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <algorithm>
#include <gsl/gsl_integration.h>


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

    double fFleming(double vib, double Cu, double mu, double vmin, double abstol = 1e-6, double reltol = 1e-6, double *err = nullptr) const;
    double gFleming(double vib, double Cu, double mu, double vmin, double abstol = 1e-6, double reltol = 1e-6, double *err =nullptr) const;

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
        std::cout << "Kaiser correction for Y00: " << y00 << std::endl;
        return -0.5 - y00 / Y10;
    }
    return -0.5;
}


// ---------------- G(v), dG, d2G ----------------
inline double rkr_procedure::G(double vib) const {
    double x = vib + 0.5;
    return we * std::pow(x, 1) + xwe * std::pow(x, 2) + ywe * std::pow(x, 3)
         + zwe * std::pow(x, 4) + awe * std::pow(x, 5) + bwe * std::pow(x, 6)
         + cwe * std::pow(x, 7) + dwe * std::pow(x, 8) + ewe * std::pow(x, 9)
         + fwe * std::pow(x,10) + gwe * std::pow(x,11) + hwe * std::pow(x,12)
         + iwe * std::pow(x,13) + jwe * std::pow(x,14) + kwe * std::pow(x,15)
         + lwe * std::pow(x,16) + Y00_0;
}

inline double rkr_procedure::dG(double vib) const {
    double x = vib + 0.5;
    return we
         + 2.0 * xwe * std::pow(x, 1)
         + 3.0 * ywe * std::pow(x, 2)
         + 4.0 * zwe * std::pow(x, 3)
         + 5.0 * awe * std::pow(x, 4)
         + 6.0 * bwe * std::pow(x, 5)
         + 7.0 * cwe * std::pow(x, 6)
         + 8.0 * dwe * std::pow(x, 7)
         + 9.0 * ewe * std::pow(x, 8)
         +10.0 * fwe * std::pow(x, 9)
         +11.0 * gwe * std::pow(x,10)
         +12.0 * hwe * std::pow(x,11)
         +13.0 * iwe * std::pow(x,12)
         +14.0 * jwe * std::pow(x,13)
         +15.0 * kwe * std::pow(x,14)
         +16.0 * lwe * std::pow(x,15);
}

inline double rkr_procedure::d2G(double vib) const {
    double x = vib + 0.5;
    return +2.0 * xwe
         + 6.0 * ywe * std::pow(x, 1)
         +12.0 * zwe * std::pow(x, 2)
         +20.0 * awe * std::pow(x, 3)
         +30.0 * bwe * std::pow(x, 4)
         +42.0 * cwe * std::pow(x, 5)
         +56.0 * dwe * std::pow(x, 6)
         +72.0 * ewe * std::pow(x, 7)
         +90.0 * fwe * std::pow(x, 8)
         +110.0 * gwe * std::pow(x, 9)
         +132.0 * hwe * std::pow(x,10)
         +156.0 * iwe * std::pow(x,11)
         +182.0 * jwe * std::pow(x,12)
         +210.0 * kwe * std::pow(x,13)
         +240.0 * lwe * std::pow(x,14);
}


// ---------------- B(v), dB ----------------
inline double rkr_procedure::B(double vib) const {
    double x = vib + 0.5;
    return be + ae * x + ye * std::pow(x,2) + _1e * std::pow(x,3) + _2e * std::pow(x,4) + _3e * std::pow(x,5) + _4e * std::pow(x,6) + _5e * std::pow(x,7) + _6e * std::pow(x,8)+ _7e * std::pow(x,9) + _8e * std::pow(x,10) + _9e * std::pow(x,11) + _10e * std::pow(x,12) + _11e * std::pow(x,13) + _12e * std::pow(x,14) +_13we * std::pow(x,15);
}

inline double rkr_procedure::dB(double vib) const {
    double x = vib + 0.5;
    return +ae
           + 2.0 * ye  * std::pow(x,1)
           + 3.0 * _1e * std::pow(x,2)
           + 4.0 * _2e * std::pow(x,3)
           + 5.0 * _3e * std::pow(x,4)
           + 6.0 * _4e * std::pow(x,5)
           + 7.0 * _5e * std::pow(x,6)
           + 8.0 * _6e * std::pow(x,7)
           + 9.0 * _7e * std::pow(x,8)
           + 10.0 * _8e * std::pow(x,9)
           + 11.0 * _9e * std::pow(x,10)
           + 12.0 * _10e * std::pow(x,11)
           + 13.0 * _11e * std::pow(x,12)
           + 14.0 * _12e * std::pow(x,13)
           + 15.0 * _13we * std::pow(x,14);
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

    std::cerr << "Warning: no extremum (dG=0) found up to v=" << vmax << "\n";
    return vmax;
}


// ---------------- Fleming integrals ----------------
struct KleinParams {
    const rkr_procedure* ve;
    double Gv;
    double dGv;
    double Bv;
};

// integrand for first fleming integral
static double integrand_f(double x, void* params) {
    KleinParams* p = static_cast<KleinParams*>(params);

    double diffE = p->Gv - p->ve->G(x);
    if (diffE <= 0.0) diffE = 0.0;               // keep real like MATLAB expects
    double denom = std::sqrt(diffE);
    if (denom == 0.0) return 0.0;

    double diff = p->dGv - p->ve->dG(x);
    return diff / denom;
}

// integrand for second fleming integral
static double integrand_g(double x, void* params) {
    KleinParams* p = static_cast<KleinParams*>(params);

    double diffE = p->Gv - p->ve->G(x);
    if (diffE <= 0.0) diffE = 0.0;
    double denom = std::sqrt(diffE);
    if (denom == 0.0) return 0.0;

    double diff = p->ve->B(x) * p->dGv - p->Bv * p->ve->dG(x);
    return diff / denom;
}

inline double rkr_procedure::fFleming(double vib, double Cu, double mu, double vmin, double abstol, double reltol, double *err) const {
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(100000);

    KleinParams params{this, G(vib), dG(vib), 0.0};

    double factor = std::sqrt(Cu / mu) * (1.0 / params.dGv);

    // MATLAB term: 2*sqrt(G(v)+Y00), where MATLAB's G already has Y00.
    double argTerm = params.Gv + Y00_0;
    if (argTerm < 0.0) argTerm = 0.0;
    double term = 2.0 * std::sqrt(argTerm);

    gsl_function F;
    F.function = &integrand_f;
    F.params   = &params;

    double result = 0.0, local_err=0.0;

    //                           AbsTol=abstol, RelTol=1e-6      
    gsl_integration_qag(&F, vmin, vib, abstol, reltol, 100000,GSL_INTEG_GAUSS61, w, &result, &local_err);
    if(err) *err =  local_err;

    gsl_integration_workspace_free(w);
    return factor * (term + result);
}

inline double rkr_procedure::gFleming(double vib, double Cu, double mu, double vmin, double abstol, double reltol, double *err) const {
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(100000);

    KleinParams params{this, G(vib), dG(vib), B(vib)};

    double factor = std::sqrt(mu / Cu) * (1.0 / params.dGv);

    // MATLAB term: 2*B(v)*sqrt(G(v)+Y00)
    double argTerm = params.Gv + Y00_0;
    if (argTerm < 0.0) argTerm = 0.0;
    double term = 2.0 * params.Bv * std::sqrt(argTerm);

    gsl_function F;
    F.function = &integrand_g;
    F.params   = &params;

    double result = 0.0, local_err=0.0;
                                   //abserr, relerr
    gsl_integration_qag(&F, vmin, vib, abstol, reltol, 100000, GSL_INTEG_GAUSS61, w, &result, &local_err);
    //gsl_integration_qags(&F, vmin, vib, abstol, reltol, 100000, w, &result, &local_err);
    if(err) *err = local_err;
    //printf("The absolute error approached by the gFleming integral is %g\n", local_err);

    gsl_integration_workspace_free(w);
    return factor * (term + result);
}

#endif
