#ifndef HELPER_HPP
#define HELPER_HPP

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <cstddef>
#include <unordered_map>
#include <functional>
#include <numeric>
#include <cctype>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <utility>
#include <limits>


//these spectroscopic constants are actually Dunham coefficients.. should rename it to DunhamCoeffs later.
struct GvParams
{
    double we{0.0};
    double xwe{0.0};
    double ywe{0.0};
    double zwe{0.0};
    double awe{0.0};
    double bwe{0.0};
    double cwe{0.0};
    double dwe{0.0};
    double ewe{0.0};
    double fwe{0.0};
    double gwe{0.0};
    double hwe{0.0};
    double iwe{0.0};
    double jwe{0.0};
    double kwe{0.0};
    double lwe{0.0};
};

struct BvParams
{
    double Be{0.0};
    double ae{0.0};
    double ye{0.0};
    double _1e{0.0};
    double _2e{0.0};
    double _3e{0.0};
    double _4e{0.0};
    double _5e{0.0};
    double _6e{0.0};
    double _7e{0.0};
    double _8e{0.0};
    double _9e{0.0};
    double _10e{0.0};
    double _11e{0.0};
    double _12e{0.0};
    double _13e{0.0};
};

struct ExtraParams {
    double m1{}, m2{}, space{},ladderspace{}, Te{}, De{}, ke{}, Re{}, Vmax{}, UseKaiser{}, NetCharge{};
};

struct GvParams_nde {
    double De{};   // dissociation energy, cm^-1
    double vD{};   // vibrational index at dissociation
    double A{};    // leading NDE coefficient
    double n{6.0}; // long-range power, e.g. 6 for C6/r^6
};

struct BvParams_nde {
    double De{};      // Dissociation energy
    double B0{};      // Leading rotational constant
    double vD{};      // Vibrational index at dissociation
    double A{};       // Leading G(v) coefficient
    double n{6.0};    // Long-range power (6 for C6/r^6)
};
struct switch_params {
    bool use_morse_fit{false};
    bool use_spline_fit{true};
    bool use_rkr_fit{true};
    double vs{0.0};  // switching vibrational index
    double ds{1.0};  // switching width
};

struct AllParams {
    std::string state_name;
    GvParams gv;
    BvParams bv;
    ExtraParams ex;
    GvParams_nde gv_nde;
    BvParams_nde bv_nde;
    switch_params sw;
};

struct MorseParams {
    double Te;   // electronic term origin (cm⁻¹)
    double De;   // dissociation energy (cm⁻¹)
    double Re;   // equilibrium bond length (Å)
    double ke;   // force constant (cm⁻¹ / Å²)

    // compute shape parameter a = √(Ke / (2 De))
    double a() const;
};

//for morse potential
inline double MorseParams::a() const {
    if (De > 0.0 && ke > 0.0)
        return std::sqrt(ke / (2.0 * De));
    return 0.0;
}

//file reader to collect spectroscopic constants.
inline std::string trim(std::string s)
{
    auto not_space = [](unsigned char c) { return !std::isspace(c); };

    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());

    return s;
}

inline bool ReadConstantsFromStream(std::istream& in, AllParams& p)
{
    std::string line;

    while (std::getline(in, line))
    {

        // Remove comments
        auto commentPos = line.find('#');
        if (commentPos != std::string::npos)
            line = line.substr(0, commentPos);

        commentPos = line.find('%');
        if (commentPos != std::string::npos)
            line = line.substr(0, commentPos);

        line = trim(line);

        if (line.empty())
            continue;

        auto eqPos = line.find('=');
        if (eqPos == std::string::npos)
            continue;

        std::string key = trim(line.substr(0, eqPos));
        std::string val = trim(line.substr(eqPos + 1));

        try
        {
            if (key == "name") {
               // p.name = val;
            }
            else if (key == "state" || key == "state_name") {
               // p.state_name = val;
            }
            else if (key == "intmethod") {
               // p.intmethod = val;
            }

            // G(v)
            else if (key == "Y1,0")  p.gv.we  = std::stod(val);
            else if (key == "Y2,0") p.gv.xwe = std::stod(val);
            else if (key == "Y3,0") p.gv.ywe = std::stod(val);
            else if (key == "Y4,0") p.gv.zwe = std::stod(val);
            else if (key == "Y5,0") p.gv.awe = std::stod(val);
            else if (key == "Y6,0") p.gv.bwe = std::stod(val);
            else if (key == "Y7,0") p.gv.cwe = std::stod(val);
            else if (key == "Y8,0") p.gv.dwe = std::stod(val);
            else if (key == "Y9,0") p.gv.ewe = std::stod(val);
            else if (key == "Y10,0") p.gv.fwe = std::stod(val);
            else if (key == "Y11,0") p.gv.gwe = std::stod(val);
            else if (key == "Y12,0") p.gv.hwe = std::stod(val);
            else if (key == "Y13,0") p.gv.iwe = std::stod(val);
            else if (key == "Y14,0") p.gv.jwe = std::stod(val);
            else if (key == "Y15,0") p.gv.kwe = std::stod(val);
            else if (key == "Y16,0") p.gv.lwe = std::stod(val);

            // B(v)
            else if (key == "Y0,1")  p.bv.Be  = std::stod(val);
            else if (key == "Y1,1")  p.bv.ae  = std::stod(val);
            else if (key == "Y2,1")  p.bv.ye  = std::stod(val);
            else if (key == "Y3,1") p.bv._1e = std::stod(val);
            else if (key == "Y4,1") p.bv._2e = std::stod(val);
            else if (key == "Y5,1") p.bv._3e = std::stod(val);
            else if (key == "Y6,1") p.bv._4e = std::stod(val);
            else if (key == "Y7,1") p.bv._5e = std::stod(val);
            else if (key == "Y8,1") p.bv._6e = std::stod(val);
            else if (key == "Y9,1") p.bv._7e = std::stod(val);
            else if (key == "Y10,1") p.bv._8e = std::stod(val);
            else if (key == "Y11,1") p.bv._9e = std::stod(val);
            else if (key == "Y12,1") p.bv._10e = std::stod(val);
            else if (key == "Y13,1") p.bv._11e = std::stod(val);
            else if (key == "Y14,1") p.bv._12e = std::stod(val);
            else if (key == "Y15,1") p.bv._13e = std::stod(val);

            // Extra
            else if (key == "m1") p.ex.m1 = std::stod(val);
            else if (key == "m2") p.ex.m2 = std::stod(val);
            //else if (key == "netcharge") p.ex.netcharge = std::stod(val);
            else if (key == "space") p.ex.space = std::stod(val);
            else if (key == "net charge") p.ex.NetCharge = std::stod(val);
            else if (key == "ladderspace") p.ex.ladderspace = std::stod(val);
            else if (key == "Vmax") p.ex.Vmax = std::stod(val);
            //else if (key == "errortol") p.ex.errortol = std::stod(val);
            else if (key == "kaiser" || key == "UseKaiser") {
                //std::cout << "UseKaiser: " << val << std::endl;
                p.ex.UseKaiser = std::stod(val);
            }
            else if (key == "Te") p.ex.Te = std::stod(val);
            else if (key == "De") p.ex.De = std::stod(val);
            else if (key == "Ke" || key == "ke") p.ex.ke = std::stod(val);
            else if (key == "re" || key == "Re") p.ex.Re = std::stod(val);

            else {
                //std::cerr << "Unknown key skipped: [" << key << "] = [" << val << "]\n";
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Failed to parse line: " << line << '\n';
            std::cerr << "Reason: " << e.what() << '\n';
            return false;
        }
    }

    return true;
}

inline bool ReadConstantsFromText(const std::string& text, AllParams& p)
{
    std::istringstream in(text);
    return ReadConstantsFromStream(in, p);
}

inline bool ReadConstantsFromFile(const std::string& filename, AllParams& p)
{
    std::ifstream in(filename);

    if (!in)
    {
        //std::cerr << "Could not open file: " << filename << '\n';
        return false;
    }

    return ReadConstantsFromStream(in, p);
}

// Reads paired arrays (r, G) from a text file into vectors.
inline void ReadRGFromFile(const std::string& filename,
                           std::vector<double>& r,
                           std::vector<double>& G)
{
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        //throw std::runtime_error("ReadRGFromFile: cannot open file '" + filename + "'");
    }

    r.clear();
    G.clear();

    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        double x = 0.0;
        double y = 0.0;

        if (!(iss >> x >> y)) {
            continue; // skip malformed lines
        }

        r.push_back(x);
        G.push_back(y);
    }

    fin.close();
}

// Writes paired arrays (r, G) to a text file.
inline void WriteRGToFile(const std::string& filename,
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

    fout << std::scientific << std::setprecision(16);
    for (size_t i = 0; i < r.size(); ++i) {
        fout << std::setw(15) << r[i]
             << "  "
             << std::setw(15) << G[i]
             << '\n';
    }

    fout.close();
}



inline void SortPaired(std::vector<double>& E, std::vector<double>& r)
{
    if (r.size() != E.size()) {
        throw std::runtime_error("SortPaired: r and E must have the same size.");
    }

    if (r.empty()) return;

    std::vector<size_t> idx(r.size());
    for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;

    std::sort(idx.begin(), idx.end(),
              [&](size_t i1, size_t i2) {
                  return r[i1] < r[i2];
              });

    std::vector<double> r_sorted, E_sorted;

    for (size_t i = 0; i < idx.size(); ++i) {
        r_sorted.push_back(r[idx[i]]);
        E_sorted.push_back(E[idx[i]]);
    }

    const double drTiny   = 1e-8;
    const double dELarge  = 1e-2;
    const double slopeMax = 1e6;

    std::vector<double> r_clean, E_clean;

    r_clean.push_back(r_sorted[0]);
    E_clean.push_back(E_sorted[0]);

    for (size_t i = 1; i < r_sorted.size(); ++i) {
        double dr = r_sorted[i] - r_clean.back();
        double dE = E_sorted[i] - E_clean.back();

        if (!std::isfinite(dr) || !std::isfinite(dE) ||
            !std::isfinite(r_sorted[i]) || !std::isfinite(E_sorted[i])) {
            continue;
        }

        if (std::abs(dr) < drTiny && std::abs(dE) > dELarge) {
            continue;
        }

        if (std::abs(dr) > 0.0) {
            double slope = std::abs(dE / dr);
            if (slope > slopeMax) {
                continue;
            }
        }

        r_clean.push_back(r_sorted[i]);
        E_clean.push_back(E_sorted[i]);
    }

    r.swap(r_clean);
    E.swap(E_clean);
}


//just detects index of given value in V.
inline size_t findClosestIndex(const std::vector<double>& V, double value)
{
    auto it = std::min_element(
        V.begin(), V.end(),
        [value](double a, double b)
        {
            return std::abs(a - value) < std::abs(b - value);
        });

    return std::distance(V.begin(), it);
}


inline std::vector<size_t> findMinimumIndices(const std::vector<double>& V,
                                              const std::vector<double>& r)
{
    if (V.size() != r.size())
        throw std::runtime_error("findMinimumIndices: V and r must have the same size.");

    if (V.empty())
        throw std::runtime_error("findMinimumIndices: empty vectors.");

    double minVal = *std::min_element(V.begin(), V.end());

    std::vector<size_t> indices;

    for (size_t i = 0; i < V.size(); ++i)
    {
        if (V[i] == minVal)
            indices.push_back(i);
    }

    return indices;
}


template<typename T>
inline T average(const std::vector<T>& v)
{
    if (v.empty())
        throw std::runtime_error("average: empty vector.");

    return std::accumulate(v.begin(), v.end(), T{}) /
           static_cast<T>(v.size());
}

//returns the i_ex_outer.
inline size_t splitandfindClosestIndex(const std::vector<double>& V,   size_t split_idx, double V_ex){

    std::vector<double> V_outer(V.begin() + split_idx, V.end());


    return findClosestIndex(V_outer,V_ex);
}



#endif