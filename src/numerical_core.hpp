#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>

struct GvParams {
    double we{}, xwe{}, ywe{}, zwe{}, awe{}, bwe{}, cwe{}, dwe{}, ewe{}, fwe{};
};

struct BvParams {
    double Be{}, ae{}, ye{}, _1e{}, _2e{};
};

struct ExtraParams {
    double m1{}, m2{}, space{},ladderspace{}, Te{}, De{}, ke{}, Re{}, Vmax{}, UseKaiser{};  // space can be an integer-like count, kept double for UI consistency
};

struct AllParams {
    std::string state_name;
    GvParams gv;
    BvParams bv;
    ExtraParams ex;
};

struct MorseParams {
    double Te;   // electronic term origin (cm⁻¹)
    double De;   // dissociation energy (cm⁻¹)
    double Re;   // equilibrium bond length (Å)
    double ke;   // force constant (cm⁻¹ / Å²)

    // compute shape parameter a = √(Ke / (2 De))
    double a() const;
};

// single-point Morse potential
double morseV(double r, const MorseParams& m);

// uniform grid version: fills r_out and V_out
void makeMorseCurve(const MorseParams& m,
                    double rmin, double rmax, double dr,
                    std::vector<double>& r_out,
                    std::vector<double>& V_out);

// evaluate Morse potential on an arbitrary r-grid
void makeMorseOnGrid(const MorseParams& m,
                     const std::vector<double>& r_in,
                     std::vector<double>& V_out);

void CubicSplineFit(const std::vector<double>& r_in,
                    const std::vector<double>& E_in,
                    double dr,
                    std::vector<double>& r_out,
                    std::vector<double>& E_out);
   

                    
void clean_up_data(const std::vector<double>& r_raw,
                    const std::vector<double>& E_raw,
                    std::vector<double>& r_clean,
                    std::vector<double>& E_clean);

                    
void LinearResample(const std::vector<double>& r_in,
                    const std::vector<double>& E_in,
                    double dr,
                    std::vector<double>& r_out,
                    std::vector<double>& E_out);

void calc(std::vector<double>& E, std::vector<double>& r,std::vector<double>& E_v, std::vector<double>& r_v, const AllParams& p);
void SortPaired(std::vector<double>& E, std::vector<double>& r);
void PlotXY(std::vector<double>&E,std::vector<double>&r);
void WriteRGToFile(const std::string& filename,
                   const std::vector<double>& r,
                   const std::vector<double>& G);


                   
// Template helper implementation
template<typename... Vecs>
std::vector<std::vector<double>> packColumns(const Vecs&... cols) {
    return { cols... };
}

// Main function implementation
template<typename... Vecs>
void WriteColumnsToFile(const std::string& filename, const Vecs&... cols)
{
    // Pack the arguments into a list of columns
    std::vector<std::vector<double>> columns = packColumns(cols...);

    if (columns.empty()) {
        throw std::runtime_error("WriteColumnsToFile: no columns provided.");
    }

    // Check all column sizes match
    const size_t N = columns[0].size();
    for (size_t i = 1; i < columns.size(); ++i) {
        if (columns[i].size() != N) {
            throw std::runtime_error("WriteColumnsToFile: column size mismatch.");
        }
    }

    std::ofstream fout(filename);
    if (!fout) {
        throw std::runtime_error("WriteColumnsToFile: cannot open file '" + filename + "'");
    }

    fout << std::scientific << std::setprecision(8);

    // Write each row
    for (size_t row = 0; row < N; ++row) {
        for (size_t col = 0; col < columns.size(); ++col) {
            fout << std::setw(15) << columns[col][row];
            if (col + 1 < columns.size()) fout << "  ";
        }
        fout << "\n";
    }

    fout.close();
}

