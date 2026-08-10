#include <chrono>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <filesystem>
#include "helper.hpp"
#include "calc_wrapper.hpp"


#include <chrono>


int main(){
// ---- External data (filled by your calc) ----
std::vector<double> r;    // X values
std::vector<double> E;    // Y values

std::vector<double> r_v;  // X values
std::vector<double> E_v;  // Y values
std::vector<double> v;    // discrete v values

CalcParam cp;
AllParams p;


// Read all parameters once from the input text
if (!ReadConstantsFromText("I2B.dat", p))
{
    printf("Error reading input parameters.\n");
    return 1;
}

// Keep original parameters unchanged
const AllParams base_p = p;

std::ofstream timing_file("cpp_timing_all.txt");
timing_file << "Step_Size N_Points Timing_ms\n";

for (int i = 100; i >= 100; --i)
{
    const double step =
        static_cast<double>(i) / 1000.0;

    // Start every run with the same original parameters
    p = base_p;

    // Change ONLY the RKR step size
    p.ex.space = step;       // use actual member name if different

    // Clear previous calculation results
    r.clear();
    E.clear();

    // If these are filled/appended by calc_wrapper,
    // clear them as well
    E_v.clear();
    r_v.clear();
    v.clear();

    // -----------------------------
    // Start benchmark
    // -----------------------------
    auto start =
        std::chrono::steady_clock::now();

    calc_wrapper(
        E,
        r,
        E_v,
        r_v,
        v,
        p,
        cp
    );

    auto end =
        std::chrono::steady_clock::now();

    // -----------------------------
    // End benchmark
    // -----------------------------

    const double elapsed_ms =
        std::chrono::duration<double, std::milli>(
            end - start
        ).count();

    const std::size_t npoints = v.size();

    timing_file
        << std::fixed
        << std::setprecision(3)
        << step << " "
        << npoints << " "
        << std::setprecision(10)
        << elapsed_ms
        << '\n';

    printf(
        "DV = %.3f   N = %zu   time = %.6f ms\n",
        step,
        npoints,
        elapsed_ms
    );
}

timing_file.close();

printf(
    "Benchmark complete. Results saved to cpp_timing_all.txt\n"
);


return 0;
}