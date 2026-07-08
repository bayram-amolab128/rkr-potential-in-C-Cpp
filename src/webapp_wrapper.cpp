#include <string>
#include <sstream>
#include <vector>
#include <exception>

#include <emscripten/bind.h>

#include "calc_wrapper.hpp"
#include "helper.hpp"

std::string ExportPotentialDAT(const std::vector<double>& r,
                               const std::vector<double>& V)
{
    std::ostringstream out;
    out << "# r(Angstrom)    V(cm^-1)\n";

    for (size_t i = 0; i < r.size(); ++i)
        out << r[i] << "    " << V[i] << "\n";

    return out.str();
}

std::string RunRKRFromText(const std::string& inputText)
{
    try
    {
        AllParams p;

        if (!ReadConstantsFromText(inputText, p))
            return "ERROR: Failed to parse input parameters.\n";

        std::vector<double> r;
        std::vector<double> V;

        std::vector<double> r_v;
        std::vector<double> V_v;
        std::vector<double> v;

        CalcParam cp;

        calc_wrapper(V, r, V_v, r_v, v, p, cp);
        calc_wrapper(V, r, V_v, r_v, v, p, cp);

        std::cout << "r size = " << r.size() << "\n";
        std::cout << "V size = " << V.size() << "\n";

        return ExportPotentialDAT(r, V);
    }
    catch (const std::exception& e)
    {
        return std::string("ERROR: ") + e.what() + "\n";
    }
}

EMSCRIPTEN_BINDINGS(rkr_module)
{
    emscripten::function("RunRKRFromText", &RunRKRFromText);
}