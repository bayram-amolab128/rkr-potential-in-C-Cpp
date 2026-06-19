#ifndef GUI_PARAMETERS_H
#define GUI_PARAMETERS_H

#include <string>

struct Parameters {
    std::string _name;
    std::string _state;
    double _kaiser;
    std::string _intmethod;
    double _vmax;
    double _space;
    double _ladderspace;
    double _errortol;
    double _mass1;
    double _mass2;
    double _netcharge;
    double _we;
    double _xwe;
    double _ywe;
    double _zwe;
    double _awe;
    double _bwe;
    double _cwe;
    double _dwe;
    double _ewe;
    double _fwe;
    double Be;
    double ae;
    double _ye;
    double _1e;
    double _2e;
    double _Te;
    double _De;
    double _ke;
    double _re;
};

static const Parameters params = {
    "Li2",  // name
    "X1Sg",  // state
    0.0,  // kaiser
    "Fleming",  // intmethod
    100.0,  // vmax
    0.01,  // space
    3.0,  // ladderspace
    1e-09,  // errortol
    7.016004,  // mass1
    7.016004,  // mass2
    0.0,  // netcharge
    351.43,  // we
    2.61,  // xwe
    0.00295,  // ywe
    0.0,  // zwe
    0.0,  // 1we
    0.0,  // 2we
    0.0,  // 3we
    0.0,  // 4we
    0.0,  // 5we
    0.0,  // 6we
    0.67264,  // Be
    0.00704,  // ae
    -4e-05,  // ye
    0.0,  // 1e
    0.0,  // 2e
    0.0,  // Te
    0.0,  // De
    0.0,  // ke
    0.0  // re
};

#endif // GUI_PARAMETERS_H
