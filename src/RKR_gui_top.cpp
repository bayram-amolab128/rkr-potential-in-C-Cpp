#include <wx/wx.h>
#include <wx/spinctrl.h>
#include <wx/dcbuffer.h>   // wxAutoBufferedPaintDC
#include <wx/filepicker.h> // wxFilePickerCtrl
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

// ---- External data (filled by your calc) ----
std::vector<double> r;    // X values
std::vector<double> E;    // Y values

std::vector<double> r_v;  // X values
std::vector<double> E_v;  // Y values
std::vector<double> v;    // discrete v values

CalcParam cp;


//molecular name and state.
std::string molecule_name = "";

// ---------- Plot panel (E vs r) ----------
class PlotPanel : public wxPanel {
public:
    PlotPanel(wxWindow* parent)
        : wxPanel(parent, wxID_ANY, wxDefaultPosition, wxSize(520, -1)) {
        SetBackgroundStyle(wxBG_STYLE_PAINT);
        Bind(wxEVT_PAINT,        &PlotPanel::OnPaint,      this);
        Bind(wxEVT_LEFT_DOWN,    &PlotPanel::OnLeftDown,   this);
        Bind(wxEVT_LEFT_UP,      &PlotPanel::OnLeftUp,     this);
        Bind(wxEVT_MOTION,       &PlotPanel::OnMouseMove,  this);
        Bind(wxEVT_LEAVE_WINDOW, &PlotPanel::OnMouseLeave, this);
        Bind(wxEVT_MOUSEWHEEL,   &PlotPanel::OnMouseWheel, this);
    }

    void SetDataXY(const std::vector<double>& x, const std::vector<double>& y) {
        x_ = x; y_ = y;

        // Initialize view to the data bounds
        if (!x_.empty() && !y_.empty() && x_.size() == y_.size()) {
            auto [xminIt, xmaxIt] = std::minmax_element(x_.begin(), x_.end());
            auto [yminIt, ymaxIt] = std::minmax_element(y_.begin(), y_.end());
            double xmin = *xminIt, xmax = *xmaxIt;
            double ymin = *yminIt, ymax = *ymaxIt;
            if (xmax - xmin < 1e-15) { xmax = xmin + 1.0; }
            if (ymax - ymin < 1e-15) { ymax = ymin + 1.0; }
            xMin_ = xmin; xMax_ = xmax;
            yMin_ = ymin; yMax_ = ymax;
        }
        Refresh();
        Update();
    }

    // Optional overlay curve (e.g., Morse)
    void SetOverlayXY(const std::vector<double>& x, const std::vector<double>& y) {
        overlayX_ = x;
        overlayY_ = y;
        Refresh();
    }
    void SetOverlayVisible(bool v) {
        showOverlay_ = v;
        Refresh();
    }

private:
    std::vector<double> x_, y_;

    // overlay
    std::vector<double> overlayX_, overlayY_;
    bool showOverlay_ = false;

    // ---- View window for panning/zooming ----
    double xMin_ = 0.0, xMax_ = 1.0;
    double yMin_ = 0.0, yMax_ = 1.0;

    // Mouse-panning state
    bool    panning_ = false;
    wxPoint lastPt_;

    struct AxisTicks {
        double tick0 = 0.0;
        double step  = 1.0;
        int    n     = 0;
        double minv  = 0.0;
        double maxv  = 1.0;
    };

    static double niceNumFixed(double range, int targetSteps = 6) {
        if (range <= 0.0) return 1.0;
        double rawStep = range / targetSteps;
        double magnitude = std::pow(10.0, std::floor(std::log10(rawStep)));
        double fraction = rawStep / magnitude;

        double niceFraction;
        if      (fraction < 1.5) niceFraction = 1.0;
        else if (fraction < 3.0) niceFraction = 2.0;
        else if (fraction < 7.0) niceFraction = 5.0;
        else                     niceFraction = 10.0;

        return std::round(niceFraction * magnitude * 1000.0) / 1000.0;
    }

    static AxisTicks makeTicks(double vmin, double vmax, int targetMajor = 6) {
        AxisTicks t;
        if (vmax < vmin) std::swap(vmin, vmax);
        if (std::abs(vmax - vmin) < 1e-15) vmax = vmin + 1.0;

        const double step  = niceNumFixed(vmax - vmin, targetMajor);
        const double graphMin = std::floor(vmin / step) * step;
        const double graphMax = std::ceil (vmax / step) * step;

        t.tick0 = graphMin;
        t.step  = step;
        t.minv  = graphMin;
        t.maxv  = graphMax;
        t.n     = (int)std::floor((graphMax - graphMin) / step + 0.5) + 1;
        if (t.n < 2) t.n = 2;
        return t;
    }

    // Helpers to check if a point is inside plot rect
    bool InPlotRect(int x, int y, int L, int T, int plotW, int plotH) const {
        return x >= L && x <= L + plotW - 1 && y >= T && y <= T + plotH - 1;
    }

    // Convert pixel (screen) to data coordinates (cursor-centric zoom)
    void PixelToData(int px, int py, int L, int T, int plotW, int plotH, double& xData, double& yData) const {
        if (plotW <= 1) plotW = 2;
        if (plotH <= 1) plotH = 2;
        double tx = (px - L) / double(plotW - 1);                 // 0..1
        double ty = (py - T) / double(plotH - 1);                 // 0..1 (downward)
        tx = std::clamp(tx, 0.0, 1.0);
        ty = std::clamp(ty, 0.0, 1.0);
        xData = xMin_ + tx * (xMax_ - xMin_);
        yData = yMax_ - ty * (yMax_ - yMin_); // inverted Y
    }

    // Zoom keeping a given (px,py) fixed in screen, i.e., centered around its data coords
    void ZoomAtPixel(int px, int py, double factor) {
        const int w = GetClientSize().GetWidth();
        const int h = GetClientSize().GetHeight();
        const int L = 70, R = 20, T = 20, B = 50;
        const int plotW = std::max(1, w - L - R);
        const int plotH = std::max(1, h - T - B);
        if (plotW < 2 || plotH < 2) return;

        double cx, cy;
        PixelToData(px, py, L, T, plotW, plotH, cx, cy);

        // factor < 1 => zoom in (shrink window); factor > 1 => zoom out
        const double minRange = 1e-9; // keep sensible floor
        double xRange = std::max(xMax_ - xMin_, minRange);
        double yRange = std::max(yMax_ - yMin_, minRange);

        double newXRange = std::max(xRange * factor, minRange);
        double newYRange = std::max(yRange * factor, minRange);

        // Re-center so that cx,cy remain under cursor
        double tx = (cx - xMin_) / xRange;
        double ty = (cy - yMin_) / yRange;

        xMin_ = cx - tx * newXRange;
        xMax_ = xMin_ + newXRange;

        yMin_ = cy - ty * newYRange;
        yMax_ = yMin_ + newYRange;

        Refresh();
    }

    void OnPaint(wxPaintEvent&) {
        wxAutoBufferedPaintDC dc(this);

        // Clear entire panel background
        dc.SetBackground(*wxWHITE_BRUSH);
        dc.Clear();

        const int w = GetClientSize().GetWidth();
        const int h = GetClientSize().GetHeight();
        const int L = 70, R = 20, T = 20, B = 50;
        const int plotW = std::max(1, w - L - R);
        const int plotH = std::max(1, h - T - B);

        if (x_.empty() || y_.empty() || x_.size() != y_.size()) {
            // Draw panel frame and message (no clipping needed)
            dc.SetPen(*wxBLACK_PEN);
            dc.DrawRectangle(L, T, plotW, plotH);
            dc.DrawText("No data or size mismatch. Press Run after your calc populates r & E.",
                        L + 10, T + 10);
            return;
        }

        // If view window isn't initialized (first paint before SetDataXY), set it now.
        if (!(xMax_ > xMin_) || !(yMax_ > yMin_)) {
            auto [xminIt, xmaxIt] = std::minmax_element(x_.begin(), x_.end());
            auto [yminIt, ymaxIt] = std::minmax_element(y_.begin(), y_.end());
            double xmin = *xminIt, xmax = *xmaxIt;
            double ymin = *yminIt, ymax = *ymaxIt;
            if (xmax - xmin < 1e-15) { xmax = xmin + 1.0; }
            if (ymax - ymin < 1e-15) { ymax = ymin + 1.0; }
            xMin_ = xmin; xMax_ = xmax;
            yMin_ = ymin; yMax_ = ymax;
        }

        // Ticks based on **current view**
        AxisTicks tx = makeTicks(xMin_, xMax_);
        AxisTicks ty = makeTicks(yMin_, yMax_);

        // NOTE: No clamping in transforms (we'll cull via clip region)
        auto X = [&](double x) {
            return L + (int)std::round(((x - xMin_) / (xMax_ - xMin_)) * (plotW - 1));
        };
        auto Y = [&](double y) {
            return T + plotH - 1 - (int)std::round(((y - yMin_) / (yMax_ - yMin_)) * (plotH - 1));
        };

        // --------- CLIP to plot rectangle to prevent artifacts ----------
        dc.SetClippingRegion(L, T, plotW, plotH);

        // Fill plot background explicitly (inside clip)
        dc.SetPen(*wxTRANSPARENT_PEN);
        dc.SetBrush(*wxWHITE_BRUSH);
        dc.DrawRectangle(L, T, plotW, plotH);

        // --- Draw background grid (inside clip) ---
        wxPen gridMajorPen(wxColour(200,200,200));
        wxPen gridMinorPen(wxColour(230,230,230)); gridMinorPen.SetStyle(wxPENSTYLE_SHORT_DASH);
        const int minorDiv = 5;

        // X major + minor grids
        for (int i = 0; i < tx.n; ++i) {
            double xv = tx.tick0 + i * tx.step;
            int xpx = X(xv);
            dc.SetPen(gridMajorPen);
            dc.DrawLine(xpx, T, xpx, T + plotH);
            if (i < tx.n - 1) {
                for (int k = 1; k < minorDiv; ++k) {
                    double xm = xv + k * (tx.step / minorDiv);
                    int xpm = X(xm);
                    dc.SetPen(gridMinorPen);
                    dc.DrawLine(xpm, T, xpm, T + plotH);
                }
            }
        }

        // Y major + minor grids
        for (int j = 0; j < ty.n; ++j) {
            double yv = ty.tick0 + j * ty.step;
            int ypx = Y(yv);
            dc.SetPen(gridMajorPen);
            dc.DrawLine(L, ypx, L + plotW, ypx);
            if (j < ty.n - 1) {
                for (int k = 1; k < minorDiv; ++k) {
                    double ym = yv + k * (ty.step / minorDiv);
                    int ypm = Y(ym);
                    dc.SetPen(gridMinorPen);
                    dc.DrawLine(L, ypm, L + plotW, ypm);
                }
            }
        }

        // --- Smooth spline-connected line (inside clip) ---
        // REPLACES your scatter-only block 1:1
        dc.SetBrush(*wxBLUE_BRUSH);
        dc.SetPen(wxPen(*wxBLUE, 2));  // thicker line for visibility

        // Build wxPoint array in pixel coordinates
        std::vector<wxPoint> pts;
        pts.reserve(x_.size());
        for (size_t i = 0; i < x_.size(); ++i) {
            pts.emplace_back(X(x_[i]), Y(y_[i]));
        }

        // Draw smooth cardinal spline through all points
        if (pts.size() >= 3) {
            dc.DrawSpline((int)pts.size(), pts.data());
        } else if (pts.size() >= 2) {
            dc.DrawLines((int)pts.size(), pts.data());  // fallback for <3 points
        }

        // Optional: draw scatter dots on top of the spline
        //dc.SetBrush(*wxBLUE_BRUSH);
        //dc.SetPen(*wxTRANSPARENT_PEN);
        //for (auto& p : pts)
        //    dc.DrawCircle(p, 2);


        // --- Optional overlay curve (red line) ---
        if (showOverlay_ && overlayX_.size() == overlayY_.size() && overlayX_.size() >= 2) {
            wxPen redPen(*wxRED, 2);
            dc.SetPen(redPen);
            for (size_t i = 1; i < overlayX_.size(); ++i) {
                int x1 = X(overlayX_[i-1]), y1 = Y(overlayY_[i-1]);
                int x2 = X(overlayX_[i]),   y2 = Y(overlayY_[i]);
                dc.DrawLine(x1, y1, x2, y2);
            }
        }

        // --- Remove clip before drawing border & labels ---
        dc.DestroyClippingRegion();

        // Border on top of content
        dc.SetPen(*wxBLACK_PEN);
        dc.SetBrush(*wxTRANSPARENT_BRUSH);
        dc.DrawRectangle(L, T, plotW, plotH);

        // --- Axis ticks and labels (outside clip, so always visible) ---
        dc.SetPen(*wxBLACK_PEN);
        dc.SetTextForeground(*wxBLACK);

        for (int i = 0; i < tx.n; ++i) {
            double xv = tx.tick0 + i * tx.step;
            int xpx = X(xv);
            dc.DrawLine(xpx, T + plotH, xpx, T + plotH + 5);
            wxString lab = wxString::Format("%.3f", xv);
            wxSize ts = dc.GetTextExtent(lab);
            dc.DrawText(lab, xpx - ts.GetWidth()/2, T + plotH + 7);
        }

        for (int j = 0; j < ty.n; ++j) {
            double yv = ty.tick0 + j * ty.step;
            int ypx = Y(yv);
            dc.DrawLine(L - 5, ypx, L, ypx);
            wxString lab = wxString::Format("%.3f", yv);
            wxSize ts = dc.GetTextExtent(lab);
            dc.DrawText(lab, L - 10 - ts.GetWidth(), ypx - ts.GetHeight()/2);
        }

        // --- Axis titles ---
        dc.DrawText(wxString::FromUTF8("Intermolecular Distance (r) in Å"),
                    L + plotW/2 - 60, T + plotH + 25);

        // Rotated Y-axis label (90° CCW)
        wxString ylabel = wxString::FromUTF8("Total Energy (E) in cm⁻¹");
        wxSize ts = dc.GetTextExtent(ylabel);
        int ylabX = L - 75;
        int ylabY = T + plotH/2 + ts.GetWidth()/2;
        dc.DrawRotatedText(ylabel, ylabX, ylabY, 90);
    }

    // ---- Panning handlers ----
    void OnLeftDown(wxMouseEvent& e) {
        const int w = GetClientSize().GetWidth();
        const int h = GetClientSize().GetHeight();
        const int L = 70, R = 20, T = 20, B = 50;
        const int plotW = std::max(1, w - L - R);
        const int plotH = std::max(1, h - T - B);

        if (InPlotRect(e.GetX(), e.GetY(), L, T, plotW, plotH)) {
            panning_ = true;
            lastPt_ = e.GetPosition();
            CaptureMouse();
            SetCursor(wxCursor(wxCURSOR_HAND));
        }
        e.Skip();
    }

    void OnLeftUp(wxMouseEvent& e) {
        if (panning_) {
            panning_ = false;
            if (HasCapture()) ReleaseMouse();
            SetCursor(wxNullCursor);
            Refresh();
        }
        e.Skip();
    }

    void OnMouseLeave(wxMouseEvent& e) {
        if (panning_) {
            panning_ = false;
            if (HasCapture()) ReleaseMouse();
            SetCursor(wxNullCursor);
            Refresh();
        }
        e.Skip();
    }

    void OnMouseMove(wxMouseEvent& e) {
        if (!panning_ || !e.Dragging() || !e.LeftIsDown()) { e.Skip(); return; }

        // Panel / plot geometry
        const int w = GetClientSize().GetWidth();
        const int h = GetClientSize().GetHeight();
        const int L = 70, R = 20, T = 20, B = 50;
        const int plotW = std::max(1, w - L - R);
        const int plotH = std::max(1, h - T - B);

        // Pixel deltas
        wxPoint pt = e.GetPosition();
        int dxPix = pt.x - lastPt_.x;
        int dyPix = pt.y - lastPt_.y;
        lastPt_ = pt;

        // Convert pixel delta to data delta
        double xRange = std::max(xMax_ - xMin_, 1e-12);
        double yRange = std::max(yMax_ - yMin_, 1e-12);

        double dxData = (double)dxPix * xRange / (double)plotW;
        double dyData = (double)dyPix * yRange / (double)plotH;

        // screen y+ downward => dragging down increases y
        xMin_ -= dxData; xMax_ -= dxData;
        yMin_ += dyData; yMax_ += dyData;

        Refresh(); // redraw with new view // glitches with old frames, try fixing it later.
        e.Skip();
    }

    // ---- Wheel zoom handler ----
    void OnMouseWheel(wxMouseEvent& e) {
        int rotation = e.GetWheelRotation();
        int delta    = e.GetWheelDelta() ? e.GetWheelDelta() : 120;
        int steps    = rotation / delta;

        if (steps == 0) { e.Skip(); return; }

        const double zoomStep = 1.2; // 20% per notch
        double factor = (steps > 0) ? std::pow(1.0/zoomStep, steps)  // zoom in
                                    : std::pow(zoomStep,   -steps);  // zoom out

        ZoomAtPixel(e.GetX(), e.GetY(), factor);
        e.Skip();
    }
};



// ---------- UI ----------
struct NumField {
    wxStaticText* label{};
    wxWindow* editor{};
    bool useText{};
};

static wxSpinCtrlDouble* makeNum(wxWindow* parent, double val=0.0, double inc=0.1, int digits=6,
                                 double lo=-1e12, double hi=+1e12, const wxSize& sz=wxSize(140, -1)) {
    auto* sc = new wxSpinCtrlDouble(parent, wxID_ANY, wxEmptyString, wxDefaultPosition, sz);
    sc->SetRange(lo, hi);
    sc->SetDigits(digits);
    sc->SetIncrement(inc);
    sc->SetValue(val);
    return sc;
}

class MyFrame : public wxFrame {
public:
    MyFrame()
    : wxFrame(nullptr, wxID_ANY, "Rydberg-Klein-Rees Potential",
              wxDefaultPosition, wxSize(1200, 800))  // normal window size
    {
        aspect_ = 1.5; // width/height fixed aspect ratio

        auto* panel = new wxPanel(this);

        // Title
        auto* titleLbl = new wxStaticText(panel, wxID_ANY, "Molecular State:");
        titleCtrl = new wxTextCtrl(panel, wxID_ANY, "", wxDefaultPosition, wxSize(300, -1));

        // G(v) (left)
        auto* gvBox = new wxStaticBox(panel, wxID_ANY, wxString::FromUTF8("G(v) parameters (in cm⁻¹)"));
        auto* gvSizer = new wxStaticBoxSizer(gvBox, wxVERTICAL);
        // Make G(v) region scrollable because there are many higher-order fields
        auto* gvScroll = new wxScrolledWindow(panel, wxID_ANY, wxDefaultPosition, wxSize(-1, 220), wxVSCROLL);
        gvScroll->SetScrollRate(5, 5);
        auto* gvGrid  = new wxFlexGridSizer(0, 2, 6, 8);
        addField(gvGrid, gvScroll, wxString::FromUTF8("Y₁,₀ "),      gv.we,    351.43,    0.01, 18,
             wxString::FromUTF8("Fundamental vibration constant Y₁,₀."), true);
        addField(gvGrid, gvScroll, wxString::FromUTF8("Y₂,₀ "),   gv.xwe,   2.61,      0.001, 18,
             wxString::FromUTF8("First anharmonicity constant Y₂,₀ for Y₁,₀."), true);
        addField(gvGrid, gvScroll, wxString::FromUTF8("Y₃,₀ "),   gv.ywe,   0.00295,   0.001, 18,
             wxString::FromUTF8("Second anharmonicity constant Y₃,₀ for Y₁,₀."), true);
        addField(gvGrid, gvScroll, wxString::FromUTF8("Y₄,₀ "),   gv.zwe,   0.0,       0.001, 18,
             wxString::FromUTF8("Third anharmonicity constant Y₄,₀ for Y₁,₀."), true);
        addField(gvGrid, gvScroll, wxString::FromUTF8("Y₅,₀ "),   gv.awe,   0.0,       0.001, 18,
             wxString::FromUTF8("Higher-order Dunham coefficient Y₅,₀ for Y₁,₀."), true);
        addField(gvGrid, gvScroll, wxString::FromUTF8("Y₆,₀ "),   gv.bwe,   0.0,       0.001, 18,
             wxString::FromUTF8("Higher-order Dunham coefficient Y₆,₀ for Y₁,₀."), true);
        addField(gvGrid, gvScroll, wxString::FromUTF8("Y₇,₀ "),   gv.cwe,   0.0,       0.001, 18,
             wxString::FromUTF8("Higher-order Dunham coefficient Y₇,₀ for Y₁,₀."), true);
        addField(gvGrid, gvScroll, wxString::FromUTF8("Y₈,₀ "),   gv.dwe,   0.0,       0.001, 18,
             wxString::FromUTF8("Higher-order Dunham coefficient Y₈,₀ for Y₁,₀."), true);
        addField(gvGrid, gvScroll, wxString::FromUTF8("Y₉,₀ "),   gv.ewe,   0.0,       0.001, 18,
             wxString::FromUTF8("Higher-order Dunham coefficient Y₉,₀ for Y₁,₀."), true);
        addField(gvGrid, gvScroll, wxString::FromUTF8("Y₁₀,₀ "),   gv.fwe,   0.0,       0.001, 18,
             wxString::FromUTF8("Higher-order Dunham coefficient Y₁₀,₀ for Y₁,₀."), true);
        addField(gvGrid, gvScroll, wxString::FromUTF8("Y₁₁,₀ "),   gv.gwe,   0.0,       1e-30, 18,
             wxString::FromUTF8("Higher-order Dunham coefficient Y₁₁,₀ for Y₁,₀."), true);
        addField(gvGrid, gvScroll, wxString::FromUTF8("Y₁₂,₀ "),   gv.hwe,   0.0,       1e-30, 18,
             wxString::FromUTF8("Higher-order Dunham coefficient Y₁₂,₀ for Y₁,₀."), true);
        addField(gvGrid, gvScroll, wxString::FromUTF8("Y₁₃,₀ "),   gv.iwe,   0.0,       1e-30, 18,
             wxString::FromUTF8("Higher-order Dunham coefficient Y₁₃,₀ for Y₁,₀."), true);
        addField(gvGrid, gvScroll, wxString::FromUTF8("Y₁₄,₀ "),   gv.jwe,   0.0,       1e-30, 18,
             wxString::FromUTF8("Higher-order Dunham coefficient Y₁₄,₀ for Y₁,₀."), true);
        addField(gvGrid, gvScroll, wxString::FromUTF8("Y₁₅,₀ "),   gv.kwe,   0.0,       1e-30, 18,
             wxString::FromUTF8("Higher-order Dunham coefficient Y₁₅,₀ for Y₁,₀."), true);
        addField(gvGrid, gvScroll, wxString::FromUTF8("Y₁₆,₀ "),   gv.lwe,   0.0,       1e-30, 18,
             wxString::FromUTF8("Higher-order Dunham coefficient Y₁₆,₀ for Y₁,₀."), true);

        gvGrid->AddGrowableCol(1, 1);
        gvScroll->SetSizer(gvGrid);
        gvScroll->FitInside();
        gvSizer->Add(gvScroll, 1, wxEXPAND | wxALL, 8);

        // B(v) (left)
        auto* bvBox = new wxStaticBox(panel, wxID_ANY, wxString::FromUTF8("B(v) parameters (in cm⁻¹)"));
        auto* bvSizer = new wxStaticBoxSizer(bvBox, wxVERTICAL);
        // Make B(v) region scrollable to keep layout compact
        auto* bvScroll = new wxScrolledWindow(panel, wxID_ANY, wxDefaultPosition, wxSize(-1, 160), wxVSCROLL);
        bvScroll->SetScrollRate(5, 5);
        auto* bvGrid  = new wxFlexGridSizer(0, 2, 6, 8);
        addField(bvGrid, bvScroll, wxString::FromUTF8("Y₀,₁ " ),   bv.Be,  0.67264,  0.001, 18,
             wxString::FromUTF8("Equilibrium rotational constant Y₀,₁."), true);
        addField(bvGrid, bvScroll, wxString::FromUTF8("Y₁,₁ " ),   bv.ae,  0.00704,  0.001, 18,
             wxString::FromUTF8("Vibration-rotation interaction constant Y₁,₁."), true);
        addField(bvGrid, bvScroll, wxString::FromUTF8("Y₂,₁ " ),   bv.ye, -0.00004,  0.00001, 18,
             wxString::FromUTF8("Dunham rotational coefficient Y₂,₁."), true);
        addField(bvGrid, bvScroll, wxString::FromUTF8("Y₃,₁ " ),  bv._1e, 0.0,      0.001, 18,
             wxString::FromUTF8("First higher-order rotational coefficient Y₃,₁."), true);
        addField(bvGrid, bvScroll, wxString::FromUTF8("Y₄,₁ " ),  bv._2e, 0.0,      0.001, 18,
             wxString::FromUTF8("Second higher-order rotational coefficient Y₄,₁."), true);
           addField(bvGrid, bvScroll, wxString::FromUTF8("Y₅,₁ " ),  bv._3e, 0.0,      1e-30, 18,
               wxString::FromUTF8("Third higher-order rotational coefficient Y₅,₁."), true);
           addField(bvGrid, bvScroll, wxString::FromUTF8("Y₆,₁ " ),  bv._4e, 0.0,      1e-30, 18,
               wxString::FromUTF8("Fourth higher-order rotational coefficient Y₆,₁."), true);
           addField(bvGrid, bvScroll, wxString::FromUTF8("Y₇,₁ " ),  bv._5e, 0.0,      1e-30, 18,
               wxString::FromUTF8("Fifth higher-order rotational coefficient Y₇,₁."), true);
           addField(bvGrid, bvScroll, wxString::FromUTF8("Y₈,₁ " ),  bv._6e, 0.0,      1e-30, 18,
               wxString::FromUTF8("Sixth higher-order rotational coefficient Y₈,₁."), true);
           addField(bvGrid, bvScroll, wxString::FromUTF8("Y₉,₁ " ),  bv._7e, 0.0,      1e-30, 18,
               wxString::FromUTF8("Seventh higher-order rotational coefficient Y₉,₁."), true);
           addField(bvGrid, bvScroll, wxString::FromUTF8("Y₁₀,₁ " ),  bv._8e, 0.0,      1e-30, 18,
               wxString::FromUTF8("Eighth higher-order rotational coefficient Y₁₀,₁."), true);
           addField(bvGrid, bvScroll, wxString::FromUTF8("Y₁₁,₁ " ),  bv._9e, 0.0,      1e-30, 18,
               wxString::FromUTF8("Ninth higher-order rotational coefficient Y₁₁,₁."), true);
           addField(bvGrid, bvScroll, wxString::FromUTF8("Y₁₂,₁ " ),  bv._10e, 0.0,     1e-30, 18,
               wxString::FromUTF8("Tenth higher-order rotational coefficient Y₁₂,₁."), true);
           addField(bvGrid, bvScroll, wxString::FromUTF8("Y₁₃,₁ " ),  bv._11e, 0.0,     1e-30, 18,
               wxString::FromUTF8("Eleventh higher-order rotational coefficient Y₁₃,₁."), true);
           addField(bvGrid, bvScroll, wxString::FromUTF8("Y₁₄,₁ " ),  bv._12e, 0.0,     1e-30, 18,
               wxString::FromUTF8("Twelfth higher-order rotational coefficient Y₁₄,₁."), true);
           addField(bvGrid, bvScroll, wxString::FromUTF8("Y₁₅,₁ " ),  bv._13e, 0.0,     1e-30, 18,
               wxString::FromUTF8("Thirteenth higher-order rotational coefficient Y₁₅,₁."), true);
        bvGrid->AddGrowableCol(1, 1);
        bvScroll->SetSizer(bvGrid);
        bvScroll->FitInside();
        bvSizer->Add(bvScroll, 0, wxEXPAND | wxALL, 8);

        // Extras (left)
        auto* extraBox = new wxStaticBox(panel, wxID_ANY, "General parameters");
        auto* extraSizer = new wxStaticBoxSizer(extraBox, wxVERTICAL);
        // 4 columns = (label, field) x 2 per row
        auto* extraGrid  = new wxFlexGridSizer(0, 4, 6, 14);
        addField(extraGrid, panel, wxString::FromUTF8("m₁ (amu) "),    extras.m1,    7.016, 0.001, 8,
                 wxString::FromUTF8("Mass of atom 1 in atomic mass units."));
        addField(extraGrid, panel, wxString::FromUTF8("m₂ (amu) "),    extras.m2,    7.016, 0.001, 8,
                 wxString::FromUTF8("Mass of atom 2 in atomic mass units."));
        //addField(extraGrid, panel, "kaiser",      extras.kaiser,1.0,   0.001, 8);
        addField(extraGrid, panel, wxString::FromUTF8("Tₑ (cm⁻¹)"),          extras.Te, 0.001, 0.0001, 8,
                 wxString::FromUTF8("Electronic term energy Tₑ."));
        addField(extraGrid, panel, wxString::FromUTF8("Vmax "),        extras.Vmax,  100.0, 1.0,   0,
                 wxString::FromUTF8("Maximum vibrational quantum number Vmax."));
        addField(extraGrid, panel, wxString::FromUTF8("step size "),   extras.space, 0.001, 0.0001, 8,
                 wxString::FromUTF8("Grid step size for the calculation."));
        // Kaiser checkbox (last row in General parameters) 
        extraGrid->Add(new wxStaticText(panel, wxID_ANY, "Use Kaiser Correction"),
               0, wxALIGN_CENTER_VERTICAL);

        kaiserCheck = new wxCheckBox(panel, wxID_ANY, "");
        kaiserCheck->SetValue(false);
        extraGrid->Add(kaiserCheck, 0, wxALIGN_CENTER_VERTICAL);


        extraGrid->AddGrowableCol(3, 1);
        extraSizer->Add(extraGrid, 0, wxEXPAND | wxALL, 8);

        // Morse parameters (right beside Extras)
        auto* morseBox = new wxStaticBox(panel, wxID_ANY, "Morse parameters");
        auto* morseSizer = new wxStaticBoxSizer(morseBox, wxVERTICAL);
        auto* morseGrid  = new wxFlexGridSizer(0, 2, 6, 8);
        addField(morseGrid, panel, wxString::FromUTF8("Dₑ (cm⁻¹) "),  morse.De,  10000.0, 1.0, 8,
                 wxString::FromUTF8("Morse potential well depth Dₑ."));
        addField(morseGrid, panel, wxString::FromUTF8("kₑ (cm⁻¹ Å⁻²) "),  morse.ke,  1.0,     0.01, 8,
                 wxString::FromUTF8("Morse force constant kₑ."));
        morseGrid->AddGrowableCol(1, 1);
        morseSizer->Add(morseGrid, 0, wxEXPAND | wxALL, 8);

        // Checkbox moved here (inside Morse parameters)
        showMorse = new wxCheckBox(panel, wxID_ANY, "Show Morse potential overlay");
        showMorse->SetValue(false);
        morseSizer->Add(showMorse, 0, wxTOP | wxLEFT | wxRIGHT, 4);


        // ---- File picker + Load + Run all in one row ----
        auto* fileRow = new wxBoxSizer(wxHORIZONTAL);

        filePicker = new wxFilePickerCtrl(
            panel, wxID_ANY, wxEmptyString,
            "Choose parameter file",
            "Text files (*.txt;*.dat;*.cfg)|*.txt;*.dat;*.cfg|All files (*.*)|*.*",
            wxDefaultPosition, wxSize(360, -1)
        );

        auto* loadBtn = new wxButton(panel, wxID_ANY, "Load Params");
        auto* runBtn  = new wxButton(panel, wxID_ANY, "Run");

        fileRow->Add(filePicker, 1, wxRIGHT | wxALIGN_CENTER_VERTICAL, 8);
        fileRow->Add(loadBtn,   0, wxRIGHT | wxALIGN_CENTER_VERTICAL, 8);
        fileRow->Add(runBtn,    0, wxALIGN_CENTER_VERTICAL);


        // Only Extras here — Morse is moved above!
        auto* extraRow = new wxBoxSizer(wxVERTICAL);
        extraRow->Add(extraSizer, 0, wxEXPAND | wxALL, 8);

        // Output (left)
        output = new wxTextCtrl(panel, wxID_ANY, "",
                                wxDefaultPosition, wxSize(-1, -1),
                                wxTE_MULTILINE | wxTE_READONLY);
        output->SetFont(wxFontInfo(10).Family(wxFONTFAMILY_TELETYPE));

        // ---- Right: Plot ----
        plot = new PlotPanel(panel);
        plot->SetMinSize(wxSize(520, -1));

        // ====== Layout ======
        auto* topH = new wxBoxSizer(wxHORIZONTAL);     // main two-column layout
        auto* leftCol  = new wxBoxSizer(wxVERTICAL);   // parameters + output
        auto* rightCol = new wxBoxSizer(wxVERTICAL);   // plot

        // Left column content
        auto* titleRow = new wxBoxSizer(wxHORIZONTAL);
        titleRow->Add(titleLbl, 0, wxALIGN_CENTER_VERTICAL | wxRIGHT, 8);
        titleRow->Add(titleCtrl, 1, wxEXPAND);

        // Right parameter column = B(v) on top, Morse underneath
        auto* rightParamCol = new wxBoxSizer(wxVERTICAL);
        rightParamCol->Add(bvSizer,    0, wxEXPAND | wxBOTTOM, 8);
        rightParamCol->Add(morseSizer, 0, wxEXPAND);
                
        // G(v) on left, (B(v)+Morse) stacked on right
        auto* twoParamCols = new wxBoxSizer(wxHORIZONTAL);
        twoParamCols->Add(gvSizer,      1, wxRIGHT | wxEXPAND, 6);
        twoParamCols->Add(rightParamCol,1, wxLEFT  | wxEXPAND, 6);



        leftCol->Add(titleRow,     0, wxALL | wxEXPAND, 12);
        leftCol->Add(twoParamCols, 0, wxLEFT | wxRIGHT | wxBOTTOM | wxEXPAND, 12);
        leftCol->Add(extraRow,     0, wxLEFT | wxRIGHT | wxBOTTOM | wxEXPAND, 12);

        leftCol->Add(new wxStaticText(panel, wxID_ANY, "Load parameters from file:"), 0, wxLEFT | wxRIGHT, 8);
        leftCol->Add(fileRow, 0, wxEXPAND | wxALL, 8);


        leftCol->Add(new wxStaticText(panel, wxID_ANY, "Output:"), 0, wxLEFT | wxRIGHT, 12);
        leftCol->Add(output,       1, wxALL | wxEXPAND, 12);   // stretch to use remaining height

        // Right column title
        auto* titleText = new wxStaticText(
            panel, wxID_ANY, "Energy Levels vs. Intermolecular Distance");
        titleText->SetFont(wxFont(16, wxFONTFAMILY_SWISS, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_BOLD));
        titleText->SetForegroundColour(*wxBLACK);

        rightCol->Add(titleText, 0, wxALIGN_CENTER_HORIZONTAL | wxTOP | wxLEFT | wxRIGHT, 12);
        rightCol->Add(plot, 1, wxLEFT | wxRIGHT | wxBOTTOM | wxEXPAND, 20);

        // Put columns side by side with 3:6 weighting (~30% : ~60%)
        topH->Add(leftCol,  3, wxEXPAND);   // proportion 3
        topH->Add(rightCol, 6, wxEXPAND);   // proportion 6

        panel->SetSizerAndFit(topH);

        // Events
        runBtn->Bind(wxEVT_BUTTON,  &MyFrame::OnRun,        this);
        loadBtn->Bind(wxEVT_BUTTON, &MyFrame::OnLoadParams, this);
        showMorse->Bind(wxEVT_CHECKBOX, &MyFrame::OnToggleMorse, this);

        // keep a reasonable min size and bind resize to maintain aspect ratio
        SetMinSize(wxSize(900, 600));
        Bind(wxEVT_SIZE, &MyFrame::OnSizedKeepAspect, this);
    }

private:
    struct { NumField we,xwe,ywe,zwe,awe,bwe,cwe,dwe,ewe,fwe, gwe, hwe, iwe, jwe, kwe, lwe; } gv;
    struct { NumField Be,ae,ye,_1e,_2e,_3e,_4e,_5e,_6e,_7e,_8e,_9e,_10e,_11e,_12e,_13e; } bv;
    struct { NumField m1,m2,space,Vmax, Te; } extras;
    struct { NumField De, ke; } morse;

    wxTextCtrl* titleCtrl{};
    wxTextCtrl* output{};
    wxFilePickerCtrl* filePicker{};
    PlotPanel*  plot{};
    wxCheckBox* showMorse{};
    wxCheckBox* kaiserCheck{};


    double aspect_{1.5}; // width/height

    // cached latest parameters (for re-building overlay on toggle)
    AllParams lastP_{};

    void addField(wxFlexGridSizer* grid, wxWindow* parent, const wxString& name,
              NumField& slot, double defVal, double inc, int digits,
              const wxString& hint = {}, bool useText = false) {
        slot.label = new wxStaticText(parent, wxID_ANY, name + ":");
        slot.useText = useText;
        if (useText) {
            auto initial = wxString::Format("%.*g", digits, defVal);
            auto* txt = new wxTextCtrl(parent, wxID_ANY, initial, wxDefaultPosition, wxSize(140, -1));
            slot.editor = txt;
        } else {
            slot.editor = makeNum(parent, defVal, inc, digits);
        }

        if (!hint.empty()) {
            slot.label->SetToolTip(hint);
            if (slot.useText) {
                static_cast<wxTextCtrl*>(slot.editor)->SetToolTip(hint);
            } else {
                static_cast<wxSpinCtrlDouble*>(slot.editor)->SetToolTip(hint);
            }
        }

        grid->Add(slot.label, 0, wxALIGN_CENTER_VERTICAL);
        grid->Add(slot.editor,  1, wxEXPAND);
    }

    // safer getter
    double val(const NumField& f, double def = 0.0) const {
        if (!f.editor) return def;
        if (f.useText) {
            auto* txt = static_cast<wxTextCtrl*>(f.editor);
            try {
                return std::stod(std::string(txt->GetValue().mb_str()));
            } catch (...) {
                return def;
            }
        }
        return static_cast<wxSpinCtrlDouble*>(f.editor)->GetValue();
    }

    // Maintain fixed aspect ratio during resize
    void OnSizedKeepAspect(wxSizeEvent& e) {
        static bool inHandler = false;
        if (inHandler) { e.Skip(); return; }
        inHandler = true;

        wxSize sz = GetSize();
        int w = sz.GetWidth();
        int h = sz.GetHeight();

        int targetH = static_cast<int>(std::round(w / aspect_));
        int targetW = static_cast<int>(std::round(h * aspect_));

        if (std::abs(targetH - h) < std::abs(targetW - w)) {
            SetSize(w, targetH);
        } else {
            SetSize(targetW, h);
        }

        inHandler = false;
        e.Skip();
    }

    

    void OnRun(wxCommandEvent&) {

        //create output file
        try {
            std::filesystem::create_directories("./output");
        } catch (const std::exception& e) {
            output->AppendText(wxString::Format("Failed to create output folder: %s\n", e.what()));
        }
        AllParams p;
        //CalcParam cp; --globally initialized.
        //p.name = std::string()
        p.state_name = std::string(titleCtrl->GetValue().mb_str());
    

        // G(v) params
        p.gv.we  = val(gv.we);   p.gv.xwe = val(gv.xwe); p.gv.ywe = val(gv.ywe);
        p.gv.zwe = val(gv.zwe);  p.gv.awe = val(gv.awe); p.gv.bwe = val(gv.bwe);
        p.gv.cwe = val(gv.cwe);  p.gv.dwe = val(gv.dwe); p.gv.ewe = val(gv.ewe);
        p.gv.fwe = val(gv.fwe);  p.gv.gwe = val(gv.gwe); p.gv.hwe = val(gv.hwe);
        p.gv.iwe = val(gv.iwe);  p.gv.jwe = val(gv.jwe); p.gv.kwe = val(gv.kwe); 
        p.gv.lwe = val(gv.lwe);

        // B(v) params
        p.bv.Be = val(bv.Be);      p.bv.ae = val(bv.ae);
        p.bv.ye = val(bv.ye);      p.bv._1e = val(bv._1e);
        p.bv._2e  = val(bv._2e);   p.bv._3e  = val(bv._3e);   p.bv._4e  = val(bv._4e);
        p.bv._5e  = val(bv._5e);   p.bv._6e  = val(bv._6e);   p.bv._7e  = val(bv._7e);
        p.bv._8e  = val(bv._8e);   p.bv._9e  = val(bv._9e);   p.bv._10e = val(bv._10e);
        p.bv._11e = val(bv._11e);  p.bv._12e = val(bv._12e);  p.bv._13e = val(bv._13e);

        // Extras
        p.ex.m1 = val(extras.m1);
        p.ex.m2 = val(extras.m2);
        p.ex.space  = val(extras.space);
        p.ex.UseKaiser = kaiserCheck->IsChecked() ? 1.0 : 0.0;

        p.ex.Vmax   = val(extras.Vmax);
        p.ex.Te     = val(extras.Te);
        //output->AppendText(wxString::Format("Value of Te %f!\n", p.ex.Te));

        // Morse (overlay) params
        p.ex.De = val(morse.De, 0.0);
        p.ex.ke = val(morse.ke, 0.0);

        std::ofstream timing_file("cpp_timing_all.txt");

        // header
        timing_file << "Step_Size N_Points Timing_ms\n";


        std::vector<double> steps;

        std::ifstream file("steps.txt");
        printf("Here is 1\n");

        if (!file.is_open())
            //throw std::runtime_error("Could not open file: " + "steps.txt");

       
         printf("Here is 2\n");
        for(int i = 0; i < 100; i++)
        {
            double step;
            file >> step;

            if (file.fail())
                break;

            steps.push_back(step);
            printf("Step %d: %f\n", i, step);
        }

        for(int i = 0; i < 100; i++)
        {
            
            //printf("Step %d: %f\n", i, steps[i]);
        }

        gsl_set_error_handler_off();        
        for(int i = 0; i < 1; i++) {  
            // Read all parameters once from the input text
        if (!ReadConstantsFromText("I2B.dat", p))
        {
            printf("Error reading input parameters.\n");
            return;
        }
        //printf("Here %d\n", i);
        // Keep original parameters unchanged
        const AllParams base_p = p;  
        // Calls heavy calc (fills E & r)
        output->AppendText(wxString::Format("Performing turning point calculations... wait!\n"));
        r.clear(); E.clear();

        p.ex.space  = steps[i];
        //auto start = std::chrono::high_resolution_clock::now();
        auto start = std::chrono::steady_clock::now();
        
        calc_wrapper(E, r, E_v, r_v, v, p, cp);   //calculates E and r values;
        

        //auto end   = std::chrono::high_resolution_clock::now();
        auto end   = std::chrono::steady_clock::now();

            // time in milliseconds
        double elapsed_ms =
        std::chrono::duration<double, std::milli>(
            end - start
        ).count();
         // number of calculated points
        std::size_t npoints = static_cast<std::size_t>(E.size() / 2.0);

        // write:
        // step_size   npoints   time_ms
        timing_file<<std::fixed << std::setprecision(6)
        << steps[i] << " "
        << npoints << " "
        << elapsed_ms
        << "\n";

        
        output->AppendText(wxString::Format("Calculation Complete, elapsed time = %f seconds.\n", elapsed_ms / 1000.0));
        }
        timing_file.close();

        //save in file for further analysis
        printf("=============\n"); 
        printf(p.state_name.c_str());
        printf("=============\n"); 

        try {
            std::filesystem::create_directories("./output");
        } catch (const std::exception& e) {
            output->AppendText(wxString::Format("Failed to create output folder: %s\n", e.what()));
        }
        WriteRGToFile("./output/Evsr.dat",r, E );
        // Plot from external vectors
        plot->SetDataXY(r, E);
        output->AppendText(wxString::Format("Plot is updated and Evsr data is saved in ./output/Evsr.dat.\n"));

        //auto [v_n, E_n] = calc_discrete(ve, p.ex.Vmax, p.ex.Te, dissociationEnergy_rel);

        // print to console
        //for (size_t i = 0; i < r_v.size(); ++i) {
           // output->AppendText(wxString::Format("v=%.1f  E=%.3f 1/cm\n", r_v[i], E_v[i]));
        //}

        WriteRGToFile("./output/Evsr_discrete.dat", r_v, E_v);
        output->AppendText(wxString::Format("Evsr_discrete.dat output file is updated with discrete energy levels.\n"));

    }

    void OnToggleMorse(wxCommandEvent&) {
        BuildOrClearMorseOverlay();
    }

    void BuildOrClearMorseOverlay() {
        if (!showMorse->GetValue()) {
            plot->SetOverlayVisible(false);
            plot->SetOverlayXY({}, {});
            return;
        }

        // Need valid parameters (De>0, Ke>0) to build a Morse curve
        MorseParams mp{ lastP_.ex.Te, lastP_.ex.De, lastP_.ex.Re, lastP_.ex.ke };
        if (!(mp.De > 0.0 && mp.ke > 0.0)) {
            output->AppendText("Morse overlay disabled: require De>0 and ke>0.\n");
            plot->SetOverlayVisible(false);
            plot->SetOverlayXY({}, {});
            return;
        }

        if (r.empty()) {
            plot->SetOverlayVisible(false);
            plot->SetOverlayXY({}, {});
            return;
        }

        auto [xminIt, xmaxIt] = std::minmax_element(r.begin(), r.end());
        double rmin = *xminIt, rmax = *xmaxIt;
        if (rmax <= rmin) { rmax = rmin + 1.0; }

        // Make a smooth 1000-point curve
        std::vector<double> rM, VM;
        const double N = 1000.0;
        const double dr = (rmax - rmin) / (N - 1.0);
        makeMorseCurve(mp, rmin, rmax, dr, rM, VM);

        plot->SetOverlayXY(rM, VM);
        plot->SetOverlayVisible(true);
    }

    // load parameters from file and update controls
    void OnLoadParams(wxCommandEvent&) {
        wxString path = filePicker->GetPath();
        if (path.empty()) {
            output->AppendText("No file selected.\n");
            return;
        }
        int updated = LoadParamsFromFile(path);
        output->AppendText(wxString::Format("Loaded %d parameter(s) from: %s\n", updated, path));
    }

    int LoadParamsFromFile(const wxString& path) {
        std::ifstream in(path.ToStdString());
        if (!in) {
            output->AppendText("Failed to open file.\n");
            return 0;
        }

        // Map parameter names -> editor widgets
        std::unordered_map<std::string, wxWindow*> m;

        // G(v)
        m["we"]  = gv.we.editor;   m["xwe"] = gv.xwe.editor;  m["ywe"] = gv.ywe.editor;
        m["zwe"] = gv.zwe.editor;  m["awe"] = gv.awe.editor;  m["bwe"] = gv.bwe.editor;
        m["cwe"] = gv.cwe.editor;  m["dwe"] = gv.dwe.editor;  m["ewe"] = gv.ewe.editor;
        m["fwe"] = gv.fwe.editor;  m["gwe"] = gv.gwe.editor;  m["hwe"] = gv.hwe.editor;
        m["iwe"] = gv.iwe.editor;  m["jwe"] = gv.jwe.editor;  m["kwe"] = gv.kwe.editor;
        m["lwe"] = gv.lwe.editor;

        // B(v)
        m["Be"]  = bv.Be.editor;   m["ae"]  = bv.ae.editor;   m["ye"]  = bv.ye.editor;
        m["_1e"] = bv._1e.editor;  m["_2e"] = bv._2e.editor;
        m["_3e"] = bv._3e.editor;  m["_4e"] = bv._4e.editor;  m["_5e"] = bv._5e.editor;
        m["_6e"] = bv._6e.editor;  m["_7e"] = bv._7e.editor;  m["_8e"] = bv._8e.editor;
        m["_9e"] = bv._9e.editor;  m["_10e"] = bv._10e.editor; m["_11e"] = bv._11e.editor;
        m["_12e"] = bv._12e.editor; m["_13e"] = bv._13e.editor;

        // Extras
        m["m1"]    = extras.m1.editor;
        m["m2"]    = extras.m2.editor;
        m["Te"]    = extras.Te.editor;
        ///m["kaiser"]= extras.kaiser.editor;
        m["Vmax"]  = extras.Vmax.editor;
        m["space"] = extras.space.editor;

        // Morse
        m["De"]    = morse.De.editor;
        m["ke"]    = morse.ke.editor;

        int count = 0;
        std::string line;

        // Clear all numeric inputs first so parameters omitted from the file do not keep old values.
        for (auto& kv : m) {
            if (!kv.second) continue;
            if (auto* sc = wxDynamicCast(kv.second, wxSpinCtrlDouble)) {
                sc->SetValue(0.0);
            } else if (auto* txt = wxDynamicCast(kv.second, wxTextCtrl)) {
                txt->SetValue("0");
            }
        }
        if (kaiserCheck) {
            kaiserCheck->SetValue(false);
        }

        auto trim = [](std::string& s){
            size_t a = s.find_first_not_of(" \t\r\n");
            size_t b = s.find_last_not_of(" \t\r\n");
            if (a == std::string::npos) { s.clear(); return; }
            s = s.substr(a, b - a + 1);
        };

        while (std::getline(in, line)) {
            if (auto pos = line.find('#'); pos != std::string::npos) line = line.substr(0, pos);
            trim(line);
            if (line.empty()) continue;

            auto eq = line.find('=');
            if (eq == std::string::npos) continue;

            std::string key = line.substr(0, eq);
            std::string val = line.substr(eq + 1);
            trim(key); trim(val);
            if (key.empty() || val.empty()) continue;

            if (key == "state") {
                titleCtrl->SetValue(val);

                ++count;
                continue;
            }

            if (key == "name") {
                printf(val.c_str());
                printf("\n");
                //titleCtrl->SetValue(val);
                ++count;
                continue;
            }


            if (key == "kaiser") {
            try {
                double d = std::stod(val);
                kaiserCheck->SetValue(d != 0.0);
                ++count;
            } catch (...) {
                output->AppendText(wxString::Format("Bad value for %s: %s\n", key, val));
            }
            continue;
            }


            auto it = m.find(key);
            if (it != m.end() && it->second) {
                try {
                    double d = std::stod(val);
                    if (auto* sc = wxDynamicCast(it->second, wxSpinCtrlDouble)) {
                        sc->SetValue(d);
                    } else if (auto* txt = wxDynamicCast(it->second, wxTextCtrl)) {
                        txt->SetValue(val);
                    }
                    ++count;
                } catch (...) {
                    output->AppendText(wxString::Format("Bad value for %s: %s\n", key, val));
                }
            } else {
                output->AppendText(wxString::Format("Unknown key: %s\n", key));
            }

        }
        return count;
    }

    void print_in_disp(){
        //print out some calculated internal values.
        output->AppendText(wxString::Format("Equillibrium Internuclear Distance: %f\n", cp.re));
    }
};

class MyApp : public wxApp {
public:
    bool OnInit() override {
        auto* frame = new MyFrame();
        frame->Show(true);
        return true;
    }
};
wxIMPLEMENT_APP(MyApp);
