/**
 *
 *    @file  RealTimeParameterPlot.cxx
 *   @brief  
 *
 *    @date  02/21/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "RealTimeParameterPlot.h"

#include "parameter_handling.h"

#include <MathUtilities.h>

#include <TApplication.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TSystem.h>

#include <algorithm>
#include <cassert>
#include <complex>
#include <iterator>
#include <memory>
#include <string>
#include <utility>
#include <vector>

TROOT gmroot("", "");
TApplication theApp("App", 0, 0);

std::unique_ptr<TGraph> amplitude_graph(const std::vector<double>& bin_low_edges,
                                        const std::vector<std::complex<double>>& ms) noexcept
{
    // Check that there are as many points in the mass shape as bins.
    assert(ms.size() == bin_low_edges.size());

    // Vector of amplitudes.
    std::vector<double> amplitudes;
    amplitudes.reserve(ms.size());
    std::transform(std::begin(ms), std::end(ms), std::back_inserter(amplitudes),
                   [](const auto& x) { return std::abs(x); });

    return std::make_unique<TGraph>(amplitudes.size(), bin_low_edges.data(), amplitudes.data());
}

std::unique_ptr<TGraph> phase_graph(const std::vector<double>& bin_low_edges,
                                    const std::vector<std::complex<double>>& ms) noexcept
{
    // Check that there are as many points in the mass shape as bins.
    assert(ms.size() == bin_low_edges.size());

    // Vector of free-amplitude phases.
    std::vector<double> phases;
    phases.reserve(ms.size());

    const auto fp = yap_to_fit_parameters(ms);
    const auto phase_shift = yap::deg(std::arg(ms.at(0))) - fp.at(1);
    std::transform(std::begin(ms), std::end(ms), std::back_inserter(phases),
                   [&](const auto& x) { return yap::deg(std::arg(x)) - phase_shift; });

    return std::make_unique<TGraph>(phases.size(), bin_low_edges.data(), phases.data());
}

// Helper function to create a TGraph whose y-values are all 0.
std::unique_ptr<TGraph> zero_graph(const std::vector<double>& bin_low_edges) noexcept
{
    const std::vector<double> zero_vector(bin_low_edges.size(), 0);
    return std::make_unique<TGraph>(zero_vector.size(), bin_low_edges.data(), zero_vector.data());
}

RealTimeParameterPlot::RealTimeParameterPlot(binned_complex_function step_mass_shape,
                                             const std::string& wave_name) noexcept :
    AmplitudeCanvasName_(wave_name + " free-amplitude amplitude"),
    AmplitudeCanvas_(std::make_unique<TCanvas>(AmplitudeCanvasName_.c_str(), AmplitudeCanvasName_.c_str())),
    FitAmplitude_(zero_graph(step_mass_shape.second)),
    GuessAmplitude_(amplitude_graph(step_mass_shape.second, step_mass_shape.first)),
    PhaseCanvasName_(wave_name + " free-amplitude phase"),
    PhaseCanvas_(std::make_unique<TCanvas>(PhaseCanvasName_.c_str(), PhaseCanvasName_.c_str())),
    FitPhase_(zero_graph(step_mass_shape.second)),
    GuessPhase_(phase_graph(step_mass_shape.second, step_mass_shape.first))
{
    assert(step_mass_shape.first.size() == step_mass_shape.second.size());

    // -- Associate plots to canvases ---
    assert(AmplitudeCanvas_);
    AmplitudeCanvas_->cd();
//    GuessAmplitude_->SetTitle("Free-amplitude amplitude"); // XXX disappears after update
    GuessAmplitude_->GetXaxis()->SetTitle("Mass [GeV/c^2]");
    GuessAmplitude_->GetYaxis()->SetTitle("Amplitude");
    GuessAmplitude_->SetLineColor(kBlue);
    GuessAmplitude_->Draw("A L*");
    FitAmplitude_->SetLineColor(kRed); 
    FitAmplitude_->SetMarkerColor(kRed);
    FitAmplitude_->Draw("SAME L*");

    assert(PhaseCanvas_);
    PhaseCanvas_->cd();
//    GuessPhase_->SetTitle("Free-amplitude phase"); // XXX disappears after update
    GuessPhase_->GetXaxis()->SetTitle("Mass [GeV/c^2]");
    GuessPhase_->GetYaxis()->SetTitle("Phase [deg]");
    GuessPhase_->SetLineColor(kBlue);
    GuessPhase_->Draw("A L*");
    FitPhase_->SetLineColor(kRed);
    FitPhase_->SetMarkerColor(kRed);
    FitPhase_->Draw("SAME L*");
    // ----------------------------------

    AmplitudeCanvas_->Modified();
    AmplitudeCanvas_->Update();

    PhaseCanvas_->Modified();
    PhaseCanvas_->Update();

    gSystem->ProcessEvents();
}

RealTimeParameterPlot::~RealTimeParameterPlot() = default;

void RealTimeParameterPlot::setCurrentParameters(const std::vector<std::complex<double>>& p) const noexcept {
    const auto x_axis     = FitAmplitude_->GetX();
    const auto n_elements = FitAmplitude_->GetN();
    const std::vector<double> mass_partition(x_axis, x_axis + n_elements);
    assert(mass_partition.size() == p.size());

    for (size_t i = 0; i < p.size(); ++ i) {
        FitAmplitude_->SetPoint(i, mass_partition[i], std::abs(p[i]));
        FitPhase_->SetPoint(i, mass_partition[i], yap::deg(std::arg(p[i])));
    }

    AmplitudeCanvas_->Modified();
    AmplitudeCanvas_->Update();

    PhaseCanvas_->Modified();
    PhaseCanvas_->Update();

    gSystem->ProcessEvents();
}
