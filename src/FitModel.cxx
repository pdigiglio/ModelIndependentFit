/**
 *
 *    @file  FitModel.cxx
 *   @brief  
 *
 *    @date  02/03/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "FitModel.h"

#include "Fit.h"
#include "MassBin.h"

#include <Exceptions.h>
#include <FreeAmplitude.h>
#include <MassAxes.h>
#include <MassRange.h>
#include <Model.h>

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iterator>
#include <memory>
#include <vector>

FitModel::FitModel(std::unique_ptr<yap::Model> m,
                   const std::string model_name) :
    Name_(model_name),
    Model_(std::move(m)),
    MassAxes_(model()->massAxes()),
    SquaredMassRanges_(yap::squared(mass_range(FitModel::Dmass(), massAxes(), model()->finalStateParticles()))),
    FreeAmplitudes_()
{
    // Check if the model is valid
    if (!model())
        throw yap::exceptions::Exception("Model nullptr", "FitModel::FitModel");

    assert(model()->locked());
    // Check for consistency
    if (!model()->consistent())
        throw yap::exceptions::Exception("Model inconsistent", "FitModel::FitModel");

    // Check if the mass axes are not empty.
    if (massAxes().empty())
        throw yap::exceptions::Exception("Mass axes empty", "FitModel::FitModel");

    // Check if the mass ranges are not empty.
    if (squaredMassRanges().empty())
        throw yap::exceptions::Exception("Mass ranges empty", "FitModel::FitModel");

    // Copy the non-fixed FreeAmplitude's in the internal storage.
    const auto model_fas = yap::free_amplitudes(*model());
    
    // Select the free amplitudes.
    yap::FreeAmplitudeVector tmp_fas;
    tmp_fas.reserve(model_fas.size());
    std::copy_if(std::begin(model_fas), std::end(model_fas), std::back_inserter(tmp_fas),
                 [](const typename yap::FreeAmplitudeVector::value_type& fa)
                 { return fa->variableStatus() != yap::VariableStatus::fixed; } );

    // Make the free amplitudes constant.
    std::transform(std::begin(tmp_fas), std::end(tmp_fas), std::back_inserter(FreeAmplitudes_),
                   [](const auto& fa){ return std::static_pointer_cast<const yap::FreeAmplitude>(fa); });
}

void FitModel::setParameters(const std::vector<std::vector<std::complex<double>>>& p) noexcept {
    // Linearize the parameter matrix.
    std::vector<std::complex<double>> lp;
    lp.reserve(std::accumulate(std::begin(p), std::end(p), size_t(0),
                               [] (const size_t s, const auto& pars) { return s + pars.size(); }));
    std::for_each(std::begin(p), std::end(p), [&](const auto& pars)
                  { std::copy(std::begin(pars), std::end(pars), std::back_inserter(lp)); });

    assert(freeAmplitudes().size() == lp.size());


    // Here I'll cast the constantness of the FreeAmplitudes away!
    for (size_t i = 0; i < freeAmplitudes().size(); ++ i)
        *std::const_pointer_cast<yap::FreeAmplitude>(freeAmplitudes()[i]) = lp[i];
}

void FitModel::fixParameters(const std::vector<std::vector<std::complex<double>>>& p) noexcept {
    assert(freeAmplitudes().size() == p.size());

    setParameters(p);
    std::for_each(std::begin(freeAmplitudes()), std::end(freeAmplitudes()), [](const auto& fa)
                  { std::const_pointer_cast<yap::FreeAmplitude>(fa)->variableStatus() = yap::VariableStatus::fixed; });
}

FitModel::~FitModel() = default;

// Helper function to find the index of a FreeAmplitude in the parameter vector.
const size_t free_amplitude_index(const std::shared_ptr<const yap::FreeAmplitude>& fa,
                                  const std::vector<std::shared_ptr<const yap::FreeAmplitude>>& fas) {
    using namespace std;
    const auto it = find(begin(fas), end(fas), fa);

    if (it == end(fas))
        throw yap::exceptions::Exception("FreeAmplitude not found", "free_amplitude_idx");

    return distance(begin(fas), it);
}

const size_t free_amplitude_index(const std::shared_ptr<const yap::FreeAmplitude>& fa,
                                  const FitModel& fit_model) {
    return free_amplitude_index(fa, fit_model.freeAmplitudes());
}

const size_t free_amplitude_index(const std::shared_ptr<const yap::FreeAmplitude>& fa,
                                  const std::shared_ptr<const FitModel>& fit_model) {
    return free_amplitude_index(fa, fit_model->freeAmplitudes());
}
