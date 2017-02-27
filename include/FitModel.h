/**
 *
 *    @file  FitModel.h
 *   @brief  
 *
 *    @date  02/03/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#ifndef  FIT_MODEL_H
#define  FIT_MODEL_H

#include <fwd/MassRange.h>
#include <fwd/FreeAmplitude.h>

#include "fwd/Fit.h"

#include <MassAxes.h>
#include <Model.h>

#include <iosfwd>
#include <memory>
#include <vector>

/// @defgroup Models
/// @brief Wrappers around `yap::Model` that provide constant access to its `yap::FreeAmplitude`'s.
/// @details The aim of these classes is to provide some utility function for the specific
///          \f$D^+ \to\pi^+\pi^-\pi^+\f$ decay. Also, they should protect the `yap::Model`
///          from erroneous changes during the fit process.

/// @ingroup Models
/// @brief Wrapper for a `yap::Model` with helper functions targeted to the \f$D^+ \to\pi^+\pi^-\pi^+\f$ decay.
class FitModel {
public:
    /// @brief Constructor.
    /// @param model      The YAP model of the decay.
    /// @param model_name The (optional) name of the model.
    explicit FitModel(std::unique_ptr<yap::Model> model,
                      const std::string model_name = "");

    /// _Default_ destructor.
    virtual ~FitModel();

    /// _Default_ move constructor.
    FitModel(FitModel&&) = default;
    /// _Default_ move assignment operator.
    FitModel& operator=(FitModel&&) = default;

    /// _Deleted_ copy constructor.
    FitModel(const FitModel&)            = delete;
    /// _Deleted_ copy assignment operator.
    FitModel& operator=(const FitModel&) = delete;

    /// Returns the model name.
    const std::string& name() const noexcept
    { return Name_; }

    /// Returns the YAP model.
    const std::unique_ptr<yap::Model>& model() const noexcept
    { return Model_; }

    /// Return the default mass axes.
    const yap::MassAxes& massAxes() const noexcept
    { return MassAxes_; }

    /// Returns the mass ranges.
    const std::vector<yap::MassRange>& squaredMassRanges() const noexcept
    { return SquaredMassRanges_; }

    /// Return the non-fixed free amplitudes.
    const std::vector<std::shared_ptr<const yap::FreeAmplitude>>& freeAmplitudes() const noexcept
    { return FreeAmplitudes_; }

    /// Use a common radial size [GeV^-1] for all the resonances.
    static constexpr const double radialSize() noexcept
    { return 3.; }

    /// Returns the mass of the D+.
    static constexpr const double Dmass() noexcept
    { return 1.86961; }

    /// Returns the mass of the \f$\pi^+\f$.
    static constexpr const double PiPlusMass() noexcept
    { return 1.3957018e-01; }

    /// Returns the mass of the \f$\pi^-\f$.
    static constexpr const double PiMinusMass() noexcept
    { return 1.3957018e-01; }

    /// @brief Sets the FreeAmplitude's in the model.
    /// @param p The value to which the amplitudes should be set.
    virtual void setParameters(const std::vector<std::complex<double>>& p) noexcept;

    /// @brief Sets the FreeAmplitude's in the model and make them fixed.
    /// @param p The value to which the amplitudes should be set.
    void fixParameters(const std::vector<std::complex<double>>& p) noexcept;

private:

    /// The (optional) name of the model.
    std::string Name_;

    /// The decay fit model.
    std::unique_ptr<yap::Model> Model_;

    /// The default mass axes.
    yap::MassAxes MassAxes_;

    /// The squared-mass ranges for the decay.
    std::vector<yap::MassRange> SquaredMassRanges_;

    /// @brief The non-fixed FreeAmplitude's of the model.
    /// @attention By chenging these, the values of the free
    /// amplitudes in the model will also change.
    std::vector<std::shared_ptr<const yap::FreeAmplitude>> FreeAmplitudes_;

protected:

    /// @brief Access the non-fixed free-amplitude vector (useful for sorting in derived classes).
    std::vector<std::shared_ptr<const yap::FreeAmplitude>>& accessFreeAmplitudes() noexcept
    { return FreeAmplitudes_; }
};

/// @brief Helper function to create a YAP model.
/// @tparam T The formalism to use for the model.
template <typename T>
inline std::unique_ptr<yap::Model> make_model()
{ return std::make_unique<yap::Model>(std::make_unique<T>()); }

/// Returns the index of a free amplitude in the FreeAmplitude vector of the fit model.
/// @param fa        The free amplitude to query for.
/// @param fit_model The model in which to look for the free amplitude.
/// @throw yap::exceptions::Exception If the free queried free amplitude is not in the vector.
const size_t free_amplitude_index(const std::shared_ptr<const yap::FreeAmplitude>& fa,
                                  const std::shared_ptr<const FitModel>& fit_model);

const size_t free_amplitude_index(const std::shared_ptr<const yap::FreeAmplitude>& fa,
                                  const FitModel& fit_model);
#endif
