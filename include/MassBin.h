/// \file

#ifndef  yap_MassBin_h
#define  yap_MassBin_h

#include "fwd/CachedValue.h"
#include "fwd/DataPartition.h"
#include "fwd/DataPoint.h"
#include "fwd/ParticleCombination.h"
#include "fwd/StatusManager.h"

#include "MassShape.h"
#include "Parameter.h"

#include <memory>

/// Class to model a mass bin \f$m_{\text{min}} \le m < m_{\text{max}}\f$.
class MassBin : public yap::MassShape
{
public:
    /// \brief Constructor.
    /// \param lower_edge The lower value of the mass.
    /// \param upper_edge The upper value of the mass.
    explicit MassBin(double lower_edge, double upper_edge);

    const std::shared_ptr<yap::PositiveRealParameter>& lowerEdge() const noexcept
    { return LowerEdge_; }

    const std::shared_ptr<yap::PositiveRealParameter>& upperEdge() const noexcept
    { return UpperEdge_; }

    using MassShape::calculate;

    /// update the calculationStatus for a DataPartition
    virtual void updateCalculationStatus(yap::StatusManager& D) const override;
    
    /// \return value for DataPoint and ParticleCombination
    /// \param d DataPoint
    /// \param pc shared_ptr to ParticleCombination
    virtual const std::complex<double> value(const yap::DataPoint& d, const std::shared_ptr<const yap::ParticleCombination>& pc) const override;

    virtual void calculate(yap::DataPartition& D,
                           const std::shared_ptr<const yap::ParticleCombination>& pc,
                           unsigned si) const override;
private:
    /// Lower edge of the bin.
    const std::shared_ptr<yap::PositiveRealParameter> LowerEdge_;
    /// Upper edge of the bin.
    const std::shared_ptr<yap::PositiveRealParameter> UpperEdge_;

    /// Access the cached dynamic amplitude.
    const std::shared_ptr<yap::RealCachedValue>& T() const noexcept
    { return T_; }

    /// Cached dynamic amplitude.
    const std::shared_ptr<yap::RealCachedValue> T_;
};

#endif
