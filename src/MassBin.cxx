#include "./MassBin.h"

#include "CachedValue.h"
#include "CalculationStatus.h"
#include "DataPartition.h"
#include "FourMomenta.h"
#include "Model.h"
#include "Parameter.h"
#include "StatusManager.h"

#include <complex>
#include <memory>

MassBin::MassBin(double lower_edge, double upper_edge) :
    yap::MassShape(),
    LowerEdge_(std::make_shared<yap::PositiveRealParameter>(lower_edge)),
    UpperEdge_(std::make_shared<yap::PositiveRealParameter>(upper_edge)),
    T_(yap::RealCachedValue::create(*this))
{
    addParameter(LowerEdge_);
    addParameter(UpperEdge_);
}

//-------------------------
void MassBin::calculate(yap::DataPartition& D,
                        const std::shared_ptr<const yap::ParticleCombination>& pc,
                        unsigned si) const
{
    auto in_range = [=](const double m) { return lowerEdge()->value() <= m && m < upperEdge()->value(); };
    for (auto& d : D) {
        const auto mass = model()->fourMomenta()->m(d, pc);
        T()->setValue(in_range(mass), d, si, D);
    }
}

//-------------------------
const std::complex<double> MassBin::value(const yap::DataPoint& d, const std::shared_ptr<const yap::ParticleCombination>& pc) const 
{ return T()->value(d, symmetrizationIndex(pc)); }

//-------------------------
void MassBin::updateCalculationStatus(yap::StatusManager& D) const
{
    if (status() == yap::VariableStatus::changed)
        D.set(*T(), yap::CalculationStatus::uncalculated);
}
