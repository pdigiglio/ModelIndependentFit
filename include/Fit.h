// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__MIFIT__H
#define __BAT__MIFIT__H

#include <BAT/BCModel.h>

#include "fwd/FitIntegrator.h"
#include "fwd/FitModel.h"
#include "fwd/RootFitData.h"

#include <fwd/FreeAmplitude.h>

#include <memory>
#include <string>

// This is a Fit header file.
// Model source code is located in file miFit/Fit.cxx

// ---------------------------------------------------------
class Fit final : public BCModel
{

public:

    /// @brief Constructor.
    /// @param fit_model  The decay model.
    /// @param fit_data   The data to fit the model to.
    /// @param integrator The model integrator.
    explicit Fit(std::shared_ptr<FitModel> fit_model,
                 std::unique_ptr<RootFitData> fit_data,
                 std::unique_ptr<FitIntegrator> integrator);

//    /// @brief Constructor.
//    explicit Fit(const std::string& name,
//                   std::unique_ptr<yap::Model> model,
//                   const unsigned integration_points,
//                   const unsigned integration_threads,
//                   const unsigned integration_batch_size) ;
//    //
//    /// Function for generating new points for integration.
//    using Generator = std::function<std::vector<yap::FourVector<double> >()>;
//    Generator& integrationPointGenerator() noexcept
//    { return IntegrationPointGenerator_; }
//
//    const unsigned numberOfIntegrationPoints() const noexcept
//    { return NumberOfIntegrationPoints_; }
//
//    const unsigned integrationBatchSize() const noexcept
//    { return IntegrationBatchSize_; }
//
//    const unsigned numberOfIntegrationThreads() const noexcept
//    { return NumberOfIntegrationThreads_; }
//
//    /// _Default_ destructor.
//    ~Fit() final = default;

    /// _Default_ destructor.
    ~Fit();

    // Overload LogLikelihood to implement model
    double LogLikelihood(const std::vector<double>& pars) final;

    // Overload LogAprioriProbability if not using built-in 1D priors
    // double LogAPrioriProbability(const std::vector<double> & pars);

    /// @brief Calculate real and imaginary parts of the FreeAmplitude's.
    /// @param pars The vector of the FreeAmplitude's being fitted.
    void CalculateObservables(const std::vector<double> & pars) final;

    /// @brief Fix a compoment of the FreeAmplitude _fa_.
    /// @param fa    The FreeAmplitude whose component one wants to fix.
    /// @param j     The FreeAmplitude component.
    /// @param value The value to fix the FreeAmplitude's component to.
    void fixComponent(const std::shared_ptr<yap::FreeAmplitude>& fa, const unsigned char j, const double value);

//    void setComponentPrior(const std::shared_ptr<yap::FreeAmplitude>& fa, const unsigned char j, BCPrior* prior);
//    void setPriors(const std::shared_ptr<yap::FreeAmplitude>& fa, BCPrior* first_component_prior, BCPrior* second_component_prior);

//    /// Returns the name of the pdf file where the plots are drawn.
//    const std::string pdfPlotsOutputName() const noexcept
//    { return OutputPath_ + GetSafeName() + "_plots.pdf"; }
//
    /// Returns the name of the log file.
    const std::string logFileName() const noexcept;
//    { return OutputPath_ + GetSafeName() + "_log.txt"; }

//    /// Access the model.
//    const std::shared_ptr<yap::Model>& model() const noexcept
//    { return Model_; }

//    yap::DataSet& fitData() noexcept
//    { return FitData_; }

//    const yap::MassAxes& axes() const noexcept
//    { return Axes_; }

//    yap::DataPartitionVector& fitPartitions() noexcept
//    { return FitPartitions_; }


    /// @brief Sets the ranges of a FreeAmplitude.
    /// @param fa        The FreeAmplitude whose range to fix.
    /// @param first_low   The low boundary of the first component.
    /// @param first_high  The high boundary of the first component.
    /// @param second_low  The low boundary of the second component.
    /// @param second_high The high boundary of the second component.
    void setRanges(const std::shared_ptr<yap::FreeAmplitude>& fa, double first_low, double first_high, double second_low, double second_high);

//    const std::vector<double>& massPartition() const noexcept
//    { return MassPartition_; }

private:

    /// The fit model for the decay.
    const std::shared_ptr<FitModel> FitModel_;

    /// The data to fit the model to.
    const std::unique_ptr<RootFitData> Data_;

    /// The model integrator.
    const std::unique_ptr<FitIntegrator> Integrator_;

//    /// Directory where the fitting output is located.
//    const std::string OutputPath_ = "output/";
//
//    /// The YAP model to fix the data.
//    const std::shared_ptr<yap::Model> Model_;
//
////    /// @brief Free amplitudes of the model.
////    /// @attention By changing these, the values of the free amplitudes
////    /// in the model will also change.
////    yap::FreeAmplitudeVector FreeAmplitudes_;
//
//    /// Data to fit the model to.
//    yap::DataSet FitData_;
//
//    /// Partitioning of the fit data.
//    yap::DataPartitionVector FitPartitions_;
//
////    /// Ranges for the masses.
////    yap::MassAxes Axes_;
//
////    /// Point generator for the integration.
////    Generator IntegrationPointGenerator_;
////
////    /// The number of integration points (?)
////    unsigned NumberOfIntegrationPoints_;
////    /// (?)
////    unsigned IntegrationBatchSize_;
////    /// (?)
////    unsigned NumberOfIntegrationThreads_;
//
////    /// Integral result.
////    yap::ModelIntegral       Integral_;
////    yap::DataSet             IntegralData_;
////    yap::DataPartitionVector IntegralDataPartitions_;
//
//    const std::vector<double> MassPartition_;
//
//    std::vector<double> OldParameters_;

    /// @brief Sets the free amplitudes to the values of _p_.
    /// @details This is the function where one can interpret the parameters
    /// as they want through the mathematical function one uses to assign them
    /// to the free amplitudes.
    /// @param p The new value of the free amplitudes.
    void setParameters(const std::vector<double>& p);

//    /// @brief Evaluates the integral of the model.
//    void integrate();
};
// ---------------------------------------------------------

///// Loads input data to fit.
//size_t load_data(yap::DataSet& data, const yap::Model& M, const yap::MassAxes& A,
//                 double initial_mass, TTree& t_mcmc, int N, int lag);

#endif
