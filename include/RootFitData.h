/**
 *
 *    @file  RootFitData.h
 *   @brief  
 *
 *    @date  02/04/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#ifndef  ROOT_FIT_DATA_H
#define  ROOT_FIT_DATA_H

#include "fwd/FitModel.h"
#include "fwd/RootFileHandler.h"

#include <DataPartition.h>
#include <DataSet.h>

#include <memory>

/// Interfaces YAP with the input ROOT data file.
class RootFitData {
public:
    /// @brief Constructor.
    /// @param input_rfh The input ROOT-file wrapper.
    /// @param fit_model The model to fit to the data.
    RootFitData(std::unique_ptr<RootFileHandler> rfh,
                std::shared_ptr<const FitModel> fit_model);

    /// Returns the fit data.
    const yap::DataSet& data() const noexcept
    { return Data_; }

    /// Returns the fit data partitions.
    const yap::DataPartitionVector& partitions() const noexcept
    { return Partitions_; }

    /// Access the fit data partitions.
    yap::DataPartitionVector& partitions() noexcept
    { return Partitions_; }

    /// Returns the path where the input data file is contained.
    const std::string path() const noexcept;

private:

    /// The input ROOT-file wrapper.
    const std::unique_ptr<RootFileHandler> RootFile_;

    /// The fit model.
    const std::shared_ptr<const FitModel> FitModel_;

    /// The data to fit the model to.
    yap::DataSet Data_;

    /// The partitioning of the fit data.
    yap::DataPartitionVector Partitions_;
};

#endif
