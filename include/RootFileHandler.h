/**
 *
 *    @file  RootFileHandler.h
 *   @brief  
 *
 *    @date  02/04/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#ifndef  ROOT_FILE_HANDLER_H
#define  ROOT_FILE_HANDLER_H

#include "fwd/TTreeHandler.h"

#include <TFile.h>

#include <memory>
#include <string>

/// @defgroup Root
/// @brief Wrappers for ROOT classes with error handling and convenience functions.

/// @ingroup Root
/// @brief Manages the input `.root` data file.
class RootFileHandler {
public:
    /// @brief Constructor.
    /// @param file_path  The path of the input .root file.
    /// @param file_name  The input .root file name (relative to the path).
    /// @param model_name The name of the model (for fecthing the branches).
    explicit RootFileHandler(const char* file_path,
                             const char* file_name,
                             const char* model_name); 

    /// Destructor.
    ~RootFileHandler();

    /// Return the name of the input `.root` file.
    const std::string& fileName() const noexcept
    { return FileName_; }

    /// Return the path where the input `.root` file is contained.
    const std::string& path() const noexcept
    { return Path_; }

    /// Return the ROOT file.
    const std::unique_ptr<TFile>& file() const noexcept
    { return File_; }

    /// Return the model name in the ROOT file.
    const std::string& modelName() const noexcept
    { return ModelName_; }

    /// Return the Markov-chain monte-carlo data tree.
    const std::unique_ptr<TTreeHandler>& mcmcTree() const noexcept
    { return MCMCTree_; }

private:
    /// The path where the input `.root` file is contained.
    std::string Path_;
    /// The input file name.
    std::string FileName_;
    /// The model name.
    std::string ModelName_;
    /// Pointer to the input `.root` data file.
    std::unique_ptr<TFile> File_;

    /// Pointer to the Monte-Carlo data generated throguh Markov chains.
    std::unique_ptr<TTreeHandler> MCMCTree_;
};

/// @brief Return the name of the Markov-chain monte-carlo data `TTree`.
/// @param rfh The root file handled to take the model name from.
inline std::string mcmc_tree_name(const RootFileHandler& rfh)
{ return rfh.modelName() + "_mcmc"; }

#endif
