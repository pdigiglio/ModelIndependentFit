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

/// Manages the input .root data file and provides convenience
/// functions to access its branches.
class RootFileHandler {
public:
    /// @brief Constructor.
    /// @param file_path  The path of the input .root file.
    /// @param file_name  The input .root file name (relative to the path).
    /// @param model_name The name one wants to give to the model.
    explicit RootFileHandler(const char* file_path,
                             const char* file_name,
                             const char* model_name = ""); 
    /// Destructor.
    ~RootFileHandler();

    /// Returns the model name.
    const std::string& modelName() const noexcept
    { return ModelName_; }

    /// Returns the name of the input .root file.
    const std::string& fileName() const noexcept
    { return FileName_; }

    const std::unique_ptr<TTreeHandler>& parameterTree() const noexcept
    { return ParameterTree_; }

    const std::unique_ptr<TTreeHandler>& mcmcTree() const noexcept
    { return MCMCTree_; }

    /// Return the name of the data TTree.
    const std::string mcmcTreeName() const
    { return "D3PI_" + ModelName_ + "_RESONANCE_mcmc"; }

    /// Return the name of the parameter TTree.
    const std::string parameterTreeName() const
    { return "D3PI_" + ModelName_ + "_RESONANCE_parameters"; }

    /// Return the ROOT file.
    const std::unique_ptr<TFile>& file() const noexcept
    { return File_; }

    /// Return the path where the input .root file is contained.
    const std::string& path() const noexcept
    { return Path_; }
private:
    /// The path where the input .root file is contained.
    const std::string Path_;
    /// The input file name.
    const std::string FileName_;

    /// The (optional) name of the model.
    const std::string ModelName_;

    /// Pointer to the input .root data file.
    std::unique_ptr<TFile> File_;

    /// Pointer to the Monte-Carlo data generated throguh Markov chains.
    std::unique_ptr<TTreeHandler> MCMCTree_;

    /// (?)
    std::unique_ptr<TTreeHandler> ParameterTree_;
};

#endif
