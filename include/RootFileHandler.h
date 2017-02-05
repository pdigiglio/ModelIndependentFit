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

#include <TFile.h>
#include <TTree.h>

#include <memory>
#include <string>

/// Manages the input .root data file and provides convenience
/// functions to access its branches.
class RootFileHandler {
public:
    /// @brief Constructor.
    /// @param file_name  The name of the input .root file.
    /// @param model_name The name one wants to give to the model.
    explicit RootFileHandler(const char* file_name, const char* model_name = "");

    /// Returns the model name.
    const std::string& modelName() const noexcept
    { return ModelName_; }

    /// Returns the name of the input .root file.
    const std::string& rootFileName() const noexcept
    { return RootFileName_; }

    const std::string fitName() const noexcept
    { return modelName() + "_fit"; }

    const std::unique_ptr<TTree>& parameterTree() const noexcept
    { return ParameterTree_; }

    const std::unique_ptr<TTree>& mcmcTree() const noexcept
    { return MCMCTree_; }

    /// Return the name of the data TTree.
    const std::string mcmcTreeName() const
    { return "D3PI_" + ModelName_ + "_RESONANCE_mcmc"; }

    /// Return the name of the parameter TTree.
    const std::string parameterTreeName() const
    { return "D3PI_" + ModelName_ + "_RESONANCE_parameters"; }

private:
    const std::string ModelName_;

    /// The input file name.
    const std::string RootFileName_;

    /// Pointer to the input .root data file.
    std::unique_ptr<TFile> RootFile_;

    /// Pointer to the Monte-Carlo data generated throguh Markov chains.
    std::unique_ptr<TTree> MCMCTree_;

    /// (?)
    std::unique_ptr<TTree> ParameterTree_;
};

#endif
