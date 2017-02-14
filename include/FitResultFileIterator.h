/**
 *
 *    @file  FitResultFileIterator.h
 *   @brief  
 *
 *    @date  02/13/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#ifndef  FIT_RESULT_FILE_ITERATOR_H
#define  FIT_RESULT_FILE_ITERATOR_H

#include "fwd/FittedFreeAmplitude.h"

#include <iosfwd>
#include <iterator>

/// Iterates over the fit-result file.
class FitResultFileIterator : public std::iterator<std::input_iterator_tag, std::string> {
public:

    using char_type    = typename std::string::value_type;
    using traits_type  = typename std::string::traits_type;
    /// The input file type.
    using istream_type = std::basic_istream<char_type, traits_type>;

private:

    /// Default constructor (private).
    explicit FitResultFileIterator() noexcept : InputFile_(nullptr) {}

public:

    /// @brief Constructor.
    /// @param input_file The (already opened) input file.
    explicit FitResultFileIterator(istream_type& input_file) :
        InputFile_(&input_file)
    { ++(*this); }

    /// Pre-increment operator (read line in).
    FitResultFileIterator& operator++();

    /// Post-increment operator (read line in).
    FitResultFileIterator operator++(int);

    /// Deference operator.
    FittedFreeAmplitude operator*() const noexcept;

    /// Arrow operator.
    FittedFreeAmplitude operator->() const noexcept;

    /// Returns the an iterator to the past-to-end line.
    static const FitResultFileIterator& end() noexcept
    { static const FitResultFileIterator END; return END; }

    /// Equality operator.
    friend const bool operator==(const FitResultFileIterator& lhs, const FitResultFileIterator& rhs) noexcept
    { return lhs.InputFile_ == rhs.InputFile_; }

private:

    /// Pointer to the input stream.
    istream_type* InputFile_;

    /// Currently loaded line.
    std::string CurrentLine_;
};

/// Unequality operator.
inline const bool operator!=(const FitResultFileIterator& lhs, const FitResultFileIterator& rhs) noexcept
{ return !(lhs == rhs); }

/// Counts the entries of the input fit-result file.
size_t count_lines(typename FitResultFileIterator::istream_type& input_file) noexcept;

#endif
