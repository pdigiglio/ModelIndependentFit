/**
 *
 *    @file  FitResultFileIterator.cxx
 *   @brief  
 *
 *    @date  02/13/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "FitResultFileIterator.h"
#include "FittedFreeAmplitude.h"

#include <numeric>
#include <sstream>
#include <string>

// Helper function to check if the line is a comment.
const bool is_comment(const std::string& line) noexcept {
    // Check if the first non-space element is a '#'.
    for (const auto& c : line)
        if (!isspace(c))
            return c == '#';

    // If no return has happened, the line is empty and should
    // be ignored as if it were a comment.
    return true;
}

FitResultFileIterator& FitResultFileIterator::operator++() {
    // Read the input lines till they're comments.
    do {
        // Check if the iterator is at the end of file.
        if (InputFile_ and !getline(*InputFile_, CurrentLine_))
            return *this = end();
    } while (is_comment(CurrentLine_));

    return *this;
}

FitResultFileIterator FitResultFileIterator::operator++(int) {
    // Store the current value.
    const auto previous = *this;

    // Read line in.
    ++(*this);

    // Return the old value.
    return previous;
}

FittedFreeAmplitude FitResultFileIterator::operator*() const noexcept {
    std::istringstream iss(CurrentLine_);

    // Bin low edge.
    double ble = 0.;
    iss >> ble;

    // Amplitude.
    double a = 0.;
    iss >> a;

    // Amplitude uncertainty.
    double au = 0.;
    iss >> au;

    // Phase difference.
    double pd = 0.;
    iss >> pd;

    // Phase-difference uncertainty.
    double pdu = 0.;
    iss >> pdu;

    // Cumulative phase.
    double cp = 0.;
    iss >> cp;

    // Cumulative-phase uncertainty.
    double cpu = 0.;
    iss >> cpu;

    return FittedFreeAmplitude(ble, a, au, pd, pdu, cp, cpu);
}

FittedFreeAmplitude FitResultFileIterator::operator->() const noexcept
{ return this->operator*(); }

size_t count_lines(typename FitResultFileIterator::istream_type& input_file) noexcept {
    return std::accumulate(FitResultFileIterator(input_file), FitResultFileIterator::end(), 0,
                           [](size_t a, const auto& b) { return ++a; });
}
