/**
 *
 *    @file  parallelism.cxx
 *   @brief  
 *
 *    @date  02/04/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "parallelism.h"
#include <thread>

const unsigned available_threads() noexcept {
    const auto threads = std::thread::hardware_concurrency();
    return ( threads ? threads : 1 );
}
