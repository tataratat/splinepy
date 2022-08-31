#pragma once

#include <cstdlib>
#include <thread>

namespace splinepy::utils {

/// N-Thread execution. Queries will be splitted into chunks and each thread
/// will execute those.
template<typename Func, typename IndexT>
void NThreadExecution(const Func& f,
                      const IndexT& total,
                      IndexT nthread /* copy */
) {
  // For any negative value, std::thread::hardware_concurrency() will be taken
  // If you are not satisfied with returned value, use positive value.
  // For more info:
  //   https://en.cppreference.com/w/cpp/thread/thread/hardware_concurrency
  // Minimal value that comes out of it will be 0.
  if (nthread < 0) {
    nthread = std::thread::hardware_concurrency();
  }

  // 0 or 1, don't create a thread.
  if (nthread <= 1) {
    f(0, total);
    return;
  }

  // get chunk size and prepare threads
  // make sure it rounds up
  const IndexT chunk_size = std::div((total + nthread - 1), nthread).quot;
  std::vector<std::thread> thread_pool;
  thread_pool.reserve(nthread);

  for (int i{0}; i < (nthread - 1); i++) {
    thread_pool.emplace_back(
        std::thread{f, i * chunk_size, (i + 1) * chunk_size});
  }
  {
    // last one
    thread_pool.emplace_back(std::thread{f, (nthread - 1) * chunk_size, total});
  }

  for (auto& t : thread_pool) {
    t.join();
  }
}

} /* namespace splinepy::utils */
