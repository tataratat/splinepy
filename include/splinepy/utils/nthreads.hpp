#pragma once

#include <cstdlib>
#include <thread>

namespace splinepy::utils {
/// N-Thread execution. Queries will be split into chunks and each thread
/// will execute those.
template<typename Func, typename IndexType>
void NThreadExecution(const Func& f,
                      const IndexType& total,
                      IndexType nthread /* copy */) {
  // For any negative value, std::thread::hardware_concurrency() will be taken
  // If you are not satisfied with returned value, use positive value.
  // For more info:
  // https://en.cppreference.com/w/cpp/thread/thread/hardware_concurrency
  // Minimal value that comes out of it will be 0.

  if (nthread < 0) {
    nthread = std::thread::hardware_concurrency();
  }

  // 0 or 1, don't create a thread.
  if (nthread <= 1) {
    f(0, total, 0);
    return;
  }

  // we don't want nthread to exceed total
  nthread = std::min(total, nthread);

  // reserve thread pool
  std::vector<std::thread> thread_pool;
  thread_pool.reserve(nthread);

  // get chunk size and prepare threads
  // make sure it rounds up
  const IndexType chunk_size = std::div((total + nthread - 1), nthread).quot;

  for (int i{}; i < (nthread - 1); i++) {
    thread_pool.emplace_back(
        std::thread{f, i * chunk_size, (i + 1) * chunk_size, i});
  }

  // last one
  thread_pool.emplace_back(
      std::thread{f, (nthread - 1) * chunk_size, total, nthread - 1});

  for (auto& t : thread_pool) {
    t.join();
  }
}
} // namespace splinepy::utils
/* namespace splinepy::utils */
