/*
MIT License

Copyright (c) 2021 Jaewook Lee

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

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
