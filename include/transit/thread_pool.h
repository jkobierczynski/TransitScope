// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <condition_variable>
#include <deque>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

namespace transit {

/// Minimal fixed-size thread pool; tasks run in submission order (FIFO).
class ThreadPool {
public:
    explicit ThreadPool(unsigned n) {
        if (n == 0) n = std::max(1u, std::thread::hardware_concurrency());
        for (unsigned i = 0; i < n; ++i) workers_.emplace_back([this] { loop(); });
    }
    ~ThreadPool() {
        { std::lock_guard<std::mutex> l(m_); stop_ = true; }
        cv_.notify_all();
        for (auto& w : workers_) w.join();
    }
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    unsigned size() const { return static_cast<unsigned>(workers_.size()); }

    template <class F>
    auto submit(F&& f) -> std::future<decltype(f())> {
        using R = decltype(f());
        auto task = std::make_shared<std::packaged_task<R()>>(std::forward<F>(f));
        auto fut = task->get_future();
        { std::lock_guard<std::mutex> l(m_); q_.emplace_back([task] { (*task)(); }); }
        cv_.notify_one();
        return fut;
    }

private:
    void loop() {
        for (;;) {
            std::function<void()> job;
            {
                std::unique_lock<std::mutex> l(m_);
                cv_.wait(l, [this] { return stop_ || !q_.empty(); });
                if (stop_ && q_.empty()) return;
                job = std::move(q_.front());
                q_.pop_front();
            }
            job();
        }
    }
    std::vector<std::thread> workers_;
    std::deque<std::function<void()>> q_;
    std::mutex m_;
    std::condition_variable cv_;
    bool stop_ = false;
};

}  // namespace transit
