#ifndef PROGRESS_HPP
#define PROGRESS_HPP

// ===========================================================================
//  ProgressReporter — one progress line at a fixed interval.
// ===========================================================================
//
//  A long step (reading a large input, finalizing many bins) writes nothing
//  for a long time, so a working task can look stuck. A background thread
//  writes one line every interval:
//
//      <prefix>progress: <text from the current line function>, <h:mm:ss> elapsed
//
//  The interval is 600 s. The environment variable KMER_GWAS_PROGRESS_SEC
//  changes it; 0 turns the lines off. The reporter only reads counters that
//  the caller updates, so it cannot change the output.
// ===========================================================================

#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>

namespace kmer {

class ProgressReporter {
public:
    // line(elapsed_seconds) returns the text after "progress: ".
    using LineFn = std::function<std::string(double)>;

    ProgressReporter(std::ostream& out, std::string prefix, LineFn line)
        : out_(out), prefix_(std::move(prefix)), line_(std::move(line)),
          start_(std::chrono::steady_clock::now())
    {
        long sec = 600;
        if (const char* e = std::getenv("KMER_GWAS_PROGRESS_SEC")) sec = std::strtol(e, nullptr, 10);
        if (sec > 0) {
            interval_ = std::chrono::seconds(sec);
            thread_   = std::thread([this] { run(); });
        }
    }

    ~ProgressReporter() { stop(); }
    ProgressReporter(const ProgressReporter&) = delete;
    ProgressReporter& operator=(const ProgressReporter&) = delete;

    // Changes the line for the next step.
    void set_line(LineFn line) {
        std::lock_guard<std::mutex> lk(m_);
        line_ = std::move(line);
    }

    void stop() {
        {
            std::lock_guard<std::mutex> lk(m_);
            stop_ = true;
        }
        cv_.notify_all();
        if (thread_.joinable()) thread_.join();
    }

    // 1234 -> "1234", 120300000 -> "120.3 M", 3400000000 -> "3.4 G".
    static std::string count_text(double n) {
        char buf[32];
        if (n >= 1e9)      std::snprintf(buf, sizeof buf, "%.1f G", n / 1e9);
        else if (n >= 1e6) std::snprintf(buf, sizeof buf, "%.1f M", n / 1e6);
        else               std::snprintf(buf, sizeof buf, "%.0f", n);
        return buf;
    }

private:
    void run() {
        std::unique_lock<std::mutex> lk(m_);
        while (!cv_.wait_for(lk, interval_, [this] { return stop_; })) {
            const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - start_).count();
            const long   t = (long)s;
            char hms[32];
            std::snprintf(hms, sizeof hms, "%ld:%02ld:%02ld", t / 3600, (t / 60) % 60, t % 60);
            out_ << prefix_ << "progress: " << line_(s) << ", " << hms << " elapsed" << std::endl;
        }
    }

    std::ostream&                         out_;
    std::string                           prefix_;
    LineFn                                line_;
    std::chrono::steady_clock::time_point start_;
    std::chrono::seconds                  interval_{0};
    std::mutex                            m_;
    std::condition_variable               cv_;
    bool                                  stop_ = false;
    std::thread                           thread_;
};

} // namespace kmer

#endif
