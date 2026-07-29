#ifndef MEM_LIMIT_HPP
#define MEM_LIMIT_HPP

// ===========================================================================
//  detect_memory_limit — how much RAM is actually enforced on this process.
// ===========================================================================
//
//  Stage 1 sizes its accumulation budget from THIS, not from a config number,
//  because the config number is wrong exactly when it matters:
//
//    * `--mem=0` / `--exclusive` — the job owns the whole node, but the config
//      still says e.g. 128 GB. The budget should track the node, not 128.
//    * a QOS/partition cap below what was requested.
//    * a bare workstation with no cgroup and no SLURM at all.
//
//  Several sources are read and the MINIMUM taken, because any one of them can
//  be the binding constraint:
//
//    1. cgroup v2 — walk /proc/self/cgroup's `0::` path from the process's own
//       cgroup up to the root, reading memory.max / memory.high at each level.
//       The leaf frequently reports "max"; the real limit lives on an ANCESTOR
//       (measured on Ibex: the job's limit sat on job_<id>/memory.max while the
//       task leaf was unlimited), so walking up and taking the min is required.
//    2. cgroup v1 — /sys/fs/cgroup/memory/memory.limit_in_bytes, ignoring the
//       near-SIZE_MAX "unlimited" sentinel.
//    3. SLURM env — SLURM_MEM_PER_NODE (MB), or SLURM_MEM_PER_CPU x
//       SLURM_CPUS_PER_TASK. Carries `--mem` even when the cluster does not
//       enforce it through a cgroup.
//    4. MemTotal from /proc/meminfo — the final floor (bare workstation).
//
//  Purely observational: it reads, never sets. Output is unaffected by the
//  budget, so this can only change how often Stage 1 spills, never the result.
// ===========================================================================

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <string>

namespace kmer {

struct MemLimit {
    size_t      bytes = 0;   // 0 => nothing could be determined
    std::string source;      // human-readable, for the startup log
};

namespace mem_detail {

// Read a single non-negative integer from a one-line sysfs file. Returns 0 for
// a missing file or the cgroup "max" (unlimited) sentinel.
inline size_t read_size_file(const std::string& path) {
    std::ifstream f(path);
    if (!f) return 0;
    std::string tok;
    if (!(f >> tok)) return 0;
    if (tok == "max") return 0;
    char* end = nullptr;
    unsigned long long v = std::strtoull(tok.c_str(), &end, 10);
    if (end == tok.c_str()) return 0;
    return (size_t)v;
}

// cgroup v1 uses a value near SIZE_MAX to mean "no limit".
inline bool v1_unlimited(size_t v) {
    return v == 0 || v > ((size_t)1 << 62);
}

// The v2 cgroup path from /proc/self/cgroup, whose sole line is "0::/<path>".
inline std::string cgroup_v2_path() {
    std::ifstream f("/proc/self/cgroup");
    std::string line;
    while (std::getline(f, line))
        if (line.rfind("0::", 0) == 0) return line.substr(3);
    return {};
}

} // namespace mem_detail

inline MemLimit detect_memory_limit() {
    size_t      best = std::numeric_limits<size_t>::max();
    std::string src  = "none";
    auto consider = [&](size_t v, const std::string& s) {
        if (v > 0 && v < best) { best = v; src = s; }
    };

    // 1. cgroup v2 — walk from the process's cgroup up to /sys/fs/cgroup.
    const std::string root = "/sys/fs/cgroup";
    std::string cg = mem_detail::cgroup_v2_path();
    if (!cg.empty()) {
        std::string dir = root + cg;
        for (;;) {
            for (const char* fn : { "/memory.max", "/memory.high" }) {
                size_t v = mem_detail::read_size_file(dir + fn);
                if (v > 0) consider(v, "cgroup2 " + dir + fn);
            }
            if (dir.size() <= root.size()) break;
            size_t slash = dir.find_last_of('/');
            dir = (slash == std::string::npos || slash < root.size())
                      ? root : dir.substr(0, slash);
        }
    }

    // 2. cgroup v1.
    {
        size_t v = mem_detail::read_size_file(root + "/memory/memory.limit_in_bytes");
        if (!mem_detail::v1_unlimited(v)) consider(v, "cgroup1 memory.limit_in_bytes");
    }

    // 3. SLURM environment.
    if (const char* m = std::getenv("SLURM_MEM_PER_NODE")) {
        size_t mb = std::strtoull(m, nullptr, 10);
        if (mb) consider(mb * 1024 * 1024, "SLURM_MEM_PER_NODE");
    } else if (const char* pc = std::getenv("SLURM_MEM_PER_CPU")) {
        size_t mb = std::strtoull(pc, nullptr, 10);
        const char* ncpu = std::getenv("SLURM_CPUS_PER_TASK");
        size_t n = ncpu ? std::strtoull(ncpu, nullptr, 10) : 1;
        if (mb && n) consider(mb * n * 1024 * 1024, "SLURM_MEM_PER_CPU x CPUS");
    }

    // 4. MemTotal — the final floor.
    {
        std::ifstream f("/proc/meminfo");
        std::string key;
        while (f >> key) {
            if (key == "MemTotal:") {
                size_t kb; std::string unit;
                if (f >> kb >> unit) consider(kb * 1024, "MemTotal");
                break;
            }
            std::string rest; std::getline(f, rest);
        }
    }

    MemLimit out;
    if (best != std::numeric_limits<size_t>::max()) { out.bytes = best; out.source = src; }
    return out;   // bytes == 0, source "none" if nothing was found
}

} // namespace kmer

#endif
