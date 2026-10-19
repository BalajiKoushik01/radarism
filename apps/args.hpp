#pragma once
#include <cstdlib>
#include <map>
#include <set>
#include <stdexcept>
#include <string>

// key=value command-line parser that rejects unknown keys.
struct Args {
    std::map<std::string, std::string> kv;
    std::set<std::string> used;
    Args(int argc, char** argv) {
        for (int i = 1; i < argc; ++i) {
            const std::string s = argv[i];
            const auto pos = s.find('=');
            if (pos == std::string::npos || pos == 0) throw std::invalid_argument("expected key=value, got: " + s);
            kv[s.substr(0, pos)] = s.substr(pos + 1);
        }
    }
    bool has(const std::string& k) const { return kv.count(k) > 0; }
    std::string str(const std::string& k, const std::string& def) {
        used.insert(k);
        auto it = kv.find(k);
        return it == kv.end() ? def : it->second;
    }
    double num(const std::string& k, double def) {
        used.insert(k);
        auto it = kv.find(k);
        if (it == kv.end()) return def;
        char* e = nullptr;
        const double v = std::strtod(it->second.c_str(), &e);
        if (it->second.empty() || *e != '\0') throw std::invalid_argument("bad number for " + k + ": " + it->second);
        return v;
    }
    void check_unused() const {
        for (const auto& [k, v] : kv)
            if (!used.count(k)) throw std::invalid_argument("unknown option: " + k + " (run with 'help')");
    }
};
