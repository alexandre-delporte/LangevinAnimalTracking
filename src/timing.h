#ifndef TIMING_H
#define TIMING_H

#include <chrono>
#include <string>
#include <map>
#include <RcppArmadillo.h>

class Timer {
private:
    std::map<std::string, double> timings;
    std::map<std::string, int> counts;
    std::chrono::high_resolution_clock::time_point start_time;
    
public:
    void start() {
        start_time = std::chrono::high_resolution_clock::now();
    }
    
    void record(const std::string& label) {
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        timings[label] += elapsed.count();
        counts[label]++;
        start_time = end_time;
    }
    
    void reset() {
        timings.clear();
        counts.clear();
    }
    
    Rcpp::List get_timings() {
        Rcpp::List result;
        for (auto& pair : timings) {
            Rcpp::List entry = Rcpp::List::create(
                Rcpp::Named("total_time") = pair.second,
                Rcpp::Named("count") = counts[pair.first],
                Rcpp::Named("mean_time") = pair.second / counts[pair.first]
            );
            result[pair.first] = entry;
        }
        return result;
    }
    
    void print_timings() {
        Rcpp::Rcout << "\n=== Profiling Results ===\n";
        Rcpp::Rcout << std::setw(30) << "Function" 
                    << std::setw(12) << "Total (s)" 
                    << std::setw(12) << "Count"
                    << std::setw(12) << "Mean (ms)\n";
        Rcpp::Rcout << std::string(66, '-') << "\n";
        
        for (auto& pair : timings) {
            Rcpp::Rcout << std::setw(30) << pair.first
                        << std::setw(12) << std::fixed << std::setprecision(3) << pair.second
                        << std::setw(12) << counts[pair.first]
                        << std::setw(12) << std::fixed << std::setprecision(3) 
                        << (pair.second / counts[pair.first] * 1000) << "\n";
        }
        Rcpp::Rcout << "========================\n\n";
    }
};

// Global timer instance
static Timer global_timer;

#endif // TIMING_H
