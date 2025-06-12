#ifndef BPF_H
#define BPF_H

#include <vector>

std::vector<double> createFrequencyAxis(double dt, int nt);
std::vector<double> createOrmsbyFilter(double f1, double f2, double f3, double f4, 
                                 const std::vector<double>& freqs);

std::vector<std::vector<double>> ormsbyBPF(
    const std::vector<std::vector<double>>& seismogram,
    double dt, double f1, double f2, double f3, double f4);

#endif // ORMSBY_FILTER_H