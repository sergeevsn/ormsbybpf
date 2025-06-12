#include "ormsbybpf.h"
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <fftw3.h>

using namespace std;

// Creating Frequency Axis
vector<double> createFrequencyAxis(double dt, int nt) {
    vector<double> freqs(nt);
    double df = 1.0 / (dt * nt);
    for (int i = 0; i < nt; ++i) {
        freqs[i] = (i <= nt / 2) ? (i * df) : ((i - nt) * df);
    }
    return freqs;
}

// Creating Ormsby Bandpass Filter
vector<double> createOrmsbyFilter(double f1, double f2, double f3, double f4, const vector<double>& freqs) {
    vector<double> filter(freqs.size(), 0.0);
    for (size_t i = 0; i < freqs.size(); ++i) {
        double f = freqs[i];
        // Positive Frequencies
        if (f > f1 && f <= f2)
            filter[i] = (f - f1) / (f2 - f1);
        else if (f > f2 && f <= f3)
            filter[i] = 1.0;
        else if (f > f3 && f <= f4)
            filter[i] = 1.0 - (f - f3) / (f4 - f3);
        // Negative Frequencies
        else if (f < -f1 && f >= -f2)
            filter[i] = (-f - f1) / (f2 - f1);
        else if (f < -f2 && f >= -f3)
            filter[i] = 1.0;
        else if (f < -f3 && f >= -f4)
            filter[i] = 1.0 - (-f - f3) / (f4 - f3);
    }
    return filter;
}

// Apply Bandpass Filtering
vector<vector<double>> ormsbyBPF(
    const vector<vector<double>>& seismogram,
    double dt, double f1, double f2, double f3, double f4) {
    
    if (seismogram.empty()) return {};
    const int ntraces = seismogram.size();
    const int nt = seismogram[0].size();

    // Check traces sizes
    for (const auto& trace : seismogram) {
        if (trace.size() != nt)
            throw runtime_error("Inconsistent trace length in seismogram.");
    }

    // Frequency axis
    vector<double> freqs = createFrequencyAxis(dt, nt);
    vector<double> filter = createOrmsbyFilter(f1, f2, f3, f4, freqs);

    // FFT Buffers
    vector<complex<double>> in(nt), out(nt);
    fftw_plan plan_forward = fftw_plan_dft_1d(nt,
        reinterpret_cast<fftw_complex*>(in.data()),
        reinterpret_cast<fftw_complex*>(out.data()),
        FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_plan plan_backward = fftw_plan_dft_1d(nt,
        reinterpret_cast<fftw_complex*>(out.data()),
        reinterpret_cast<fftw_complex*>(in.data()),
        FFTW_BACKWARD, FFTW_ESTIMATE);
    
    vector<vector<double>> result(ntraces, vector<double>(nt));

    for (int trace = 0; trace < ntraces; ++trace) {
        // Forward FFT
        std::copy(seismogram[trace].begin(), seismogram[trace].end(), in.begin());
        fftw_execute(plan_forward);

        // Apply filter 
        for (int i = 0; i < nt; ++i) {
            out[i] *= filter[i];
        }

        // Backward FFT
        fftw_execute(plan_backward);

        // Normalizing and saving result
        for (int i = 0; i < nt; ++i) {
            result[trace][i] = in[i].real() / nt;
        }
    }
    
    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);

    return result;
}
