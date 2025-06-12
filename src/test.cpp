#include "ormsbybpf.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>

using namespace std;

constexpr const char* input_file = "../data/input_data.bin";
constexpr const char* output_file = "../data/output_data.bin";
constexpr const char* filter_file = "../data/filter.txt";
constexpr double dt = 0.004;
constexpr int n_traces = 216;
constexpr int n_samples = 501;
constexpr double f1 = 10.0;
constexpr double f2 = 18.0;
constexpr double f3 = 40.0;
constexpr double f4 = 50.0;

// Function to read binary seismogram file
vector<vector<double>> readSeismogram(const string& filename, 
                                     int n_traces, int n_samples) {
    ifstream file(filename, ios::binary);
    if (!file) {
        throw runtime_error("Cannot open file: " + filename);
    }

    vector<vector<double>> seismogram(n_traces, vector<double>(n_samples));
    
    for (int trace = 0; trace < n_traces; ++trace) {
        file.read(reinterpret_cast<char*>(seismogram[trace].data()), 
                n_samples * sizeof(double));
        if (!file) {
            throw runtime_error("Error reading data from file");
        }
    }
    
    return seismogram;
}

// Function to write filter response to text file
void writeFilter(const string& filename, const vector<double>& freqs, const vector<double>& amps) {  
    ofstream file(filename);

   if (!file.is_open()) {
        throw runtime_error("Error opening file: " + filename);  
    }
    
    file << "Freq,Amp\n";
    
    for (size_t i = 0; i < freqs.size(); ++i) {
        file << freqs[i] << "," << amps[i] << "\n";
    }       
    file.close();
}        

// Function to write binary seismogram file
void writeSeismogram(const string& filename, 
                    const vector<vector<double>>& seismogram) {
    ofstream file(filename, ios::binary);
    if (!file) {
        throw runtime_error("Cannot create file: " + filename);
    }

    for (const auto& trace : seismogram) {
        file.write(reinterpret_cast<const char*>(trace.data()), 
                 trace.size() * sizeof(double));
        if (!file) {
            throw runtime_error("Error writing data to file");
        }
    }
}

int main(int argc, char* argv[]) {    

    try {        
        
        // Read input seismogram
        cout << "Reading seismogram from " << input_file << "...\n";
        auto seismogram = readSeismogram(input_file, n_traces, n_samples);

        // Save Ormsby filter to text file
        vector<double> freqs = createFrequencyAxis(dt, n_samples);
        vector<double> filter = createOrmsbyFilter(f1, f2, f3, f4, freqs);
        writeFilter(filter_file, freqs, filter);   

        // Apply bandpass filter
        cout << "Applying Ormsby bandpass filter (" 
             << f1 << ", " << f2 << ", " << f3 << ", " << f4 << " Hz)...\n";
        auto filtered = ormsbyBPF(seismogram, dt, f1, f2, f3, f4);

        // Write filtered seismogram
        cout << "Writing filtered seismogram to " << output_file << "...\n";
        writeSeismogram(output_file, filtered);

        cout << "Processing completed successfully.\n";
    } 
    catch (const exception& e) {
        cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}