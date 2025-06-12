import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq

input_file = 'data/input_data.bin'
output_file = 'data/filtered_data.bin'
n_traces = 216
n_samples = 501
dt = 0.004
filter_file = 'data/filter.txt'
comparison_pic_file = 'pic/seismic_comparison.png'
filter_pic_file = 'pic/filter.png'

def read_binary_seismogram(filename, n_traces, n_samples):
   
    with open(filename, 'rb') as f:
        data = np.fromfile(f, dtype=np.float64)  
    return data.reshape(n_traces, n_samples)

def plot_filter(filter_file, pic_file):

    freqs = []
    amps = []
    with open(filter_file) as f:
        for line in f.readlines()[1:]:
            freq, amp = [float(s) for s in line.split(',')]
            if freq >=0:
                freqs.append(freq)
                amps.append(amp)
    plt.figure(figsize=(10,5))            
    plt.title('Ormsby Filter Response')
    plt.plot(freqs, amps)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude')
    plt.grid()
    plt.savefig(pic_file)
    print(f"Filter response plot is saved to {filter_pic_file}")



def plot_comparison(input_file, output_file, n_traces, n_samples, dt):
  
    input_data = read_binary_seismogram(input_file, n_traces, n_samples)
    output_data = read_binary_seismogram(output_file, n_traces, n_samples)
    
    time = np.arange(n_samples) * dt
    freq = fftfreq(n_samples, dt)[:n_samples//2]
    
    input_fft = fft(input_data, axis=1)[:, :n_samples//2]
    output_fft = fft(output_data, axis=1)[:, :n_samples//2]
    input_spectrum = np.mean(np.abs(input_fft), axis=0)
    output_spectrum = np.mean(np.abs(output_fft), axis=0)

    vmin = np.percentile(input_data, 2)
    vmax = np.percentile(input_data, 98)
    
    # Создание фигуры
    fig, axs = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle('Before and After Ormsby Bandpass Filter', fontsize=16)
    
    # Отрисовка сейсмограмм
    axs[0, 0].imshow(input_data.T, aspect='auto', cmap='seismic', 
                    extent=[0, n_traces, time[-1], 0], vmin=vmin, vmax=vmax)
    axs[0, 0].set_title('Initial Seismogram')
    axs[0, 0].set_xlabel('Time (s)')
    axs[0, 0].set_ylabel('Trace Number')
    
    axs[0, 1].imshow(output_data.T, aspect='auto', cmap='seismic', 
                    extent=[0, n_traces, time[-1], 0], vmin=vmin, vmax=vmax)
    axs[0, 1].set_title('Filtered Seismogram')
    axs[0, 1].set_xlabel('Time (s)')
    axs[0, 1].set_ylabel('Trace Number')
    
    # Отрисовка спектров
    axs[1, 0].plot(freq, input_spectrum)  
    axs[1, 0].set_title('Initial Spectrum')
    axs[1, 0].set_xlabel('Frequency (Hz)')
    axs[1, 0].set_ylabel('Amplitude')

    axs[1, 0].grid(True)

    axs[1, 1].plot(freq, output_spectrum)    
    axs[1, 1].set_title('Filtered Spectrum')
    axs[1, 1].set_xlabel('Frequency (Hz)')
    axs[1, 1].set_ylabel('Amplitude')

    axs[1, 1].grid(True)
    
    
   
    plt.tight_layout()
    plt.savefig(comparison_pic_file, dpi=300)
    plt.close()
    print(f"Comparison picture is saved to {filter_pic_file}")

if __name__ == "__main__":
  
    plot_filter(filter_file, filter_pic_file)
    plot_comparison(input_file, output_file, 
                   n_traces, n_samples, dt)