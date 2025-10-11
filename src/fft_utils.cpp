#include "fft_utils.h"
#include <cmath>
#include <complex>
#include <algorithm>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void FFTProcessor::fft(std::vector<std::complex<double>>& input) {
    size_t N = input.size();
    
    if (N <= 1) return;
    
    // Check if N is power of 2
    if (!isPowerOfTwo(N)) {
        std::cerr << "Warning: FFT input size " << N << " is not a power of 2. Consider padding." << std::endl;
        return;
    }
    
    // Bit-reversal permutation
    bitReverse(input);
    
    // Cooley-Tukey FFT
    for (size_t len = 2; len <= N; len *= 2) {
        double angle = -2.0 * M_PI / len;
        std::complex<double> wlen(cos(angle), sin(angle));
        
        for (size_t i = 0; i < N; i += len) {
            std::complex<double> w(1.0, 0.0);
            
            for (size_t j = 0; j < len / 2; j++) {
                std::complex<double> u = input[i + j];
                std::complex<double> v = input[i + j + len / 2] * w;
                
                input[i + j] = u + v;
                input[i + j + len / 2] = u - v;
                
                w *= wlen;
            }
        }
    }
}

void FFTProcessor::ifft(std::vector<std::complex<double>>& input) {
    size_t N = input.size();
    
    // Conjugate input
    for (auto& x : input) {
        x = std::conj(x);
    }
    
    // Forward FFT
    fft(input);
    
    // Conjugate output and scale
    for (auto& x : input) {
        x = std::conj(x) / static_cast<double>(N);
    }
}

std::vector<std::complex<double>> FFTProcessor::fftCorrelation(
    const std::vector<std::complex<double>>& signal,
    const std::vector<std::complex<double>>& reference) {
    
    // Ensure both vectors have the same size
    size_t N = std::max(signal.size(), reference.size());
    
    // Pad to next power of 2 for efficiency
    size_t paddedSize = nextPowerOfTwo(N);
    
    // Copy and pad signals
    std::vector<std::complex<double>> fft_signal(paddedSize, 0.0);
    std::vector<std::complex<double>> fft_reference(paddedSize, 0.0);
    
    std::copy(signal.begin(), signal.end(), fft_signal.begin());
    std::copy(reference.begin(), reference.end(), fft_reference.begin());
    
    // Forward FFT of both signals
    fft(fft_signal);
    fft(fft_reference);
    
    // Cross-correlation: fft_signal .* conj(fft_reference)
    std::vector<std::complex<double>> product(paddedSize);
    for (size_t i = 0; i < paddedSize; i++) {
        product[i] = fft_signal[i] * std::conj(fft_reference[i]);
    }
    
    // Inverse FFT
    ifft(product);
    
    // For cross-correlation, the result needs to be properly shifted
    // The peak should appear at the correct lag position
    std::vector<std::complex<double>> result(N);
    
    // Simply copy the first N elements - the correlation is already correctly positioned
    std::copy(product.begin(), product.begin() + N, result.begin());
    
    return result;
}

std::vector<std::complex<double>> FFTProcessor::fftConvolution(
    const std::vector<std::complex<double>>& signal,
    const std::vector<std::complex<double>>& reference) {
    
    // Ensure both vectors have the same size
    size_t N = std::max(signal.size(), reference.size());
    
    // Pad to next power of 2 for efficiency
    size_t paddedSize = nextPowerOfTwo(N);
    
    // Copy and pad signals
    std::vector<std::complex<double>> fft_signal(paddedSize, 0.0);
    std::vector<std::complex<double>> fft_reference(paddedSize, 0.0);
    
    std::copy(signal.begin(), signal.end(), fft_signal.begin());
    std::copy(reference.begin(), reference.end(), fft_reference.begin());
    
    // Forward FFT of both signals
    fft(fft_signal);
    fft(fft_reference);
    
    // Convolution: fft_signal .* fft_reference (NO conjugate)
    std::vector<std::complex<double>> product(paddedSize);
    for (size_t i = 0; i < paddedSize; i++) {
        product[i] = fft_signal[i] * fft_reference[i];  // No conjugate for convolution
    }
    
    // Inverse FFT
    ifft(product);
    
    // Copy the first N elements
    std::vector<std::complex<double>> result(N);
    std::copy(product.begin(), product.begin() + N, result.begin());
    
    return result;
}

std::vector<std::complex<double>> FFTProcessor::multiply(
    const std::vector<std::complex<double>>& a,
    const std::vector<std::complex<double>>& b) {
    
    size_t N = std::min(a.size(), b.size());
    std::vector<std::complex<double>> result(N);
    
    for (size_t i = 0; i < N; i++) {
        result[i] = a[i] * b[i];
    }
    
    return result;
}

std::vector<std::complex<double>> FFTProcessor::padToPowerOfTwo(
    const std::vector<std::complex<double>>& input) {
    
    size_t newSize = nextPowerOfTwo(input.size());
    std::vector<std::complex<double>> padded(newSize, 0.0);
    std::copy(input.begin(), input.end(), padded.begin());
    
    return padded;
}

void FFTProcessor::bitReverse(std::vector<std::complex<double>>& input) {
    size_t N = input.size();
    size_t j = 0;
    
    for (size_t i = 1; i < N; i++) {
        size_t bit = N >> 1;
        
        while (j & bit) {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        
        if (i < j) {
            std::swap(input[i], input[j]);
        }
    }
}

bool FFTProcessor::isPowerOfTwo(size_t n) {
    return n > 0 && (n & (n - 1)) == 0;
}

size_t FFTProcessor::nextPowerOfTwo(size_t n) {
    if (n <= 1) return 1;
    
    size_t power = 1;
    while (power < n) {
        power *= 2;
    }
    
    return power;
}