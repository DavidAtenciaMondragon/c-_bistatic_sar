#ifndef FFT_UTILS_H
#define FFT_UTILS_H

#include <vector>
#include <complex>

/**
 * @brief Efficient FFT implementation using Cooley-Tukey algorithm
 * 
 * This is a self-contained FFT implementation that doesn't require external libraries.
 * Optimized for radar signal processing applications.
 */
class FFTProcessor {
public:
    /**
     * @brief Compute Forward FFT
     * @param input Input signal (will be modified in-place)
     */
    static void fft(std::vector<std::complex<double>>& input);
    
    /**
     * @brief Compute Inverse FFT
     * @param input Input signal (will be modified in-place)
     */
    static void ifft(std::vector<std::complex<double>>& input);
    
    /**
     * @brief FFT-based correlation equivalent to MATLAB: ifft(fft(signal) .* conj(fft(reference)))
     * @param signal Input signal
     * @param reference Reference signal
     * @return Correlated signal
     */
    static std::vector<std::complex<double>> fftCorrelation(
        const std::vector<std::complex<double>>& signal,
        const std::vector<std::complex<double>>& reference);
    
    /**
     * @brief FFT-based convolution equivalent to MATLAB: ifft(fft(signal) .* fft(reference))
     * @param signal Input signal
     * @param reference Reference signal
     * @return Convolved signal
     */
    static std::vector<std::complex<double>> fftConvolution(
        const std::vector<std::complex<double>>& signal,
        const std::vector<std::complex<double>>& reference);
    
    /**
     * @brief Element-wise multiplication of two complex vectors
     * @param a First vector
     * @param b Second vector
     * @return Element-wise product
     */
    static std::vector<std::complex<double>> multiply(
        const std::vector<std::complex<double>>& a,
        const std::vector<std::complex<double>>& b);
    
    /**
     * @brief Pad vector to next power of 2 for efficient FFT
     * @param input Input vector
     * @return Padded vector
     */
    static std::vector<std::complex<double>> padToPowerOfTwo(
        const std::vector<std::complex<double>>& input);
    
private:
    /**
     * @brief Bit-reversal permutation for FFT
     * @param input Input vector
     */
    static void bitReverse(std::vector<std::complex<double>>& input);
    
    /**
     * @brief Check if number is power of 2
     * @param n Number to check
     * @return True if power of 2
     */
    static bool isPowerOfTwo(size_t n);
    
    /**
     * @brief Find next power of 2 greater than or equal to n
     * @param n Input number
     * @return Next power of 2
     */
    static size_t nextPowerOfTwo(size_t n);
};

#endif // FFT_UTILS_H