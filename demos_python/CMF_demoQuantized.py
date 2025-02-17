import numpy as np
import matplotlib.pyplot as plt
import time
import pycmf  # Your binding module that exposes medfilt_circ2d_quant

def deg_to_rad(x):
    """
    Convert degrees to radians in the range (-pi, pi].
    Equivalent to MATLAB: (x/360 - 0.5) * 2*pi.
    """
    return (x / 360 - 0.5) * 2 * np.pi

def main():
    # Load wind direction data from the text file.
    # Adjust the path as needed.
    filename = '../data/CMF_windDir.txt'
    try:
        windDir = np.loadtxt(filename)
    except Exception as e:
        print("Error loading text file:", e)
        return

    # Ensure the data is 1D (or squeeze out extra dimensions).
    windDir = np.squeeze(windDir)
    
    # Convert from degrees to radians in (-pi, pi]
    windDir2Pi = deg_to_rad(windDir)
    
    # Create quantization levels for full degree resolution (0 to 359 degrees).
    quant = deg_to_rad(np.arange(0, 360))
    
    # Set filter size:
    # Data is recorded every 10 minutes; one day has 24*6 = 144 samples.
    # MATLAB uses R = 24*6 + 1.
    T = 24 * 6 + 1
    R = 1  # Since the data is 1D
    
    # Convert the data to Fortran order (column-major) for compatibility with the C++ code.
    windDir2Pi = np.asfortranarray(windDir2Pi)

    # the data has to be in the second dimension for compatibility with the C++ code
    windDir2Pi_reshaped = windDir2Pi.reshape((1, windDir2Pi.size))

    # Perform filtering using your quantized median filter.
    start_time = time.time()
    circleMedian = pycmf.medfilt_circ2d_quant(windDir2Pi_reshaped, R, T, quant)
    elapsed_time = time.time() - start_time
    print("Filtering took {:.3f} seconds".format(elapsed_time))
    
    # restore the original shape    
    circleMedian = np.squeeze(circleMedian)

    # Plot the original and filtered signals.
    plt.figure(figsize=(10, 8))
    
    plt.subplot(2, 1, 1)
    plt.plot(windDir2Pi, '.')
    plt.title('Data (converted to radians)')
    plt.axis('tight')
    
    plt.subplot(2, 1, 2)
    plt.plot(circleMedian, '.')
    plt.title('Arc distance median filter')
    plt.axis('tight')
    
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()