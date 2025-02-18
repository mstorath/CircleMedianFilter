import numpy as np
import pycirclemedianfilter as cmf
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import time
from PIL import Image
import os

def phase_to_rgb(x):
    """
    Convert a phase image (with values in [-pi, pi]) to an RGB image using HSV.
    The hue is computed as: hue = (x/(2*pi)) + 0.5, with saturation and value fixed to 1.
    """
    # Compute hue in [0, 1]. (For x in [-pi, pi], x/(2*pi) is in [-0.5,0.5], then add 0.5.)
    hue = x / (2 * np.pi) + 0.5
    # Create saturation and value arrays (all ones, same shape as hue)
    sat = np.ones_like(hue)
    val = np.ones_like(hue)
    # Stack into an HSV image of shape (..., 3)
    hsv = np.stack((hue, sat, val), axis=-1)
    # Convert HSV to RGB
    rgb = colors.hsv_to_rgb(hsv)
    return rgb

def main():
    # Load the InSAR image. Replace 'CMF_imgInSAR.png' with the correct path if needed.
    path = os.path.abspath(os.path.dirname(__file__))
    filename = os.path.join(path, '../data/CMF_imgInSAR.png')

    img = Image.open(filename)
    
    
    img = np.array(img, dtype=np.float64)
    
    # If the image was normalized to [0,1] (common with some imread functions),
    # scale it to [0,255] first.
    if img.max() <= 1.0:
        img = img * 255.0

    # Scale from [0, 255] to [-pi, pi]
    img2pi = (img / 255.0 - 0.5) * 2.0 * np.pi
    
    # Set filter size (R is used for both dimensions, so here T=R)
    R = 3
    
    # Convert the image to Fortran order (column-major) for compatibility with the C++ code
    img2pi = np.asfortranarray(img2pi)

    # Apply the circle-median filter and measure execution time.
    start_time = time.time()
    circle_median = cmf.medfilt_circ2d(img2pi, R, R)
    elapsed_time = time.time() - start_time
    print("Filtering time: {:.3f} seconds".format(elapsed_time))
    
    # Plot the original InSAR data and the filtered result side by side.
    fig, axs = plt.subplots(1, 3, figsize=(6, 6))

    axs[0].imshow(phase_to_rgb(img2pi))
    axs[0].set_title('InSAR Data')
    axs[0].axis('off')
    
    axs[1].imshow(phase_to_rgb(circle_median))
    axs[1].set_title('Arc distance median filter (3x3)')
    axs[1].axis('off')

    axs[2].imshow(img2pi - circle_median)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()