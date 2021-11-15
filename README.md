# frc

`frc` is a Python library for calculating Fourier Ring Correlation (FRC) curves and their associated resolution.

This is particularly useful for determining the resolution of images taken with super-resolution microscopy techniques.

Its most computationally intensive functions are all implemented either in NumPy or Rust, making this library very fast.

### Fourier Ring Correlation

Introduced in 1982, Fourier Ring Correlation compares two 2D images, which are assumed to be noise-independent. In Fourier space, mage structure dominates the lower spatial frequencies (in the Fourier domain), while noise dominates at the higher frequencies.

In Fourier space, there are rings of constant spatial frequency around the origin. By calculating the correlation between the rings in the two images, you get the FRC curve:

![formula](https://render.githubusercontent.com/render/math?math=%5Ccolor%7Bgray%7D%5Ctext%7BFRC%7D%28q%29%20%3D%20%5Cfrac%7B%5Csum_%7B%5Cvec%7Bq%7D%20%5Cin%20%5Ctext%7Bcircle%7D%7D%20%5Cwidehat%7Bf_1%7D%28%5Cvec%7Bq%7D%29%20%5Cwidehat%7Bf_2%7D%28%5Cvec%7Bq%7D%29%5E%7B%5Ctextbf%7B%2A%7D%7D%7D%7B%5Csqrt%7B%5Csum_%7B%5Cvec%7Bq%7D%20%5Cin%20%5Ctext%7Bcircle%7D%7D%20%5Clvert%5Cwidehat%7Bf_1%7D%28%5Cvec%7Bq%7D%29%5Crvert%5E2%7D%20%5Csqrt%7B%5Csum_%7B%5Cvec%7Bq%7D%20%5Cin%20%5Ctext%7Bcircle%7D%7D%20%5Clvert%5Cwidehat%7Bf_2%7D%28%5Cvec%7Bq%7D%29%5Crvert%5E2%7D%7D%0A)

See the accompanying `FRC.pdf` for additional details.

At some spatial frequency, the signal cannot be separated from the noise. What spatial frequency exactly depends on what threshold function is used. The standard 0.143 and 1/2-bit thresholds are both available, as well as an algorithm to compute the intersection and resulting resolution.



### Installation

You can download this library from PyPI:

```shell
pip install frc
```

This library depends on [tmtenbrink/rustfrc](https://www.github.com/tmtenbrink/rustfrc), [DIPlib](https://github.com/DIPlib/diplib), [scipy](https://scipy.org/) and of course, [numpy](https://numpy.org/).

### Usage

The code snippet below (which for illustration purposes assumes you have a numpy array or DIP image representing your input, its associated scale and also matplotlib to plot the result) will calculate the FRC curve and the associated resolution using the standard 1/7 threshold.

```python
import frc
import numpy as np
import matplotlib.pyplot as plt

... # get an image and scale

img = np.array(img)
# Input can be a numpy array or DIP image
img = frc.util.square_image(img, add_padding=False)
img = frc.util.apply_tukey(img)
# Apply 1FRC technique
frc_curve = frc.one_frc(img)

img_size = img.shape[0]
xs_pix = np.arange(len(frc_curve)) / img_size
# scale has units [pixels <length unit>^-1] corresponding to original image
xs_nm_freq = xs_pix * scale
frc_res, res_y, thres = frc.frc_res(xs_nm_freq, frc_curve, img_size)
plt.plot(xs_nm_freq, thres(xs_nm_freq))
plt.plot(xs_nm_freq, frc_curve)
plt.show()
```