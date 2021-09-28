# Differential Emission Line Filter (DELF)

Apply a Differential Emission Line Filter (DELF) to an astronomical spectrum data cube.

## Background

The usual approach to find point sources such as planetary nebulæ (PNe) in astronomical observations has been to observe the object region using imaging techniques. In that approach, the region is observed both on-band and off-band using narrow bandpass filters; a comparison between the two images reveals objects such as PNe. Such an approach can work with PNe as they emit nearly all their intensity in a few emission lines; where the forbidden emission line of oxygen, [OIII]5007, is typically the strongest one.

Astro-d11 (DELF) presents an alternative approach where a data cube based on integral-field spectroscopy observations provides means to use two very narrow bandpasses near the emission line when subracting the background signal. In comparison to the imaging approach, the narrow "filters" represented by the bandpasses should make it possible to detect fainter objects!

The algorithm is first described in the paper mentioned in the Links section below.


## Method

Two narrow bandpasses, a blue and a red bandpass, are offset from the current wavelength (layer) towards bluer (lower) and redder (higher) pixels, beginning at an initial offset (`offset`). The total width of the red and blue bandpasses is set using the parameter `cwidth`; either bandpass is skipped for the bluest (lowest) and reddest (highest) pixels. The initially offset bandpasses are thereafter shifted away from the layer as needed in such a way that telluric and [optionally also] emission lines are avoided. Additionally, the subtracted continuum value is normalized with a reference spectrum of a pre-selected aperture with few emission-line features, using the same bandpasses. The location and size of the reference aperture must be set using the parameters `x`, `y`, and `apr`.

The reference spectrum (`rspec`) and its continuum bandpasses (`rspec_blue` and `rspec_red`) are defined with `n_blue` and `n_red` layers in the blue and red bandpasses, respectively. Likewise, the flux and continuum bandpasses of each spatial element are defined with `img`, `img_blue`, and `img_red`, respectively; using the same bandpasses as the reference spectrum! The continuum is then subtracted from the input data cube for the current layer `i` using the following equation:

```
corr_factor[i] = rspec[i] /
     ((sum(rspec_blue[i]) / n_blue[i] + sum(rspec_red[i]) / n_red[i]) / 2)
out[i] = img[i] - corr_factor[i] *
      (sum(img_blue[i]) / n_blue[i] + sum(img_red[i]) / n_red[i]) / 2
```


### Telluric lines

The list of telluric lines is specified using the parameter `telluriclines`, which needs to be set to the name of a plain-text file where each line contains the wavelength of a telluric line in the first column (the unit is Angstrom, Å); the default line list file is *telluric_lines_hires.dat*, which is available in the *data* directory. The bandpass width can be adjusted using the parameter `bwidth` [Angstrom], where the default value is 3.0 Å.


### Emission lines

The list of emission lines is specified using the parameter `emissionlines`, which needs to be set to the name of a plain-text file where each line contains the wavelength of an emission line in the first column (the unit is Angstrom, Å); a default line list file is provided in *emission_lines-ground_based-noFe.dat*, which is also available in the *data* directory.

The procedure is to create a spatially dependent emission-line mask by looping through all spatial elements and emission-line entries. For this purpose, and to save execution time, the data can be binned on the spatial axes to create spectra with higher signal-to-noise before the fitting. See the parameter `bin`.

The emission line redshift can be set using the parameter `vel_z` (unit km/s; default is 0 km/s), and an additional permitted offset is specified using the parameter `dwl` (unit Angstrom; default is 1.0 Å). For each spatial element and emission line, a section of the object spectrum is fitted using the tool __mpfit.py__ (see link below). A fitted line results in the bandpass centered on the wavelength to be masked. The emission line bandpass width is set using the parameter `bwidth` [Angstrom], where the default width is 3.0 Å.

Please Note! The fitting procedure of individual emission lines is slow. So it might be a wise idea to begin with a small number of emission lines in the list to see that everything works properly before increasing the number.


### Resulting Image

The filtered image is written to a file, adding a set of header keywords that indicate waht argument values were used (`d11_x`, `d11_y`, `d11_apr`, and `d11_cwid`) for the parameters `x, y, apr, cwidth`. The output filename can be set explicitly using the parameter `ofilename`), otherwise the input filename is used with the added suffix *_d11*.


## Links

The filter is described in the paper _Toward Precision Cosmology with Improved PNLF Distances Using VLT-MUSE I. Methodology and Tests_, Martin M. Roth, George H. Jacoby, Robin Ciardullo, Brian D. Davis, Owen Chase, and Peter M. Weilbacher 2021, [The Astrophysical Journal](https://iopscience.iop.org/journal/0004-637X), [916, 21, 44 pp. (PDF)](https://ui.adsabs.harvard.edu/link_gateway/2021ApJ...916...21R/PUB_PDF) [[*ApJ* abstract page](https://www.doi.org/10.3847/1538-4357/ac02ca), [NASA ADS](https://ui.adsabs.harvard.edu/abs/2021ApJ...916...21R/abstract)].

This tool is also available in the integral-field spectroscopy data-reduction package __p3d__, which is available at https://p3d.sourceforge.io, where the tool is named __p3d_d11__. While __p3d__ is written using the Interactive Data Language (IDL), it can be used without a license using the IDL Virtual Machine.


## Installation

The code is available in the [python package index](https://pypi.org/project/astro-delf) and can be installed using `pip`
```
pip install astro-d11
```

The tool __mpfit.py__ for Python 3 needs to be retrieved separately from [GitHub/astrolibpy/mpfit](https://github.com/segasai/astrolibpy/blob/master/mpfit/mpfit.py).
