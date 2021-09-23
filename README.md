# Differential Emission Line Filter (DELF)

Applies a Differential Emission Line Filter (DELF) to an astronomical spectrum data cube.

## Background

The usual approach to find point sources such as planetary nebulae in astronomical observations has been to observe the object region using imaging techniques. In that approach, the region is observed both on-band and off-band using narrow bandpass filters. Such an approach can work with planetary nebulae as they emit nearly all their intensity in a few emission lines; where the forbidden emission line of oxygen, [OIII]5007, is typically the strongest one.

Astro-d11 (DELF) presents an alternative approach where a datacube based on integral-field spectroscopy observations provides the means to use two very narrow bandpasses near the emission line when subracting the background signal.


## Method

The blue and red bandpasses are offset from the current wavelength (layer) towards bluer and redder pixels, beginning at an initial offset (`offset`). The total width of the red and blue bandpass is set using the parameter `cwidth`; either bandpass is skipped for the bluest (lowest) and reddest (highest) pixels. The initially offset bandpasses are thereafter shifted away from the layer as needed in such a way that telluric and [optionally also] emission lines are excluded. Additionally, the subtracted continuum value is normalized with a spectrum selected in a pre-selected aperture with few emission-line features, using the same bandpasses. The aperture is defined using the parameters `x`, `y`, and `apr`.


### Telluric lines

The list of telluric lines is specified using the parameter `telluriclines`, which needs to be set to the name of a plain-text file where each line contains the wavelength of a telluric line in the first column (the unit is Angstrom, Å); the default line list file is `telluric_lines_hires.dat`, which is available in the `data` directory. The telluric-line bandpass width can be adjusted using the parameter `twidth` [Angstrom], where the default value is 3.0 Å.


### Emission lines

The list of emission lines is specified using the parameter `emissionlines`, which needs to be set to the name of a plain-text file where each line contains the wavelength of an emission line in the first column (the unit is Angstrom, Å); a default line list file is provided in `emission_lines-ground_based-noFe.dat`, which is also available in the `data` directory.

The procedure is to create a spatially dependent emission-line mask by looping through all spatial elements and emission-line entries. For this purpose, the emission line redshift can be set using the parameter `vel_z` (unit km/s; default is 0 km/s), and an additional permitted offset is specified using the parameter `dwl` (unit Angstrom; default is 1.0 Å). For each spatial element and emission line, a section of the object spectrum is fitted using the tool `mpfit.py` (see link below). A fitted line results in the bandpass centered on the wavelength to be masked. The emission line bandpass width is set using the parameter `bwidth` [Angstrom], where the default width is 3.0 Å.

Please Note! The fitting procedure of individual emission lines is slow.


## Links

The filter is described in the paper *Toward Precision Cosmology with Improved PNLF Distances Using VLT-MUSE I. Methodology and Tests*, Martin M. Roth, George H. Jacoby, Robin Ciardullo, Brian D. Davis, Owen Chase, and Peter M. Weilbacher 2021, [The Astrophysical Journal](https://iopscience.iop.org/journal/0004-637X), [916, 21, 44 pp. (PDF)](https://ui.adsabs.harvard.edu/link_gateway/2021ApJ...916...21R/PUB_PDF) [[*ApJ* abstract page](https://www.doi.org/10.3847/1538-4357/ac02ca), [NASA ADS](https://ui.adsabs.harvard.edu/abs/2021ApJ...916...21R/abstract)].

This tool is also available in the integral-field spectroscopy data-reduction package `p3d`, which is available at https://p3d.sourceforge.io, where the tool is called `p3d_d11`. While `p3d` is written using the Interactive Data Language (IDL), it can be used without a license using the IDL Virtual Machine.


## Installation

The code is available in the [python package index](https://pypi.org/project/astro-delf) and can be installed using `pip`
```
pip install astro-d11
```

The tool mpfit for Python 3 needs to be downloaded separately from [GitHub/astrolibpy/mpfit](https://github.com/segasai/astrolibpy/blob/master/mpfit/mpfit.py).
