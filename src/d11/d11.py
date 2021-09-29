#!/usr/bin/env python
##
## d11.py
##
## --
##
## Author: Christer Sandin
##         Sandin Advanced Visualization (SAV)
##         Stockholm, SWEDEN
##
## d11: Differential Emission Line Filter (DELF)
##
## Copyright 2021
##           Leibniz Institute for Astrophysics Potsdam (AIP)
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##
## 1. Redistributions of source code must retain the above copyright notice,
##    this list of conditions and the following disclaimer.
##
## 2. Redistributions in binary form must reproduce the above copyright notice,
##    this list of conditions and the following disclaimer in the documentation
##    and/or other materials provided with the distribution.
##
## 3. Neither the name of the copyright holder nor the names of its
##    contributors may be used to endorse or promote products derived from this
##    software without specific prior written permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
## AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
## IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
## ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
## LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
## CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
## SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
## INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
## CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
## ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
## POSSIBILITY OF SUCH DAMAGE.

import numpy as np
import scipy
import logging

def d11_spec_sec(i_0, i_1, tmask, emask=None, pos=True):

    N = len(tmask)
    i_str = ""
    i_type = 0
    use_emask = False
    if emask is not None: use_emask = True

    # Step over masked pixels:

    if pos:
        if i_0 > N - 1: i_0 = N - 1
        if i_1 > N: i_1 = N
        if i_0 == i_1 - 1: return (i_0, i_1, "", 0)

        if use_emask:
            mtmask = max(tmask[i_0 : i_1])
            memask = max(emask[i_0 : i_1])

            while mtmask > 0 or memask > 0:
                it__str = ""
                if mtmask > 0: it_str = str(tmask[i_0]) + " (red)"
                ie_str = ""
                if memask > 0: ie_str = "er:" + str(memask) + "; "
                i_str = ie_str + it_str

                i_type = 0
                if mtmask > 0.0: i_type = 1
                if memask > 0: i_type += 2
                i_0 = i_0 + 1
                i_1 = i_1 + 1

                if i_0 > N - 1: i_0 = N - 1
                if i_1 > N: i_1 = N
                if i_0 == i_1 - 1: return (i_0, i_1, "", 0)
        else:
            while max(tmask[i_0 : i_1]) > 0:
                i_str = str(tmask[i_0]) + " (red)"
                i_type = 1
                i_0 = i_0 + 1
                i_1 = i_1 + 1

                if i_0 > N - 1: i_0 = N - 1
                if i_1 > N: i_1 = N
                if i_0 == i_1 - 1: return (i_0, i_1, "", 0)
    else:
        if i_0 < 0: i_0 = 0
        if i_1 < 1: i_1 = 1
        if i_0 == i_1 - 1: return (i_0, i_1, "", 0)

        if use_emask:
            mtmask = max(tmask[i_0 : i_1])
            memask = max(emask[i_0 : i_1])

            while mtmask > 0 or memask > 0:
                it__str = ""
                if mtmask > 0: it_str = str(tmask[i_1]) + " (blue)"
                ie_str = ""
                if memask > 0: ie_str = "eb:" + str(memask) + "; "
                i_str = ie_str + it_str

                i_type = 0
                if mtmask > 0.0: i_type = 1
                if memask > 0: i_type += 2
                i_0 = i_0 - 1
                i_1 = i_1 - 1

                if i_0 < 0: i_0 = 0
                if i_1 < 1: i_1 = 1
                if i_0 == i_1 - 1: return (i_0, i_1, "", 0)
        else:
            while max(tmask[i_0 : i_1]) > 0:
                i_str = str(tmask[i_1 - 1]) + " (blue)"
                i_type = 1
                i_0 = i_0 - 1
                i_1 = i_1 - 1

                if i_0 < 0: i_0 = 0
                if i_1 < 1: i_1 = 1
                if i_0 == i_1 - 1: return (i_0, i_1, "", 0)

    return (i_0, i_1, i_str, i_type)


def d11_mpfit(w_init, dwl, cdisp, x=None, y=None, w_too=None,
              fit_intensity_limit=0.0, fit_flux_continuum_fraction=0.0,
              xstr="", verbose=None, debug=False, contall=False):

    from mpfit import mpfit
    from d11_mpfit_f import d11_mpfit_f

    error = 0
    ok_fit = 0

    # Fit the data with a Gaussian.

    norm = sum(y) / y.size
    y /= norm

    n_fits = 2 + 3
    n_tied = 0 if w_too is None else w_too.size
    if n_tied > 0: n_fits += n_tied


    # Configure the fit.

    parinfo = [{'value':0., 'fixed':0, 'limited':[0, 0], 'limits':[0., 0.]}
	       for i in range(n_fits)]

    # line center.
    parinfo[2]['value'] = w_init
    parinfo[2]['limited'] = [1, 1]
    parinfo[2]['limits'][0] = w_init - dwl / cdisp
    parinfo[2]['limits'][1] = w_init + dwl / cdisp

    # sigma.
    sigma_min = 0.8 * 2.0 / np.sqrt(8.0 * np.log(2.0))
    sigma_max = 1.2 * 2.0 / np.sqrt(8.0 * np.log(2.0))
    parinfo[3]['value'] = 2.0 / np.sqrt(8.0 * np.log(2.0))
    parinfo[3]['limited'] = [1, 1]
    parinfo[3]['limits'][0] = sigma_min
    parinfo[3]['limits'][1] = sigma_max

    # intensity.
    parinfo[4]['value'] = 1e3
    parinfo[4]['limited'] = [1, 0]

    # Intensities of additional lines.
    if n_tied > 0:
        coff = np.zeros(n_tied)
        for i in range(0, n_tied):
            ii = 5 + i

            # line center.
            coff[i] = w_too[i] - w_init

            # intensity.
            parinfo[ii]['value'] = 1e3
            parinfo[ii]['limited'] = [1, 0]
    else:
        coff = np.zeros(1)

    p0 = np.zeros(n_fits)
    for i in range(0, n_fits): p0[i] = parinfo[i]['value']

    dy = np.sqrt(y)
    fctargs = {'x':x, 'y':y, 'dy':dy, 'coff':coff, 'n_tied':n_tied}

    # Perform the fit.
    quiet = 1 if verbose < 4 else 0
    m = mpfit(d11_mpfit_f, p0, functkw=fctargs, parinfo=parinfo, \
              maxiter=100, quiet=quiet)

    y *= norm
    dy *= norm
    m.params[0 : 2] *= norm  # Constant offset and slope
    m.params[4 :] *= norm  # Intensities

    ## From the MPFIT documentation:
    #   *If* you can assume that the true reduced chi-squared value is
    #   unity - meaning that the fit is implicitly assumed to be of good
    #   quality - then the estimated parameter uncertainties can be
    #   computed by scaling PERROR by the measured chi-squared value.
    #   sigpar *= sqrt(chisq / dof)

    yfit = d11_mpfit_f(m.params, x=x, y=y, dy=dy, \
                       coff=coff, n_tied=n_tied, funceval=True)

    yfit_p2 = d11_mpfit_f(m.params, x=m.params[2], y=0.0, dy=0.0, \
                       coff=coff, n_tied=n_tied, funceval=True)

    # Need to determine if the fit is good...this could be improved
    yfit_bg = (m.params[0] + m.params[1]*m.params[2])
    ok_fit = m.nfev > 1 and m.status > 0 and \
        m.params[3] > sigma_min and \
        m.params[3] < sigma_max and \
        m.params[4] > fit_intensity_limit and \
        (yfit_p2-yfit_bg)/yfit_bg > fit_flux_continuum_fraction


    # Debugging: plot diagnostic properties.

    if debug and not contall:
        pass
        # Could include this...or not

    return (m.params[2], ok_fit, m.params[4], error)


def d11_filter(i, offset, dwave, spec, data, axis_s=1, ix=None, iy=None,
               mask=None, emask=None, inmsg="", nwidth=1, dwidth=1,
               verbose=0, error=0, debug=False):

    # Calculate the contribution to the flux on the feature blue side.

    sb_0 = i - (offset+dwave)
    sb_1 = i - offset + 1
    (sb_0, sb_1, sb_str, sb_type) = d11_spec_sec(sb_0, sb_1, mask,
                                                 emask=emask, pos=False)

    nb = sb_1 - sb_0 if sb_1 > sb_0 else 0
    if nb > 0: I_blue = spec[sb_0 : sb_1].sum()


    # Calculate the contribution to the flux on the feature red side.

    sr_0 = i + offset
    sr_1 = i + (offset+dwave) + 1

    (sr_0, sr_1, sr_str, sr_type) = d11_spec_sec(sr_0, sr_1, mask, emask=emask)

    nr = sr_1 - sr_0 if sr_1 > sr_0 else 0
    if nr > 0: I_red = spec[sr_0 : sr_1].sum()


    # Sum up all flux in the blue and red bands around the feature.

    if ix is None and iy is None:

        if axis_s == 1:
            if nb > 0: img_blue = np.sum(data[:, :, sb_0 : sb_1], axis=2)
            if nr > 0: img_red  = np.sum(data[:, :, sr_0 : sr_1], axis=2)
            img_i = data[:, :, i]
        elif axis_s == 2:
            if nb > 0: img_blue = np.sum(data[:, sb_0 : sb_1, :], axis=1)
            if nr > 0: img_red  = np.sum(data[:, sr_0 : sr_1, :], axis=1)
            img_i = data[:, i, :]
        else:
            if nb > 0: img_blue = np.sum(data[sb_0 : sb_1, :, :], axis=0)
            if nr > 0: img_red  = np.sum(data[sr_0 : sr_1, :, :], axis=0)
            img_i = data[i, :, :]

    else:

        if axis_s == 1:
            if nb > 0: img_blue = np.sum(data[iy, ix, sb_0 : sb_1], axis=2)
            if nr > 0: img_red  = np.sum(data[iy, ix, sr_0 : sr_1], axis=2)
            img_i = data[iy, ix, i]
        elif axis_s == 2:
            if nb > 0: img_blue = np.sum(data[iy, sb_0 : sb_1, ix], axis=1)
            if nr > 0: img_red  = np.sum(data[iy, sr_0 : sr_1, ix], axis=1)
            img_i = data[iy, ix, :]
        else:
            if nb > 0: img_blue = np.sum(data[sb_0 : sb_1, iy, ix], axis=0)
            if nr > 0: img_red  = np.sum(data[sr_0 : sr_1, iy, ix], axis=0)
            img_i = data[i, iy, ix]


    # Correct the layer i flux by subtracting blue and red contrib.

    if nb == 0 and nr == 0:

        img = img_i  # No correction possible
        msg1 = "not corrected."
        msg2 = "red band: -, blue band: -."
        msg3 = ""

    elif nb == 0:

        if np.isnan(img_red).all():
            img = img_i
        else:
            corr_fac = spec[i]*nr/I_red
            img = img_i - corr_fac*img_red/nr
        msg1 = "corrected using red band."
        msg2 = "blue band: {0:{width}}:{1:{width}} [px] " \
            "(n = {2:{dwidth}}), red band: {3:{width}}:{4:{width}}" \
            " [px] (n = {5:{dwidth}})".format("-", "-", 0, \
                                              sr_0 + 1, sr_1, sr_1 - sr_0, \
                                              width=nwidth, dwidth=dwidth)
        msg3 = "."
        if sr_str != "":
            msg3 = " :: Line in bandpass: {0}.".format(sr_str)
            if sr_type & 1 == 1: msg3 += " {tl.}"
            if sr_type & 2 == 2: msg3 += " {em.}"
            msg3 += "."

    elif nr == 0:

        if np.isnan(img_blue).all():
            img = img_i
        else:
            corr_fac = spec[i]*nb/I_blue
            img = img_i - corr_fac*img_blue/nb
        msg1 = "corrected using blue band."
        msg2 = "blue band: {0:{width}}:{1:{width}} [px] " \
            "(n = {2:{dwidth}}), red band: {3:{width}}:{4:{width}}" \
            " [px] (n = {5:{dwidth}})".format(sb_0 + 1, sb_1, sb_1 - sb_0, \
                                              "-", "-", 0, \
                                              width=nwidth, dwidth=dwidth)
        msg3 = "."
        if sb_str != "":
            msg3 = " :: Line in bandpass: {0}.".format(sb_str)
            if sb_type & 1 == 1: msg3 += " {tl.}"
            if sb_type & 2 == 2: msg3 += " {em.}"
            msg3 += "."

    else:

        corr_fac = spec[i]/(I_blue/nb + I_red/nr)
        img = img_i - corr_fac*(img_blue/nb + img_red/nr)
        msg1 = "corrected using blue and red bands."
        msg2 = "blue band: {0:{width}}:{1:{width}} [px] " \
            "(n = {2:{dwidth}}), red band: {3:{width}}:{4:{width}}" \
            " [px] (n = {5:{dwidth}})".format(sb_0 + 1, sb_1, sb_1 - sb_0, \
                                              sr_0 + 1, sr_1, sr_1 - sr_0, \
                                              width=nwidth, dwidth=dwidth)
        msg3 = ""
        msgp = ""
        if sb_str != "" or sr_str != "": msg3 = " :: Line in bandpass: "

        msg3b = ""
        if sb_str != "":
            msg3b = "{0}".format(sb_str)
            if sb_type & 1 == 1: msg3b += " {tl.}"
            if sb_type & 2 == 2: msg3b += " {em.}"
            msgp = "."

        msg3r = ""
        if sr_str != "":
            if msg3b != "": msg3b = msg3b + ", "
            msg3r = "{0}".format(sr_str)
            if sr_type & 1 == 1: msg3r += " {tl.}"
            if sr_type & 2 == 2: msg3r += " {em.}"
            msgp = "."

        msg3 = msg3 + msg3b + msg3r + msgp

    if verbose >= 1:
        if verbose == 1:
            log_str = inmsg + msg1
        else:
            log_str = inmsg + msg2 + msg3
        print(log_str)
        logging.info(log_str)

    return img


def d11(filename, x, y, apr, cwidth, ofilename=None, offset=5, wave=None,
        spec=None, emissionlines=None, noemissionlines=None, dwl=1.0,
        vel_z=0.0, fit_intensity_limit=0.0,
        fit_flux_continuum_fraction=0.0, bin=1, limit=0.0,
        telluriclines=None, bwidth=3.0, commentslines=None,
        overwrite=False, verbose=0, debug=False):
    """astro-d11: astronomical spectrum data cube continuum subtraction
    filter

    Apply a Differential Emission Line Filter (DELF) to an astronomical
    spectrum data cube.

    *Background*

    The usual approach to find point sources such as planetary nebulæ
    (PNe) in astronomical observations has been to observe the object
    region using imaging techniques. In that approach, the region is
    observed both on-band and off-band using narrow bandpass filters; a
    comparison between the two images reveals objects such as PNe. Such
    an approach can work with PNe as they emit nearly all their
    intensity in a few emission lines; where the forbidden emission line
    of oxygen, [OIII]5007, is typically the strongest one.

    Astro-d11 (DELF) presents an alternative approach where a data cube
    based on integral-field spectroscopy observations provides means to
    use two very narrow bandpasses near the emission line when
    subracting the background signal. In comparison to the imaging
    approach, the narrow "filters" represented by the bandpasses should
    make it possible to detect fainter objects!

    The algorithm is first described in the paper mentioned in the Links
     section below.


    *Method*

    Two narrow bandpasses, a blue and a red bandpass, are offset from
    the current wavelength (layer) towards bluer (lower) and redder
    (higher) pixels, beginning at an initial offset (//offset//). The
    total width of the red and blue bandpasses is set using the
    parameter //cwidth//; either bandpass is skipped for the bluest
    (lowest) and reddest (highest) pixels. The initially offset
    bandpasses are thereafter shifted away from the layer as needed in
    such a way that telluric and [optionally also] emission lines are
    avoided. Additionally, the subtracted continuum value is normalized
    with a reference spectrum of a pre-selected aperture with few
    emission-line features, using the same bandpasses. The location and
    size of the reference aperture must be set using the parameters
    //x//, //y//, and //apr//.

    The reference spectrum (/rspec/) and its continuum bandpasses
    (/rspec_blue/ and /rspec_red/) are defined with /n_blue/ and /n_red/
    layers in the blue and red bandpasses, respectively. Likewise, the
    flux and continuum bandpasses of each spatial element are defined
    with /img/, /img_blue/, and /img_red/, respectively; using the same
    bandpasses as the reference spectrum! The continuum is then
    subtracted from the input data cube for the current layer /i/ using
    the following equation:

    corr_factor[i] = rspec[i] / \
     ((sum(rspec_blue[i])/n_blue[i] + sum(rspec_red[i])/n_red[i])/2)
    out[i] = img[i] - corr_factor[i] * \
      (sum(img_blue[i])/n_blue[i] + sum(img_red[i])/n_red[i])/2


    *Telluric lines*

    The list of telluric lines is specified using the parameter
    //telluriclines//, which needs to be set to the name of a plain-text
    file where each line contains the wavelength of a telluric line in
    the first column (the unit is Angstrom, Å); the default line list
    file is 'telluric_lines_hires.dat', which is available in the 'data'
    directory. The bandpass width can be adjusted using the parameter
    //bwidth// [Angstrom], where the default value is 3.0 Å.


    *Emission lines*

    The list of emission lines is specified using the parameter
    //emissionlines//, which needs to be set to the name of a plain-text
    file where each line contains the wavelength of an emission line in
    the first column (the unit is Angstrom, Å); a default line list file
    is provided in 'emission_lines-ground_based-noFe.dat', which is also
    available in the 'data' directory.

    The procedure is to create a spatially dependent emission-line mask
    by looping through all spatial elements and emission-line entries.
    For this purpose, and to save execution time, the data can be binned
    on the spatial axes to create spectra with higher signal-to-noise
    before the fitting. See the parameter //bin//.

    The emission line redshift can be set using the parameter //vel_z//
    (unit km/s; default is 0 km/s), and an additional permitted offset
    is specified using the parameter //dwl// (unit Angstrom; default is
    1.0 Å). For each spatial element and emission line, a section of the
    object spectrum is fitted using the tool __mpfit.py__ (see link
    below). A fitted line results in the bandpass centered on the
    wavelength to be masked. The emission line bandpass width is set
    using the parameter //bwidth// [Angstrom], where the default width
    is 3.0 Å.

    Please Note! The fitting procedure of individual emission lines is
    slow. So it might be a wise idea to begin with a small number of
    emission lines in the list to see that everything works properly
    before increasing the number.


    *Resulting image*

    The filtered image is written to a file, adding a set of header
    keywords that indicate waht argument values were used (//d11_x//,
    //d11_y//, //d11_apr//, and //d11_cwid//) for the parameters
    //x, y, apr, cwidth//. The output filename can be set explicitly
    using  the parameter //ofilename), otherwise the input filename is
    used with the added suffix '_d11'.


    *Links*

    The filter is described in the paper "Toward Precision Cosmology
    with Improved PNLF Distances Using VLT-MUSE I. Methodology and
    Tests", Martin M. Roth, George H. Jacoby, Robin Ciardullo,
    Brian D. Davis, Owen Chase, and Peter M. Weilbacher 2021,
    The Astrophysical Journal (ApJ), 916, 21, 44 pp.

      https://iopscience.iop.org/journal/0004-637X

    PDF file:
      https://ui.adsabs.harvard.edu/link_gateway/2021ApJ...916...21R/PUB_PDF

    ApJ abstract page:
      https://www.doi.org/10.3847/1538-4357/ac02ca

    NASA ADS:
      https://ui.adsabs.harvard.edu/abs/2021ApJ...916...21R/abstract


    This tool is also available in the integral-field spectroscopy
    data-reduction package p3d, which is available at

      https://p3d.sourceforge.io

    where the DELF-tool is named *p3d_d11*. While p3d is written
    using the Interactive Data Language (IDL), it can be used without a
    license using the IDL Virtual Machine.


    *Calling sequence*

    The program is used with the following keywords and options:

    d11.py <file> x y apr cwidth [-f] [-e <file>] [-d <value>] \
        [-z <value>] [-t <file>] [-b <value>] [-u <char>] \
        [-o <file>] [-w] [-v <int>]

    <file>:
      The name of the data cube file. The file needs to be stored using
      the FITS format. An attempt is made at locating the dispersion
      axis in the data cube using the CTYPEx header keywords (x is an
      integer in the range 1-3), which needs to be set to either AWAV
      or WAVE.

    x:
      The x-coordinate of the region that is used to create a reference
      spectrum. The value is specified in pixel units. There is no
      default as this value has to be chosen by identifying a region in
      the data cube where there is little change.

    y:
      The y-coordinate of the region that is used to create a reference
      spectrum. The value is specified in pixel units. There is no
      default as this value has to be chosen by identifying a region in
      the data cube where there is little change.

    apr:
      The reference spectrum region half-width. The value is specified
      in pixel units. There is no default as this value has to be chosen
      by identifying a region in the data cube where there is little
      change.

    cwidth:
      The full (band)width of the continuum region that includes both
      the region towards lower and higher pixels away from the current
      layer on the dispersion axis. The value is specified in wavelength
      units (Angstrom). There is no default as this value has to be
      chosen depending on the data.

    offset <value>:
      The initial offset towards lower and higher pixels when defining
      the continuum image is set using this keyword. The unit is pixels,
      and the default value is 5 pixels.

    emissionlines <string>:
      The name of a plain-text file that lists wavelengths of
      possibly redshifted emission lines that should be excluded in the
      calculation of the continuum regions. Red and blue shifted
      emission lines are identified by fitting a Gaussian profile to a
      potentially existing emission line, assuming a redshift set using
      the parameter //vel_z// while allowing an offset from that value
      (//dwl//). The wavelength unit is Angstrom. The default value is:
      'data/emission_lines-ground_based-noFe.dat'.

    noemissionlines:
      Do not use the default emission-line file and do not fit any
      emission lines.

    dwl <value>:
      A value that defines a maximum allowed deviation from the provided
      center wavelengths of emission lines. The unit is Ångström [Å].

    vel_z <value>:
      A scalar value that specifies the redshift of all regular emission
      lines as a velocity. The unit assumed is km/s. The redshift is
      recovered using the equation (where c is the light speed):
        z = sqrt((1 + vel_z / c) / (1 - vel_z / c)) - 1

    fit_intensity_limit:
      A scalar decimal value that defines a lower limit value on the
      fitted emission line intensity for the fit to be considered OK.

    fit_flux_continuum_fraction:
      A scalar decimal value that defines a lower limit value on the
      ratio between the emission line flux at the line center and the
      continuum for the fit to be considered OK.

    bin:
      A scalar integer that is used to bin the input data on the spatial
      axes before fitting emission lines. For example, if bin is set to
      10, the spectra of ten spatial elements are summed together on
      both spatial axes to form a binned spectrum out of 100 unbinned
      spectra. In the case that the number of spatial elements on either
      axis is not evenly divisible with the bin number, an additional
      bin is added that contains as many spatial elements, counting from
      the back.

    telluriclines <string>:
      The name of a plain-text file that lists wavelengths of telluric
      lines that should be excluded in the calculation of the
      continuum regions. The wavelength unit is Angstrom. The default
      value is: 'data/telluric_lines_hires.dat'.

    bwidth <value>:
      The bandwidth of bandpasses to ignore centered on telluric and
      emission lines. The value is specified in wavelength units
      (Angstrom). The default value is 3.0 Å.

    commentslines <string>:
      Specify a character that identifies lines with comments in the
      telluric line-list file. The default value is '#'.

      Note! There appears to be a bug in numpy.loadtxt, which is why
        this parameter cannot be used. The only comment character
        accepted is '#'.

    ofilename <string>:
      The name of the resulting filtered file is usually the same as the
      input file, with the added suffix '_d11'. Use this keyword to
      provide an own filename.

    overwrite:
      Any existing file with the same name as the output file will not
      be overwritten unless this keyword is used to overwrite the file.

    verbose <value>:
      An integer that specifies the verbosity of the filter processing.
      The default is to show no information (0). Set the verbosity to
      1 (some information) or 2 (all information).

    help
      Show this information and exit.


    Here is an *example* of how this tool is launched from the shell or
    the console:
    $ d11.py -o datacube_d11.fits -v 1 datacube.fits 12.0 25.0 50.0 20.0
    """

    import os
    import sys
    from pathlib import Path
    import inspect
    import math
    from astropy.io import fits
    import time

    screxe = os.path.basename(__file__)
    screxe = screxe[0:screxe.find(".")] + ": "

    # Check that the filename argument contains the name of an existing
    # FITS file with a data cube.

    if not isinstance(filename, str):
        msg = screxe + "<filename>,  the first argument, must be a " \
            "string; " + str(type(filename)) + "."
        raise RuntimeError(msg)

    if not os.path.isfile(filename):
        msg = screxe + "<filename> must contain the name of an existing FITS" \
            "-type file with a data cube."
        raise RuntimeError(msg)

    if ofilename is None:
        if filename.endswith(".fits"):
            idx = filename.rfind(".")
            ofilename = filename[0:idx] + "_d11.fits"
        else:
            ofilename = filename + "_d11.fits"

    if not isinstance(ofilename, str):
        msg = screxe + "<ofilename> must be a string; " + \
            str(type(ofilename)) + "."
        raise RuntimeError(msg)

    if os.path.exists(ofilename) and not overwrite:
        msg = screxe + "The output file already exists. Remove it first or u" \
            "se the --overwrite option."
        print(screxe + "  ofilename=\"" + ofilename + "\"")
        raise RuntimeError(msg)

    comments = "#"
    if commentslines is not None:
        if not isinstance(commentslines, str):
            msg = screxe + "<commentslines> must be a one-character string; " \
                + str(type(commentslines)) + "."
            raise RuntimeError(msg)
        if len(commentslines) != 1:
            msg = screxe + "<commentslines> must be a one-character string; " \
                + str(type(commentslines)) + "."
            raise RuntimeError(msg)
        comments = commentslines

    use_emissionlines = False
    skip_emissionlines = False
    if noemissionlines is not None:
        if noemissionlines: skip_emissionlines = True

    if not skip_emissionlines:
        if emissionlines is not None and isinstance(emissionlines, str):
            if not os.path.isfile(emissionlines):
                msg = screxe + "<emissionlines> must contain the name of an " \
                    "existing plain-text file listing emission lines (Angstr" \
                    "om)."
                raise RuntimeError(msg)

            use_emissionlines = True
        else:
            exefile = inspect.getabsfile(inspect.currentframe())
            path = Path(exefile)
            del exefile
            path = path.parent.parent.parent.absolute()
            emissionlines = os.path.join(path, "data",
                                     "emission_lines-ground_based-noFe.dat")
            use_emissionlines = True
            del path
    del skip_emissionlines

    if use_emissionlines:
        elines = np.loadtxt(emissionlines, comments=comments, usecols=(0))
        print(elines)
        if not isinstance(dwl, float):
            msg = screxe + "<dwl> must be set to a decimal value (Angstrom)."
            raise RuntimeError(msg)

        if not isinstance(vel_z, float):
            msg = screxe + "<vel_z> must be set to a decimal value (km/s)."
            raise RuntimeError(msg)

        clight = 2.99792458e10
        vel_z *= 1e5  # km/s => cm/s
        z = np.sqrt((1.0 + vel_z / clight) / (1.0 - vel_z / clight)) - 1.0

        if not isinstance(fit_intensity_limit, float):
            msg = screxe + "<fit_intensity_limit> must be set to a decimal v" \
                "alue; fit_intensity_limit >= 0."
            raise RuntimeError(msg)

        if fit_intensity_limit < 0.0:
            msg = screxe + "<fit_intensity_limit> must be set to a decimal v" \
                "alue; fit_intensity_limit >= 0."
            raise RuntimeError(msg)

        if not isinstance(fit_flux_continuum_fraction, float):
            msg = screxe + "<fit_flux_continuum_fraction> must be set to a d" \
                "ecimal value; fit_flux_continuum_fraction >= 0."
            raise RuntimeError(msg)

        if fit_flux_continuum_fraction < 0.0:
            msg = screxe + "<fit_flux_continuum_fraction> must be set to a d" \
                "ecimal value; fit_flux_continuum_fraction >= 0."
            raise RuntimeError(msg)

        if not isinstance(bin, int):
            msg = screxe + "<bin> must be set to an integer; bin >= 1."
            raise RuntimeError(msg)

        if bin <= 0:
            msg = screxe + "<bin> must be set to an integer; bin >= 1."
            raise RuntimeError(msg)

    use_telluriclines = False
    if telluriclines is not None and isinstance(telluriclines, str):
        if not os.path.isfile(telluriclines):
            msg = screxe + "<telluriclines> must contain the name of an exis" \
                "ting plain-text file listing telluric lines (Angstrom)."
            raise RuntimeError(msg)

        use_telluriclines = True
    else:
        exefile = inspect.getabsfile(inspect.currentframe())
        path = Path(exefile)
        del exefile
        path = path.parent.parent.parent.absolute()
        telluriclines = os.path.join(path, "data", "telluric_lines_hires.dat")
        use_telluriclines = True
        del path

    if use_telluriclines:
        tlines = np.loadtxt(telluriclines, comments=comments, usecols=(0))


    if verbose >= 1:
        idx = filename.rfind(".")
        logfile = filename[0:idx] + "_d11.log"

        if os.path.exists(logfile) and not overwrite:
            msg = screxe + "The log file already exists. New information wil" \
                "l be appended at the end of the file unless the " \
                "--overwrite option is used."
            print(screxe + "  logfile=\"" + logfile + "\"")
            raise RuntimeError(msg)
        if os.path.exists(ofilename) and overwrite: os.remove(logfile)

        logging.basicConfig(filename=logfile, format="%(asctime)s %(message)s",
                            level=logging.DEBUG)
        log_str = [screxe + "Apply a differential emission line filter (DELF" \
                   ") on an astronomical data cube.", \
                   screxe + "The tool was run using the following options:", \
                   screxe + "      filename = \"" + filename + \
                   "\"", \
                   screxe + "             x = " + str(x), \
                   screxe + "             y = " + str(y), \
                   screxe + "           apr = " + str(apr), \
                   screxe + "        cwidth = " + str(cwidth), \
                   screxe + "        bwidth = " + str(bwidth), \
                   screxe + "        offset = " + str(offset)]
        if use_telluriclines or use_emissionlines:
            log_str.append(screxe + " commentslines = \"" + \
                           comments + "\"")
        if use_emissionlines:
            log_str = [log_str, \
                       screxe + " emissionlines = \"" + emissionlines + \
                       "\"", \
                       screxe + "           dwl = " + str(dwl), \
                       screxe + "         vel_z = " + str(vel_z), \
                       screxe + "         fit_intensity_limit = " + \
                       str(fit_intensity_limit), \
                       screxe + " fit_flux_continuum_fraction = " + \
                       str(fit_flux_continuum_fraction), \
                       screxe + "           bin = " + str(bin)]
        if use_telluriclines:
            log_str.append(screxe + \
                           " telluriclines = \"" + telluriclines + "\"")
        log_str.append(screxe + "     ofilename = \"" + ofilename + "\"")
        for log_str_i in log_str:
            print(log_str_i)
            logging.info(log_str_i)

    try:
        with fits.open(filename) as hdul:
            pass
    except OSError:
        log_str = screxe + "<filename> must contain the name of an existing " \
            "FITS-type file with a data cube."
        print(log_str)
        if verbose >= 1: logging.error(log_str)
        raise RuntimeError(log_str)

    with fits.open(filename) as hdul:
        hdr0 = hdul[0].header
        hdr1 = hdul[1].header


        # Retrieve required values from the header.

        if hdr1["naxis"] != 3:
            log_str = screxe + "<filename> must contain the name of an exist" \
                "ing FITS-type file with a data cube."
            print(log_str)
            if verbose >= 1: logging.error(log_str)
            raise RuntimeError(log_str)

        ctype1 = hdr1["ctype1"]
        ctype2 = hdr1["ctype2"]
        ctype3 = hdr1["ctype3"]
        if ctype1 == "AWAV" or ctype1 == "WAVE":
            axis_s = 1
            axis_x = 2
            axis_y = 3
            cdisp = hdr1["cd1_1"]
        elif ctype2 == "AWAV" or ctype2 == "WAVE":
            axis_x = 1
            axis_s = 2
            axis_y = 3
            cdisp = hdr1["cd2_2"]
        else:
            axis_x = 1
            axis_y = 2
            axis_s = 3
            cdisp = hdr1["cd3_3"]

        axis_ss = str(axis_s)    
        if verbose >= 2:
            logging.info(screxe + "Dispersion axis is " + axis_ss)

        nwave = hdr1["naxis" + axis_ss]
        crval = hdr1["crval" + axis_ss]
        crpix = hdr1["crpix" + axis_ss]

        if verbose >= 2:
            logging.info(screxe + "NAXIS" + axis_ss + "=" + str(nwave) +
                         ", CRVAL" + axis_ss + "=" + str(crval) +
                         ", CDELT" + axis_ss + "=" + str(cdisp) +
                         ", CRPIX" + axis_ss + "=" + str(crpix))

        xsize = hdr1["naxis" + str(axis_x)]
        ysize = hdr1["naxis" + str(axis_y)]

        bxsize = xsize // bin
        if xsize % bin != 0: bxsize = bxsize + 1
        bysize = ysize // bin
        if ysize % bin != 0: bysize = bysize + 1

        n1width = max([len(str(xsize)), len(str(ysize))])
        nwidth = max([len(str(bxsize)), len(str(bysize))])


        # Setup a wavelength array.

        wave = cdisp*(np.arange(nwave, dtype=float) - crpix + 1.0) + crval
        idx = np.arange(nwave)

        dwave = round(0.5*cwidth/cdisp)


        # Create a line mask that optionally accounts for telluric lines.

        mask = np.zeros(nwave)
        if use_telluriclines:
            twave = round(0.5*bwidth/cdisp)
            for wave_tell in tlines:
                w_i = (wave_tell-crval)/cdisp - 1.0 + crpix

                w_i__low = math.floor(w_i - bwidth)
                if w_i__low < 0: w_i__low = 0
                if w_i__low > nwave - 1: w_i__low = nwave - 1

                w_i__high = math.ceil(w_i + bwidth) + 1
                if w_i__high < 1: w_i__high = 1
                if w_i__high > nwave: w_i__high = nwave

                if w_i__high - w_i__low > 1:
                    mask[w_i__low : w_i__high] = wave_tell

        #bin_start = dwave + offset
        #bin_final = nwave - (dwave+offset)


        # Read the data; this could take quite some time.

        if verbose >= 2:
            log_str = screxe + "Read the FITS file data block: start"
            print(log_str)
            logging.info(log_str)

        data = hdul[1].data

        if verbose >= 2:
            log_str = screxe + "Read the FITS file data block: done"
            print(log_str)
            logging.info(log_str)

        odata = np.empty_like(data)


        # Set the aperture and check that it is consistent.

        x_0 = round(x - apr)
        x_1 = round(x + apr) + 1
        if x_0 < 0: x_0 = 0
        if x_1 > xsize: x_1 = xsize
        nx = 0
        if x_1 > x_0: nx = x_1 - x_0

        y_0 = round(y - apr)
        y_1 = round(y + apr) + 1
        if y_0 > 0: y_0 = 0
        if y_1 > ysize: y_1 = ysize
        ny = 0
        if y_1 > y_0: ny = y_1 - y_0

        if verbose >= 1:
            log_str = screxe + "Aperture: [" + str(x_0 + 1) + ":" + \
                str(x_1) + ", " + str(y_0 + 1) + ":" + str(y_1) + \
                "] [px, px] :: red and blue bandwidth: " + str(dwave + 1) + \
                " [px]."
            print(log_str)
            logging.info(log_str)
            del log_str


        #========================================------------------------------
        # Optionally, create a line mask that accounts for redshifted emission
        # lines.

        contall = 0
        if use_emissionlines:

            xarr = np.arange(nwave, dtype=float)
            e_mask = np.zeros((nwave, xsize, ysize))

            twave = round(0.5*bwidth/cdisp)

            # Loop over all spatial elements.
            for ixy in range(0, bxsize * bysize):

                bix = ixy % bxsize
                biy = ixy // bxsize

                if bin == 1:

                    # Data are not binned.

                    if axis_s == 1:
                        xy_spec = data[biy, bix, :]
                    elif axis_s == 2:
                        xy_spec = data[biy, :, bix]
                    else:
                        xy_spec = data[:, biy, bix]

                else:

                    # Data are binned.

                    if bix == bxsize - 1:
                        ix_1 = xsize
                        ix_0 = ix_1 - bin
                    else:
                        ix_0 = bix * bin
                        ix_1 = ix_0 + bin

                    if biy == bysize - 1:
                        iy_1 = ysize
                        iy_0 = iy_1 - bin
                    else:
                        iy_0 = biy * bin
                        iy_1 = iy_0 + bin

                    nq = (ix_1 - ix_0) * (iy_1 - iy_0)

                    if axis_s == 1:
                        xy_spec = np.sum(data[iy_0 : iy_1, ix_0 : ix_1, :],
                                         axis=(0, 1))
                    elif axis_s == 2:
                        xy_spec = np.sum(data[iy_0 : iy_1, :, ix_0 : ix_1],
                                         axis=(0, 2))
                    else:
                        xy_spec = np.sum(data[:, iy_0 : iy_1, ix_0 : ix_1],
                                         axis=(1, 2))

                # Check for NaN-element-only spectra and skip them.

                count = np.argwhere(~np.isnan(xy_spec)).size/2

                if count == 0:
                    log_str = screxe + "Spectrum bin [" + \
                        str(bix + 1).rjust(nwidth) + ", " + \
                        str(biy + 1).rjust(nwidth) + "] / [" + \
                                 str(bxsize).rjust(nwidth) + ", " + \
                                 str(bysize).rjust(nwidth) + \
                        "] :: There were no finite pixels at all in the " \
                        "spectrum - skip."
                    print(log_str)
                    logging.info(log_str)
                    del log_str

                    continue


                # Loop over all emission lines.

                for i in range(0, len(elines)):

                    # Initially, use a constant redshift across the field.

                    w_init = (elines[i]*(1.0 + z) - crval)/cdisp - 1.0 + crpix

                    if w_init < - 2.0 or w_init > nwave + 2.0: continue


                    # Locate the possiby redshifted emission line.

                    # Collect all emission lines in the interval:.
                    x__low = math.floor(w_init - 2*bwidth)
                    if x__low < (- 2): x__low = - 2
                    if x__low > nwave + 1: x__low = nwave + 1

                    x__high = math.ceil(w_init + 2*bwidth) + 1
                    if x__high < (- 1): x__high = - 1
                    if x__high > nwave + 2: x__high = nwave + 2

                    w_init_oo = (elines*(1.0+z) - crval)/cdisp - 1.0 + crpix
                    tmp = np.asarray((w_init_oo.astype(int) >= x__low) & \
                                     (w_init_oo.astype(int) <= x__high) & \
                                     (w_init_oo != w_init)).nonzero()
                    if np.array(tmp).size > 0:
                        e_count = np.array(tmp).size
                        w_too = tmp[0][:]
                        del tmp
                    else:
                        e_count = 0
                        w_too = None


                    # The elements next to the emission line must be finite.
                    xi__low = math.floor(w_init - 2)
                    if xi__low < 0: xi__low = 0
                    if xi__low > nwave - 1: xi__low = nwave - 1

                    xi__high = math.ceil(w_init + 2) + 1
                    if xi__high < 1: xi__high = 1
                    if xi__high > nwave: xi__high = nwave

                    x__low = math.floor(w_init - 2*bwidth - 2)
                    if x__low < 0: x__low = 0
                    if x__low > nwave - 1: x__low = nwave - 1

                    x__high = math.ceil(w_init + 2*bwidth + 2) + 1
                    if x__high < 1: x__high = 1
                    if x__high > nwave: x__high = nwave

                    if x__low + 1 >= x__high: continue

                    # Extract the spectrum part around the wavelength.
                    x_sec = xarr[x__low : x__high]
                    spec_sec = xy_spec[x__low : x__high]


                    #==============================--------------------
                    # Check for NaN-elements:

                    count = np.argwhere(~np.isnan(spec_sec)).size//2
                    nan_count = np.argwhere(np.isnan(
                        xy_spec[xi__low : xi__high])).size/2
                    if nan_count > 0 or count < 7:
                        log_st = "no"
                        if count >  0: log_st = "only " + str(count)
                        log_str = screxe + "Spectrum bin [" + \
                            str(bix + 1).rjust(nwidth) + ", " + \
                            str(biy + 1).rjust(nwidth) + "] / [" + \
                            str(bxsize).rjust(nwidth) + ", " + \
                            str(bysize).rjust(nwidth) + "] :: There were " + \
                            log_st + " finite pixels in the spectrum - skip."
                        print(log_str)
                        logging.info(log_str)
                        del log_str

                        continue


                    #if ~z_scalar then begin
                    #   z_use += z_value[kl]
                    #   z_use_n ++
                    #endif


                    #==============================--------------------
                    # Fit the emission line.

                    xstr = str(bix + 1) + ", " + str(biy + 1) + ", " + \
                        str(elines[i])

                    (w_i, ok_fit, ok_intensity, error) = \
                        d11_mpfit(w_init, dwl, cdisp, x=x_sec, y=spec_sec, \
                                  w_too=w_too, \
                                  fit_intensity_limit=fit_intensity_limit, \
                                  fit_flux_continuum_fraction=\
                                  fit_flux_continuum_fraction, \
                                  xstr=xstr, verbose=verbose, \
                                  debug=debug, contall=contall)
                    if error != 0: return

                    if e_count > 0: del w_too

                    if verbose >=3:
                        ok_str = "yes, intensity = " + str(ok_intensity) \
                            if ok_fit else "no"
                        log_str = screxe + "Spectrum bin [" + \
                            str(bix + 1).rjust(nwidth) + ", " + \
                            str(biy + 1).rjust(nwidth) + "] / [" + \
                            str(bxsize).rjust(nwidth) + ", " + \
                            str(bysize).rjust(nwidth) + "] :: Line fit was o" \
                            "k: " + ok_str
                        print(log_str)
                        logging.info(log_str)
                        del log_str

                    if ok_fit:
                        w_i__low = math.floor(w_i - twave)
                        if w_i__low < 0: w_i__low = 0
                        if w_i__low > nwave - 1: w_i__low = nwave - 1

                        w_i__high = math.ceil(w_i + twave) + 1
                        if w_i__high < 1: w_i__high = 1
                        if w_i__high > nwave: w_i__high = nwave

                        if w_i__high - w_i__low > 1:
                            if bin == 1:
                                e_mask[w_i__low : w_i__high, bix, biy] = i
                            else:
                                e_mask[w_i__low : w_i__high, \
                                       ix_0 : ix_1, iy_0 : iy_1] = i



        # Sum the flux in the selected aperture for all layers.

        if axis_s == 1:
            spec = np.nansum(data[y_0 : y_1, x_0 : x_1, :], axis=(0, 1))
        elif axis_s == 2:
            spec = np.nansum(data[y_0 : y_1, :, x_0 : x_1], axis=(0, 2))
        else:
            spec = np.nansum(data[:, y_0 : y_1, x_0 : x_1], axis=(1, 2))


        # Step through the cube on the dispersion axis and subtract scaled
        # blue and red continuum bands from the data cube bins.

        nwidth = len(str(nwave))
        dwidth = len(str(dwave))

        for i in range(0, nwave):

            log_str = screxe + "Layer " + str(i + 1).rjust(nwidth) + " / " \
                + str(nwave) + " :: "
            #log_str = print(tmp.format(i + 1, width=nwidth), end="")


            if use_emissionlines:

                #==============================================================
                #==============================================================
                #==============================================================
                # Loop through each spatial element separately to create an
                # image that depends on both emission lines and telluric lines.
                #==============================================================
                #==============================================================
                #==============================================================

                for ixy in range(0, xsize * ysize):

                    ix = ixy % xsize
                    iy = ixy // xsize

                    log_str__i = log_str + " [" + str(ix+1).rjust(n1width) + \
                        ", " + str(iy+1).rjust(n1width) + "] :: "

                    img = d11_filter(i, offset, dwave, spec, data,
                                     axis_s=axis_s, ix=ix, iy=iy, mask=mask,
                                     emask=e_mask[:, ix, iy], inmsg=log_str__i,
                                     nwidth=nwidth, dwidth=dwidth,
                                     verbose=verbose, debug=debug)

                    if axis_s == 1:
                        odata[iy, ix, i] = img
                    elif axis_s == 2:
                        odata[iy, i, ix] = img
                    else:
                        odata[i, iy, ix] = img


            else:

                #==============================================================
                #==============================================================
                #==============================================================
                # Only consider telluric lines.
                #==============================================================
                #==============================================================
                #==============================================================

                img = d11_filter(i, offset, dwave, spec, data, axis_s=axis_s,
                                 mask=mask, inmsg=log_str, nwidth=nwidth,
                                 dwidth=dwidth, verbose=verbose, debug=debug)

                if axis_s == 1:
                    odata[:, :, i] = img
                elif axis_s == 2:
                    odata[:, i, :] = img
                else:
                    odata[i, :, :] = img






        # Add data processing header entries and write a file with the
        # resulting data.

        hdr1["history"] = "d11: data cube processed with d11 to subtract the" \
            " continuum"
        hdr1["history"] = "d11: time of processing: " + time.asctime()
        hdr1["d11_x"] = (x, "d11: ap. center x coord. for continuum template")
        hdr1["d11_y"] = (y, "d11: ap. center y coord. for continuum template")
        hdr1["d11_apr"] = (apr, "d11: aperture radius to sample the continuu" \
                            "m")
        hdr1["d11_cwid"] = (cwidth, "d11: bandwidth of offband continuum [An" \
                            "gstrom]")

        empty_primary = fits.PrimaryHDU(header=hdr0)
        image_hdu     = fits.ImageHDU(odata, header=hdr1)

        hduw = fits.HDUList([empty_primary, image_hdu])
        if overwrite and os.path.exists(ofilename): os.remove(ofilename)
        hduw.writeto(ofilename)

        if verbose >= 1:
            log_str = screxe + "Wrote resulting data to the file " + ofilename
            print(log_str)
            logging.info(log_str)
            del log_str

        
# Create a launcher in case the program is launched from system shell:
if __name__ == "__main__":
    from argparse import ArgumentParser, RawDescriptionHelpFormatter
 
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                            description=__doc__)

    ###########################################################################
    # Parsing command-line arguments and options:
    ###########################################################################

    parser.add_argument("filename", help="The name of a FITS file with a dat" \
                        "a cube; with two spatial and one spectral dimension.")
    parser.add_argument("x", help="Reference region aperture x coordinate " \
                        "[pixel].", type=float)
    parser.add_argument("y", help="Reference region aperture y coordinate " \
                        "[pixel].", type=float)
    parser.add_argument("apr", help="Reference region aperture size [pixel].",
                        type=float)
    parser.add_argument("cwidth", help="Total (blue + red) continuum bandwid" \
                        "th for subtraction [Angstrom].", type=float)

    parser.add_argument("-f", "--offset", action="store", type=int, \
                        help="Specifies the (initial) offset of the red and " \
                        "blue continuum regions away from the current layer " \
                        "(wavelength) [pixel].")
    parser.add_argument("-e", "--emissionlines", action="store", type=str, \
                        help="Specifies the name of a plain-text file listin" \
                        "g [possibly] redshifted emission lines [Angstrom].")
    parser.add_argument("-n", "--noemissionlines", action="store_true", \
                        help="Do not fit any emission lines.")
    parser.add_argument("-d", "--dwl", action="store", type=float, \
                        help="A scalar value that specifies the allowed devi" \
                        "ation of each fitted line from specified line cente" \
                        "r wavelengths [Angstrom].")
    parser.add_argument("-z", "--vel_z", action="store", type=float, \
                        help="A scalar value that specifies the redshift of " \
                        "emission lines [km/s].")
    parser.add_argument("-l", "--fit_intensity_limit", action="store", \
                        type=float, help="A scalar value that specifies the " \
                        "lower limit of an acceptable fitted emission line i" \
                        "nensity.")
    parser.add_argument("-r", "--fit_flux_continuum_fraction", \
                        action="store", type=float, help="A scalar value tha" \
                        "t defines a lower limit on the ratio between the em" \
                        "ission line flux at the line center and the continu" \
                        "um for the fit to be considered OK")
    parser.add_argument("-g", "--bin", action="store", type=int, \
                        help="A scalar integer that bins this many elements " \
                        "on both spatial axes in one bin before fitting emis" \
                        "sion lines.")
    parser.add_argument("-t", "--telluriclines", action="store", type=str, \
                        help="Specifies the name of a plain-text file listin" \
                        "g telluric lines [Angstrom].")
    parser.add_argument("-q", "--bwidth", action="store", type=float, \
                        help="The width of telluric and emission line bandpa" \
                        "sses [Angstrom].")
    parser.add_argument("-u", "--commentslines", action="store", type=str, \
                        help="Specifies a comment character to use with the " \
                        "plain-text telluric and emission lines file [defaul" \
                        "t: '#'].")
    parser.add_argument("-o", "--ofilename", action="store", type=str, \
                        help="The output file name.")
    parser.add_argument("-w", "--overwrite", action="store_true", \
                        help="Overwrite existing output files.")
    parser.add_argument("-v", "--verbose", action="store", type=int, \
                        help="Be verbose on what is done; valid values are: " \
                        "1, 2, 3, and 4.")
    parser.add_argument("--debug", action="store_true", help="Debugging mode.")

    args = parser.parse_args()

    if args.offset is not None:
        offset = args.offset
    else:
        offset = 5

    if args.dwl is not None:
        dwl = args.dwl
    else:
        dwl = 1.0

    if args.vel_z is not None:
        vel_z = args.vel_z
    else:
        vel_z = 0.0

    if args.fit_intensity_limit is not None:
        fit_intensity_limit = args.fit_intensity_limit
    else:
        fit_intensity_limit = 0.0

    if args.fit_flux_continuum_fraction is not None:
        fit_flux_continuum_fraction = args.fit_flux_continuum_fraction
    else:
        fit_flux_continuum_fraction = 0.0

    if args.bin is not None:
        bin = args.bin
    else:
        bin = 1

    if args.bwidth is not None:
        bwidth = args.bwidth
    else:
        bwidth = 3.0

    if args.verbose is not None:
        verbose = args.verbose
    else:
        verbose = 0

    d11(args.filename, args.x, args.y, args.apr, args.cwidth,
        offset=offset, emissionlines=args.emissionlines,
        noemissionlines=args.noemissionlines, dwl=dwl, vel_z=vel_z,
        fit_intensity_limit=fit_intensity_limit,
        fit_flux_continuum_fraction=fit_flux_continuum_fraction,
        bin=bin, telluriclines=args.telluriclines, bwidth=bwidth,
        commentslines=args.commentslines, ofilename=args.ofilename,
        overwrite=args.overwrite, verbose=verbose, debug=args.debug)
