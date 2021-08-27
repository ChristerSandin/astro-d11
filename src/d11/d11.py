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

def d11_spec_sec(i_0, i_1, mask, pos=True):

    N = len(mask)
    i_str = ""

    # Step over masked pixels:

    if pos:
        if i_0 > N: i_0 = N
        if i_1 > N: i_1 = N
        if i_0 == i_1:
            return (i_0, i_1, "")

        while max(mask[i_0 : i_1]) > 0:
            i_str = str(mask[i_0]) + " (red)"
            i_0 = i_0 + 1
            i_1 = i_1 + 1

            if i_0 > N: i_0 = N
            if i_1 > N: i_1 = N
            if i_0 == i_1:
                return (i_0, i_1, "")
    else:
        if i_0 < 0: i_0 = 0
        if i_1 < 0: i_1 = 0
        if i_0 == i_1:
            return (i_0, i_1, "")

        while max(mask[i_0 : i_1]) > 0:
            i_str = str(mask[i_1 - 1]) + " (blue)"
            i_0 = i_0 - 1
            i_1 = i_1 - 1

            if i_0 < 0: i_0 = 0
            if i_1 < 0: i_1 = 0
            if i_0 == i_1:
                return (i_0, i_1, "")

    return (i_0, i_1, i_str)


def d11(filename, x, y, apr, cwidth, ofilename='DEF',
        offset=5, wave=None, spec=None, telluric=None, twidth=3.0,
        commentstelluric=None, overwrite=False, verbose=0):
    """astro-d11: astronomical spectrum datacube continuum subtraction
    filter

    The filter works as follows.

    A datacube with astronomical spectra is assumed to be provided in a
    file that is saved in the FITS format. Where two axes are spatial
    and a third is spectral.

    The user specifies a quadratic region where there is little
    disturbance and few features; this is done using the first two
    arguments (x and y [pixel] coordinates) and the half-side (apr, also
    a pixel value). The flux in the region is summed up to create a
    reference spectrum using all layers (pixels on the wavelength axis)
    in the datacube.

    A new datacube is constructed layer by layer, where a continuum
    image is created for each layer individually based on two regions
    offset (--offset) towards lower and higher pixels of the current
    layer. The total width of the continuum regions is set using the
    fourth argument (cwidth, the width is specified in Angstrom). The
    continuum image is at first normalized (corrected) with the
    corresponding continuum regions of the reference spectrum.
    Thereafter, the normalized continuum image is subtracted from the
    current layer image.

    Creating the continuum image and continuum reference spectrum for
    each layer, the two continuum regions are offset towards either
    lower or higher pixels until no telluric line falls within the
    respective region. Telluric line wavelengths are read from a
    plain-text file (--telluric) that provides a wavelength value per
    line. The full width of any individual telluric line region can be
    set (--twidth).

    The filtered image is written to a file, adding a set of header
    keywords that indicate waht argument values were used ("d11_x",
    "d11_y", "d11_apr", and "d11_cwid") for (x, y, apr, cwidth). The
    output filename can be set explicitly (using --ofilename), but
    otherwise the input filename is used with the added suffix "_d11".

    The program is used with the following keywords and options:

    d11.py <file> x y apr cwidth [-d] [-t <file>] [-q <value>] \
        [-u <char>] [-o <file>] [-w] [-v <int>]

    <file>:
      The name of the datacube file. The file needs to be stored using
      the FITS format. An attempt is made at locating the dispersion
      axis in the datacube using the CTYPEx header keywords (x is an
      integer in the range 1-3), which needs to be set to either AWAV
      or WAVE.

    x:
      The x-coordinate of the region that is used to create a reference
      spectrum. The value is specified in pixel units. There is no
      default as this value has to be chosen by identifying a region in
      the datacube where there is little change.

    y:
      The y-coordinate of the region that is used to create a reference
      spectrum. The value is specified in pixel units. There is no
      default as this value has to be chosen by identifying a region in
      the datacube where there is little change.

    apr:
      The reference spectrum region half-width. The value is specified
      in pixel units. There is no default as this value has to be chosen
      by identifying a region in the datacube where there is little
      change.

    cwidth:
      The full (band)width of the continuum region that includes both
      the region towards lower and higher pixels away from the current
      layer on the dispersion axis. The value is specified in wavelength
      units (Angstrom). There is no default as this value has to be
      chosen depending on the data.

    [-d] or [--offset]:
      The initial offset towards lower and higher pixels when defining
      the continuum image is set using this keyword. The unit is pixels,
      and the default value is 5 pixels.

    [-t] or [--telluric]:
      The name of a plain-text file that lists the wavelength of
      telluric lines that should be excluded in the calculation of the
      continuum regions. The wavelength unit is Angstrom. The default
      value is: "data/telluric_lines_hires.dat".

    [-q] or [--twidth]:
      The bandwidth of bandpasses to ignore centered on telluric lines.
      The value is specified in wavelength units (Angstrom). The default
      value is 3.0 Å.

    [-u] or [--commentstelluric]:
      Specify a character that identifies lines with comments in the
      telluric line-list file. The default value is "#".

    [-o] or [--ofilename]:
      The name of the resulting filtered file is usually the same as the
      input file, with the added suffix "_d11". Use this keyword to
      provide an own filename.

    [-w] or [--overwrite]:
      Any existing file with the same name as the output file will not
      be overwritten unless this keyword is used to overwrite the file.

    [-v] or [--verbose]:
      An integer that specifies the verbosity of the filter processing.
      The default is to show no information (0). Set the verbosity to
      1 (some information) or 2 (all information).

    [-h] or [--help]:
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
    import numpy as np
    from astropy.io import fits
    import time
    import logging

    screxe = os.path.basename(__file__)
    screxe = screxe[0:screxe.find(".")] + ": "

    # Check that the filename argument contains the name of an existing
    # FITS file with a datacube.

    if not isinstance(filename, str):
        msg = screxe + "<filename>,  the first argument, must be a " \
            "string; " + str(type(filename)) + "."
        raise RuntimeError(msg)

    if not os.path.isfile(filename):
        msg = screxe + "<filename> must contain the name of an existing FITS" \
            "-type file with a datacube."
        raise RuntimeError(msg)

    if ofilename == "DEF" or ofilename is None:
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

    use_telluric = False
    if telluric is not None and isinstance(telluric, str):
        if not os.path.isfile(telluric):
            msg = screxe + "<telluric> must contain the name of an existing " \
                "plain-text file listing telluric lines (Angstrom)."
            raise RuntimeError(msg)

        use_telluric = True
    else:
        exefile = inspect.getabsfile(inspect.currentframe())
        path = Path(exefile)
        path = path.parent.parent.parent.absolute()
        telluric = os.path.join(path, "data", "telluric_lines_hires.dat")
        use_telluric = True

    comments = "#"
    if commentstelluric is not None:
        if not isinstance(commentstelluric, str):
            msg = screxe + "<commentstelluric> must be a one-character s" \
                "tring; " + str(type(commentstelluric)) + "."
            raise RuntimeError(msg)
        if len(commentstelluric) != 1:
            msg = screxe + "<commentstelluric> must be a one-character s" \
                "tring; " + str(type(commentstelluric)) + "."
            raise RuntimeError(msg)
        comments = commentstelluric

    if use_telluric:
        tlines = np.loadtxt(telluric, comments=comments)


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

        logging.basicConfig(filename=logfile,
                            format="%(asctime)s %(message)s",
                            level=logging.DEBUG)
        log_str = [screxe + "Apply a differential emission line filter (DELF" \
                   ") on an astronomical datacube.", \
                   screxe + "The tool was run using the following options:", \
                   screxe + "  filename = \"" + filename + "\"", \
                   screxe + "  x = " + str(x), \
                   screxe + "  y = " + str(y), \
                   screxe + "  apr = " + str(apr), \
                   screxe + "  cwidth = " + str(cwidth), \
                   screxe + "  twidth = " + str(twidth), \
                   screxe + "  offset = " + str(offset)]
        if use_telluric:
            log_str = [log_str, screxe + "  commentstelluric = \"" +
                       comments + "\"", \
                       screxe + "  telluric = \"" + telluric + "\""]
        log_str = [log_str, screxe + "  ofilename = \"" + ofilename + "\""]
        for log_str_i in log_str: logging.info(log_str_i)

    try:
        with fits.open(filename) as hdul:
            pass
    except OSError:
        log_str = screxe + "<filename> must contain the name of an existing " \
            "FITS-type file with a datacube."
        print(log_str)
        if verbose >= 1: logging.error(log_str)
        raise RuntimeError(log_str)

    with fits.open(filename) as hdul:
        hdr0 = hdul[0].header
        hdr1 = hdul[1].header


        # Retrieve required values from the header.

        if hdr1["naxis"] != 3:
            log_str = screxe + "<filename> must contain the name of an exist" \
                "ing FITS-type file with a datacube."
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
                         ", CRPIX" + axis_ss + "=" + str(crpix))

        xsize = hdr1["naxis" + str(axis_x)]
        ysize = hdr1["naxis" + str(axis_y)]


        # Setup a wavelength array.

        wave = cdisp*(np.arange(nwave, dtype=float) - crpix + 1.0) + crval
        idx = np.arange(nwave)

        dwave = round(0.5*cwidth/cdisp)


        # Create a telluric line mask.

        mask = np.zeros(nwave)
        if use_telluric:
            twave = round(0.5*twidth/cdisp)
            for wave_tell in tlines:
                w_i = (wave_tell-crval)/cdisp - 1.0 + crpix

                w_i__low = math.floor(w_i - twidth)
                if w_i__low < 0: w_i__low = 0
                if w_i__low > nwave - 1: w_i__low = nwave - 1

                w_i__high = math.ceil(w_i + twidth) + 1
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

        # Sum the flux in the selected aperture for all layers.

        if axis_s == 1:
            spec = np.nansum(data[y_0 : y_1, x_0 : x_1, :], axis=(0, 1))
        elif axis_s == 2:
            spec = np.nansum(data[y_0 : y_1, :, x_0 : x_1], axis=(0, 2))
        else:
            spec = np.nansum(data[:, y_0 : y_1, x_0 : x_1], axis=(1, 2))


        # Step through the cube on the dispersion axis and subtract scaled
        # blue and red continuum bands from the datacube bins.

        nwidth = len(str(nwave))
        dwidth = len(str(dwave))

        for i in range(0, nwave):

            if verbose >= 1:
                log_str = screxe + "Layer {0:{width}} / " + str(nwave) + " :: "
                print(log_str.format(i + 1, width=nwidth), end="")
                logging.info(log_str)

            # Calculate the contribution to the flux on the feature blue side.

            sb_0 = i - (offset+dwave)
            sb_1 = i - offset + 1
            (sb_0, sb_1, sb_str) = d11_spec_sec(sb_0, sb_1, mask, pos=False)

            nb = 0
            if sb_1 > sb_0:
                nb = sb_1 - sb_0
                I_blue = spec[sb_0 : sb_1].sum()

            # Calculate the contribution to the flux on the feature red side.

            sr_0 = i + offset
            sr_1 = i + (offset+dwave) + 1

            (sr_0, sr_1, sr_str) = d11_spec_sec(sr_0, sr_1, mask)

            nr = 0
            if sr_1 > sr_0:
                nr = sr_1 - sr_0
                I_red  = spec[sr_0 : sr_1].sum()


            # Sum up all flux in the blue and red bands around the feature.

            if axis_s == 1:
                if nb > 0:
                    img_blue = np.sum(data[:, :, sb_0 : sb_1], axis=2)
                if nr > 0:
                    img_red  = np.sum(data[:, :, sr_0 : sr_1], axis=2)
                img_i = data[:, :, i]
            elif axis_s == 2:
                if nb > 0:
                    img_blue = np.sum(data[:, sb_0 : sb_1, :], axis=1)
                if nr > 0:
                    img_red  = np.sum(data[:, sr_0 : sr_1, :], axis=1)
                img_i = data[:, i, :]
            else:
                if nb > 0:
                    img_blue = np.sum(data[sb_0 : sb_1, :, :], axis=0)
                if nr > 0:
                    img_red  = np.sum(data[sr_0 : sr_1, :, :], axis=0)
                img_i = data[i, :, :]


            # Correct the layer i flux by subtracting blue and red contrib.

            if nb == 0 and nr == 0:
                img = img_i  # No correction possible
                msg1 = "not corrected."
                msg2 = "red band: -, blue band: -."
            elif nb == 0:
                if np.isnan(img_red).all():
                    img = img_i
                else:
                    corr_fac = spec[i]*nr/I_red
                    img = img_i - corr_fac*nr/img_red
                msg1 = "corrected using red band."
                msg2 = "blue band: {0:{width}}:{1:{width}} [px] " \
                    "(n = {2:{dwidth}}), red band: {3:{width}}:{4:{width}}" \
                    " [px] (n = {5:{dwidth}})".format("-", "-", \
                        0, sr_0 + 1, sr_1, sr_1 - sr_0, \
                        width=nwidth, dwidth=dwidth)
                msg3 = "."
                if sr_str != "":
                    msg3 = " :: Telluric line at continuum bandpass: " \
                        "{0}.".format(sr_str)
            elif nr == 0:
                if np.isnan(img_blue).all():
                    img = img_i
                else:
                    corr_fac = spec[i]*nb/I_blue
                    img = img_i - corr_fac*nb/img_blue
                msg1 = "corrected using blue band."
                msg2 = "blue band: {0:{width}}:{1:{width}} [px] " \
                    "(n = {2:{dwidth}}), red band: {3:{width}}:{4:{width}}" \
                    " [px] (n = {5:{dwidth}})".format(sb_0 + 1, \
                        sb_1, sb_1 - sb_0, "-", "-", \
                        0, width=nwidth, dwidth=dwidth)
                msg3 = "."
                if sb_str != "":
                    msg3 = " :: Telluric line at continuum bandpass: " \
                        "{0}.".format(sb_str)
            else:
                corr_fac = spec[i]/(I_blue/nb + I_red/nr)
                img = img_i - corr_fac*(img_blue/nb + img_red/nr)
                msg1 = "corrected using blue and red bands."
                msg2 = "blue band: {0:{width}}:{1:{width}} [px] " \
                    "(n = {2:{dwidth}}), red band: {3:{width}}:{4:{width}}" \
                    " [px] (n = {5:{dwidth}})".format(sb_0 + 1, \
                        sb_1, sb_1 - sb_0, sr_0 + 1, sr_1, \
                        sr_1 - sr_0, width=nwidth, dwidth=dwidth)
                msg3 = ""
                msgp = ""
                if sb_str != "" or sr_str != "":
                    msg3 = " :: Telluric line at continuum bandpass: "

                msg3b = ""
                if sb_str != "":
                    msg3b = "{0}".format(sb_str)
                    msgp = "."

                msg3r = ""
                if sr_str != "":
                    if msg3b != "": msg3b = msg3b + ", "
                    msg3r = "{0}".format(sr_str)
                    msgp = "."

                msg3 = msg3 + msg3b + msg3r + msgp

            if verbose >= 1:
                if verbose == 1:
                    log_str = msg1
                else:
                    log_str = msg2 + msg3
                print(log_str)
                logging.info(log_str)

            if axis_s == 1:
                odata[:, :, i] = img
            elif axis_s == 2:
                odata[:, i, :] = img
            else:
                odata[i, :, :] = img


        # Add data processing header entries and write a file with the
        # resulting data.

        hdr1["history"] = "d11: datacube processed with d11 to subtract the " \
            "continuum"
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

        
# Create a launcher in case the program is launched from system shell:
if __name__ == "__main__":
    from argparse import ArgumentParser, RawDescriptionHelpFormatter
 
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                            description=___doc___)

    ###########################################################################
    # Parsing command-line arguments and options:
    ###########################################################################

    parser.add_argument("filename", help="The name of a FITS file with a dat" \
                        "acube.")
    parser.add_argument("x", help="Aperture x coordinate [pixel].", type=float)
    parser.add_argument("y", help="Aperture y coordinate [pixel].", type=float)
    parser.add_argument("apr", help="Aperture radius [pixel].", type=float)
    parser.add_argument("cwidth", help="Bandwidth of continuum band for subt" \
                        "raction.", type=float)

    parser.add_argument("-d", "--offset", action="store", type=int, \
                        help="Specify a name of a plain-text file with tellu" \
                        "ric lines.")
    parser.add_argument("-t", "--telluric", action="store", type=str, \
                        help="Specify a name of a plain-text file with tellu" \
                        "ric lines.")
    parser.add_argument("-q", "--twidth", action="store", type=float, \
                        help="Telluric line bandwidth [Angstrom].")
    parser.add_argument("-u", "--commentstelluric", action="store", type=str, \
                        help="Specify a comment character to use with the pl" \
                        "ain-text telluric lines file [default: '#'].")
    parser.add_argument("-o", "--ofilename", action="store", type=str, \
                        help="Specify a name of the output file.")
    parser.add_argument("-w", "--overwrite", action="store_true", \
                        help="Overwrite existing output files.")
    parser.add_argument("-v", "--verbose", action="store", type=int, \
                        help="Be verbose on what the tool does to the data.")

    args = parser.parse_args()

    if args.offset is not None:
        offset = args.offset
    else:
        offset = 5

    if args.twidth is not None:
        twidth = args.twidth
    else:
        twidth = 3.0

    if args.verbose is not None:
        verbose = args.verbose
    else:
        verbose = 0

    d11(args.filename, args.x, args.y, args.apr, args.cwidth,
        offset=offset, telluric=args.telluric, twidth=twidth,
        commentstelluric=args.commentstelluric,
        ofilename=args.ofilename, overwrite=args.overwrite,
        verbose=verbose)
