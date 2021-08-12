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

def d11(filename, x, y, apr, cwidth, ofilename='DEF',
        offset=5, wave=None, spec=None, overwrite=False, verbose=0):

    import os
    import sys
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

    if verbose is None: verbose = 0
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
                   screxe + "  offset = " + str(offset), \
                   screxe + "  ofilename = \"" + ofilename + "\""]
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

        bin_start = dwave + offset
        bin_final = nwave - (dwave+offset)


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


        # Sum the flux in the selected aperture for all layers.

        if axis_s == 1:
            spec = np.nansum(data[y_0 : y_1, x_0 : x_1, :], axis=(0, 1))
        elif axis_s == 2:
            spec = np.nansum(data[y_0 : y_1, :, x_0 : x_1], axis=(0, 2))
        else:
            spec = np.nansum(data[:, y_0 : y_1, x_0 : x_1], axis=(1, 2))


        # Step through the cube on the dispersion axis and subtract scaled
        # blue and red continuum bands from the datacube bins.

        for i in range(bin_start, bin_final):

            if verbose >= 1:
                log_str = screxe + "Layer " + str(i) + "/" + str(nwave) + "."
                print(log_str)
                logging.info(log_str)

            # Calculate the contribution to the flux on the feature blue side.

            sb_0 = i - (offset+dwave)
            sb_1 = i - offset + 1
            if sb_0 < 0: sb_0 = 0
            if sb_1 < 0: sb_1 = 0

            if sb_1 > sb_0:
                nb = sb_1 - sb_0
                I_blue = spec[sb_0 : sb_1].sum()
            else:
                nb = 0

            # Calculate the contribution to the flux on the feature red side.

            sr_0 = i + offset
            sr_1 = i + (offset+dwave) + 1

            if sr_0 > nwave: sr_0 = nwave
            if sr_1 > nwave: sr_1 = nwave

            if sr_1 > sr_0:
                nr = sr_1 - sr_0
                I_red  = spec[sr_0 : sr_1].sum()
            else:
                nr = 0


            # Sum up all flux in the blue and red bands around the feature.

            if axis_s == 1:
                if nb > 0:
                    img_blue = np.nansum(data[:, :, sb_0 : sb_1], axis=2)
                if nr > 0:
                    img_red  = np.nansum(data[:, :, sr_0 : sr_1], axis=2)
                img_i = data[:, :, i]
            elif axis_s == 2:
                if nb > 0:
                    img_blue = np.nansum(data[:, sb_0 : sb_1, :], axis=1)
                if nr > 0:
                    img_red  = np.nansum(data[:, sr_0 : sr_1, :], axis=1)
                img_i = data[:, i, :]
            else:
                if nb > 0:
                    img_blue = np.nansum(data[sb_0 : sb_1, :, :], axis=0)
                if nr > 0:
                    img_red  = np.nansum(data[sr_0 : sr_1, :, :], axis=0)
                img_i = data[i, :, :]


            # Correct the layer i flux by subtracting blue and red contrib.

            if nb == 0 and nr == 0:
                img = img_i  # No correction possible
                msg = "not corrected."
            elif nb == 0:
                corr_fac = spec[i]*nr/I_red
                img = img_i - corr_fac*nr/img_red
                msg = "corrected using red band."
            elif nr == 0:
                corr_fac = spec[i]*nb/I_blue
                img = img_i - corr_fac*nb/img_blue
                msg = "corrected using blue band."
            else:
                corr_fac = spec[i]/(I_blue/nb + I_red/nr)
                img = img_i - corr_fac*(img_blue/nb + img_red/nr)
                msg = "corrected using blue and red bands."

            if verbose >= 1:
                log_str = screxe + "Flux " + msg
                if verbose >= 2: print(log_str)
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
                            description="""
Here is an example of how this tool is launched from the shell or the console
  $ d11.py -o datacube_d11.fits -v 1 datacube.fits 12.0 25.0 2.0 3.0
       """)

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

    parser.add_argument("-o", "--ofilename", action="store", type=str, \
                        help="Specify a name of the output file.")
    parser.add_argument("-w", "--overwrite", action="store_true", \
                        help="Overwrite existing output files.")
    parser.add_argument("-v", "--verbose", action="store", type=int, \
                        help="Be verbose on what the tool does to the data.")

    args = parser.parse_args()

    d11(args.filename, args.x, args.y, args.apr, args.cwidth,
        ofilename=args.ofilename, overwrite=args.overwrite,
        verbose=args.verbose)
