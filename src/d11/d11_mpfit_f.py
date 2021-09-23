## d11_mpfit_f.py
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

def d11_mpfit_f(p, fjac=None, x=None, y=None, dy=None, coff=None,
                n_tied=None, funceval=False):
    # Parameter values are passed in "p"
    # If fjac==None then partial derivatives should not be
    # computed.  It will always be None if MPFIT is called with default
    # flag.

    import numpy as np

    smax = 26.0
    f = p[0] + p[1] * x

    p_3 = 1e-20
    if p[3] > p_3: p_3 = p[3]

    u = ((x - p[2])/p_3)**2
    mask = np.zeros(x.size, dtype=int)
    mask[u < smax ** 2] = 1
    f += p[4] / (np.sqrt(2.0 * np.pi) * p[3]) * mask * np.exp(- 0.5 * u * mask)
    del mask

    # Add the tied lines.
    for i in range(0, n_tied):
        ii = 5 + i
        u = ((x - (p[0]+coff[i]))/p_3)**2
        mask = np.zeros(x.size, dtype=int)
        mask[u < smax ** 2] = 1
        f += p[ii] / (np.sqrt(2.0 * np.pi) * p[3]) * mask \
            * np.exp(- 0.5 * u * mask)
        del mask

    # Non-negative status value means MPFIT should continue, negative means
    # stop the calculation.
    status = 0
    if funceval:
        return (f)
    else:
        return (status, (y-f)/dy)
