"""
Overhaul of cosmouline's star module, for alipy2.
This module contains stuff for geometric matching algorithms.
"""

import sys
import os
import math
import numpy as np
import operator  # For sorting
import copy
import scipy.linalg
import scipy.spatial


class Star:
    """
    Simple class to represent a single source (usually stars, but not necessarily).
    In this module we often manipulate lists of such Star objects.
    """

    def __init__(self, x=0.0, y=0.0, name="untitled", flux=-1.0, props={}, fwhm=-1.0, elon=-1.0):
        """
        flux : Some "default" or "automatic" flux, might be a just good guess. Used for sorting etc.
        If you have several fluxes, colours, store them in the props dict.
        props : A placeholder dict to contain other properties of your choice (not required nor used by the methods).
        """
        self.x = float(x)
        self.y = float(y)
        self.name = str(name)
        self.flux = float(flux)
        self.props = props
        self.fwhm = float(fwhm)
        self.elon = float(elon)

    def copy(self):
        return copy.deepcopy(self)

    def __getitem__(self, key):
        """
        Used for sorting list of stars.
        """
        if key == 'flux':
            return self.flux
        if key == 'fwhm':
            return self.fwhm
        if key == 'elon':
            return self.elon

    def __setitem__(self, key, value):
        """
        Don't use it!
        """
        if key == 'flux':
            self.flux = value
        if key == 'fwhm':
            self.fwhm = value
        if key == 'elon':
            self.elon = value

    def __str__(self):
        """
        A string representation of a source.
        """
        return "%10s : (%8.2f,%8.2f) | %12.2f | %5.2f %5.2f" % (self.name, self.x, self.y, self.flux,
                                                                self.fwhm, self.elon)

    def coords(self, full=False):
        """
        Returns the coords in form of an array.

        :param full: If True, I also include flux, fwhm, elon
        :type full: boolean

        """
        if full:
            return np.array([self.x, self.y, self.flux, self.fwhm, self.elon])
        else:
            return np.array([self.x, self.y])

    def distance(self, otherstar):
        """
        Returns the distance between the two stars.
        """
        return math.sqrt(np.sum((self.coords() - otherstar.coords())**2))

    def distanceandsort(self, otherstarlist):
        """
        Returns a list of dicts(star, dist, origpos), sorted by distance to self.
        The 0th star is the closest.

        otherstarlist is not modified.
        """

        returnlist = []
        for i, star in enumerate(otherstarlist):
            dist = self.distance(star)
            returnlist.append({'star': star, 'dist': dist, 'origpos': i})
        returnlist = sorted(returnlist, key=operator.itemgetter('dist'))  # sort stars according to dist

        return returnlist

# And now some functions to manipulate list of such stars ###


def listtoarray(starlist, full=False):
    """
    Transforms the starlist into a 2D numpy array for fast manipulations.
    First index is star, second index is x or y

    :param full: If True, I also include flux, fwhm, elon
    :type full: boolean

    """
    return np.array([star.coords(full=full) for star in starlist])


def area(starlist, border=0.01):
    """
    Returns the area covered by the stars.
    Border is relative to max-min
    """
    if len(starlist) == 0:
        return np.array([0, 1, 0, 1])

    if len(starlist) == 1:
        star = starlist[0]
        return np.array([star.x - 0.5, star.x + 0.5, star.y - 0.5, star.y + 0.5])

    a = listtoarray(starlist)
    (xmin, xmax) = (np.min(a[:, 0]), np.max(a[:, 0]))
    (ymin, ymax) = (np.min(a[:, 1]), np.max(a[:, 1]))
    xw = xmax - xmin
    yw = ymax - ymin
    xmin = xmin - border*xw
    xmax = xmax + border*xw
    ymin = ymin - border*yw
    ymax = ymax + border*yw
    return np.array([xmin, xmax, ymin, ymax])


def readsexcat(sexcat, hdu=0, verbose=True, maxflag=3, posflux=True, minfwhm=1.0, propfields=[]):
    """
    sexcat is either a string (path to a file), or directly an asciidata catalog object as returned by pysex

    :param hdu: The hdu containing the science data from which I should build the catalog.
    0 will select the only available extension. If multihdu, 1 is usually science.

    We read a sextractor catalog with astroasciidata and return a list of stars.
    Minimal fields that must be present in the catalog :

        * NUMBER
        * X_IMAGE
        * Y_IMAGE
        * FWHM_IMAGE
        * ELONGATION
        * FLUX_AUTO
        * FLAGS

    maxflag : maximum value of the FLAGS that you still want to keep. Sources with higher values will be skipped.
        * FLAGS == 0 : all is fine
        * FLAGS == 2 : the flux is blended with another one; further info in the sextractor manual.
        * FLAGS == 4    At least one pixel of the object is saturated (or very close to)
        * FLAGS == 8    The object is truncated (too close to an image boundary)
        * FLAGS is the sum of these ...

    posflux : if True, only stars with positive FLUX_AUTO are included.

    propfields : list of FIELD NAMES to be added to the props of the stars.

    I will always add FLAGS as a propfield by default.

    """
    returnlist = []

    if isinstance(sexcat, str):

        import asciidata
        if not os.path.isfile(sexcat):
            print("Sextractor catalog does not exist :")
            print(sexcat)
            sys.exit(1)

        if verbose:
            print("Reading %s " % (os.path.split(sexcat)[1]))
        mycat = asciidata.open(sexcat)

    else:  # then it's already a asciidata object
        mycat = sexcat

    # We check for the presence of required fields :
    minimalfields = ["NUMBER", "X_IMAGE", "Y_IMAGE", "FWHM_IMAGE", "ELONGATION", "FLUX_AUTO", "FLAGS"]
    minimalfields.extend(propfields)
    availablefields = [col.colname for col in mycat]
    for field in minimalfields:
        if field not in availablefields:
            print("Field %s not available in your catalog file !" % (field))
            sys.exit(1)

    if verbose:
        print("Number of sources in catalog : %i" % (mycat.nrows))

    propfields.append("FLAGS")
    propfields = list(set(propfields))

    if mycat.nrows == 0:
        if verbose:
            print("No stars in the catalog :-(")
    else:
        for i, num in enumerate(mycat['NUMBER']):
            if float(mycat['FLAGS'][i]) > maxflag:
                continue
            flux = float(mycat['FLUX_AUTO'][i])
            if posflux and (flux < 0.0):
                continue
            fwhm = float(mycat['FWHM_IMAGE'][i])
            if float(fwhm) <= minfwhm:
                continue

            props = dict([[propfield, mycat[propfield][i]] for propfield in propfields])

            newstar = Star(x=mycat['X_IMAGE'][i], y=mycat['Y_IMAGE'][i], name=str(num), flux=flux,
                           props=props, fwhm=mycat['FWHM_IMAGE'][i], elon=mycat['ELONGATION'][i])

            returnlist.append(newstar)

    if verbose:
        print("I've selected %i sources" % (len(returnlist)))

    return returnlist


def sortstarlistbyflux(starlist):
    """
    We sort starlist according to flux : highest flux first !
    """
    sortedstarlist = sorted(starlist, key=operator.itemgetter('flux'))
    sortedstarlist.reverse()
    return sortedstarlist


class SimpleTransform:
    """
    Represents an affine transformation consisting of rotation, isotropic scaling, and shift.
    [x', y'] = [[a -b], [b a]] * [x, y] + [c d]
    """

    def __init__(self, v=(1, 0, 0, 0)):
        """
        v = (a, b, c, d)
        """
        self.v = np.asarray(v)

    def getscaling(self):
        return math.sqrt(self.v[0]*self.v[0] + self.v[1]*self.v[1])

    def getrotation(self):
        """
        The CCW rotation angle, in degrees
        """
        return math.atan2(self.v[1], self.v[0]) * (180.0/math.pi)  # % 360.0
    
    def get_shift(self):
        return self.v[2:4]

    def __str__(self):
        return "Rotation %+11.6f [deg], scale %8.6f" % (self.getrotation(), self.getscaling())

    def inverse(self):
        """
        Returns the inverse transform !
        """

        # To represent affine transformations with matrices, we can use homogeneous coordinates.
        homo = np.array([
            [self.v[0], -self.v[1], self.v[2]],
            [self.v[1],  self.v[0], self.v[3]],
            [0.0, 0.0, 1.0]
        ])

        inv = scipy.linalg.inv(homo)

        return SimpleTransform((inv[0, 0], inv[1, 0], inv[0, 2], inv[1, 2]))

    def matrixform(self):
        """
        Special output for scipy.ndimage.interpolation.affine_transform
        Returns (matrix, offset)
        """

        return (np.array([[self.v[0], -self.v[1]], [self.v[1], self.v[0]]]), self.v[2:4])

    def apply(self, x, y):
        """
        Applies the transform to a point (x, y)
        """
        xn = self.v[0]*x - self.v[1]*y + self.v[2]
        yn = self.v[1]*x + self.v[0]*y + self.v[3]
        return (xn, yn)

    def applystar(self, star):
        transstar = star.copy()
        (transstar.x, transstar.y) = self.apply(transstar.x, transstar.y)
        return transstar

    def applystarlist(self, starlist):
        return [self.applystar(star) for star in starlist]


def fitstars(uknstars, refstars, verbose=True):
    """
    I return the transform that puts the unknown stars (uknstars) onto the refstars.
    If you supply only two stars, this is using linalg.solve() -- perfect solution.
    If you supply more stars, we use linear least squares, i.e. minimize the 2D error.

    Formalism inspired by :
    http://math.stackexchange.com/questions/77462/
    """

    assert len(uknstars) == len(refstars)
    if len(uknstars) < 2:
        if verbose:
            print("Sorry I cannot fit a transform on less than 2 stars.")
        return None

    # ukn * x = ref
    # x is the transform (a, b, c, d)

    ref = np.hstack(listtoarray(refstars))  # a 1D vector of lenth 2n

    uknlist = []
    for star in uknstars:
        uknlist.append([star.x, -star.y, 1, 0])
        uknlist.append([star.y, star.x, 0, 1])
    ukn = np.vstack(np.array(uknlist))  # a matrix

    if len(uknstars) == 2:
        trans = scipy.linalg.solve(ukn, ref)
    else:
        trans = scipy.linalg.lstsq(ukn, ref)[0]

    return SimpleTransform(np.asarray(trans))


#     def teststars(self, uknstars, refstars, r=5.0, verbose=True):
#         """
#         We apply the trans to the uknstarlist, and check for correspondance with the refstarlist.
#         Returns the number of uknstars that could be matched to refstars within r [unit : reference image pixels !].
#         """
#
#         transuknstars = self.applystarlist(uknstars)
#
#         transukn = listtoarray(transuknstars)
#         ref = listtoarray(refstars)
#         #print "Unknown stars   : ", transukn.shape[0]
#         #print "Reference stars : ", ref.shape[0]
#
#         mindists = np.min(scipy.spatial.distance.cdist(ref, transukn), axis=0)
#         print " must become smarter !!!!"
#         nbmatch = np.sum(mindists < r)
#         if verbose:
#             print "Tested match on %4i references : OK for %4i/%4i unknown stars (r = %.1f)." % (len(refstars),
#                                                                                     nbmatch, len(uknstars), r)
#
#         return nbmatch
#
#
#     def refinestars(self, uknstars, refstars, r=5.0, verbose=True):
#         """
#         I refit myself to all matching stars.
#         """
#
#         transuknstars = self.applystarlist(uknstars)
#         transukn = listtoarray(transuknstars)
#         ref = listtoarray(refstars)
#
#         # Brute force...
#         dists = scipy.spatial.distance.cdist(ref, transukn)
#         uknmindistindexes = np.argmin(dists, axis=0) # For each ukn, the index of the closest ref
#         uknmindist = np.min(dists, axis=0) # The corresponding distances
#         uknkeepers = uknmindist < r
#
#         matchuknstars = []
#         matchrefstars = []
#         for i in range(len(uknkeepers)):
#             if uknkeepers[i] == True:
#                 matchuknstars.append(uknstars[i])
#                 matchrefstars.append(refstars[uknmindistindexes[i]])
#         if verbose:
#             print "Refining (before / after) :"
#             print self
#         self.fitstars(matchuknstars, matchrefstars)
#         if verbose:
#             print self
#


def identify(uknstars, refstars, trans=None, r=5.0, verbose=True, getstars=False):
    """
    Allows to:
     * get the number or matches, i.e. evaluate the quality of the trans
     * get corresponding stars from both lists (without the transform applied)

    :param getstars: If True, I return two lists of corresponding stars, instead of just the number of matching stars
    :type getstars: boolean

    Inspired by the "formpairs" of alipy 1.0 ...

    """

    if trans is not None:
        ukn = listtoarray(trans.applystarlist(uknstars))
    else:
        ukn = listtoarray(uknstars)
    ref = listtoarray(refstars)

    dists = scipy.spatial.distance.cdist(ukn, ref)  # Big table of distances between ukn and ref
    mindists = np.min(dists, axis=1)  # For each ukn, the minimal distance
    minok = mindists <= r  # booleans for each ukn
    minokindexes = np.argwhere(minok).flatten()  # indexes of uknstars with matches

    if verbose:
        print("%i/%i stars with distance < r = %.1f (mean %.1f, median %.1f, std %.1f)" % (np.sum(minok),
                                                                                           len(uknstars), r,
                                                                                           np.mean(mindists[minok]),
                                                                                           np.median(mindists[minok]),
                                                                                           np.std(mindists[minok])))
    matchuknstars = []
    matchrefstars = []

    for i in minokindexes:  # we look for the second nearest ...
        sortedrefs = np.argsort(dists[i, :])
        firstdist = dists[i, sortedrefs[0]]
        seconddist = dists[i, sortedrefs[1]]
        if seconddist > 2.0*firstdist:  # Then the situation is clear, we keep it.
            matchuknstars.append(uknstars[i])
            matchrefstars.append(refstars[sortedrefs[0]])
        else:
            pass  # Then there is a companion, we skip it.

    if verbose:
        print("Filtered for companions, keeping %i/%i matches" % (len(matchuknstars), np.sum(minok)))

    if getstars is True:
        return (matchuknstars, matchrefstars)
    else:
        return len(matchuknstars)
