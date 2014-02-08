import numpy as np
import os, pyfits, pickle
from glob import glob
from astropy import coordinates as coord
from astropy import units as u
import reddening

class photometry(object):
    """A class to store photometry data and artificial stars as output
    by DOLPHOT. This class is modified to preserve the full raw data
    catalogs without any culling and without extinction corrections.

    Parameters
    ----------
    name : string
        The absolute path of the DOLPHOT photometry file in this case
        can also use the shortened phot2 format output by the GHOSTS
        pipeline.
    reference : string
        The absolute path of the drizzled reference image
    filter1_name : string
        Name of the first HST filter
    filter2_name : string
        Name of the second HST filter

    Attributes
    ----------
    filename : string
        The absolute path of the photometry file
    reference : string
        The absolute path of the drizzled reference image
    coords : astropy coordinates object
        The coordinates of the aperture specified as an astropy object
    exposures1 : float array
        Array of the exposure lengths in seconds from the first filter
    exposures2 : float array
        Array of the exposure lengths in seconds from the second filter

    x_in : float array
        Input x positions of the detections on the reference image if
        applicable
    y_in : float array
        Input y positions of the detections on the reference image if
        applicable
    mag1_in : float array
        Input magnitude of the stars in the first filter if applicable
    mag2_in : float array
        Input magnitude of the stars in the second filter if applicable

    x : float array
        x positions of the detections on the reference image
    y : float array
        y positions of the detections on the reference image
    type : unsigned int array
        Object types for the detections (1 = good star, 2 = faint star
        witout PSF, 3 = elongated object, 4 = object too sharp, 5 =
        extended object)
    mag1 : float array
        Extinction corrected magnitudes in the first filter
    mag1_err : float array
        Errors on the magnitudes from the first filter
    chi1 : float array
        Chi value of fit to detection from first filter
    sn1 : float array
        Signal to noise ration of detection from first filter
    sharp1 : float array
        Sharpness of detection from first filter (positive if too
        sharp, negative if too broad)
    round1 : float array
        Roundness of detection from first filter (used to
        differentiate extended objects and diffraction spikes)
    crowd1 : float array
        Crowding of detection from first filter. Measures how much
        brighter in magnitudes the star would be if nearby stars were
        not simultaneously fit.
    flag1 : unsigned int array
        Error flag for detection from first filter (values added to
        error flag are: +1 = photometry aperture extends off chip, +2
        = too many bad or saturated pixels, +4 = saturated at center
        of star, +8 = extreme case of previous errors)
    mag2 : float array
        as before, but for second filter
    mag2_err : float array
        as before, but for second filter
    chi2 : float array
        as before, but for second filter
    sn2 : float array
        as before, but for second filter
    sharp2 : float array
        as before, but for second filter
    round2 : float array
        as before, but for second filter
    crowd2 : float array
        as before, but for second filter
    flag2 : unsigned int array
        as before, but for second filter

    """
    def __init__(self, name, reference, filter1_name="F606W", filter2_name="F814W"):
        # Set basic information for photometry run
        self.name = name
        self.reference = reference
        self.get_ra_dec()
        self.get_info()
        # Initialise the photometry data
        self.x_in     = np.asfarray([])
        self.y_in     = np.asfarray([])
        self.mag1_in  = np.asfarray([])
        self.mag2_in  = np.asfarray([])
        self.x        = np.asfarray([])
        self.y        = np.asfarray([])
        self.type     = np.asarray([], np.uint8)
        self.back1    = np.asfarray([])
        self.mag1     = np.asfarray([])
        self.mag1_err = np.asfarray([])
        self.chi1     = np.asfarray([])
        self.sn1      = np.asfarray([])
        self.sharp1   = np.asfarray([])
        self.round1   = np.asfarray([])
        self.crowd1   = np.asfarray([])
        self.flag1    = np.asarray([], np.uint8)
        self.back2    = np.asfarray([])
        self.mag2     = np.asfarray([])
        self.mag2_err = np.asfarray([])
        self.chi2     = np.asfarray([])
        self.sn2      = np.asfarray([])
        self.sharp2   = np.asfarray([])
        self.round2   = np.asfarray([])
        self.crowd2   = np.asfarray([])
        self.flag2    = np.asarray([], np.uint8)

        # Load in the photometry file
        if os.path.exists(name):
            if (name[-7:]=='phot.gz'):
                self.loadDOLPHOT()
                self.extinction_mag1, self.extinction_mag2 = reddening.find_ra_dec(self.ra,self.dec)
            elif (name.find('fake')>0):
                self.loadDOLPHOTartificial()
                self.extinction_mag1, self.extinction_mag2 = (0,0)
            else:
                print "Photometry file format not recognized"
        else:
            print "Photometry files", name, "not found"
    
    def get_ra_dec(self):
        """A method that returns the RA and Dec from the reference
        image
        
        Parameters & Attributes
        -----------------------
        As defined in base class
        """
        header = pyfits.getheader(self.reference,0)
        ra = float(header['RA_APER'])
        dec = float(header['DEC_APER'])
        self.coords = coord.ICRSCoordinates(ra, dec, unit=(u.degree, u.degree))

    def get_info(self):
        """Set some rudimentary information from the *.phot.info file
        if available

        Parameters & Attributes
        -----------------------
        As defined in base class
        """

        self.exposures1=np.asfarray([])
        self.exposures2=np.asfarray([])
        file_dir = self.name[:self.name.rfind('/')]
        base_name = os.path.basename(self.name)
        phot_file_name = "{}/{}.phot.info".format(file_dir, base_name[:base_name.find('.')])
        with open(phot_file_name, "r") as filestream:
            for line in filestream:
                elements = line.split()
                # If there are 6 elements in the file then we are
                # looking at the exposure information
                if (len(elements) == 6):
                    # See if this line corresponds to the first filter
                    # or second and append the exposure time to the
                    # relevant list
                    if (line.find(self.filter1_name.upper())>=0):
                        self.exposures1=np.append(self.exposures1, float(elements[-1]))
                    elif (line.find(self.filter2_name.upper())>=0):
                        self.exposures2=np.append(self.exposures2, float(elements[-1]))


    def loadDOLPHOT(self):
        """Loads data from the original DOLPHOT photometry file
        
        Parameters & Attributes
        -----------------------
        As defined in base class

        Notes
        -----
        No filtering of the raw data, and no extinction correction
        """
        print self.name
        phot = np.loadtxt(self.name)
        
        self.x        = np.asfarray(phot[:,3])
        self.y        = np.asfarray(phot[:,4])
        self.type     = np.asarray(phot[:,11], np.uint8)
        # Apply extinction correction
        self.back1    = np.asfarray(phot[:,13])
        self.mag1     = np.asfarray(phot[:,16])
        self.mag1_err = np.asfarray(phot[:,18])
        self.chi1     = np.asfarray(phot[:,19])
        self.sn1      = np.asfarray(phot[:,20])
        self.sharp1   = np.asfarray(phot[:,21])
        self.round1   = np.asfarray(phot[:,22])
        self.crowd1   = np.asfarray(phot[:,23])
        self.flag1    = np.asarray(phot[:,24], np.uint8)
        # Apply extinction correction
        self.back2    = np.asfarray(phot[:,26])
        self.mag2     = np.asfarray(phot[:,29])
        self.mag2_err = np.asfarray(phot[:,31])
        self.chi2     = np.asfarray(phot[:,32])
        self.sn2      = np.asfarray(phot[:,33])
        self.sharp2   = np.asfarray(phot[:,34])
        self.round2   = np.asfarray(phot[:,35])
        self.crowd2   = np.asfarray(phot[:,36])
        self.flag2    = np.asarray(phot[:,37], np.uint8)


    def loadDOLPHOTartificial(self):
        """Loads data from the artificial star tests output by DOLPHOT
        
        Parameters & Attributes
        -----------------------
        As defined in base class

        Notes
        -----
        No filtering of the raw data and no extinction correction
        """

        file_dir = self.name[:self.name.rfind('/')]
        base_name = os.path.basename(self.name)
        phot_file_name = "{}/{}.phot.info".format(file_dir, base_name[:base_name.find('.')])
        posF606 = False
        posF814 = False
        nframe=0
        with open(phot_file_name, "r") as filestream:
            # Find the number of exposures in the combined photometry
            nframe = int(filestream.readline().split()[0])         
            # Find the filter order 
            i = 1
            while(not (bool(posF606) & bool(posF814))):
              try:
                line = filestream.readline().split()[3]
                if ((line == 'F606W') & (not posF606)):
                    posF606 = i
                if ((line == 'F814W') & (not posF814)):
                    posF814 = i
                i+=1
              except:
                i+=1
                continue
            posmin = min(posF606,posF814)
            posF606 -= posmin
            posF814 -= posmin
            print posF606,posF814
        
        # Calculate the offset to the combined photometry sets for each filter
        offset = 3 + nframe*2
        
        # Find, load and append all the individual artificial star tests for the field
        files = glob(self.name[:self.name.rfind('_')])
        for fn in files:
            phot = np.loadtxt(fn)
            self.x_in     = np.append(self.x, np.asfarray(phot[:,2]))
            self.y_in     = np.append(self.y, np.asfarray(phot[:,3]))
            self.mag1_in  = np.append(self.mag1_in, np.asfarray(phot[:,5+posF606*2]))
            self.mag2_in  = np.append(self.mag2_in, np.asfarray(phot[:,5+posF814*2]))
            self.x        = np.append(self.x, np.asfarray(phot[:,3 + offset]))
            self.y        = np.append(self.y, np.asfarray(phot[:,4 + offset]))
            self.type     = np.append(self.type, np.asarray(phot[:,11 + offset], np.uint8))
            self.back1    = np.append(self.back1, np.asfarray(phot[:,13 + offset]))
            self.mag1     = np.append(self.mag1, np.asfarray(phot[:,16 + offset]))
            self.mag1_err = np.append(self.mag1_err, np.asfarray(phot[:,18 + offset]))
            self.chi1     = np.append(self.chi1, np.asfarray(phot[:,19 + offset]))
            self.sn1      = np.append(self.sn1, np.asfarray(phot[:,20 + offset]))
            self.sharp1   = np.append(self.sharp1, np.asfarray(phot[:,21 + offset]))
            self.round1   = np.append(self.round1, np.asfarray(phot[:,22 + offset]))
            self.crowd1   = np.append(self.crowd1, np.asfarray(phot[:,23 + offset]))
            self.flag1    = np.append(self.flag1, np.asarray(phot[:,24 + offset], np.uint8))
            self.back2    = np.append(self.back2, np.asfarray(phot[:,26 + offset]))
            self.mag2     = np.append(self.mag2, np.asfarray(phot[:,29 + offset]))
            self.mag2_err = np.append(self.mag2_err, np.asfarray(phot[:,31 + offset]))
            self.chi2     = np.append(self.chi2, np.asfarray(phot[:,32 + offset]))
            self.sn2      = np.append(self.sn2, np.asfarray(phot[:,33 + offset]))
            self.sharp2   = np.append(self.sharp2, np.asfarray(phot[:,34 + offset]))
            self.round2   = np.append(self.round2, np.asfarray(phot[:,35 + offset]))
            self.crowd2   = np.append(self.crowd2, np.asfarray(phot[:,36 + offset]))
            self.flag2    = np.append(self.flag2, np.asarray(phot[:,37 + offset], np.uint8))

def load_from_dir(dir,drz_dir="/Volumes/data/processed/drc",outfile="photometry_full.dat"):
    """Loads photometry files into a list of photometry classes and
    saves as a binary file
        
    Parameters
    ----------
    dir : string
        The directory where the photometry files are located
    drz_dir : string
        The directory where the correspoding drizzled reference files
        are located. If not specified, defaults here to
        /Volumes/data/processed/drc
    outfile : string
        The name of the file to save the photometry classes to as a
        pickled object. Will be orverwritten if already exists, and if
        not specified defaults to "photometry.dat" in the current
        directory
    """
    # Load in the files
    files = glob("{}/*phot.gz".format(dir))
    if (os.path.exists(outfile)):
        os.system('rm {}'.format(outfile))
    photometry_list = []
    counter=0
    with open(outfile,'wb') as out:
        for file_name in files:
            print counter
            absolute_file_name = os.path.abspath(file_name)
            base_name = os.path.basename(absolute_file_name)
            drc_file_name = "{}/{}_drc.chip1.fits".format(drz_dir,base_name[:base_name.find('.')])
            photometry_list.extend([photometry(absolute_file_name, drc_file_name)])
            counter+=1
        # Using pickle to save custom class object
        pickle.dump(photometry_list, out)
        print 'wrote',len(photometry_list),'files'

if __name__=="__main__":
    # If this file is called from the command line, will load
    # photometry data from the directory "/Volumes/data/HStars/dolphot
    load_from_dir('/Volumes/data/empty/')
