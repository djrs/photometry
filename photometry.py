import numpy as np
import os, pyfits, pickle
from glob import glob
from astropy import coordinates as coord
from astropy import units as u
import reddening

class photometry(object):
    """A class to store photometry data output from DOLPHOT

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
        self.filter1_name = filter1_name
        self.filter2_name = filter2_name
        self.get_ra_dec()
        self.get_info()
        self.extinction_mag1, self.extinction_mag2 = reddening.find_ra_dec(self.ra,self.dec)

        # Load in the photometry file
        if os.path.exists(name):
            if (name[-7:]=='phot.gz'):
                self.loadDOLPHOT()
            elif (name[-16:]=='phot_abridged.gz'):
                self.loadDOLPHOTabridged()
            else:
                print "Photometry file format not recognized"
        else:
            print "Photometry files", name, "not found"
            self.x        = np.asfarray([0])
            self.y        = np.asfarray([0])
            self.type     = np.asarray([0], np.uint8)
            self.mag1     = np.asfarray([0] - self.extinction_mag1)
            self.mag1_err = np.asfarray([0])
            self.chi1     = np.asfarray([0])
            self.sn1      = np.asfarray([0])
            self.sharp1   = np.asfarray([0])
            self.round1   = np.asfarray([0])
            self.crowd1   = np.asfarray([0])
            self.flag1    = np.asarray([0], np.uint8)
            self.mag2     = np.asfarray([0] - self.extinction_mag2)
            self.mag2_err = np.asfarray([0])
            self.chi2     = np.asfarray([0])
            self.sn2      = np.asfarray([0])
            self.sharp2   = np.asfarray([0])
            self.round2   = np.asfarray([0])
            self.crowd2   = np.asfarray([0])
            self.flag2    = np.asarray([0], np.uint8)
    
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
        """Loads data from the original DOLPHOT photometry file and
        applies some broad selection cuts
        
        Parameters & Attributes
        -----------------------
        As defined in base class

        Notes
        -----
        Filters the selection based on cuts on object type, signal to
        noise ration, and sharpness
        """
        print self.name
        phot = np.loadtxt(self.name)
        
        # Select objects of type 1 with error flags = 1 or 2
        sel = (phot[:,24]<=2) & (phot[:,37]<=2) & (phot[:,11]==1) 
        # Select objects with a signal/noise ratio >=5
        sel = (sel) & (phot[:,20]>=5.0) & (phot[:,33]>=5.0)
        # Select objects with a combined sharpness > -0.06
        sel = (sel & (phot[:,21]+phot[:,34] > -0.06))
        # Select objects with a combined sharpness < 1.30
        sel = (sel & (phot[:,21]+phot[:,34] < 1.30))

        print "Culled data from", len(sel), "detections to", np.sum(sel)

        self.x        = np.asfarray(phot[sel,3])
        self.y        = np.asfarray(phot[sel,4])
        self.type     = np.asarray(phot[sel,11], np.uint8)
        # Apply extinction correction
        self.mag1     = np.asfarray(phot[sel,16] - self.extinction_mag1)
        self.mag1_err = np.asfarray(phot[sel,18])
        self.chi1     = np.asfarray(phot[sel,19])
        self.sn1      = np.asfarray(phot[sel,20])
        self.sharp1   = np.asfarray(phot[sel,21])
        self.round1   = np.asfarray(phot[sel,22])
        self.crowd1   = np.asfarray(phot[sel,23])
        self.flag1    = np.asarray(phot[sel,24], np.uint8)
        # Apply extinction correction
        self.mag2     = np.asfarray(phot[sel,29] - self.extinction_mag2)
        self.mag2_err = np.asfarray(phot[sel,31])
        self.chi2     = np.asfarray(phot[sel,32])
        self.sn2      = np.asfarray(phot[sel,33])
        self.sharp2   = np.asfarray(phot[sel,34])
        self.round2   = np.asfarray(phot[sel,35])
        self.crowd2   = np.asfarray(phot[sel,36])
        self.flag2    = np.asarray(phot[sel,37], np.uint8)

            
    def loadDOLPHOTabridged(self):
        """Loads data from the GHOSTS photometry file and applies some
        broad selection cuts
        
        Parameters & Attributes
        -----------------------
        As defined in base class

        Notes
        -----
        Filters the selection based on cuts on object type, signal to
        noise ration, and sharpness
        """
        print self.name
        phot = np.loadtxt(self.name)
        
        # Select objects of type 1 with error flags = 1 or 2
        sel = (phot[:,12]<=2) & (phot[:,22]<=2) & (phot[:,2]==1) 
        # Select objects with a signal/noise ratio >=5
        sel = (sel) & (phot[:,8]>=5.0) & (phot[:,18]>=5.0)
        # Select objects with a combined sharpness > -0.06
        sel = (sel & (phot[:,9]+phot[:,19] > -0.06))
        # Select objects with a combined sharpness < 1.30
        sel = (sel & (phot[:,9]+phot[:,19] < 1.30))

        print "Culled data from", len(sel), "detections to", np.sum(sel)

        self.x        = np.asfarray(phot[sel,0])
        self.y        = np.asfarray(phot[sel,1])
        self.type     = np.asarray(phot[sel,2], np.uint8)
        # Apply extinction correction
        self.mag1     = np.asfarray(phot[sel,4] - self.extinction_mag1)
        self.mag1_err = np.asfarray(phot[sel,6])
        self.chi1     = np.asfarray(phot[sel,7])
        self.sn1      = np.asfarray(phot[sel,8])
        self.sharp1   = np.asfarray(phot[sel,9])
        self.round1   = np.asfarray(phot[sel,10])
        self.crowd1   = np.asfarray(phot[sel,11])
        self.flag1    = np.asarray(phot[sel,12], np.uint8)
        # Apply extinction correction
        self.mag2     = np.asfarray(phot[sel,14] - self.extinction_mag2)
        self.mag2_err = np.asfarray(phot[sel,16])
        self.chi2     = np.asfarray(phot[sel,17])
        self.sn2      = np.asfarray(phot[sel,18])
        self.sharp2   = np.asfarray(phot[sel,19])
        self.round2   = np.asfarray(phot[sel,20])
        self.crowd2   = np.asfarray(phot[sel,21])
        self.flag2    = np.asarray(phot[sel,22], np.uint8)

def load_from_dir(dir,drz_dir="/Volumes/data/processed/drc",outfile="photometry.dat"):
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
    # Load in the files, try the abridged version first and if they
    # don't exist look for the raw DOLPHOT output
    files = glob("{}/*phot_abridged.gz".format(dir))
    if (len(files)<1):
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
    load_from_dir('/Volumes/data/HStars/dolphot')
