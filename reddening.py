import mechanize
from mechanize import ParseResponse, urlopen, urljoin

def find_ra_dec(ra,dec):
    """
    Return the foreground extinction in a patch of sky specified by
    equatorial coordinates.

    Parameters
    ----------
    ra : float
        The Right Accesion in J2000.0 equatorial coordinates
    dec : float
        The declination in J2000.0 equatorial coordinates

    Returns
    -------
    A_F606W : float
        The foreground extinction in the ACS/F606W band pass as
        measured in magnitudes

    A_F814W : float
        The foreground extinction in the ACS/F814W band pass as
        measured in magnitudes

    Notes
    -----
    This code queries the NASA Extragalactic Database (NED), using the
    mechanize library to find the extinction. It matches the nearest
    object with 10 arcmin and queries the galactic reddening, E(B-V),
    for that object. This is then converted to an extinction in ACS
    bandpasses using the conversions from Sirianni et al. (2005). If
    anything goes wrong with the web query, both return values are set
    to zero.

    Examples
    --------
    >>> extinction_f606w, extinction_f814w = find_ra_dec(120.2, 32.6)
    (0.143948, 0.09518799999999998)

    """
    
    if (ra<-180) | (ra>360) | (dec<-90) | (dec>90):
        print 'incorrect ra and dec specified'
        return 0,0
    try:
        # Read in the forms from the NED near position query page
        response = urlopen("http://ned.ipac.caltech.edu/forms/nearposn.html")
        forms = ParseResponse(response, backwards_compat=False)
        form = forms[0]
        form["lon"] = "{:9.5f}d".format(ra)
        form["lat"] = "{:9.5f}".format(dec)
        form["radius"] = "10.0"
        opener = mechanize.OpenerFactory(mechanize.SeekableResponseOpener).build_opener()

        # Submit the page and read back the data
        data = opener.open(form.click())
        forms = ParseResponse(data, backwards_compat=False)
        linedata = data.get_data().split('\n')

        # Initialize some parameters used to parse the returned website
        searchFlag,reddeningFlag=True,False
        reddening = 0
        url = ''
        for line in linedata:
            # Read in URL to the first object results page if a list of objects is returned
            if (line[:16]=="</strong><A HREF") & searchFlag:
                url=line[line.find('=')+1:line.find('TARGET')]
                searchFlag=False
            # Find if a single object was returned (with corresponding extinction data)
            if (line[:29]=='<h3><font color="#CCCCCCC">1 '):
                reddeningFlag=True
            # If it is a single object then read in the E(B) and E(V) values
            if (reddeningFlag) & (line[:21]=="<tr><td>A<sub>&lambda") & (reddening==0):
                elements=line.split('<td>')
                b_red=elements[3]
                v_red=elements[4]
                reddening=float(b_red[:b_red.rfind('<')])-float(v_red[:v_red.rfind('<')])

        # If a single object was not returned, then follow the URL to the nearest object
        if (reddeningFlag==False) & (url[1:8]=='cgi-bin'):
            response = urlopen('http://ned.ipac.caltech.edu/'+url)
            for line in response.read().split('\n'):
                # Find and read the E(B) and E(V) values
                if (line[:21]=="<tr><td>A<sub>&lambda") & (reddening==0):
                    elements=line.split('<td>')
                    b_red=elements[3]
                    v_red=elements[4]
                    reddening=float(b_red[:b_red.rfind('<')])-float(v_red[:v_red.rfind('<')])

    except:
        # Catching all errors as there are numerous ways in which the web query can fail.
        # Usually these can be remedied by just running the query again, so just printing
        # a warning for now.
        print "WARNING: Extinction query failed, returning (0, 0)"
        reddening=0

    # Convert the E(B-V) value to extinctions in the F606W and F814W bands
    return reddening*2.716, reddening*1.796
        
