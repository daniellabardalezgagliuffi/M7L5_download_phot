# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 15:13:52 2016

@author: daniella
"""
def download_vizier(inputfile,outputfile):

    import splat
    import time
    import numpy
    import numpy as np
    import pandas as pd
    from astroquery.vizier import Vizier
    from astropy import coordinates
    from astropy.coordinates import SkyCoord
    from astropy.coordinates import Angle
    from astropy import units as u
    from astropy.table import Column, Table, join
    
    #Your input file needs these columns: INDEX, RA (deg), Dec (deg)    
    #full = pd.read_csv('/Users/daniella/Python/Thesis/targetlist/Input_Lists/all_lists_combined.csv')
    #folder = '/Users/daniella/Python/Thesis/Targetlist/Input_Lists/'
    folder = '/Users/daniella/Research/M7L5Sample/' 
    
#   full = pd.read_csv(folder+inputfile)
#    print(full)
#    full['INDEX'] = np.arange(len(full))
#    full.columns
#    print(full['ra_deg'])
    fullfull = pd.read_csv(folder+inputfile)
    full = fullfull.ix[1000:]
    #try these catalogs:
    #2MASS: JHK
    #SDSS: riz
    #SSSPM: RI
    #WISE: W1W2W3
    #ULAS: YJHK
    #PPMXL: PMRA, PMDEC
        #GAIA: G, parallax
        
#    VIZIERdf = pd.DataFrame(index=np.arange(len(full)),columns=['INDEX','2MASS_RA','2MASS_DEC','2MASS_NAME','2MASS_J',\
    VIZIERdf = pd.DataFrame(index=np.arange(len(full))+1000,columns=['INDEX','2MASS_RA','2MASS_DEC','2MASS_NAME','2MASS_J',\
    '2MASS_J_E','2MASS_H','2MASS_H_E','2MASS_KS','2MASS_KS_E','2MASS_QFLG','2MASS_CNTR','2MASS_OBSDATE',\
    '2MASS_OBSJD','2MASS_Rmag','2MASS_MATCHES','SDSS_RA','SDSS_DEC','SDSS_NAME','SDSS_cl','SDSS_OBSDATE',\
    'SDSS_QFLG','SDSS_r','SDSS_r_E','SDSS_i','SDSS_i_E','SDSS_z','SDSS_z_E','SDSS_ObjID','SDSS_SPSHAPE','SDSS_SPCLASS',\
    'SDSS_SPT','SDSS_pmRA','SDSS_pmRA_E','SDSS_pmDE','SDSS_pmDE_E','SDSS_MATCHES','WISE_RA','WISE_DEC',\
    'WISE_NAME','WISE_W1','WISE_W1_E','WISE_W2','WISE_W2_E','WISE_W3','WISE_W3_E','WISE_2M_JMAG',\
    'WISE_2M_JERR','WISE_2M_HMAG','WISE_2M_HERR','WISE_2M_KMAG','WISE_2M_KERR','WISE_pmRA','WISE_pmRA_E',
    'WISE_pmDE','WISE_pmDE_E','WISE_2M','WISE_2M_sep','WISE_MATCHES','ULAS_RA','ULAS_DEC','ULAS_NAME',\
    'ULAS_Y','ULAS_Y_E','ULAS_J','ULAS_J_E','ULAS_H','ULAS_H_E','ULAS_K','ULAS_K_E','ULAS_EPOCH',\
    'ULAS_cl','ULAS_pmRA','ULAS_pmRA_E','ULAS_pmDE','ULAS_pmDE_E','ULAS_MATCHES','PPMXL_RA','PPMXL_DEC',
    'PPMXL_pmRA','PPMXL_pmRA_E','PPMXL_pmDE','PPMXL_pmDE_E','PPMXL_epochRA','PPMXL_epochDE','PPMXL_J',
    'PPMXL_J_E','PPMXL_H','PPMXL_H_E','PPMXL_K','PPMXL_K_E','PPMXL_USNO_R2','PPMXL_USNO_I','PPMXL_MATCHES',
    'GAIA_G','GAIA_plx','GAIA_plx_e','GAIA_pmRA','GAIA_pmDE','GAIA_pmRA_e','GAIA_pmDE_e'])
    
    VIZIERdf['INDEX'] = full['INDEX'].values
    
    # now do a VIZIER search for sources based on coordinates
    print('\nVIZIER search')
    v = Vizier()
    v.TIMEOUT = 5000
        
    vizier_radius = 15*u.arcsec
    
    start = time.time()    
    
#    for i,des in enumerate(range(len(full))):
    for i in np.arange(len(full))+1000:
        des = i
        c = SkyCoord(full['RA (deg)'][i],full['Dec (deg)'][i],unit=(u.deg,u.deg),frame='icrs')
  #      c = SkyCoord(full['ra_head'][i],full['dec_head'][i],unit=(u.hour,u.deg),frame='icrs')
        result = Vizier.query_region(c,radius=vizier_radius, catalog=['SDSS9','2MASS','AllWISE','ULAS9','PPMXL','GAIA'])    
        if 'II/246/out' in result.keys():
            tmass = result['II/246/out']
        else:
            tmass = np.nan
        if 'V/139/sdss9' in result.keys():
            sdss = result['V/139/sdss9']
        else:
            sdss = np.nan
        if 'II/328/allwise' in result.keys():
            wise = result['II/328/allwise']
        else:
            wise = np.nan
        if 'II/319/las9' in result.keys():
            ulas = result['II/319/las9']
        else:
            ulas = np.nan
        if 'I/317/sample' in result.keys():
            ppmxl = result['I/317/sample']
        else:
            ppmxl = np.nan
        if 'I/337/gaia' in result.keys():
            gaia = result['I/337/gaia']
        else:
            gaia = np.nan            
 
#        VIZIERdf['INDEX'][i] = full['INDEX'][i]

        if isinstance(tmass,Table):
            print('\nSource {} Designation = {} {} match(es) in 2MASS'.format(i+1,des,len(tmass)))
            #many sources found
            n_tmass = len(tmass)
            if len(tmass) > 1:      # take the closest position
                sep = [c.separation(SkyCoord(tmass['RAJ2000'][lp],tmass['DEJ2000'][lp],unit=(u.deg,u.deg))).arcsecond for lp in numpy.arange(len(tmass))]
                tmass['sep'] = sep
                tmass.sort('sep')
                while len(tmass)>1:
                    tmass.remove_row(1) 
            #one source found
            else:
                tmass['sep'] = [c.separation(SkyCoord(tmass['RAJ2000'][0],tmass['DEJ2000'][0],unit=(u.deg,u.deg))).arcsecond]
            print(tmass)
            VIZIERdf['2MASS_RA'][i] = tmass['RAJ2000'][0]
            VIZIERdf['2MASS_DEC'][i] = tmass['DEJ2000'][0]
            VIZIERdf['2MASS_NAME'][i] = 'J{}'.format(tmass['_2MASS'][0].decode())
            VIZIERdf['2MASS_J'][i] = tmass['Jmag'][0]
            VIZIERdf['2MASS_J_E'][i] = tmass['e_Jmag'][0]
            VIZIERdf['2MASS_H'][i] = tmass['Hmag'][0]
            VIZIERdf['2MASS_H_E'][i] = tmass['e_Hmag'][0]
            VIZIERdf['2MASS_KS'][i] = tmass['Kmag'][0]
            VIZIERdf['2MASS_KS_E'][i] = tmass['e_Kmag'][0]
            VIZIERdf['2MASS_QFLG'][i] = tmass['Qflg'][0].decode()
            VIZIERdf['2MASS_CNTR'][i] = tmass['Cntr'][0]
            VIZIERdf['2MASS_OBSDATE'][i] = tmass['Date'][0].decode()
            VIZIERdf['2MASS_OBSJD'][i] = tmass['JD'][0]
            VIZIERdf['2MASS_Rmag'][i] = tmass['Rmag'][0]
            VIZIERdf['2MASS_MATCHES'][i] = n_tmass
        
        if isinstance(sdss,Table):
            print('\nSource {} Designation = {} {} match(es) in SDSS9'.format(i+1,des,len(sdss)))
        #many sources found
            n_sdss = len(sdss)
            if len(sdss) > 1:      # take the closest position
                sep = [c.separation(SkyCoord(sdss['RAJ2000'][lp],sdss['DEJ2000'][lp],unit=(u.deg,u.deg))).arcsecond for lp in numpy.arange(len(sdss))]
                sdss['sep'] = sep
                sdss.sort('sep')
                while len(sdss)>1:
                    sdss.remove_row(1) 
        #one source found
            else:
                sdss['sep'] = [c.separation(SkyCoord(sdss['RAJ2000'][0],sdss['DEJ2000'][0],unit=(u.deg,u.deg))).arcsecond]
            print(sdss)
            VIZIERdf['SDSS_RA'][i] = sdss['RAJ2000'][0]
            VIZIERdf['SDSS_DEC'][i] = sdss['DEJ2000'][0]
            VIZIERdf['SDSS_NAME'][i] = sdss['SDSS9'][0].decode()
            VIZIERdf['SDSS_cl'][i] = sdss['cl'][0]         #classification by shape: 3=galaxy, 6=star
            VIZIERdf['SDSS_OBSDATE'][i] = sdss['ObsDate'][0]
            VIZIERdf['SDSS_QFLG'][i] = sdss['Q'][0]         #(0=unknown): 1=bad 2=acceptable 3=good 4=missing 5=hole
            VIZIERdf['SDSS_r'][i] = sdss['rpmag'][0]    #PSF magnitude
            VIZIERdf['SDSS_r_E'][i] = sdss['e_rpmag'][0]
            VIZIERdf['SDSS_i'][i] = sdss['ipmag'][0]
            VIZIERdf['SDSS_i_E'][i] = sdss['e_ipmag'][0]
            VIZIERdf['SDSS_z'][i] = sdss['zpmag'][0]
            VIZIERdf['SDSS_z_E'][i] = sdss['e_zpmag'][0]
            VIZIERdf['SDSS_ObjID'][i] = sdss['objID'][0]
            VIZIERdf['SDSS_SPSHAPE'][i] = sdss['spType'][0].decode()
            VIZIERdf['SDSS_SPCLASS'][i] = sdss['spCl'][0].decode()
            VIZIERdf['SDSS_SPT'][i] = sdss['subClass'][0].decode()
            VIZIERdf['SDSS_pmRA'][i] = sdss['pmRA'][0]
            VIZIERdf['SDSS_pmRA_E'][i] = sdss['e_pmRA'][0]
            VIZIERdf['SDSS_pmDE'][i] = sdss['pmDE'][0]
            VIZIERdf['SDSS_pmDE_E'][i] = sdss['e_pmDE'][0]
            VIZIERdf['SDSS_MATCHES'][i] = n_sdss
        
        if isinstance(wise,Table):
            print('\nSource {} Designation = {} {} match(es) in AllWISE'.format(i+1,des,len(wise)))
        #many sources found
            n_wise = len(wise)
            if len(wise) > 1:      # take the closest position
                sep = [c.separation(SkyCoord(wise['RAJ2000'][lp],wise['DEJ2000'][lp],unit=(u.deg,u.deg))).arcsecond for lp in numpy.arange(len(wise))]
                wise['sep'] = sep
                wise.sort('sep')
                while len(wise)>1:
                    wise.remove_row(1) 
        #one source found
            else:
                wise['sep'] = [c.separation(SkyCoord(wise['RAJ2000'][0],wise['DEJ2000'][0],unit=(u.deg,u.deg))).arcsecond]
            print(wise)
            VIZIERdf['WISE_RA'][i] = wise['RAJ2000'][0]
            VIZIERdf['WISE_DEC'][i] = wise['DEJ2000'][0]
            VIZIERdf['WISE_NAME'][i] = wise['AllWISE'][0].decode()
            VIZIERdf['WISE_W1'][i] = wise['W1mag'][0]         
            VIZIERdf['WISE_W1_E'][i] = wise['e_W1mag'][0]
            VIZIERdf['WISE_W2'][i] = wise['W2mag'][0]       
            VIZIERdf['WISE_W2_E'][i] = wise['e_W2mag'][0]
            VIZIERdf['WISE_W3'][i] = wise['W3mag'][0]
            VIZIERdf['WISE_W3_E'][i] = wise['e_W3mag'][0]
            VIZIERdf['WISE_2M_JMAG'][i] = wise['Jmag'][0]
            VIZIERdf['WISE_2M_JERR'][i] = wise['e_Jmag'][0]        
            VIZIERdf['WISE_2M_HMAG'][i] = wise['Hmag'][0]
            VIZIERdf['WISE_2M_HERR'][i] = wise['e_Hmag'][0] 
            VIZIERdf['WISE_2M_KMAG'][i] = wise['Kmag'][0]
            VIZIERdf['WISE_2M_KERR'][i] = wise['e_Kmag'][0] 
            VIZIERdf['WISE_pmRA'][i] = wise['pmRA'][0]
            VIZIERdf['WISE_pmRA_E'][i] = wise['e_pmRA'][0]
            VIZIERdf['WISE_pmDE'][i] = wise['pmDE'][0]
            VIZIERdf['WISE_pmDE_E'][i] = wise['e_pmDE'][0]
            VIZIERdf['WISE_2M'][i] = wise['_2Mkey'][0]
            VIZIERdf['WISE_2M_sep'][i] = wise['d2M'][0]     #arcsec
            VIZIERdf['WISE_MATCHES'][i] = n_wise
        
        if isinstance(ulas,Table):
            print('\nSource {} Designation = {} {} match(es) in ULAS9'.format(i+1,des,len(ulas)))
        #many sources found
            n_ulas = len(ulas)
            if len(ulas) > 1:      # take the closest position
                sep = [c.separation(SkyCoord(ulas['RAJ2000'][lp],ulas['DEJ2000'][lp],unit=(u.deg,u.deg))).arcsecond for lp in numpy.arange(len(ulas))]
                ulas['sep'] = sep
                ulas.sort('sep')
                while len(ulas)>1:
                    ulas.remove_row(1) 
        #one source found
            else:
                ulas['sep'] = [c.separation(SkyCoord(ulas['RAJ2000'][0],ulas['DEJ2000'][0],unit=(u.deg,u.deg))).arcsecond]
            print(ulas)
            VIZIERdf['ULAS_RA'][i] = ulas['RAJ2000'][0]
            VIZIERdf['ULAS_DEC'][i] = ulas['DEJ2000'][0]
            VIZIERdf['ULAS_NAME'][i] = ulas['ULAS'][0].decode()
            VIZIERdf['ULAS_Y'][i] = ulas['Ymag'][0]        
            VIZIERdf['ULAS_Y_E'][i] = ulas['e_Ymag'][0]
            VIZIERdf['ULAS_J'][i] = ulas['Jmag2'][0]         
            VIZIERdf['ULAS_J_E'][i] = ulas['e_Jmag2'][0]
            VIZIERdf['ULAS_H'][i] = ulas['Hmag'][0]
            VIZIERdf['ULAS_H_E'][i] = ulas['e_Hmag'][0]
            VIZIERdf['ULAS_K'][i] = ulas['Kmag'][0]
            VIZIERdf['ULAS_K_E'][i] = ulas['e_Kmag'][0]        
            VIZIERdf['ULAS_EPOCH'][i] = ulas['Epoch'][0]
            VIZIERdf['ULAS_cl'][i] = ulas['cl'][0]         #-2=galaxy,-1=star
            VIZIERdf['ULAS_pmRA'][i] = ulas['pmRA'][0]
            VIZIERdf['ULAS_pmRA_E'][i] = ulas['e_pmRA'][0]
            VIZIERdf['ULAS_pmDE'][i] = ulas['pmDE'][0]
            VIZIERdf['ULAS_pmDE_E'][i] = ulas['e_pmDE'][0]
            VIZIERdf['ULAS_MATCHES'][i] = n_ulas
        
        if isinstance(ppmxl,Table):
            print('\nSource {} Designation = {} {} match(es) in PPMXL'.format(i+1,des,len(ppmxl)))
        #many sources found
            n_ppmxl = len(ppmxl)
            if len(ppmxl) > 1:      # take the closest position
                sep = [c.separation(SkyCoord(ppmxl['RAJ2000'][lp],ppmxl['DEJ2000'][lp],unit=(u.deg,u.deg))).arcsecond for lp in numpy.arange(len(ppmxl))]
                ppmxl['sep'] = sep
                ppmxl.sort('sep')
                while len(ppmxl)>1:
                    ppmxl.remove_row(1) 
        #one source found
            else:
                ppmxl['sep'] = [c.separation(SkyCoord(ppmxl['RAJ2000'][0],ppmxl['DEJ2000'][0],unit=(u.deg,u.deg))).arcsecond]
            print(ppmxl)
            VIZIERdf['PPMXL_RA'][i] = ppmxl['RAJ2000'][0]
            VIZIERdf['PPMXL_DEC'][i] = ppmxl['DEJ2000'][0]
            VIZIERdf['PPMXL_pmRA'][i] = ppmxl['pmRA'][0]
            VIZIERdf['PPMXL_pmRA_E'][i] = ppmxl['e_pmRA'][0]
            VIZIERdf['PPMXL_pmDE'][i] = ppmxl['pmDE'][0]
            VIZIERdf['PPMXL_pmDE_E'][i] = ppmxl['e_pmDE'][0]        
            VIZIERdf['PPMXL_epochRA'][i] = ppmxl['epRA'][0]         
            VIZIERdf['PPMXL_epochDE'][i] = ppmxl['epDE'][0]
            VIZIERdf['PPMXL_J'][i] = ppmxl['Jmag'][0]        
            VIZIERdf['PPMXL_J_E'][i] = ppmxl['e_Jmag'][0]
            VIZIERdf['PPMXL_H'][i] = ppmxl['Hmag'][0]
            VIZIERdf['PPMXL_H_E'][i] = ppmxl['e_Hmag'][0]
            VIZIERdf['PPMXL_K'][i] = ppmxl['Kmag'][0]
            VIZIERdf['PPMXL_K_E'][i] = ppmxl['e_Kmag'][0]        
            VIZIERdf['PPMXL_USNO_R2'][i] = ppmxl['r2mag'][0]
            VIZIERdf['PPMXL_USNO_I'][i] = ppmxl['imag'][0] 
            VIZIERdf['PPMXL_MATCHES'][i] = n_ppmxl

        if isinstance(gaia,Table):
            print('\nSource {} Designation = {} {} match(es) in GAIA'.format(i+1,des,len(gaia)))
            #many sources found
            n_gaia = len(gaia)
            if len(gaia) > 1:      # take the closest position
                sep = [c.separation(SkyCoord(gaia['_RAJ2000'][lp],gaia['_DEJ2000'][lp],unit=(u.deg,u.deg))).arcsecond for lp in np.arange(len(gaia))]
                gaia['sep'] = sep
                gaia.sort('sep')
                while len(gaia)>1:
                    gaia.remove_row(1) 
            #one source found
            else:
                gaia['sep'] = [c.separation(SkyCoord(gaia['_RAJ2000'][0],gaia['_DEJ2000'][0],unit=(u.deg,u.deg))).arcsecond]
            print(gaia)
            VIZIERdf['GAIA_G'][i] = gaia['__Gmag_'][0]
            VIZIERdf['GAIA_plx'][i] = gaia['Plx'][0]
            VIZIERdf['GAIA_plx_e'][i] = gaia['e_Plx'][0]
            VIZIERdf['GAIA_pmRA'][i] = gaia['pmRA'][0]
            VIZIERdf['GAIA_pmDE'][i] = gaia['pmDE'][0]
            VIZIERdf['GAIA_pmRA_e'][i] = gaia['e_pmRA'][0]
            VIZIERdf['GAIA_pmDE_e'][i] = gaia['e_pmDE'][0]
            
        end = time.time()
        print(end-start)
        
    VIZIERdf.to_csv('VIZIERinfo_'+outputfile)
    
    print('Open file at '+folder+'VIZIERinfo_'+outputfile)

