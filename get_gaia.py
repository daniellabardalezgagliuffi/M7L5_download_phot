#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 10:31:29 2018

@author: daniella
"""

def get_gaia(df):

    import time
    import numpy as np
    import pandas as pd
    from astroquery.vizier import Vizier
    from astropy.coordinates import SkyCoord
    from astropy.coordinates import Angle
    from astropy import units as u
    from astropy.table import Column, Table, join
    from astroquery.gaia import Gaia
    Gaia.login(user='dbardale',password='DBGgaia2018!')
    
    qry = """
        SELECT TOP 10 g.*, t.*
        FROM gaiadr1.tmass_original_valid AS t
        JOIN gaiadr2.tmass_neighbourhood AS xt ON xt.tmass_oid=t.tmass_oid
        JOIN gaiadr2.gaia_source AS g ON g.source_id=xt.source_id
        WHERE g.phot_g_mean_mag IS NOT NULL
        """

    bkg=Gaia.launch_job_async(qry).get_results().to_pandas()
    bkg['abs_g']=bkg.phot_g_mean_mag-5*np.log10(1000./bkg.parallax)+5.
    columns=bkg.columns.tolist()
    
    GAIAdf=pd.DataFrame(index=df.index,columns=columns)
    
    for k,row in df.iterrows():
        qry="""
            SELECT g.*, t.*
            FROM gaiadr1.tmass_original_valid AS t
            LEFT OUTER JOIN gaiadr2.tmass_neighbourhood AS xt ON xt.tmass_oid = t.tmass_oid
            LEFT OUTER JOIN gaiadr2.gaia_source AS g ON xt.source_id = g.source_id
            where 1=CONTAINS(POINT('ICRS', t.ra, t.dec),CIRCLE('ICRS', {}, {}, 5./3600))
            """.format(row['RA (deg)'],row['Dec (deg)'])
        
        data=Gaia.launch_job_async(qry).get_results().to_pandas()
        
        data['INDEX']=k
        print(k,data.shape)
        GAIAdf=GAIAdf.append(data).dropna(how='all')
        Gaia.logout()
    
    return GAIAdf
   
