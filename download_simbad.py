# -*- coding: utf-8 -*-
"""
Created on Sun Jun 19 09:38:41 2016

@author: daniella
"""

def download_simbad(inputfile,outputfile):

    import splat
    import time
    import numpy
    import numpy as np
    import pandas as pd
    from astroquery.simbad import Simbad
    from astropy import coordinates
    from astropy.coordinates import SkyCoord
    from astropy.coordinates import Angle
    from astropy import units as u
    from astropy.table import Column, Table, join
    
    #full = pd.read_csv('/Users/daniella/Python/Thesis/targetlist/Input_Lists/all_lists_combined.csv')
    folder = '/Users/daniella/Research/M7L5Sample/'
    #folder = '/Users/daniella/Python/Thesis/Targetlist/Input_Lists/'
    
    full = pd.read_csv(folder+inputfile)
    print(full)
    SIMBADdf = pd.DataFrame(index=np.arange(len(full)),columns=['INDEX','SIMBAD_RA','SIMBAD_DEC','SIMBAD_NAME',
                            'DESIGNATION','LIT_TYPE','LIT_TYPE_REF','OBJECT_TYPE','LUMINOSITY_CLASS',
                            'METALLICITY_CLASS','SIMBAD_OTYPE','SIMBAD_SPT','SIMBAD_SPT_REF','PARALLAX',
                            'PARALLAX_E','PARALLAX_REF','MU_RA','MU_DEC','MU','MU_E','MU_REF','RV',
                            'RV_E','RV_REF','VSINI','VSINI_E','VSINI_REF','J_2MASS','J_2MASS_E',
                            'H_2MASS','H_2MASS_E','KS_2MASS','KS_2MASS_E','r_SDSS','r_SDSS_E','i_SDSS',
                            'i_SDSS_E','z_SDSS','z_SDSS_E','V_JOHN','V_JOHN_E','R_JOHN','R_JOHN_E',
                            'I_JOHN','I_JOHN_E'])

    SIMBADdf['INDEX'] = full['INDEX']
    # now do a SIMBAD search for sources based on coordinates
    print('\nSIMBAD search')
    sb = Simbad()
    votfields = ['otype','parallax','sptype','propermotions','rot','rvz_radvel','rvz_error',\
    'rvz_bibcode','fluxdata(B)','fluxdata(V)','fluxdata(R)','fluxdata(I)','fluxdata(r)',\
    'fluxdata(i)','fluxdata(z)','fluxdata(J)','fluxdata(H)','fluxdata(K)']
    for v in votfields:
        sb.add_votable_fields(v)
        
    simbad_radius = 15*u.arcsec
    
    start = time.time()
    
    for i,des in enumerate(range(len(full))):
        print(i)
    #    t_sim = sb.query_region(full['Designation'][i],radius=simbad_radius)
        c = SkyCoord(full['RA (deg)'][i],full['Dec (deg)'][i],unit=(u.deg,u.deg),frame='icrs')
#        c = SkyCoord(full['ra_head'][i],full['dec_head'][i],unit=(u.hour,u.deg),frame='icrs')

        t_sim = sb.query_region(c,radius=simbad_radius)
        if isinstance(t_sim,Table):
            print('\nSource {} Designation = {} {} match(es)'.format(i+1,des,len(t_sim)))
    #        c = splat.designationToCoordinate(full['Designation'][i].split(' ')[1])
        #many sources found
            if len(t_sim) > 1:      # take the closest position
                sep = [c.separation(SkyCoord(t_sim['RA'][lp],t_sim['DEC'][lp],unit=(u.deg,u.deg))).arcsecond for lp in numpy.arange(len(t_sim))]
                t_sim['sep'] = sep
                t_sim.sort('sep')
                while len(t_sim)>1:
                    t_sim.remove_row(1) 
        #one source found
            else:
                t_sim['sep'] = [c.separation(SkyCoord(t_sim['RA'][0],t_sim['DEC'][0],unit=(u.hourangle,u.degree))).arcsecond]
            print(t_sim)
            SIMBADdf['SIMBAD_RA'][i] = c.ra.value
            SIMBADdf['SIMBAD_DEC'][i] = c.dec.value
            SIMBADdf['SIMBAD_NAME'][i] = t_sim['MAIN_ID'][0]
            SIMBADdf['DESIGNATION'][i] = 'J{}{}'.format(t_sim['RA'][0],t_sim['DEC'][0]).replace(' ','').replace('.','')          
            SIMBADdf['LIT_TYPE'][i] = t_sim['SP_TYPE'][0]
            SIMBADdf['LIT_TYPE_REF'][i] = t_sim['SP_BIBCODE'][0]
            if b'III' in t_sim['SP_TYPE'][0]:
                SIMBADdf['LUMINOSITY_CLASS'][i] = 'III'
                SIMBADdf['OBJECT_TYPE'][i] = 'GIANT'
            if b'IV' in t_sim['SP_TYPE'][0]:
                SIMBADdf['LUMINOSITY_CLASS'][i] = 'IV'
                SIMBADdf['OBJECT_TYPE'][i] = 'SUBGIANT'
            if b'VI' in t_sim['SP_TYPE'][0] or b'sd' in t_sim['SP_TYPE'][0]:
                SIMBADdf['METALLICITY_CLASS'][i] = 'sd'
            if b'VII' in t_sim['SP_TYPE'][0]:
                SIMBADdf['LUMINOSITY_CLASS'][i] = 'VII'
                SIMBADdf['OBJECT_TYPE'][i] = 'WD'
            SIMBADdf['SIMBAD_NAME'][i] = t_sim['MAIN_ID'][0]
            SIMBADdf['SIMBAD_OTYPE'][i] = t_sim['OTYPE'][0]
            SIMBADdf['SIMBAD_SPT'][i] = t_sim['SP_TYPE'][0]
            SIMBADdf['SIMBAD_SPT_REF'][i] = t_sim['SP_BIBCODE'][0]
            SIMBADdf['PARALLAX'][i] = t_sim['PLX_VALUE'][0]
            SIMBADdf['PARALLAX_E'][i] = t_sim['PLX_ERROR'][0]
            SIMBADdf['PARALLAX_REF'][i] = t_sim['PLX_BIBCODE'][0]
            SIMBADdf['MU_RA'][i] = t_sim['PMRA'][0]
            SIMBADdf['MU_DEC'][i] = t_sim['PMDEC'][0]
            SIMBADdf['MU'][i] = (t_sim['PMRA'][0]**2+t_sim['PMDEC'][0]**2)**0.5
            SIMBADdf['MU_E'][i] = t_sim['PM_ERR_MAJA'][0]
            SIMBADdf['MU_REF'][i] = t_sim['PM_BIBCODE'][0]
            SIMBADdf['RV'][i] = t_sim['RVZ_RADVEL'][0]
            SIMBADdf['RV_E'][i] = t_sim['RVZ_ERROR'][0]
            SIMBADdf['RV_REF'][i] = t_sim['RVZ_BIBCODE'][0]
            SIMBADdf['VSINI'][i] = t_sim['ROT_Vsini'][0]
            SIMBADdf['VSINI_E'][i] = t_sim['ROT_err'][0]
            SIMBADdf['VSINI_REF'][i] = t_sim['ROT_bibcode'][0]
            SIMBADdf['r_SDSS'][i] = t_sim['FLUX_r'][0]
            SIMBADdf['r_SDSS_E'][i] = t_sim['FLUX_ERROR_r'][0]
            SIMBADdf['i_SDSS'][i] = t_sim['FLUX_i'][0]
            SIMBADdf['i_SDSS_E'][i] = t_sim['FLUX_ERROR_i'][0]
            SIMBADdf['z_SDSS'][i] = t_sim['FLUX_z'][0]
            SIMBADdf['z_SDSS_E'][i] = t_sim['FLUX_ERROR_z'][0]
            SIMBADdf['J_2MASS'][i] = t_sim['FLUX_J'][0]
            SIMBADdf['J_2MASS_E'][i] = t_sim['FLUX_ERROR_J'][0]
            SIMBADdf['H_2MASS'][i] = t_sim['FLUX_H'][0]
            SIMBADdf['H_2MASS_E'][i] = t_sim['FLUX_ERROR_H'][0]
            SIMBADdf['KS_2MASS'][i] = t_sim['FLUX_K'][0]
            SIMBADdf['KS_2MASS_E'][i] = t_sim['FLUX_ERROR_K'][0]
            SIMBADdf['V_JOHN'][i] = t_sim['FLUX_V'][0]
            SIMBADdf['V_JOHN_E'][i] = t_sim['FLUX_ERROR_V'][0]
            SIMBADdf['R_JOHN'][i] = t_sim['FLUX_R'][0]
            SIMBADdf['R_JOHN_E'][i] = t_sim['FLUX_ERROR_R'][0]
            SIMBADdf['I_JOHN'][i] = t_sim['FLUX_I'][0]
            SIMBADdf['I_JOHN_E'][i] = t_sim['FLUX_ERROR_I'][0]
    end = time.time()
    print(end-start)
    
    SIMBADdf['RV_REF'] = SIMBADdf['RV_REF'].map(lambda x: x.decode() if pd.notnull(x) else np.nan)
    SIMBADdf['MU_REF'] = SIMBADdf['MU_REF'].map(lambda x: x.decode() if pd.notnull(x) else np.nan)
    SIMBADdf['PARALLAX_REF'] = SIMBADdf['PARALLAX_REF'].map(lambda x: x.decode() if pd.notnull(x) else np.nan)
    SIMBADdf['SIMBAD_SPT_REF'] = SIMBADdf['SIMBAD_SPT_REF'].map(lambda x: x.decode() if pd.notnull(x) else np.nan)
    SIMBADdf['SIMBAD_SPT'] = SIMBADdf['SIMBAD_SPT'].map(lambda x: x.decode() if pd.notnull(x) else np.nan)
    SIMBADdf['SIMBAD_OTYPE'] = SIMBADdf['SIMBAD_OTYPE'].map(lambda x: x.decode() if pd.notnull(x) else np.nan)
    SIMBADdf['LIT_TYPE_REF'] = SIMBADdf['LIT_TYPE_REF'].map(lambda x: x.decode() if pd.notnull(x) else np.nan)
    SIMBADdf['LIT_TYPE'] = SIMBADdf['LIT_TYPE'].map(lambda x: x.decode() if pd.notnull(x) else np.nan)
    SIMBADdf['SIMBAD_NAME'] = SIMBADdf['SIMBAD_NAME'].map(lambda x: x.decode() if pd.notnull(x) else np.nan)
    SIMBADdf.to_csv('SIMBADinfo_'+outputfile)
    
    print('Open file at '+folder+'SIMBADinfo_'+outputfile)


#for i,des in enumerate(full['Designation']):
#    t_sim = sb.query_region(des,radius=simbad_radius)
# source found in query
#        if isinstance(t_sim,Table):
#            print('\nSource {} Designation = {} {} match(es)'.format(i+1,des,len(t_sim)))
#            c = splat.designationToCoordinate(des)
# many sources found
#            if len(t_sim) > 1:      # take the closest position
#                sep = [c.separation(SkyCoord(str(t_sim['RA'][lp]),str(t_sim['DEC'][lp]),unit=(u.hourangle,u.degree))).arcsecond for lp in numpy.arange(len(t_sim))]
#                t_sim['sep'] = sep
#                t_sim.sort('sep')
#                while len(t_sim)>1:
#                    t_sim.remove_row(1) 
# one source found
#            else:
#                t_sim['sep'] = [c.separation(SkyCoord(str(t_sim['RA'][0]),str(t_sim['DEC'][0]),unit=(u.hourangle,u.degree))).arcsecond]
#            print(t_sim)
# update coordinates
#            c = splat.properCoordinates('{} {}'.format(t_sim['RA'][0],t_sim['DEC'][0]))
#            t_src['DESIGNATION'][i] = splat.coordinateToDesignation(c)
#            t_src['RA'][i] = c.ra.value
#            t_src['DEC'][i] = c.dec.value

# check if source is in the library already; if so, fill in source info
#            if t_sim['MAIN_ID'][0] in splat.DB_SOURCES['SIMBAD_NAME']:
#                for c in t_src.keys():
#                    t_src[c][i] = splat.DB_SOURCES[c][numpy.where(splat.DB_SOURCES['SIMBAD_NAME'] == t_sim['MAIN_ID'][0])][0]
#                    t_spec['SOURCE_KEY'][i] = t_src['SOURCE_KEY'][i]

# grab library spectra and see if any were taken on the same date (possible redundancy)
#                matchlib = splat.searchLibrary(idkey=t_src['SOURCE_KEY'][i])
#                if t_spec['OBSERVATION_DATE'][i] in matchlib['OBSERVATION_DATE']:
# previous observation on this date found - retain in case this is a better spectrum
#                    mkey = matchlib['DATA_KEY'][numpy.where(matchlib['OBSERVATION_DATE'] == t_spec['OBSERVATION_DATE'][i])]
#                    print('Previous spectrum found in library for data key {}'.format(mkey))
#                    compdict[str(t_spec['DATA_KEY'][i])]['comparison'] = splat.Spectrum(mkey)
#                    compdict[str(t_spec['DATA_KEY'][i])]['comparison_type'] = 'repeat spectrum: {}'.format(mkey)
# no previous observation on this date - retain the spectrum with the highest S/N
#                else:
#                    if len(matchlib) > 1:
#                        matchlib.sort('MEDIAN_SNR')
#                        matchlib.reverse()
#                    compdict[str(t_spec['DATA_KEY'][i])]['comparison'] = splat.Spectrum(matchlib['DATA_KEY'][0])
#                    compdict[str(t_spec['DATA_KEY'][i])]['comparison_type'] = 'alternate spectrum: {} taken on {}'.format(matchlib['DATA_KEY'][0],matchlib['OBSERVATION_DATE'][0])
# SIMBAD source is not in the library - fill in source information
#            else:
#                full['simbad_name'][i] = t_sim['MAIN_ID'][0]
#                t_src['DESIGNATION'][i] = 'J{}{}'.format(t_sim['RA'][0],t_sim['DEC'][0]).replace(' ','').replace('.','')
#                coord = splat.properCoordinates(t_src['DESIGNATION'][i])
#                t_src['RA'][i] = coord.ra.value
#                t_src['DEC'][i] = coord.dec.value
#                t_src['LIT_TYPE'][i] = t_sim['SP_TYPE'][0]
#                t_src['LIT_TYPE_REF'][i] = t_sim['SP_BIBCODE'][0]
#                t_src['OBJECT_TYPE'][i] = 'VLM'
#                if 'I' in t_sim['SP_TYPE'][0] and 'V' not in t_sim['SP_TYPE'][0]:
#                    t_src['LUMINOSITY_CLASS'][i] = 'I{}'.format(t_sim['SP_TYPE'][0].split('I',1)[1])
#                    t_src['OBJECT_TYPE'][i] = 'GIANT'
#                if 'VI' in t_sim['SP_TYPE'][0] or 'sd' in t_sim['SP_TYPE'][0]:
#                    t_src['METALLICITY_CLASS'][i] = 'sd'
#                t_src['SIMBAD_NAME'][i] = t_sim['MAIN_ID'][0]
#                t_src['SIMBAD_OTYPE'][i] = t_sim['OTYPE'][0]
#                t_src['SIMBAD_SPT'][i] = t_sim['SP_TYPE'][0]
#                t_src['SIMBAD_SPT_REF'][i] = t_sim['SP_BIBCODE'][0]
#                t_src['PARALLAX'][i] = t_sim['PLX_VALUE'][0]
#                t_src['PARALLAX_E'][i] = t_sim['PLX_ERROR'][0]
#                t_src['PARALLEX_REF'][i] = t_sim['PLX_BIBCODE'][0]
#                t_src['MU_RA'][i] = t_sim['PMRA'][0]
#                t_src['MU_DEC'][i] = t_sim['PMDEC'][0]
#                try:            # this is in case MU is not present
#                t_src['MU'][i] = (t_sim['PMRA'][0]**2+t_sim['PMDEC'][0]**2)**0.5
#                t_src['MU_E'][i] = t_sim['PM_ERR_MAJA'][0]
#                except:
#                    pass
#                t_src['MU_REF'][i] = t_sim['PM_BIBCODE'][0]
#                t_src['RV'][i] = t_sim['RVZ_RADVEL'][0]
#                t_src['RV_E'][i] = t_sim['RVZ_ERROR'][0]
#                t_src['RV_REF'][i] = t_sim['RVZ_BIBCODE'][0]
#                t_src['VSINI'][i] = t_sim['ROT_Vsini'][0]
#                t_src['VSINI_E'][i] = t_sim['ROT_err'][0]
#                t_src['VSINI_REF'][i] = t_sim['ROT_bibcode'][0]
#                t_src['J_2MASS'][i] = t_sim['FLUX_J'][0]
#                t_src['J_2MASS_E'][i] = t_sim['FLUX_ERROR_J'][0]
#                t_src['H_2MASS'][i] = t_sim['FLUX_H'][0]
#                t_src['H_2MASS_E'][i] = t_sim['FLUX_ERROR_H'][0]
#                t_src['KS_2MASS'][i] = t_sim['FLUX_K'][0]
#                t_src['KS_2MASS_E'][i] = t_sim['FLUX_ERROR_K'][0]

# no source found in SIMBAD: just check library - TO BE DONE LATER
#        else:
#            if t_src['SHORTNAME'][i] in splat.DB_SOURCES['SHORTNAME']:
#                for c in t_src.keys():
#                    print(c, t_src[c][i])
#                    print(splat.DB_SOURCES[c][numpy.where(splat.DB_SOURCES['SIMBAD_NAME'] == t_sim['MAIN_ID'])][0])
#                    t_src[c][i] = splat.DB_SOURCES[c][numpy.where(splat.DB_SOURCES['SHORTNAME'] == t_src['SHORTNAME'][i])][0]
#                    t_spec['SOURCE_KEY'][i] = t_src['SOURCE_KEY'][i]

# grab library spectra and see if any were taken on the same date (possible redundancy)
#                print(t_src['SOURCE_KEY'][i])
#                matchlib = splat.searchLibrary(idkey=t_src['SOURCE_KEY'][i])
#                if t_spec['OBSERVATION_DATE'][i] in matchlib['OBSERVATION_DATE']:
# previous observation on this date found - retain in case this is a better spectrum
#                    mkey = matchlib['DATA_KEY'][numpy.where(matchlib['OBSERVATION_DATE'] == t_spec['OBSERVATION_DATE'][i])]
#                    print('Previous spectrum found in library for data key {}'.format(mkey))
#                    compdict[str(t_spec['DATA_KEY'][i])]['comparison'] = splat.Spectrum(mkey)
#                    compdict[str(t_spec['DATA_KEY'][i])]['comparison_type'] = 'repeat spectrum: {}'.format(mkey)
# no previous observation on this date - retain the spectrum with the highest S/N
#                else:
#                    if len(matchlib) > 1:
#                        matchlib.sort('MEDIAN_SNR')
#                        matchlib.reverse()
#                    compdict[str(t_spec['DATA_KEY'][i])]['comparison'] = splat.Spectrum(matchlib['DATA_KEY'][0])
#                    compdict[str(t_spec['DATA_KEY'][i])]['comparison_type'] = 'alternate spectrum: {} taken on {}'.format(matchlib['DATA_KEY'][0],matchlib['OBSERVATION_DATE'][0])
#            else:
#                t_src['NAME'][i] = t_src['DESIGNATION'][i]

# generate check plots
#    junk = [sp.normalize() for sp in splist]
#    junk = [compdict[c]['comparison'].normalize() for c in compdict.keys()]
#    plotlist = []
#    legends = []
#    clist = []
#    for c in compdict.keys():
#        plotlist.append([compdict[c]['observed'],compdict[c]['comparison']])
#        legends.extend([compdict[c]['observed'].name,'{} {}'.format(compdict[c]['comparison'].name,compdict[c]['comparison_type'])])
#        clist.extend(['k','r'])
#    splat.plotSpectrum(plotlist,multiplot=True,layout=[2,2],multipage=True,legends=legends,output=review_folder+'review_plots.pdf',colors=clist,fontscale=0.5)

# output database updates
#    t_src.remove_column('SHORTNAME')
#    t_src.remove_column('SELECT')
#    t_spec.remove_column('SELECT')   
#    t_spec.remove_column('SOURCE_SELECT')
#    for i in numpy.arange(len(t_spec['NOTE'])):
#        t_spec['NOTE'][i] = compdict[str(t_spec['DATA_KEY'][i])]['comparison_type']
#    t_src.write(review_folder+'source_update.csv',format='ascii.csv')
#    t_spec.write(review_folder+'spectrum_update.csv',format='ascii.csv')
#    for c in compdict.keys():
#        print(c,compdict[c]['observed'],compdict[c]['comparison'],compdict[c]['comparison_type'])
#    print('\n')


# open up windows to review spreadsheets
# NOTE: WOULD LIKE TO MAKE THIS AUTOMATICALLY OPEN FILE
#    app = QtGui.QApplication(sys.argv)
#    window = Window(10, 5)
#    window.resize(640, 480)
#    window.show()
#    app.exec_()

#    print('\nSpectral plots and update speadsheets now available in {}'.format(review_folder))
#    response = raw_input('Please review and edit, and press any key when you are finished...')

# NEXT STEP - MOVE FILES TO APPROPRIATE PLACES, UPDATE MAIN DATABASES


#'Designation','RA (deg)', 'DEC (deg)', 'RA (hh mm ss.ss)','DEC (dd mm ss.ss)','Opt SpT', 'NIR SpT', 'simbad_spt','SpTn','sdss_rmag','sdss_rmag_e','sdss_imag','sdss_imag_e','sdss_zmag','sdss_zmag_e','jmag','jmag_e','hmag','hmag_e','kmag','kmag_e','wise_w1','wise_w1_e','wise_w2','wise_w2_e'
# STOPPED HERE
