import interp_to_P
import matplotlib
import sys
import numpy as np
from metpy.units import units

#varlist = list(['QCLOUD','QRAIN','QSNOW','QGRAUP'])
#varlist = list(['TEMPERATURE','QVAPOR','QICE','QSNOW','QGRAUP','QCLOUD','QRAIN'])
varlist = list(['W','TEMPERATURE','QVAPOR','QNICE','QICE','QSNOW','QGRAUP','QCLOUD','QRAIN',\
                'PRI_INU','PRI_IDE','PRS_IAU','PRI_IHM','PRS_SCW'])

period = sys.argv[1]
#fix = 'Trn' 
fix = 'Stlshmn'
plevs = np.arange(100000,10000,-3000) * units.Pa
years = list([2009,2012,2014,2015,2019])
months = list(['12','01','02'])

#for y in years:
#    for m in months:
#        for v in varlist:
#            print(v)
#            interp_to_P.interp_p(plevs,y,m,period,v)


interp_to_P.interp_p_varlist(plevs,2015,'02',period,varlist,fix)
