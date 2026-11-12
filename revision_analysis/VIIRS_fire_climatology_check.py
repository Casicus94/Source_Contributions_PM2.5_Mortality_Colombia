# VIIRS fire-climatology check for Reviewer 3, major comment 4.
# Step 1 (shell, on the full FIRMS file) - South America totals by year, vegetation fires (type 0), confidence n/h:
#   awk -F, 'NR>1 && $15==0 && $10!="l" {y=substr($6,1,4); m=substr($6,6,2)+0; f=$13; A[y]+=f; if($1>=0){NS[y]+=f}; if(m<=4){J[y]+=f}; if($1>=0 && m<=4){NJ[y]+=f}} END{for(y in A) print y","A[y]","J[y]","NS[y]","NJ[y]}' fire_archive_SV-C2_508148.csv
# Step 2 (shell) - crop to the inner WRF domain:
#   awk -F, 'NR==1{print "lat,lon,date,frp,conf,daynight,type";next} $1>=-2 && $1<=13 && $2>=-80 && $2<=-66 {print $1","$2","$6","$13","$10","$14","$15}' fire_archive_SV-C2_508148.csv > VIIRS_2013_2023_D02.csv
# Step 3 (this script) - Colombia and department totals; needs AWB_vs_FINN_check.py in the same folder for the masks.
import pandas as pd, numpy as np, importlib.util, shapefile
from shapely.geometry import shape
spec=importlib.util.spec_from_file_location('m','AWB_vs_FINN_check.py'); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
W='../World_Countries/'
d=pd.read_csv('../Fires/VIIRS_2013_2023_D02.csv')
print('rows',len(d), d.type.value_counts().to_dict(), d.conf.value_counts().to_dict())
d=d[(d.type==0)&(d.conf!='l')].copy()
d['year']=d.date.str[:4].astype(int)
col=m.load_regions(W)['Colombia']
d['col']=m.inside(col,d.lon.values,d.lat.values)
sf=shapefile.Reader(W+'GDBIDEAM_departamentos.shp',encoding='latin-1')
d['dept']=''
dc=d[d.col]
for r in sf.shapeRecords():
    g=shape(r.shape.__geo_interface__); name=m.norm(r.record[2])
    idx=dc.index[m.inside(g,dc.lon.values,dc.lat.values)]
    d.loc[idx,'dept']=name
d.to_pickle('viirs_tagged.pkl')
c=d[d.col]
t=pd.DataFrame({'domain_FRP':d.groupby('year').frp.sum(),'COL_FRP':c.groupby('year').frp.sum(),'COL_count':c.groupby('year').size()})
t['COL_FRP_rel2019']=t.COL_FRP/t.COL_FRP.loc[2019]
print((t/1).round(2).to_string())
clim=t.drop(2019).mean(); print('\nR (mean 2013-2023 excl 2019 / 2019): FRP COL',round(clim.COL_FRP/t.COL_FRP.loc[2019],3),' count COL',round(clim.COL_count/t.COL_count.loc[2019],3),' FRP domain',round(clim.domain_FRP/t.domain_FRP.loc[2019],3))
print('rank of 2019 in COL FRP:', int(t.COL_FRP.rank(ascending=False).loc[2019]),'of',len(t))
