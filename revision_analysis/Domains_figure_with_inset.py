# Figure 1 (domains + city markers) with Aburrá Valley inset. Revision for R3 minor comment 4.
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.ticker as mticker, matplotlib.colors
from netCDF4 import Dataset
import cartopy.crs as ccrs, cartopy.io.shapereader as shpreader
from cartopy.feature import ShapelyFeature
from shapely.geometry.polygon import LinearRing
from matplotlib.patches import Rectangle
np.bool = np.bool_
P='../'  # run from Scripts_submit/
crs=ccrs.PlateCarree()
v=Dataset(P+'GMRTv3_6_20190709topo.grd')
x=np.arange(v['x_range'][0],v['x_range'][1]+0.5*v['spacing'][0],v['spacing'][0])
y=np.arange(v['y_range'][0],v['y_range'][1]+0.5*v['spacing'][1],v['spacing'][1])
X,Y=np.meshgrid(x,y); z=np.reshape(v['z'][:],X.shape)[::-1,:]
W=P+'World_Countries/'
shps=[W+'World_Countries.shp',W+'GDBIDEAM_departamentos.shp',W+'COL_no_dptos.shp']; lws=[0.1,0.15,0.3]
cities=[('Itagui',6.181957,-75.613805),('Caldas',6.085799,-75.635587),('Medellin',6.241542,-75.584248),
('Cali',3.419928,-76.510929),('Envigado',6.160965,-75.548687),('Bogota',4.649475,-74.108089),
('Bucaramanga',7.114798,-73.117709),('Barbosa',6.426635,-75.360970),('Copacabana',6.343290,-75.517163)]
markers=['o','*','s','v','^','>','<','p','P']
terrain=matplotlib.colors.LinearSegmentedColormap.from_list('t',plt.cm.terrain(np.linspace(0.3,1,200)))
water='#97b6e1'
def base(ax,ext,lw_scale=1,sub=1):
    m=(X>=ext[0]-1)&(X<=ext[1]+1)&(Y>=ext[2]-1)&(Y<=ext[3]+1)
    r=np.where(m.any(0))[0]; c=np.where(m.any(1))[0]
    xs=slice(r[0],r[-1]+1,sub); ys=slice(c[0],c[-1]+1,sub)
    cf=ax.contourf(X[ys,xs],Y[ys,xs],z[ys,xs],levels=np.arange(0,5006,100),cmap=terrain,extend='neither')
    ax.contourf(X[ys,xs],Y[ys,xs],z[ys,xs],levels=[-20000,0],colors=[water])
    ax.set_extent(ext,crs=crs)
    for k,s in enumerate(shps):
        ax.add_feature(ShapelyFeature(shpreader.Reader(s).geometries(),crs,facecolor='none',edgecolor='k',linewidth=lws[k]*lw_scale))
    return cf
fig=plt.figure(figsize=(10,8))
ax=fig.add_axes([0.05,0.08,0.58,0.84],projection=crs)
cf=base(ax,[-83,-64,-5,14],sub=2)
cax=fig.add_axes([0.645,0.18,0.018,0.64]); cb=fig.colorbar(cf,cax=cax); cb.set_label('Elevation (m)')
ring=LinearRing(list(zip([-66,-80,-80,-66],[-2,-2,13,13])))
ax.add_geometries([ring],crs,facecolor='none',edgecolor='k',linewidth=2)
gl=ax.gridlines(crs=crs,draw_labels=True,linewidth=0.75,color='k',alpha=0.5,linestyle='--')
gl.xlocator=mticker.FixedLocator([-68,-72,-76,-80]); gl.ylocator=mticker.FixedLocator([-5,0,5,10]); gl.right_labels=False
ax.text(-82.75,13.25,'D01',fontweight='bold'); ax.text(-79.75,12.25,'D02',fontweight='bold')
for (n,la,lo),mk in zip(cities,markers):
    ax.scatter(lo,la,marker=mk,s=60,edgecolor='red',facecolors='none',linewidth=1.5,label=n,transform=crs,zorder=5)
iext=[-75.85,-75.2,5.95,6.6]
ax.add_patch(Rectangle((iext[0],iext[2]),iext[1]-iext[0],iext[3]-iext[2],fill=False,ec='k',lw=1.2,transform=crs,zorder=6))
h,l=ax.get_legend_handles_labels()
fig.legend(h,l,frameon=False,ncol=2,loc='lower center',bbox_to_anchor=(0.85,0.615),fontsize=10.5,columnspacing=0.8,handletextpad=0.3,labelspacing=0.6)
# inset
ia=fig.add_axes([0.73,0.20,0.24,0.36],projection=crs)
base(ia,iext,lw_scale=4)
for (n,la,lo),mk in zip(cities,markers):
    if iext[0]<lo<iext[1] and iext[2]<la<iext[3]:
        ia.scatter(lo,la,marker=mk,s=80,edgecolor='red',facecolors='none',linewidth=1.8,transform=crs,zorder=5)
        dx,ha=(-0.03,'right') if n in ('Itagui','Caldas') else (0.03,'left')
        ia.text(lo+dx,la+0.012,n,fontsize=8,ha=ha,transform=crs,zorder=6)
g2=ia.gridlines(crs=crs,draw_labels=True,linewidth=0.5,color='k',alpha=0.4,linestyle='--')
g2.xlocator=mticker.FixedLocator([-75.7,-75.4]); g2.ylocator=mticker.FixedLocator([6.1,6.4])
g2.top_labels=False; g2.left_labels=False; g2.right_labels=True; g2.xlabel_style=g2.ylabel_style={'size':8}
ia.set_title('Aburrá Valley',fontsize=10,fontweight='bold')
fig.savefig('../Manuscript/Figures/Domains_MPI-C_WRF.jpg',dpi=300,bbox_inches='tight')
fig.savefig('../Manuscript/Figures/Domains_MPI-C_WRF.pdf',dpi=300,bbox_inches='tight')
