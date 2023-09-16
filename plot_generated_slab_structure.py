#### create figure
import geopandas as gpd
import pygmt
import numpy as np
import pandas as pd
new_ages=[57, 47, 37]
lat=[]
lon=[]
lat1=[]
lon1=[]
lat2=[]
lon2=[]

data = np.loadtxt("Pheonix_depth392.dat", dtype=float)
orid_data = np.loadtxt("Pheonix.txt", dtype=float)
print(data)
print(orid_data)
#data1 = np.loadtxt("Pheonix_depth214.dat", dtype=float)
#print(data)
#print(data[0][1])

#exit()

#### create figure
import geopandas as gpd
fig=pygmt.Figure()
#pygmt.makecpt(cmap="seis", series=[corrected_ages.min(), corrected_ages.max(), time_step])
#pygmt.makecpt(cmap="seis", series=[0, 120, 5])
#fig.basemap(region = [220, 300, 10, 60], projection="M8c", frame=True)
fig.basemap(region="g", projection="W12/12c", frame=True)
fig.coast(shorelines="0.8p,white")

##### loop over the data points
for i in range(len(data)):
    lat.append(data[i][0])
    lon.append(data[i][1])
    lat1.append(orid_data[i][0])
    lon1.append(orid_data[i][1])
    #lat2.append(data1[i][0])
    #lon2.append(data1[i][1])
#print(lat)
#print(lon)

fig.plot(x=lon, y=lat, style="c0.07c", fill="red", pen="black")
fig.plot(x=lon1, y=lat1, style="c0.04c", fill="green", pen="black")
#fig.plot(x=lon2, y=lat2, style="c0.05c", fill="blue", pen="black")
#for i in range(len(data)):
#   traced_age=seed_point_ages[i]
#    new_ages=np.array([i + traced_age for i in corrected_ages])
#    pygmt.makecpt(cmap="seis", series=[0, 80, 10])
#    fig.plot(x=points2[i][0], y=points2[i][1],
#    fill=new_ages,
#    style="c0.05c",
#    cmap=True,
#    pen="black",    
#    )
    #exit()
    #print(new_ages)
print("I am here")
rgi = gpd.read_file('reconstructed_40.00Ma.shp')
fig.plot(data=rgi, pen="0.5p,black")
fig.colorbar(frame='af+l"Age (Ma)"')
fig.show()


