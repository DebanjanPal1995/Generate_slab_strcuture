#@ Written by Debanjan Pal (2023), this code is used to generate slab structre, this code incorprates the track-point-through-topologies.ipynb module 
# uploaded in Github by Lauren Ilano (https://github.com/GPlates/pygplates-tutorials/blob/master/notebooks/track-point-through-topologies.ipynb) 

import pygplates
import numpy as np
import matplotlib.pyplot as plt
import pandasas pd
#from sphere_to_cartesian import *
import project_slab_to_depth as slab_geom
import sys
from tabulate import tabulate


#File_data1 = np.loadtxt("Line1.txt", dtype=float)
#File_data2 = np.loadtxt("Line2.txt", dtype=float)

#a=get_bearing_dip_depth(File_data1, File_data2)


##### the function traces the point through time ############3#
def track_point(rotation_filename, input_topology_filename, time_step, oldest_seed_time, start_time, end_time, File_data, PlateID):    
    
    topology_features = pygplates.FeatureCollection(input_topology_filename)
    rotation_model = pygplates.RotationModel(rotation_filename)
    
    # Empty array for storing Long/Lat of 
    coordinates = []
    time_list = np.arange(oldest_seed_time,time_step,-time_step)

    for seed_time in time_list:

        # Location of seed point for Kerguelenc
        #seed_geometry = pygplates.PointOnSphere(-50, 80)
    
        # Seed point for Hawaii
        seed_geometry = pygplates.MultiPointOnSphere([(lat,lon) for lat, lon in File_data])

    
        for time in np.arange(seed_time,end_time,-time_step):
        #print max_time, time

            # Get the plate polygons for this time
            resolved_topologies = []
            pygplates.resolve_topologies(topology_features, rotation_model, resolved_topologies, time)

            # make plate partitioner from polygons
            plate_partitioner = pygplates.PlatePartitioner(resolved_topologies, rotation_model)

            # Find the plate id of the polygon that contains the point
            partitioned_inside_geometries = []
            plate_partitioner.partition_geometry(seed_geometry, partitioned_inside_geometries)
            #PlateID = partitioned_inside_geometries[0][0].get_feature().get_reconstruction_plate_id()
            #PlateID = 550

            #print PlateID

            # Get the stage rotation that will move the point from where it is at the current time
            # to its location at the next time step, based on the plate id that contains the point at the 
            # current time
            stage_rotation = rotation_model.get_rotation(time-time_step, PlateID, time, anchor_plate_id=)

            # use the stage rotation to reconstruct the tracked point from position at current time 
            # to position at the next time step
            seed_geometry = stage_rotation * seed_geometry
    
        print('seed time = %d, plume is within plate %i' % (seed_time, PlateID))
        #point_longitude.append(seed_geometry.to_lat_lon_list().get_longitude())
        coordinates.append(seed_geometry.to_lat_lon_list())

    a=np.array(coordinates)
    return(a)



##### define function to isolate lat & lon, this function is meant for plotting
def lat_lon(coord, no_of_time_steps,no_of_points):

    for x in range(no_of_points):
        globals()['lons%s' %x] = [coord[i][x][1] for i in range(no_of_time_steps)]
        globals()['lats%s' %x] = [coord[i][x][0] for i in range(no_of_time_steps)]

    coordinates=[[eval("lons" + str(x)), eval("lats" + str(x))] for x in range (no_of_points)]
    return(coordinates)

## Input Parmeters
rotation_filename = ''
input_topology_filename = ''
time_step = 
oldest_seed_time = .
start_time=
end_time= 
#### load the gmt file which has the points
file = np.loadtxt("Pheonix.txt", dtype=float) ##### This I am initializng to make the dimenisional similarity for file and File_data2, since file will be
# used to track ages from nc file (has format lon, lat) and for file (has format lat lon).
File_data2 = np.loadtxt("Pheonix.txt", dtype=float) 

#@ Call track_point function to get the lats and lons through time
# Izanagi: 926, Pheonix:919, Farallon: 902, Meso-Tethys:530
track_points2 = track_point(rotation_filename,input_topology_filename, time_step, oldest_seed_time, start_time, end_time, File_data2, PlateID=) 

#@ Age calculation is here
#@ sample age from the age file
import pygmt
import netCDF4 as nc    
import math
total_time_steps=int((start_time-end_time)/time_step + 1)
ages=np.arange(start=end_time,stop=start_time+time_step,step=time_step,dtype=int)
ages1=np.flip(ages)
avg_age=[]
#print(ages1)
file[:,[0,1]]=File_data2[:,[1,0]] ## interchaning columns to extract ages 

#@ This is written to track the ages at every points
#for i in range(141, 166, 1): ## this will calculate the age of the subducting slabs from 141 to 165, 140 is not required since it is at the surface
    #print(i)
    #trace_ages = np.array(pygmt.grdtrack(points=file, grid='./Age/Y19_seafloor_age_mask_'+str(i)+'.0Ma.nc')) ## extract the ages at the coordinates 
    #print("Information", pygmt.grdinfo(grid='./Age/Y19_seafloor_age_mask_'+str(i)+'.0Ma.nc'))
    #seed_point_ages=trace_ages[:,2]
    #seed_point_ages = [0 if math.isnan(x) else x for x in seed_point_ages] ## to replace NaN values with 40 Ma
    #avg_age1=sum(seed_point_ages)/len(seed_point_ages)
    #avg_age.append((sum(seed_point_ages)/len(seed_point_ages))) ## average age at reconstruction end time
    #print(" The average subducting age at time", avg_age1, i )
    #corrected_ages = [x - 40 for x in ages1] ## 40 is subtracted since I am calculating the slab location at 40 Ma
    #new_ages=np.array([i + avg_age for i in corrected_ages]) ## these represents the deeper the slab is older the slab (Simpified assumption)

#@ This part of the code calculates the duration of subduction
corrected_ages = [x - end_time for x in ages1]
duration_of_subduction=np.flip(corrected_ages)


#age_array=np.array(avg_age)

points2=lat_lon(track_points2, total_time_steps,len(File_data2))



#### create figure
import geopandas as gpd
fig=pygmt.Figure()
#pygmt.makecpt(cmap="seis", series=[new_ages.min(), new_ages.max(), 1])
pygmt.makecpt(cmap="seis", series=[140, 170, 1])
#fig.basemap(region = [220, 300, 10, 60], projection="M8c", frame=True)
fig.basemap(region="g", projection="W12/12c", frame=True)
fig.coast(shorelines="0.8p,white")

##### loop over the data points
#for i in range(len(File_data2)):
    #traced_age=seed_point_ages[i] ##
    #new_ages=np.array([i + traced_age for i in corrected_ages]) ##
#    fig.plot(x=points2[i][0], y=points2[i][1],
#    fill=ages1,
#    style="c0.08c",
#    cmap=True,
#    pen="black",    
#   )
#for i in range(len(File_data2)):
 #   traced_age=seed_point_ages[i]
 #   new_ages=np.array([i + traced_age for i in corrected_ages])
  #  pygmt.makecpt(cmap="seis", series=[0, 80, 10])
  #  fig.plot(x=points2[i][0], y=points2[i][1],
  #  fill=new_ages,
  # style="c0.05c",
  #  cmap=True,
  #  pen="black",    
  #  )
    #exit()
    #print(new_ages)
#rgi = gpd.read_file('reconstructed_140.00Ma.shp')
#fig.plot(data=rgi, pen="0.5p,black")
#fig.colorbar(frame='af+l"Age (Ma)"')
#fig.show()

#exit()


#exit()
#@ This is to call the present-day trech location
#trench_loc=np.asarray([track_points2[total_time_steps-1],track_points2[total_time_steps-1]])
coord1 = np.asarray((track_points2[total_time_steps-1])) 

#@ call project_slab_to_depth to calculate new Lat/long when projected to depth
newp=[]
extractp=[]

for i in reversed(range(total_time_steps-1)): ### this loop runs from youngest age to oldest, i.e. from trench to depth
    #slab_loc=np.asarray([points2[i][1],points2[i][0]])
    coord2=np.asarray((track_points2[i])) # coord2 calls the list of co-ordinates for a particular age
    #current_age=avg_age[i]
    clenewp=np.asarray(slab_geom.get_bearing_dip_depth(coord1,coord2,corrected_ages[i]))
    get_depth=int(slab_geom.get_depth(coord1,coord2,corrected_ages[i]))
    #@ extract the projected points
    #print("Trench location", coord1)
    #np.savetxt('Trench_location', coord1)
    #print("Original coordiantes", coord2)
    #np.savetxt('Farallon_surface'+str(get_depth)+".dat", coord2)
    #np.savetxt('new location', clenewp.transpose())
    #print("Slab loc", clenewp.transpose())
    #exit()
    np.savetxt("Pheonix_depth"+str(i)+".dat", clenewp.transpose())
    #exit()


    

exit()



#### create figure
import geopandas as gpd
fig=pygmt.Figure()
#pygmt.makecpt(cmap="seis", series=[corrected_ages.min(), corrected_ages.max(), time_step])
#pygmt.makecpt(cmap="seis", series=[0, 120, 5])
#fig.basemap(region = [220, 300, 10, 60], projection="M8c", frame=True)
fig.basemap(region="g", projection="W12/12c", frame=True)
fig.coast(shorelines="0.8p,white")

##### loop over the data points
#for i in range(len(File_data1)):
    #fig.plot(x=points1[i][0], y=points1[i][1],
    #fill=ages1,
    #style="c0.08c",
    #cmap=True,
    #pen="black",    
    #)
for i in range(len(File_data2)):
    traced_age=seed_point_ages[i]
    new_ages=np.array([i + traced_age for i in corrected_ages])
    pygmt.makecpt(cmap="seis", series=[0, 80, 10])
    fig.plot(x=points2[i][0], y=points2[i][1],
    fill=new_ages,
    style="c0.05c",
    cmap=True,
    pen="black",    
    )
    #exit()
    #print(new_ages)
rgi = gpd.read_file('reconstructed_40.00Ma.shp')
fig.plot(data=rgi, pen="0.5p,black")
fig.colorbar(frame='af+l"Age (Ma)"')
fig.show()




