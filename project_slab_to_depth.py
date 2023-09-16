## Written by Debanjan Pal,2023. This code outputs new Lat/Long given the slab length and dip.
#  This codes first computes the shortest distance between two points over the sphere. One simplified assumption is that if Lat/Long diff between points is less than
# 70 degrees only creates only 4 degree difference in calcualting distance over sphere or straight line joining theses points.
# Therefore to take the vertical and horizontal projection of the slab I apply trigonometric relations in cartesian co-oridnates since the error is only 4 degree if
# the difference is 70 degrees.

import pygplates
import numpy as np
import math
import haversine as hs
from haversine import Unit
import pyproj
import time
from math import asin, atan2, cos, degrees, radians, sin
import sys

def get_bearing_dip_depth(coord1, coord2, ages):
    slab_length=0
    proj_hori_slab_length=0
    proj_vert_slab_depth=0
    lat3=[]
    lon3=[]
    lat2=[]
    lon2=[]
    add_depth=0
    slab_depth=[]
    store_depths=[]
    store_depths1=[]
    depth=0
    j=0
    count=0
    R=6371
   

    for i in range(len(coord1)):
        geodesic = pyproj.Geod(ellps='WGS84')
        fwd_azimuth,back_azimuth,distance = geodesic.inv(coord1[i][1], coord1[i][0], coord2[i][1], coord2[i][0])

        if (fwd_azimuth > 0 and fwd_azimuth < 180):
            azimuth=fwd_azimuth
    
        else:
            azimuth=360+fwd_azimuth
        #print(azimuth)    
        #@ calculate distace in km between two points 
        point1 = pygplates.PointOnSphere(coord1[i][0], coord1[i][1])
        point2 = pygplates.PointOnSphere(coord2[i][0], coord2[i][1])
        slab_length=pygplates.GeometryOnSphere.distance(point1,point2)*6371
        #a=pygplates.GeometryOnSphere.distance(point1,point2)*6371
        #print("First", a)
        #print("slab length",slab_length)
        #earth_mean_radius_in_kms = pygplates.Earth.mean_radius_in_kms
        #@ calculate dip and depth from age of the slabs (Hu et al., 2022)
        #@ workflow 1) Find the average age of the slab predicted at every time_step
        #           2) αs = 42 − 0.067 × Asub − 7.2 × OPN (shallow dip between 0 to 125 km) 
        #           3) αd = 73 − 0.07 × Asub − 7.1 × OPN  (steep dip between 125 to 410 km)
        
        alpha_sh = 42-0.067*ages-7.2
        depth=math.sin(math.radians(alpha_sh))*slab_length
        #print(alpha_sh, depth, slab_length,ages)
        #print("depth", depth, slab_length)
        #alpha_st = 73-0.07*ages-7.1
        #proj_vert_slab_length2=math.sin(math.radians(alpha_sh))*slab_length[i]

        if (depth < 125):
         #alpha_sh = 42-0.067*ages-7.2
         #@ calculate horizontal projected distance
         #proj_hori_slab_length.append(math.cos(math.radians(alpha_sh))*slab_length[i])
         proj_hori_slab_length=math.cos(math.radians(alpha_sh))*slab_length
         #@ calculate vertical projected distance
         #proj_vert_slab_depth.append(math.sin(math.radians(alpha_sh))*slab_length[i])
         proj_vert_slab_depth=math.sin(math.radians(alpha_sh))*slab_length
         #print("The dip is shallow", alpha_sh, slab_length)
         store_depths.append(proj_vert_slab_depth) ### this will be used finally to average out the ages
         print("The dip is shallow", i, proj_vert_slab_depth)
            
                         
        if (depth > 125 and depth < 400):
            alpha_st = 73-0.07*ages-7.1
            #@ calculate horizontal projected distance
            #proj_hori_slab_length.append(math.cos(math.radians(alpha_sh))*slab_length[i])
            proj_hori_slab_length=math.cos(math.radians(alpha_st))*slab_length
            #@ calculate vertical projected distance
            #proj_vert_slab_depth.append(math.sin(math.radians(alpha_sh))*slab_length[i])
            proj_vert_slab_depth=math.sin(math.radians(alpha_st))*slab_length
            #print("The dip is shallow", alpha_sh, depth)
            #depth = math.sin(math.radians(alpha_sh))*slab_length[i]
            store_depths.append(proj_vert_slab_depth)
            #print("This is for deep depths", store_depths1)
            #add_depth=add_depth+proj_vert_slab_depth
            #count=count+1
            #depth.append(math.sin(math.radians(alpha_sh))*slab_length[i])
            print("The dip is steep", i, proj_vert_slab_depth)  
            
            



        #@ This function computes the new points on a sphere at a given distance
        #for i in range(len(coord1)):   
        R=6371-proj_vert_slab_depth
        lat1 = radians(coord1[i][0])
        lon1 = radians(coord1[i][1])
        #print("latt", lat1, lon1)
        a = radians(azimuth)

        #print("It is here", proj_hori_slab_length_array[i])
        #exit()
        #current_slab_length=proj_hori_slab_length_array[i]
        #lat: initial latitude, in degrees
        #lon: initial longitude, in degrees
        #d: target distance from initial
        #bearing: (true) heading in degrees
        #R: optional radius of sphere, defaults to mean radius of earth
        #Returns new lat/lon coordinate {d}km from initial, in degrees
        #R=6371
        lat2.append(asin(sin(lat1) * cos(proj_hori_slab_length/R) + cos(lat1) * sin(proj_hori_slab_length/R) * cos(a)))
        lon2.append(lon1 + atan2(
            sin(a) * sin(proj_hori_slab_length/R) * cos(lat1),
            cos(proj_hori_slab_length/R) - sin(lat1) * sin(lat2[i])
        ))
            #return (degrees(lat2), degrees(lon2),)    
        #print("Cal checkl", degrees(lat2[i]), degrees(lon2[i]))

    #original_stdout = sys.stdout # Save a reference to the original standard output
    #with open("NA"+str(ages) + ".dat", 'w', encoding='utf-8') as f:
    #    sys.stdout = f # Change the standard output to the file we created.
    #    print(np.degrees(lat2), np.degrees(lon2))
    #    sys.stdout = original_stdout # Reset the standard output to its original value  
    if(len(store_depths) > 0):  
        #final=np.concatenate((store_depths,store_depths1))
        print(store_depths, len(store_depths))
        print("slab depths", sum(store_depths)/len(store_depths))
        #print(add_depth/count, add_depth, count)
    else:
        exit()
    return(np.degrees(lat2),np.degrees(lon2), store_depths)
    #return(lon2,lat2,proj_vert_slab_length,slab_depth)
    #return (math.degrees(lon2), math.degrees(lat2))



#@ I have written this fuction to call the depths for generating files, need to modify later to make it simple
#but works for now
#def get_depth(coord1, coord2, ages):
    slab_length=[]
    proj_hori_slab_length=0
    proj_vert_slab_depth=0
    lat3=[]
    lon3=[]
    lat2=[]
    lon2=[]
    slab_depth=[]
    store_depths=[]
    depth=0
    j=0
    count=0
    R=6371
   

    for i in range(len(coord1)):
      
        #@ calculate distace in km between two points 
        point1 = pygplates.PointOnSphere(coord1[i][0], coord1[i][1])
        point2 = pygplates.PointOnSphere(coord2[i][0], coord2[i][1])
        slab_length=pygplates.GeometryOnSphere.distance(point1,point2)*6371
        #a=pygplates.GeometryOnSphere.distance(point1,point2)*6371
        #print("First", a)
        #print("slab length",slab_length)
        #earth_mean_radius_in_kms = pygplates.Earth.mean_radius_in_kms
        #@ calculate dip and depth from age of the slabs (Hu et al., 2022)
        #@ workflow 1) Find the average age of the slab predicted at every time_step
        #           2) αs = 42 − 0.067 × Asub − 7.2 × OPN (shallow dip between 0 to 125 km) 
        #           3) αd = 73 − 0.07 × Asub − 7.1 × OPN  (steep dip between 125 to 410 km)
        
        alpha_sh = 42-0.067*ages-7.2
        depth=math.sin(math.radians(alpha_sh))*slab_length
        #print("depth", depth, slab_length)
        #alpha_st = 73-0.07*ages-7.1
        #proj_vert_slab_length2=math.sin(math.radians(alpha_sh))*slab_length[i]
    
        if (depth < 125):
            #alpha_sh = 42-0.067*ages-7.2
            #@ calculate horizontal projected distance
            #proj_hori_slab_length.append(math.cos(math.radians(alpha_sh))*slab_length[i])
            proj_hori_slab_length=math.cos(math.radians(alpha_sh))*slab_length
            #@ calculate vertical projected distance
            #proj_vert_slab_depth.append(math.sin(math.radians(alpha_sh))*slab_length[i])
            proj_vert_slab_depth=math.sin(math.radians(alpha_sh))*slab_length
            #print("The dip is shallow", alpha_sh, slab_length)
            store_depths.append(math.sin(math.radians(alpha_sh))*slab_length)
            #print("I am here", slab_depth)
            #print("The dip is shallow", depth)        
                         
        if (depth > 125 and depth < 400):
            alpha_st = 73-0.07*ages-7.1
            #@ calculate horizontal projected distance
            #proj_hori_slab_length.append(math.cos(math.radians(alpha_sh))*slab_length[i])
            proj_hori_slab_length=math.cos(math.radians(alpha_st))*slab_length
            #@ calculate vertical projected distance
            #proj_vert_slab_depth.append(math.sin(math.radians(alpha_sh))*slab_length[i])
            proj_vert_slab_depth=math.sin(math.radians(alpha_st))*slab_length
            #print("The dip is shallow", alpha_sh, depth)
            #depth = math.sin(math.radians(alpha_sh))*slab_length[i]
            store_depths.append(math.sin(math.radians(alpha_st))*slab_length)
            #depth.append(math.sin(math.radians(alpha_sh))*slab_length[i])
            #print("The dip is steep", depth)
  
    if(len(store_depths) > 0):  
        slab_depth= sum(store_depths)/len(store_depths)
        #print(slab_depth)
        return(slab_depth)
    else:
        exit()
