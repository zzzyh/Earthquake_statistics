# -*- coding: utf-8 -*-
"""
Created on Sun Sep 24 00:58:02 2017

@author: HuY
"""

import numpy as np
import geopy
from geopy.distance import VincentyDistance
from scipy.spatial import ConvexHull
#import matplotlib.pyplot as plt

R0E=6371 #km
#*************** Sepherical Spatial to Cartision *********************
class Location:
    R0=R0E
    def __init__(self,point):
        self.lat=point[0]
        self.long=point[1]
        self.depth=point[2]
    
    def latlong(self):
        self.ll=np.array([self.lat,self.long,self.depth])
        return self.ll
    
    def cartitian(self):
        if self.lat>0:
            self.theta=(90-self.lat)*np.pi/180
        else:
            self.theta=(self.lat+90)*np.pi/180              
        if self.long>0:
            self.phi=self.long*np.pi/180
        else:
            self.phi=(self.long+180)*np.pi/180        
            
        self.r=self.R0-self.depth 
            
        x=self.r*np.sin(self.theta)*np.cos(self.phi)
        y=self.r*np.sin(self.theta)*np.sin(self.phi)
        z=self.r*np.cos(self.theta) 
            
        self.xyz=np.array([x, y, z])
        return self.xyz
    
    def nextpointsurf(self,azimuth,distance):
        # given: lat, lon, azimuth = azimuth in degrees, d = distance in kilometers
        origin = geopy.Point(self.lat, self.long)
        #print(origin)
        nxtpoint=geopy.point.Point((VincentyDistance(kilometers=distance).destination(origin, azimuth)))
        #print(nxtpoint)
        nxtpointdc = nxtpoint.format_decimal()
        [latnxt,longnxt]=nxtpointdc.split(',')
        nxtpoint2 = Location([float(latnxt),float(longnxt),self.depth])
        return nxtpoint2

    def nextpointbelow(self,depth):
        lat=self.lat
        long=self.long
        newdepth=self.depth+depth
        nxtpoint=Location([float(lat),float(long),float(newdepth)])
        return nxtpoint    
    
#************ functions **********************
def volume(points):
    vol=ConvexHull(points).volume
    return vol
    
def determine(vol1,vol2,torr):
    if abs(vol1-vol2)/vol1 <=torr:
        return True
    else:
        return False
    
#************* seismic data structure *************
class earthquake:
    def __init__(self,lat,long,depth,date,time,M):
        self.lat=lat
        self.long=long
        self.depth=depth
        self.date=date
        self.time=time
        self.M=M
        
    def loc_cartitian(self):
         ll=np.array([self.lat,self.long,self.depth])
         ll2=Location(ll)
         self.cartitian=ll2.cartitian()
         return self.cartitian
     
#************* select box *******************
def box1(start,azimuth_v,dist,azimuth_trench,length,depth):
    pt1=Location(start)
    pt1s=pt1.latlong()
    pt1c=pt1.cartitian()
    boxs=pt1s
    boxc=pt1c
    
    pt2=pt1.nextpointsurf(azimuth_v,dist)
    pt2s=pt2.latlong()
    pt2c=pt2.cartitian()
    boxs=np.row_stack((boxs,pt2s))
    boxc=np.row_stack((boxc,pt2c))

    pt3=pt1.nextpointsurf(azimuth_trench,length)
    pt3s=pt3.latlong()
    pt3c=pt3.cartitian()
    boxs=np.row_stack((boxs,pt3s))
    boxc=np.row_stack((boxc,pt3c))

    pt4=pt3.nextpointsurf(azimuth_v,dist)
    pt4s=pt4.latlong()
    pt4c=pt4.cartitian()
    boxs=np.row_stack((boxs,pt4s))
    boxc=np.row_stack((boxc,pt4c))

    pt5=pt1.nextpointbelow(depth)
    pt5s=pt5.latlong()
    pt5c=pt5.cartitian()
    boxs=np.row_stack((boxs,pt5s))
    boxc=np.row_stack((boxc,pt5c))
    
    pt6=pt2.nextpointbelow(depth)
    pt6s=pt6.latlong()
    pt6c=pt6.cartitian()
    boxs=np.row_stack((boxs,pt6s))
    boxc=np.row_stack((boxc,pt6c))

    pt7=pt3.nextpointbelow(depth)
    pt7s=pt7.latlong()
    pt7c=pt7.cartitian()
    boxs=np.row_stack((boxs,pt7s))
    boxc=np.row_stack((boxc,pt7c))

    pt8=pt4.nextpointbelow(depth)
    pt8s=pt8.latlong()
    pt8c=pt8.cartitian()
    boxs=np.row_stack((boxs,pt8s))
    boxc=np.row_stack((boxc,pt8c))
    
    return boxc

def box2(start, end, azimuth_v, dist, depth):
    pt1=Location(start)
    pt1s=pt1.latlong()
    pt1c=pt1.cartitian()
    boxs=pt1s
    boxc=pt1c
    
    pt2=pt1.nextpointsurf(azimuth_v,dist)
    pt2s=pt2.latlong()
    pt2c=pt2.cartitian()
    boxs=np.row_stack((boxs,pt2s))
    boxc=np.row_stack((boxc,pt2c))
    
    pt3=Location(end)
    pt3s=pt3.latlong()
    pt3c=pt3.cartitian()
    boxs=np.row_stack((boxs,pt3s))
    boxc=np.row_stack((boxc,pt3c))

    pt4=pt3.nextpointsurf(azimuth_v,dist)
    pt4s=pt4.latlong()
    pt4c=pt4.cartitian()
    boxs=np.row_stack((boxs,pt4s))
    boxc=np.row_stack((boxc,pt4c))

    pt5=pt1.nextpointbelow(depth)
    pt5s=pt5.latlong()
    pt5c=pt5.cartitian()
    boxs=np.row_stack((boxs,pt5s))
    boxc=np.row_stack((boxc,pt5c))
    
    pt6=pt2.nextpointbelow(depth)
    pt6s=pt6.latlong()
    pt6c=pt6.cartitian()
    boxs=np.row_stack((boxs,pt6s))
    boxc=np.row_stack((boxc,pt6c))

    pt7=pt3.nextpointbelow(depth)
    pt7s=pt7.latlong()
    pt7c=pt7.cartitian()
    boxs=np.row_stack((boxs,pt7s))
    boxc=np.row_stack((boxc,pt7c))

    pt8=pt4.nextpointbelow(depth)
    pt8s=pt8.latlong()
    pt8c=pt8.cartitian()
    boxs=np.row_stack((boxs,pt8s))
    boxc=np.row_stack((boxc,pt8c))
    
    return boxc
