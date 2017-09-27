# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 14:02:24 2017

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
def box(start,azimuth_v,dist,azimuth_trench,length,depth):
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

######################################################################
#********************* Read file **************************
input=open("Plate_boundary_CA2.csv","r")
input.readline()
inputtable=input.readlines()
#print(inputtable)
num=0
for line in inputtable:
    num=num+1
    name="Plate_boundary_seismic_CA"+str(num)+".txt"
    fout=open(name,"w")    
    #print("1")
    print(line)
    line2=line.strip()
    boundary=line2.split(',')
    latpt=float(boundary[5])
    longpt=float(boundary[4])
    depthpt=float(boundary[14])
    depthpt=-0.001*depthpt
    azimuth_v=float(boundary[11])
    azimuth_trench=float(boundary[9])
    length=float(boundary[8])

#******************** Constants *************************

# for intermediate-depth earthquakes
#1943	:CO\NA	:SUB	5	-93.656	14.016	-94.055	14.222	48.8	298	71.9	33	-71.6	-6.5	-6048	149
#latpt=14.016
#longpt=-93.656
#depthpt=6.048
    start=np.array([latpt,longpt,depthpt])  #[lat,long,depth] deg,deg,km
    Maxd_MEQ=410  #km
    dip_angle=30
    depth=350  #km
#azimuth_v=33
#azimuth_trench=298
#length=48.8
    dip_slab=dip_angle*np.pi/180
    dist=Maxd_MEQ/np.tan(dip_slab)
#print(dist)
    azimuth_oppv=azimuth_v+180
    lengthouter=100
    torr=0.1
# for outer-rise earthquakes
# select earthquake events **************
    f=open("Central_America_ISC_45_InterMD_all.txt", "r")
#    fout=open("CA_selected1.txt","w")
    print('Lat','\t','Long','\t','Depth_km','\t','date','\t','time','\t','M',file=fout)
    f.readline()
    table=f.readlines()
#print(table)
#selected=[]
    boxc=box(start,azimuth_v,dist,azimuth_trench,length,depth)
    volarea=volume(boxc)
    #print(volarea)
    n=0
#print(volarea)
    for line in table:
        line2=line.strip()
        event=line2.split('\t')
        lat=float(event[0])
        long=float(event[1])
        D=float(event[2])
        date=event[3]
        time=event[4]
        M=float(event[5])
        event2=earthquake(lat,long,depth,date,time,M)
    #event2.loc_cartitian()
    #print(event2.loc_cartitian())
        boxnew=np.row_stack((boxc,event2.loc_cartitian()))
        #print(boxnew)
        volseis=volume(boxnew)
    #print(volseis)
        dvol=abs(volarea-volseis)/volarea
        #print(dvol)
        if determine(volarea,volseis,torr):
            #print('yes')
        #selected.append(earthquake(lat,long,D,date,time,M)
            n=n+1
            print(lat,'\t',long,'\t',D,'\t',date,'\t',time,'\t',M,file=fout)               
            #print(n)
        f.close()
        
    rate=n/length
    print('rate=',rate)                
        
    fout.close()


    




##************ Plot ************************
#
#xs=boxc[:,0]
#ys=boxc[:,1]
#zs=boxc[:,2]
##print(xs)
##print(ys)
##print(zs)
##print(type(xs))
##print(boxc)
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(xs, ys, zs, zdir='z')
#ax.set_xlabel('X Label')
#ax.set_ylabel('Y Label')
#ax.set_zlabel('Z Label')
#
#plt.show()


#pt2=pt1.nextpointsurf(azimuth_v,dist)
#pt2=Location(pt2,R0E)

#pt4=pt2.nextpointsurf(azimuth_v,dist)
#pt4=GenNxtPt(pt4)
#pt5=pt1.nextpointbelow(depth)
#pt6=pt2.nextpointbelow(depth)
#pt7=pt3.nextpointbelow(depth)
#pt8=pt4.nextpointbelow(depth)


#box=[]
#pt11=Location.cartitian(pt1)
#box.append(pt11)
#pt21=Location.cartitian(pt2)
#box.append(pt21)
#pt31=Location.cartitian(pt3)
#box.append(pt31)
#pt41=Location.cartitian(pt4)
#box.append(pt41)
#pt51=Location.cartitian(pt5)
#box.append(pt51)
#pt61=Location.cartitian(pt6)
#box.append(pt61)
#pt71=Location.cartitian(pt7)
#box.append(pt71)
#pt81=Location.cartitian(pt8)
#box.append(pt81)
#
#vol=volume(box)
#print(vol)
## Test         
##*************** Generate Array *************************
##test
##origin
#lat1=-20.0
#long1=-10.0
#depth1=30.0
#point1=Location(lat1,long1,depth1,R0E)
#print(point1.cartitian())
#point11=Location(lat1,long1,depth1+20,R0E)
#print(point11.cartitian())
#
##per to origin
#azmth=33
#dist=100
#pt2=point1.nextpoint(33,100)
#point2=Location(pt2[0],pt2[1],pt2[2],R0E)
#print(point2.cartitian())
#point21=Location(pt2[0],pt2[1],pt2[2]+20,R0E)
#print(point21.cartitian())
#
##along boundary
#azmth=301
#dist=44
#pt3=point1.nextpoint(301,44)
#point3=Location(pt3[0],pt3[1],pt3[2],R0E)
#print(point3.cartitian())
#point31=Location(pt3[0],pt3[1],pt3[2]+20,R0E)
#print(point31.cartitian())
#
##last point
#azmth=33
#dist=100
#pt4=point3.nextpoint(33,100)
#point4=Location(pt4[0],pt4[1],pt4[2],R0E)
#print(point4.cartitian())
#point41=Location(pt4[0],pt4[1],pt4[2]+20,R0E)
#print(point41.cartitian())
#
#
#
##nxt=point1.nextpoint(azmth,dist)
##print(point1.nextpoint(azmth,dist))
##[lat,long]=nxt.split(',')
##print(lat)
##print(type(float(long)))
##print(type(Decimal(point1.nextpoint(azmth,dist))))
##print(geopy.point.Point.format_decimal(point1.nextpoint(azmth,dist)))
#
#
#
##************** Determine if inside *********************
##points=np.array([[0,0,0],[2,0,0],[0,2,0],[0,0,2],[2,2,0],[2,0,2],[0,2,2],[2,2,2],[1,1,1]])                 
##volume=ConvexHull(points).volume
##print(volume)
#
#points=np.array([point1.cartitian(),point2.cartitian(),point3.cartitian(),point4.cartitian(),point11.cartitian(),point21.cartitian(),point31.cartitian(),point41.cartitian()])
#print(points)
#volume=ConvexHull(points).volume
#print(volume)
#
#
##************ Plot ************************
#xs=points[:,0]
#ys=points[:,1]
#zs=points[:,2]
#print(type(xs))
#print(xs)
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(xs, ys, zs, zdir='z')
#ax.set_xlabel('X Label')
#ax.set_ylabel('Y Label')
#ax.set_zlabel('Z Label')
#
#plt.show()

##***************** Generate next points ******************
#class GenNxtPt:
#    def __init__(self,start):
#        self.lat=start[0]
#        self.long=start[1]
#        self.depth=start[2]
#    
#    def nextpointsurf(self,azimuth,distance):
#        # given: lat, lon, azimuth = azimuth in degrees, d = distance in kilometers
#        origin = geopy.Point(self.lat, self.long)
#        nxtpoint = geopy.point.Point.format_decimal((VincentyDistance(kilometers=distance).destination(origin, azimuth)))
#        [latnxt,longnxt]=nxtpoint.split(',')
#        nxtpoint2 = np.array([float(latnxt),float(longnxt),self.depth])
#        return nxtpoint2
#
#    def nextpointbelow(self,depth):
#        lat=self.lat
#        long=self.long
#        depth=self.depth+depth
#        nxtpoint=np.array([lat,long,depth])
#        return nxtpoint