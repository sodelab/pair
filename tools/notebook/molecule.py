#!/usr/bin/env python 

import os,sys
import cirpy
import element
from math import *
from copy import *

class molecule:
    def __init__(self,name,nmols):
        print("Intialized new molecule class: %s"%(name))
        self.name = name
        self.nmols = nmols
        self.coord = []
        self.mass = []
        self.atom = []
        self.frag = []

        return

    def resolve(self):
        res=cirpy.resolve(self.name, 'xyz')
        print("Resolved chemical identification of %s"%(self.name))
        r=res.split("\n")
        self.natoms=int(r[0])
        self.frag=[0]*self.natoms

        for i,line in enumerate(r):
            if i<=1:
                continue
            elif i<self.natoms+2:
                tmp=[float(line.split()[1]),
                     float(line.split()[2]),
                     float(line.split()[3])]
                self.coord.append(tmp)
                self.atom.append(line.split()[0])
                self.mass.append(element.ELEMENTS[self.atom[-1]].mass)
            else:
                continue
        return
    
    def init_dimer(self,m1,m2):
        
        
        self.name=m1.name+"-"+m2.name
        self.nmols=2
        self.coord=deepcopy(m1.coord)+deepcopy(m2.coord)
        self.mass=deepcopy(m1.mass)+deepcopy(m2.mass)
        self.natoms=m1.natoms+m2.natoms
        self.atom=deepcopy(m1.atom)+deepcopy(m2.atom)
        self.frag=[0]*m1.natoms+[1]*m2.natoms

        return
    
    def xyz(self):
        
        string="%d\n\n"%(self.natoms)
        for c,coord in enumerate(self.coord):
            string+="%-2s "%(self.atom[c])
            string+="%15.10f %15.10f %15.10f\n"%(coord[0],
                                                 coord[1],
                                                 coord[2])
            
        return string
        
    def print_xyz(self):
        
        string="%d\n\n"%(self.natoms)
        for c,coord in enumerate(self.coord):
            string+="%-2s "%(self.atom[c])
            string+="%15.10f %15.10f %15.10f\n"%(coord[0],
                                                 coord[1],
                                                 coord[2])
        print(string)
        return
    
    def calculate_com(self):
    
        if self.nmols==1:
            x=0.0
            y=0.0
            z=0.0

            totalmass=0.0

            for a,atom in enumerate(self.coord):

                x+=atom[0]*self.mass[a]
                y+=atom[1]*self.mass[a]
                z+=atom[2]*self.mass[a]
                totalmass+=self.mass[a]

            x=x/totalmass
            y=y/totalmass
            z=z/totalmass

            self.com=[x,y,z]

        elif self.nmols==2:

            x1=x2=0.0
            y1=y2=0.0
            z1=z2=0.0

            totalmass1=totalmass2=0.0

            for c,coord in enumerate(self.coord):

                if self.frag[c]==0:
                    x1+=coord[0]*self.mass[c]
                    y1+=coord[1]*self.mass[c]
                    z1+=coord[2]*self.mass[c]
                    totalmass1+=self.mass[c]
                elif self.frag[c]==1:
                    x2+=coord[0]*self.mass[c]
                    y2+=coord[1]*self.mass[c]
                    z2+=coord[2]*self.mass[c]
                    totalmass2+=self.mass[c]

            x1=x1/totalmass1
            y1=y1/totalmass1
            z1=z1/totalmass1
                   
            x2=x2/totalmass2
            y2=y2/totalmass2
            z2=z2/totalmass2

            self.com=[x1,y1,z1,x2,y2,z2]
                   
        return
   
    def calculate_vector(self):
        """Function to calculate normalized vector."""

        self.calculate_com()
        
        for l in range(3):
            x=self.com[3]-self.com[0]
            y=self.com[4]-self.com[1]
            z=self.com[5]-self.com[2]
            
        norm=sqrt(x*x + y*y + z*z)
        self.vector=[x/norm,y/norm,z/norm]
        self.radius=norm
        
        return

    
    def center(self):

        self.calculate_com()
        f=self.frag
                   
        for c,coord in enumerate(self.coord):
            coord[0]-=self.com[0+f[c]*3]
            coord[1]-=self.com[1+f[c]*3]
            coord[2]-=self.com[2+f[c]*3]         

        return

    def rotate_theta(self,theta,mID):

        if mID!=0 and self.nmols==1:
            print("This is a single molecule. The molecule ID should be zero.")  
            
        for c,coord in enumerate(self.coord):
            
            if self.frag[c]==mID:
                x=coord[0]
                y=coord[1]
                z=coord[2]

                r=sqrt(x*x + y*y + z*z)
                if r<1.0e-10:
                    continue
                    
                t=acos(z/r)
                print(t)
                p=atan2(y,x)
                #t-=radians(theta)

                coord[0]=r*sin(t)*cos(p)
                coord[1]=r*sin(t)*sin(p)
                coord[2]=r*cos(t)

        return
    
    def rotate_phi(self,phi,mID):
        
        if mID!=0 and self.nmols==1:
            print("This is a single molecule. The molecule ID should be zero.")  
            
        for c,coord in enumerate(self.coord):
            
            if self.frag[c]==mID:
                x=coord[0]
                y=coord[1]
                z=coord[2]

                r=sqrt(x*x + y*y + z*z)
                if r<1.0e-10:
                    continue
                    
                t=acos(z/r)
                p=atan2(y,x)
                p+=radians(phi)

                coord[0]=r*sin(t)*cos(p)
                coord[1]=r*sin(t)*sin(p)
                coord[2]=r*cos(t)

        return        

    def position(self,R,theta,phi):
        
        if self.nmols==1:
            print("This is a single molecule. The `position_molecule` function is for dimers.")        

        x=R*sin(radians(theta))*cos(radians(phi))
        y=R*sin(radians(theta))*sin(radians(phi))
        z=R*cos(radians(theta))

        for c,coord in enumerate(self.coord):
            
            if self.frag[c]==1:
                coord[0]+=x
                coord[1]+=y
                coord[2]+=z

        return    
    
    def setup_dimer(self):
        """This function that each monomer is initially at the origin, before
        transformation."""
        
        self.calculate_com()
        #self.R=0.0 # initial dimer distance
        #self.tv=[0.0]*3 # initial dimer vector
        self.center()
        return
        
    def transform_dimer(self,value,transform,mID):
        """This function prepares the dimer according to the type and magnitude 
        of transformation. Since translations shift the dimer off the origin, this
        operation is always perfomed last in the `set dimer` function."""
        
        if transform=="translation": # make translational vector
            self.radius=value
            
        elif transform=="phi" or transform=="theta": # rotate molecules
            
            for c,coord in enumerate(self.coord):
            
                if self.frag[c]==mID:
                    x=coord[0]
                    y=coord[1]
                    z=coord[2]

                    r=sqrt(x*x + y*y + z*z)
                    if r<1.0e-10:
                        continue

                    if transform=="phi":
                        p=atan2(y,x)
                        p+=radians(value)
                        
                        coord[0]=r*sin(0)*cos(p)
                        coord[1]=r*sin(0)*sin(p)
                        coord[2]=r*cos(0)
                        
                    elif transform=="theta":
                        p=atan2(y,x)
                        p+=radians(value)
                        
                        coord[0]=r*sin(0)*cos(p)
                        coord[1]=r*sin(0)*sin(p)
                        coord[2]=r*cos(0)
            
        else:
            raise IOError("Unknown transformation.")
            
        return
            
    def set_dimer(self):
        """"""
        
        for c,coord in enumerate(self.coord):
            if self.frag[c]==1:
                coord[0]+=self.vector[0]*self.radius
                coord[1]+=self.vector[1]*self.radius
                coord[2]+=self.vector[2]*self.radius
                
        return

if __name__=="__main__":
    pass
