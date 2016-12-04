#!/usr/bin/env python

import numpy as np
from ROOT import TVector3
from ROOT import TRotation
import sys

class CloverType:
    Yale, Tohoku, IMP = range(3)

class clover(object):
    def __init__(self,r,theta,phi,clovertype=CloverType.Yale):
        self.xdispl = 0
        self.ydispl = 0
        self.placement = TVector3(1.,1.,1.)
        self.placement.SetMag(r)
        self.placement.SetTheta(theta*np.pi/180)
        self.placement.SetPhi(phi*np.pi/180)
        self.type = clovertype

        if clovertype == CloverType.Yale:
            displ = 2.2 #cm
            self.xdispl = displ
            self.ydispl = displ
            self.covergap = 1.0 #cm
            rcrystal = 2.5 #cm
            self.segdispl = 4*rcrystal/(3*np.pi)

        if clovertype == CloverType.Tohoku:
            self.xdispl = self.ydispl = 0
            self.covergap = 0
            rcrystal = 0
            self.segdispl = 4*rcrystal/(3*np.pi)

        if clovertype == CloverType.IMP:
            displ = 2.2 # needs verification
            self.xdispl = displ
            self.ydispl = displ
            self.covergap = 0.5 # needs verification
            self.xsegdispl = displ/2
            self.ysegdispl = displ/2

    def print_clover_placement(self,idxclover,coord=None):
        if coord == "Spherical":
             print("slot.{0:02d}.X.0.vec_sph: ".format(idxclover)
                  + str(np.around(self.placement.Mag(),10)) + " "
                  + str(np.around(self.placement.Theta()*180/np.pi,10)) + " "
                  + str(np.around(self.placement.Phi()*180/np.pi,10)))
        else:
            print("slot.{0:02d}.X.0.vec: ".format(idxclover)
                  + str(np.around(self.placement.X(),10)) + " "
                  + str(np.around(self.placement.Y(),10)) + " "
                  + str(np.around(self.placement.Z(),10)))
        print ("")
    def print_crystal_placement(self,idxclover,hemisphere='Southern',coord=None):
        theta_rotation = TRotation()
        phi_rotation = TRotation()
        theta_rotation.RotateY(self.placement.Theta())
        phi_rotation.RotateZ(self.placement.Phi())
        if coord == "Spherical":
             print("slot.{0:02d}.X.0.vec_sph: ".format(idxclover)
                  + str(np.around(self.placement.Mag(),10)) + " "
                  + str(np.around(self.placement.Theta()*180/np.pi,10)) + " "
                  + str(np.around(self.placement.Phi()*180/np.pi,10)))
        else:
            print("slot.{0:02d}.X.0.vec: ".format(idxclover)
                  + str(np.around(self.placement.X(),10)) + " "
                  + str(np.around(self.placement.Y(),10)) + " "
                  + str(np.around(self.placement.Z(),10)))
        if hemisphere == 'Northern':
            crystal_offsets = [(-1,+1),(+1,+1),(+1,-1),(-1,-1)]
        else:
            crystal_offsets = [(+1,+1),(-1,+1),(-1,-1),(+1,-1)]
        for k,(xsign,ysign) in enumerate( crystal_offsets ):
            crystal_offset = TVector3(xsign*self.xdispl,ysign*self.ydispl,self.covergap)
            crystal_offset = phi_rotation*theta_rotation*crystal_offset
            crystal = TVector3(self.placement.X()+crystal_offset.X(),
                               self.placement.Y()+crystal_offset.Y(),
                               self.placement.Z()+crystal_offset.Z())
            if coord == "Spherical":
                print("slot.{0:02d}.{1}.{2}.vec_sph: ".format(idxclover,chr(0x40+k+1),0)
                      + str(np.around(crystal.Mag(),10)) + " "
                      + str(np.around(crystal.Theta()*180/np.pi,10)) + " "
                      + str(np.around(crystal.Phi()*180/np.pi,10)))
            else:
                print("slot.{0:02d}.{1}.{2}.vec: ".format(idxclover,chr(0x40+k+1),0)
                      + str(np.around(crystal.X(),10)) + " "
                      + str(np.around(crystal.Y(),10)) + " "
                      + str(np.around(crystal.Z(),10)))
        print("")

    def print_segment_placement(self,idxclover,hemisphere='Southern',coord=None):
        theta_rotation = TRotation()
        phi_rotation = TRotation()
        theta_rotation.RotateY(self.placement.Theta())
        phi_rotation.RotateZ(self.placement.Phi())
        if coord == "Spherical":
             print("slot.{0:02d}.X.0.vec_sph: ".format(idxclover)
                  + str(np.around(self.placement.Mag(),10)) + " "
                  + str(np.around(self.placement.Theta()*180/np.pi,10)) + " "
                  + str(np.around(self.placement.Phi()*180/np.pi,10)))
        else:
            print("slot.{0:02d}.X.0.vec: ".format(idxclover)
                  + str(np.around(self.placement.X(),10)) + " "
                  + str(np.around(self.placement.Y(),10)) + " "
                  + str(np.around(self.placement.Z(),10)))
        if hemisphere == 'Northern':
            crystal_offsets = [(-1,+1),(+1,+1),(+1,-1),(-1,-1)]
        else:
            crystal_offsets = [(+1,+1),(-1,+1),(-1,-1),(+1,-1)]
        for k,(xsign,ysign) in enumerate( crystal_offsets ):
            crystal_offset = TVector3(xsign*self.xdispl,ysign*self.ydispl,self.covergap)
            crystal_offset = phi_rotation*theta_rotation*crystal_offset
            crystal = TVector3(self.placement.X()+crystal_offset.X(),
                               self.placement.Y()+crystal_offset.Y(),
                               self.placement.Z()+crystal_offset.Z())
            if coord == "Spherical":
                print("slot.{0:02d}.{1}.{2}.vec_sph: ".format(idxclover,chr(0x40+k+1),0)
                      + str(np.around(crystal.Mag(),10)) + " "
                      + str(np.around(crystal.Theta()*180/np.pi,10)) + " "
                      + str(np.around(crystal.Phi()*180/np.pi,10)))
            else:
                print("slot.{0:02d}.{1}.{2}.vec: ".format(idxclover,chr(0x40+k+1),0)
                      + str(np.around(crystal.X(),10)) + " "
                      + str(np.around(crystal.Y(),10)) + " "
                      + str(np.around(crystal.Z(),10)))
            if self.type == CloverType.Yale:
                for s,segsign in enumerate([+1, -1]):
                    segment_offset = TVector3(xsign*self.xdispl+segsign*self.segdispl,ysign*self.ydispl,self.covergap)
                    segment_offset = phi_rotation*theta_rotation*segment_offset
                    segment = TVector3(self.placement.X()+segment_offset.X(),
                                       self.placement.Y()+segment_offset.Y(),
                                       self.placement.Z()+segment_offset.Z())
                    if coord == "Spherical":
                        print("slot.{0:02d}.{1}.{2}.vec_sph: ".format(idxclover,chr(0x40+k+1),s+1)
                          + str(np.around(segment.Mag(),10)) + " "
                          + str(np.around(segment.Theta()*180/np.pi,10)) + " "
                          + str(np.around(segment.Phi()*180/np.pi,10)))
                    else:
                        print("slot.{0:02d}.{1}.{2}.vec: ".format(idxclover,chr(0x40+k+1),s+1)
                              + str(np.around(segment.X(),10)) + " "
                              + str(np.around(segment.Y(),10)) + " "
                              + str(np.around(segment.Z(),10)))
            elif self.type == CloverType.IMP:
                if xsign == +1 and ysign == +1: # A
                    segment_offsets = [(+1,+1),(-1,+1),(-1,-1),(+1,-1)]
                if xsign == -1 and ysign == +1: # B
                    segment_offsets = [(-1,+1),(-1,-1),(+1,-1),(+1,+1)]
                if xsign == -1 and ysign == -1: # C
                    segment_offsets = [(-1,-1),(+1,-1),(+1,+1),(-1,+1)]
                if xsign == +1 and ysign == -1: # D
                    segment_offsets = [(+1,-1),(+1,+1),(-1,+1),(-1,-1)]

                for s,(xsegsign,ysegsign) in enumerate( segment_offsets ):
                    segment_offset = TVector3(xsign*self.xdispl+xsegsign*self.xsegdispl,
                                              ysign*self.ydispl+ysegsign*self.ysegdispl,
                                              self.covergap)
                    segment_offset = phi_rotation*theta_rotation*segment_offset
                    segment = TVector3(self.placement.X()+segment_offset.X(),
                                       self.placement.Y()+segment_offset.Y(),
                                       self.placement.Z()+segment_offset.Z())
                    if coord == "Spherical":
                        print("slot.{0:02d}.{1}.{2}.vec_sph: ".format(idxclover,chr(0x40+k+1),s+1)
                          + str(np.around(segment.Mag(),10)) + " "
                          + str(np.around(segment.Theta()*180/np.pi,10)) + " "
                          + str(np.around(segment.Phi()*180/np.pi,10)))
                    else:
                        print("slot.{0:02d}.{1}.{2}.vec: ".format(idxclover,chr(0x40+k+1),s+1)
                              + str(np.around(segment.X(),10)) + " "
                              + str(np.around(segment.Y(),10)) + " "
                              + str(np.around(segment.Z(),10)))




        print("")


    def set_placement(self, _placement):
        self.placement = _placement

if __name__=='__main__':
    coord=''
    if len(sys.argv) > 1:
        coord = sys.argv[1]

    # 45 degree LaBr3 detectors
    for i in range(0,4):
        r = 20.8
        theta = 45.0
        phi = 45-90.0*i
        if phi > 180:
            phi -= 360
        if phi < -180:
            phi += 360
        placed = clover(r,theta,phi,CloverType.Yale)
        placed.print_clover_placement(i+1,coord=coord)

    # 90 degree detectors (yale)
    for i in range(0,8):
        r = 20.8
        theta = 90.0
        phi = 3*45./2 - 45*i
        if phi > 180:
            phi -= 360
        if phi < -180:
            phi += 360
        placed = clover(r,theta,phi,CloverType.Yale)
        index = (i+4)+1
        if index < 9:
            placed.print_segment_placement(index,coord=coord)
        else:
            placed.print_segment_placement(index,hemisphere='Northern',coord=coord)

    # 135 degree detectors (IMP+ANL)
    for i in range(0,4):
        r = 20.8
        theta = 135.0
        phi = 45-90.0*i
        if phi > 180:
            phi -= 360
        if phi < -180:
            phi += 360
        place_imp = clover(r,theta,phi,CloverType.IMP)
        place_yale = clover(r,theta,phi,CloverType.Yale)
        index = (i+12)+1
        if index > 14:
            place_imp.print_segment_placement(index,coord=coord)
        else:
            #place_yale.print_segment_placement(index,coord=coord)
            place_yale.print_crystal_placement(index,coord=coord)
