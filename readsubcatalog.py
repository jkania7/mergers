"""
This script finds the pairs that meet the criteria, returns the pairs and information about them
created by Annth and modified by jwk
"""

import numpy
from readsubhalo import *

"""
Snapshots : 
z = 6.0 -> SnapID : '026'
z = 5.5 -> SnapID : '027'
z = 5.0 -> SnapID : '029'
z = 4.5 -> SnapID : '031'
z = 4.0 -> SnapID : '034'
z = 3.5 -> SnapID : '037' 
z = 3.0 -> SnapID : '039' 
z = 2.5 -> SnapID : '050'
z = 2.0 -> SnapID : '058'
"""



def loaddata_sg(snapid,mstarmin,z):
    dirname2 =  './'
    class1 = SnapDir(snapid,dirname2)
    sg = class1.readsubhalo()
    chk4 = numpy.array(class1.load('subhalo','type'))
    chkno = (chk4 == 0)
    sg = sg[chkno]
    h = 0.702
    sg_pos = numpy.array(sg['pos']) # in h^{-1}kpc
    sg_vel = numpy.array(sg['vel']) # in km/sec
    sg_logstrmass = numpy.log10(numpy.array(sg['massbytype'][:,4])*(10.0**10)/h)
    """sg_logstrmass = numpy.log10(numpy.array(sg[\'massbytype\'][:,4])*(10.0**10)/h will through RuntimeWarning: divide by zero encountered in log10,
    this causes the log to be negative infinity, this point gets thrown out at chk1 below, so we don't have to worry."""
    boxlen = 100000.0
    chk1 = (sg_logstrmass > mstarmin)
    pos = sg_pos[chk1]
    vel = [x/numpy.sqrt(1+z) for x in sg_vel[chk1]] #sg_vel[chk1]
    logstrmass = sg_logstrmass[chk1]
    return pos, vel, logstrmass, boxlen, h

def getdis_periodic(pos1,pos2,boxlen):

    permax = boxlen/2.0
    Rdiff = pos1 - pos2
    Rmask = 1 - (numpy.abs(Rdiff/permax)).astype(numpy.int32)
    Roffset = numpy.ma.masked_array((-1)*boxlen*((Rdiff/permax).astype(numpy.int32)),Rmask)
    Rdiff = Rdiff + Roffset.filled(0)     

    return Rdiff

def getmergerpairs(pos,vel,logstrmass,boxlen,rsepmin,rsepmax, vsepmin, vsepmax, mrmin, mrmax):
    fpairs = []
    fmergedata = []
    ntot = len(pos)
    for i in range(0,ntot):
             rdiff1 = getdis_periodic(pos,pos[i],boxlen)
             rprojdis = pow((rdiff1[:,0])**2 + (rdiff1[:,1])**2,0.5)
             vdiff1 = vel - vel[i]
             deltav = pow((vdiff1[:,0])**2 + (vdiff1[:,1])**2,0.5)
             strmratio = 10**(logstrmass - logstrmass[i])
             cmbmass = 10**logstrmass + 10**logstrmass[i]#added jwk
             chksep = ((rprojdis > rsepmin) & (rprojdis < rsepmax) & (deltav > vsepmin) & (deltav < vsepmax) & (strmratio > mrmin) & (strmratio < mrmax))
             cpairs = numpy.arange(0,ntot)[chksep]
             pairs = numpy.zeros((len(cpairs),2),dtype = numpy.int32)
             pairs[:,0] = i
             pairs[:,1] = cpairs  
             mergedata = numpy.zeros((len(cpairs),4),dtype = numpy.float32)#Changed from 3 to 4 by jwk
             mergedata[:,0] = rprojdis[chksep]
             mergedata[:,1] = deltav[chksep]
             mergedata[:,2] = strmratio[chksep]
             mergedata[:,3] = cmbmass[chksep]#added jwk for merger rate determination
             #print i, ntot #used to make sure it loops over all pos
             fpairs.append(pairs)
             fmergedata.append(mergedata)
    fpairs = numpy.concatenate(fpairs)
    fmergedata = numpy.concatenate(fmergedata)
    return fpairs, fmergedata

def getmergerdata(snapid,lstrmmin, rpmin, rpmax, deltavmin, deltavmax, strratiomin, strratiomax, z):
 
    pos, vel, logstrmass, boxlen, h = loaddata_sg(snapid,lstrmmin, z)
    fpairs, fmergedata = getmergerpairs(pos,vel,logstrmass,boxlen,rpmin,rpmax, deltavmin, deltavmax, strratiomin, strratiomax)
    return fpairs, fmergedata

"""
snapid = '058'
lstrmin = 9.5
rpmin = 5.0
rpmax = 25.0
deltavmin = 0
deltavmax = 500.0
strratiomin = 0.25
strrationmax = 1.0
fpairs, fmergedata = getmergerdata(snapid,lstrmmin, rpmin, rpmax, deltavmin, deltavmax, strratiomin, strratiomax)
"""
