import readsubcatalog as r
import numpy as np
import os, itertools
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


snapidList = ['026', '027', '029', '031','034','037', '039','050','058']
mstarmin = 9.5 #asked for in the email
runs = len(snapidList)
dir = "./imgs/" #directory to put the plots
pos = [None]*runs ; vel = [None]*runs; logstrmass = [None]*runs; boxlen = [None]*runs; h = [None]*runs
Rdiff = [None]*runs

for i, snapid in enumerate(snapidList):
    print("Running interation {0} with snapid {1}".format(i, snapid))

    [pos[i], vel[i], logstrmass[i], boxlen[i], h[i]] = r.loaddata_sg(snapid, mstarmin)
    #gets the data points

    #Plots the position of the points 
    if not os.path.exists(dir):
        os.makedirs(dir)
    fig = plt.figure(figsize=(12,8),dpi=80)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    ax = fig.add_subplot(111, projection='3d')
    x=[];y=[];z=[];
    for j in pos[i]:
        x.append(j[0]); y.append(j[1]); z.append(j[2])
        plt.title("Plotting snapid {0} (contains {1} subhalos, mass $\geq$ {2})".format(snapid,len(x),mstarmin), y=1.08)
    ax.scatter(x, y, z)
    plt.xlabel(r'h^{-1}kpc')
    plt.ylabel(r'h^{-1}kpc')
    ax.set_zlabel(r"h^{-1}kpc")
    plt.savefig(dir + "snapid{}.pdf".format(snapid))
    #plt.show()

    #for j, pair in enumerate(itertools.combinations(pos[i],2)):
    #    Rdiff[i].append( getdis_periodic(pair[0],pair[1]), h[i])


    #print r.getmergerdata(snapid, 9.5 , 5.0, 25.0, 0, 500, 0.25, 1.0)
