import readsubcatalog as r
import numpy as np
import os, itertools, sys
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

snapidDic = {'026':6.0, '027':5.5, '029':5.0, '031':4.5,'034':4.0,'037':3.5, '039':3.0,'050':2.5,'058':2.0}

snapidList = ['026', '027', '029', '031','034','037', '039','050','058']
mstarmin = 9.5 #asked for in the email
runs = len(snapidList)
dir = "./imgs/" #directory to put the plots
pos = [None]*runs ; vel = [None]*runs; logstrmass = [None]*runs; boxlen = [None]*runs; h = [None]*runs
Rdiff = [None]*runs; fpairs = [None]*runs; fmergedata = [None]*runs;
plotty = raw_input('\033[33mWould you like to make plots? (Y/N):\033[0m ')
#plotty = 'n' #for testing

for i, snapid in enumerate(snapidList):
    print("Running interation {0} with snapid {1}".format(i, snapid))
    
    [pos[i], vel[i], logstrmass[i], boxlen[i], h[i]] = r.loaddata_sg(snapid, mstarmin, snapidDic[snapid])
    #gets the data points
    if plotty.upper() == 'Y':
        #Plots the position of the points 
        if not os.path.exists(dir):
            os.makedirs(dir)
        fig = plt.figure(figsize=(18,9.3),dpi=150)
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        ax = fig.add_subplot(111, projection='3d')
        x=[];y=[];z=[];
        for j in pos[i]:
            x.append(j[0]); y.append(j[1]); z.append(j[2])
        plt.title("z = {0} (contains {1} subhalos, mass $\geq$ {2})".format(snapidDic[snapid],len(x),mstarmin), y=1.08)
        ax.scatter(x, y, z)
        plt.xlabel(r'h^{-1}kpc')
        plt.ylabel(r'h^{-1}kpc')
        ax.set_zlabel(r"h^{-1}kpc")
        plt.savefig(dir + "z{0}.pdf".format(snapidDic[snapid]))
        #plt.show()

        xvel = []; yvel = []; zvel = []
        for q in vel[i]:
            xvel.append(q[0]); yvel.append(q[1]); zvel.append(q[2])

        fig = plt.figure(figsize=(18,9.3),dpi=150)
        ax = fig.gca(projection='3d')
        plt.xlabel(r'h^{-1}kpc')
        plt.ylabel(r'h^{-1}kpc')
        plt.title("Velocity Vectors z = {0} (contains {1} subhalos, mass $\geq$ {2})".format(snapidDic[snapid],len(x),mstarmin), y=1.08)
        ax.quiver(x, y, z, xvel, yvel, zvel, length=1500)
        plt.savefig(dir + "z{0}_vec.pdf".format(snapidDic[snapid]))
        #plt.show()



        #fig = plt.figure(2)
        #ax = fig.gca(projection='3d')
        #ax.quiver(x, y, z, xvel, yvel, zvel)
        #plt.show()
    
        #for j, pair in enumerate(itertools.combinations(pos[i],2)):
        #    Rdiff[i].append( getdis_periodic(pair[0],pair[1]), h[i])
    
    #get the merger data
    [fpairs[i], fmergedata[i]] =  r.getmergerdata(snapid, 9.5 , 5.0, 25.0, 0, 500, 0.25, 1.0, snapidDic[snapid])

    #for j in fpairs[i][0]:
    #    if j in fpairs[i][1]:
    #        print "Double counted"
    #        sys.exit()

#Plots total mergers at a given z    
z=[None]*len(snapidList)
totalmergers = [None]*len(snapidList); fracmergers = [None]*len(snapidList)
for i,j in enumerate(snapidList):
    z[i] = snapidDic[j]
    totalmergers[i] = len(fpairs[i])
    fracmergers[i] = len(fpairs[i])/len(pos[i])
    
f, ax = plt.subplots(2, sharex=True) #plt.figure(dpi=100)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
ax[0].plot(z, totalmergers, 'ro')
ax[0].set_title('Mergers at Redshift z')
ax[0].set_xlim([snapidDic[snapidList[-1]]-0.25, snapidDic[snapidList[0]]+0.25])
ax[0].set_ylim([min(totalmergers)-5, max(totalmergers)+5])
ax[0].grid(True)
#ax[0].set_yscale('log')
ax[0].set_xlabel('z')
ax[0].set_ylabel('Number of Mergers')

ax[1].plot(z, fracmergers, 'g^')
ax[1].set_xlabel('z')
ax[1].grid(True)
ax[1].set_ylim([min(fracmergers)-0.01, max(fracmergers)+0.01])
ax[1].set_ylabel('Fraction of subhalos that merge')
plt.savefig(dir + "MergersAtZ.pdf")
plt.show()
