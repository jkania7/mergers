"""
Main script, calls the scripts and creates the plots. Created by J. Kania as part of part of undergraduate research at CMU
"""
import readsubcatalog as r
import numpy as np
import os, itertools, sys, math, timeit
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.integrate import quad
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

def friedmann(a, OmegaM, OmegaL, H0):#, OmegaR) #simplified friedmann equations
    Omega0 = OmegaM + OmegaL #+OmegaR
    #return 1/H0*(1/(OmegaR/a**2 + OmegaM/a+ OmegaL*a**2+ (1-Omega0))^(-1/2)                 
    return 1.0/H0*(OmegaM/a+ OmegaL*a**2)**(-1.0/2.0)

def mergerfreq(z, Ms, h0):
    #z = (1.0/a)-1.0 #Not needed, everything is done already in z
    Ms = Ms/(3.0*10.0**10.0/h0)#stellar mass of the pair in units of 3*10^10*H^-1*Msun as defined in K&W p1489
    #for vp < 300 km/s and rp < 30kpc/h from K&W p1497
    T0 = 2038.0*h0
    f1 = -165.0*h0/(10.0**5.0)
    f2 = 690.0*h0/(10.0**5.0)
    return (T0**(-0.5)+f1*z+f2*(np.log10(Ms)-10.0))**2.0 #K&W (9) on p1496 [Units of Myr]

def mergerfit(x, a, b, c):
    return a + b*(1.0+x)**c
                 
if __name__ == "__main__":
    start = timeit.default_timer() #used to find the runtime of the program
    
    snapidDic = {'026':6.0, '027':5.5, '029':5.0, '031':4.5,'034':4.0,'037':3.5, '039':3.0,'050':2.5,'058':2.0,'063':1.5,'068':1.0,'073':0.6,'085':0.06}
    snapidList = ['026', '027', '029', '031','034','037', '039','050','058','063','068','073','085']#=list(snapidDic) #this does not keep order in z
    """
    snapidDic = {'026':6.0, '029':5.0,'034':4.0,'039':3.0,'058':2.0,'068':1.0} #a sparse list for testing
    snapidList = ['026', '029','034','039','058','068']#=list(snapidDic) #this does not keep order in z
    """
    OmegaL = 0.725
    OmegaM = 0.275
    h0 = 0.701
    H0 = h0*100.0/(3.086e19)*3.154e7 #the terms following the 100.0 are MPc to km and year to seconds respectivley 
    mstarmin = 9.5 #asked for in the email
    dir = "./imgs/" #directory to put the plots
    #End user defined variables

    if not os.path.exists(dir):
            os.makedirs(dir)
    
    runs = len(snapidDic)
    pos = [None]*runs ; vel = [None]*runs; logstrmass = [None]*runs; boxlen = [None]*runs; h = [None]*runs
    Rdiff = [None]*runs; fpairs = [None]*runs; fmergedata = [None]*runs;
    plotty = raw_input('\033[33mWould you like to make plots? (Y/N):\033[0m ')
    #plotty = 'n' #for testing

    for i, snapid in enumerate(snapidList):
        print("\033[1;34mRunning interation {0} with snapid {1} (z = {2})\033[0m".format(i, snapid, snapidDic[snapid]))
    
        [pos[i], vel[i], logstrmass[i], boxlen[i], h[i]] = r.loaddata_sg(snapid, mstarmin, snapidDic[snapid])
        #print("The len(pos) = {0}, len(pos[i] = {1}, len(pos[i][1] =  {2}".format(len(pos),len(pos[i]), len(pos[i][1])))
        #gets the data points
                
        if plotty.upper() == 'Y':
            #Plots the position of the points 
            fig = plt.figure()#(figsize=(18,9.3),dpi=150)
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

            fig = plt.figure()#(figsize=(18,9.3),dpi=150)
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
    totalmergers = [None]*len(snapidList); fracmergers = [None]*len(snapidList); numHalos = [None]*len(snapidList)
    mergerRate = [None]*len(snapidList); #deltaT = [None]*len(snapidList); error = [None]*len(snapidList);
    fracMergerRate = [None]*len(snapidList)
    print("fpairs[i] = {0}".format(fpairs[i]))
    print("fmergedata[i] = {0}".format(fmergedata[i]))
    for i,j in enumerate(snapidList):
        z[i] = float(snapidDic[j])
        #z0 = z[i] + 0.5 #used for first guess of merger rate
        #deltaT[i], error[i] = quad(friedmann,1.0/(1.0+z0), 1.0/(1.0+z[i]),args=(OmegaM, OmegaL, H0)) #to calculate times between z intervals, used as inital guess for merger rate
        totalmergers[i] = float(len(fpairs[i]))
        fracmergers[i] = float(len(fpairs[i]))/float(len(pos[i]))
        numHalos[i] = float(len(pos[i]))
        #mergerRate[i] = float(len(fpairs[i]))/float(deltaT[i])*10e9 #first guess at merger rates
        #fracMergerRate[i] =  fracmergers[i]/float(deltaT[i])*10e9
        totalMergT = 0.0#total merger time
        for l in range(len(fpairs[i])):
            print("totalMergT = {0}".format(totalMergT))
            totalMergT = totalMergT + mergerfreq(z[i], fmergedata[i][l][3], h0)
        #mergerRate[i] = float(len(fpairs[i]))*float(totalMergT)*10.0**3.0#10^3 to go from /Myr to /Gyr
        mergerRate[i] = float(totalMergT)*10.0**3.0#10^3
        fracMergerRate[i] =  fracmergers[i]*float(totalMergT)*10.0**3.0
    
    stop = timeit.default_timer()#stops the timer before plotting
    #print mergerRate
    #print fracMergerRate
    
    f, ax = plt.subplots(3, sharex=True) #plt.figure(dpi=100)
    ax[0].plot(z, totalmergers, 'ro')
    ax[0].set_title(r'Mergers at $z$')
    ax[0].set_xlim([snapidDic[snapidList[-1]]-0.15, snapidDic[snapidList[0]]+0.15])
    ax[0].set_ylim([min(totalmergers)-3e1, max(totalmergers)+3e1])
    ax[0].grid(True)
    #ax[0].set_yscale('log')
    ax[0].set_ylabel(r'Number of Mergers')

    ax[1].plot(z, fracmergers, 'g^')
    ax[1].grid(True)
    ax[1].set_ylim([min(fracmergers)*0.6, max(fracmergers)*1.1])
    ax[1].set_ylabel(r'Merger Fraction')
    
    ax[2].plot(z, numHalos, 'bo')
    ax[2].set_xlabel(r'$z$')
    ax[2].grid(True)
    ax[2].set_ylim([min(numHalos)-2e3, max(numHalos)+2e3])
    ax[2].set_ylabel(r'Number of subhalos')
    
    plt.savefig(dir + "Merger.pdf")
    plt.show()


    f2, ax2 = plt.subplots(3, sharex=True) #plt.figure(dpi=100)
    ax2[0].plot(z, mergerRate, 'ro')
    ax2[0].set_title(r'Merger Rate at $z$')
    ax2[0].set_xlim([snapidDic[snapidList[-1]]-0.15, snapidDic[snapidList[0]]+0.15])
    ax2[0].set_ylim([min(mergerRate)-1e1, max(mergerRate)+1e1])
    ax2[0].grid(True)
    #ax2[0].set_yscale('log')
    ax2[0].set_ylabel(r'Mergers/Gyr')

    ax2[1].plot(z, fracMergerRate, 'g^')
    ax2[1].grid(True)
    ax2[1].set_ylim([min(fracMergerRate)-1e-1, max(fracMergerRate)+1e-1])
    ax2[1].set_ylabel(r'Mergere/(Gyr*subhalo)')
    
    ax2[2].plot(z, numHalos, 'bo')
    ax2[2].set_xlabel(r'$z$')
    ax2[2].grid(True)
    ax2[2].set_ylim([min(numHalos)-2e3, max(numHalos)+2e3])
    ax2[2].set_ylabel(r'Number of subhalos')
    
    plt.savefig(dir + "MergerRate.pdf")
    #mng = plt.get_current_fig_manager()
    #mng.frame.Maximize(True)
    plt.show()

    
    zdense = np.linspace(0, snapidDic[snapidList[0]]+0.25, num=150) 
    #z14 = (1.0+zdense)**4.0 #first guess at fitted curve 

    pars1, pcoc1 = curve_fit(mergerfit, z, fracMergerRate, p0=[-1,1,4], bounds=((-4, 0, 3),(0, 4.5, 4.5)))
    print("For a+b*(1+z)^c; [a,b,c] = {0}".format(pars1))
    
    fig10 = plt.plot(z, fracMergerRate, 'g^', label='MBii')
    fig11 = plt.plot(zdense, mergerfit(zdense, *pars1), 'r--', label='(1+z)^4')
    #plt.legend([fig10, fig11], ["MBii", "$(1+z)^4$"])
    #plt.legend(loc='upper left')
    plt.xlabel(r'$z$')
    plt.title(r'Log(Merger Rate)')
    plt.ylabel(r'LOG(Merger Fraction/Gyr)')
    plt.yscale('log')
    plt.xlim([-0.3 , snapidDic[snapidList[0]]+0.3])
    plt.grid(True)
    plt.savefig(dir + "MergerRateLog.pdf")
    plt.show()

    pars2, pcoc2 = curve_fit(mergerfit, z, fracmergers, p0=[-1,1,4], bounds=([-4, 0, 3], [0, 4.5, 4.5]) )
    print("For a+b*(1+z)^c; [a,b,c] = {0}".format(pars2))

    fig12 = plt.plot(z, fracmergers, 'bs')
    plt.plot(zdense, mergerfit(zdense, *pars2), 'r--', label='(1+z)^4')
    plt.xlabel(r'$z$')
    plt.title(r'Log(Merger Fraction)')
    plt.ylabel(r'Merger Fraction')
    plt.yscale('log')
    plt.xlim([-0.3 , snapidDic[snapidList[0]]+0.3])
    plt.grid(True)
    plt.savefig(dir + "MergerFracLog.pdf")
    plt.show()
    
    #print("deltaT = {0}".format(deltaT))

    
    print("Runtime = {0:.1f} min".format((stop - start)/60.0))
