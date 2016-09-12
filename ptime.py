#!/usr/bin/python
#### ptime.py

####  The bwwserver_nwchem.py and ptime.py programs were written by
####
####     Eric J. Bylaska 
####     Environmental Molecular Sciences Laboratory
####     Pacific Northwest National Laboratory
####     Richland, Washington 99354
####
####  Please consider citing the the following reference when using this program
####
####     Eric J. Bylaska, Jonathan Q. Weare, and John H. Weare (2012) "Parallel in 
####     time simulations using high level quantum chemistry methods and complex 
####     potentials" JCTC, submitted
####

import sys, os, time, socket, pickle, math

#### atomic symbols ####
def_symbols = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K',
               'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
               'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
               'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
               'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra',
               'Ac', 'Th', 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db',
               'Sg', 'Bh', 'Hs', 'Mt']
#### atomic masses ####
def_masses = [1.007825, 4.00260, 7.0160, 9.01218, 11.00931, 12.0, 14.00307, 15.99491, 18.9984, 19.99244, 22.9898,
              23.98504, 26.98154, 27.97693, 30.97376, 31.97207, 34.96885, 39.9624, 38.96371, 39.96259, 44.95592, 45.948,
              50.9440, 51.9405, 54.9381, 55.9349, 58.9332, 57.9353, 62.9298, 63.9291, 68.9257, 73.9219, 74.9216,
              79.9165, 78.9183, 83.912, 84.9117, 87.9056, 88.9054, 89.9043, 92.9060, 97.9055, 97.9072, 101.9037,
              102.9048, 105.9032, 106.90509, 113.9036, 114.9041, 117.9018, 120.9038, 129.9067, 126.9004, 131.9042,
              132.9051, 137.9050, 138.9061, 139.9053, 140.9074, 143.9099, 144.9128, 151.9195, 152.9209, 157.9241,
              159.9250, 163.9288, 164.9303, 165.9304, 168.9344, 173.9390, 174.9409, 179.9468, 180.948, 183.9510,
              186.9560, 189.9586, 192.9633, 194.9648, 196.9666, 201.9706, 204.9745, 207.9766, 208.9804, 209.9829,
              210.9875, 222.0175, 223.0198, 226.0254, 227.0278, 232.0382, 231.0359, 238.0508, 237.0482, 244.0642,
              243.0614, 247.0704, 247.0703, 251.0796, 252.0829, 257.0950, 258.0986, 259.1009, 262.1100, 261.1087,
              262.1138, 266.1219, 262.1229, 267.1318, 268.1388]

#############################################
# #
# my_sockconnect                #
#                                           #
#############################################
#
# Returns connected socket to hostport or None if failed
# where hostport = "ipaddress:port" (e.g. hostport="localhost:50001")
#
def my_sockconnect(hostport):
    l = hostport.index(':');
    r = hostport.rindex(':')
    if (r == l):
        hp = (hostport[:l], eval(hostport[l + 1:]))
    else:
        hp = (hostport[:l], eval(hostport[l + 1:r]))
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.settimeout(10.0)
    try:
        s.connect(hp)
        s.settimeout(None)
        return s
    except:
        s.close()
        s = None
    return s


#######################################################
#                                                     #
#                   runcalcs                          #
#                                                     #
#######################################################
#
# This routine calculates the energy-gradient using the dictionary job
# of the  list of xs geometries.
#
def runcalcs(hostsports, job, xs, mstart=None):
    if (mstart == None): mstart = 0
    egrads = []
    jobs = []
    for i in range(mstart, len(xs)):
        jobs.append(job.copy())
        jobs[i - mstart]['xyz'] = xs[i]
    njobs = len(jobs);
    nhp = len(hostsports)

    ## set up connections ###
    ss = [];
    sscpus = [];
    ncpus = 0
    for hp in hostsports[:njobs]:
        s = my_sockconnect(hp)
        if (s != None):
            ss.append(s)
            l = hp.index(':');
            r = hp.rindex(':')
            if (l == r):
                cpus = 1
            else:
                cpus = eval(hp[r + 1:])
            sscpus.append(cpus);
            ncpus += cpus
    if (len(ss) > 0):
        ns = len(ss)
        jobcount = [0] * ns
        for i in range(ns):                  jobcount[i] = int(njobs * sscpus[i] / float(ncpus))
        for i in range(njobs - sum(jobcount)): jobcount[i % ns] += 1
        ## send jobs to connections ##
        i = -1;
        istart = 0
        for s in ss:
            iend = istart + jobcount[i + 1]
            s.sendall(pickle.dumps(jobs[istart:iend]) + 'ericdone')  #send data
            istart = iend
            i += 1
        for s in ss:
            msg = ''
            done = False
            while not done:
                data = s.recv(1024)  #receive data
                msg = msg + data
                done = msg[len(msg) - 8:len(msg)] == 'ericdone'
            msg = msg[:len(msg) - 8]
            data = pickle.loads(msg)
            egrads += data
            s.close()

    return egrads


def util_date():
    tt = time.localtime()
    dd = "%d-%d-%d-%d:%d" % (tt[0], tt[1], tt[2], tt[3], tt[4])
    return dd


#######################################################
#                                                     #
#                   QHess2 class                      #
#                                                     #
#######################################################
#
# This class implements a Broyden quasi-Newton algorithm
#
class QHess2:
    qsize = 0
    qx = []
    qa = []
    fk = []
    qrho = []
    qindx = []

    def __init__(self, n0, qmax0=5, qalpha0=0.0):
        self.n = n0
        self.qmax = qmax0
        self.qalpha = qalpha0
        self.qsize = 0
        self.qx = [0] * self.qmax * self.n
        self.qa = [0] * self.qmax * self.n
        self.fk = [0] * self.n
        self.qrho = [0] * self.qmax
        self.qindx = range(self.qmax)

    def size(self):
        return self.qsize

    def AddHistory(self, x, a):
        rho = 0.0

        while ((self.qsize > 0) and (abs(rho) < 1.0e-9)):
            shift = self.qindx[self.qsize - 1] * self.n
            rho = 0.0
            for i in range(self.n):
                rho += (x[i] - self.qx[shift + i]) * (x[i] - self.qx[shift + i])

            if (abs(rho) < 1.0e-9):
                self.qsize -= 1;

        if (self.qsize >= self.qmax):
            itmp = self.qindx[0]
            for i in range(1, self.qmax):
                self.qindx[i - 1] = self.qindx[i]
            self.qindx[self.qmax - 1] = itmp
            self.qsize -= 1

        shift = self.qindx[self.qsize] * self.n
        for i in range(self.n):
            self.qx[shift + i] = x[i]
            self.qa[shift + i] = a[i]
        self.qrho[self.qindx[self.qsize]] = rho
        self.qsize += 1

    def Broyden(self, f):
        Bf = [0] * self.n
        for i in range(self.n): self.fk[i] = f[i]

        for m in range(self.qsize - 1, 0, -1):
            shift1 = self.qindx[m] * self.n
            shift0 = self.qindx[m - 1] * self.n
            gamma = 0.0

            for i in range(self.n):
                gamma += (self.qx[shift1 + i] - self.qx[shift0 + i]) * self.fk[i]

            gamma /= self.qrho[self.qindx[m]]
            for i in range(self.n):
                Bf[i] += (self.qa[shift1 + i] - self.qa[shift0 + i]) * gamma
                self.fk[i] -= (self.qx[shift1 + i] - self.qx[shift0 + i]) * gamma

        for i in range(self.n):
            Bf[i] += self.qalpha * self.fk[i]
        return Bf


##################################################
#                                                #
#               mdverlet_serial                  #
#                                                #
##################################################
#
#  Standard velocity verlet integration
#
#  Entry - hostports: list of hostport strings
#          job: python dictionary for the energy-force calculation
#          M: length time segment
#          timestep: velocity verlet time step 
#          Perror: not used
#          maxresiderr: not used
#          mass: array of atomic masses in au
#          x0,v0,a0: arrays of initial positions, velocities,and accelerations in au
#          u0: initial energy, not used
#  Exit - returns (x1,v1,a1,e1,ke1,u1,iterations) tuple where
#          x1,v1,a1: arrays of final positions, velocities,and accelerations in au
#          e1,ke1,u1: final total energy, kinetic energy, and potential energy
#          iterations: number of quasi-Newton steps taken
#
def mdverlet_serial(hostsports, job, M, timestep, Perror, maxresiderr, mass, x0, v0, a0, u0):
    n = len(x0);
    nion = n / 3
    eout = 0.0;
    keout = 0.0;
    uout = 0.0
    x1 = [0] * n;
    v1 = [0] * n;
    a1 = [0] * n
    iterations = M

    for it in range(M):
        for i in range(n):
            x1[i] = x0[i] + v0[i] * timestep + 0.5 * a0[i] * timestep * timestep
        result = runcalcs(hostsports, job, [x1])
        u1 = result[0][0]
        for ii in range(nion):
            a1[3 * ii] = result[0][1][3 * ii] / mass[ii]
            a1[3 * ii + 1] = result[0][1][3 * ii + 1] / mass[ii]
            a1[3 * ii + 2] = result[0][1][3 * ii + 2] / mass[ii]
        for i in range(n):
            v1[i] = v0[i] + 0.5 * (a0[i] + a1[i]) * timestep
        for i in range(n): x0[i] = x1[i]
        for i in range(n): v0[i] = v1[i]
        for i in range(n): a0[i] = a1[i]

    uout = u1
    keout = 0.0
    for ii in range(nion):
        keout += 0.5 * mass[ii] * (
            v1[3 * ii] * v1[3 * ii] + v1[3 * ii + 1] * v1[3 * ii + 1] + v1[3 * ii + 2] * v1[3 * ii + 2])
    eout = uout + keout

    return (x1, v1, a1, eout, keout, uout, iterations)


##################################################
#                                                #
#               mdverlet_parallel                #
#                                                #
##################################################
#
#  Quasi-Newton parallel in time velocity verlet integration.
#
#  Entry - hostports: list of hostport strings
#          job: python dictionary for the energy-force calculation
#          M: length of parallel in time segment
#          timestep: velocity verlet time step 
#          Perror: not used
#          maxresiderr: maximum path residual error
#          mass: array of atomic masses in au
#          x0,v0,a0: arrays of initial positions, velocities,and accelerations in au
#          u0: initial energy, not used
#  Exit - returns (x1,v1,a1,e1,ke1,u1,iterations) tuple where
#          x1,v1,a1: arrays of final positions, velocities,and accelerations in au
#          e1,ke1,u1: final total energy, kinetic energy, and potential energy
#          iterations: number of quasi-Newton steps taken
#
def mdverlet_parallel(hostsports, job, M, timestep, Perror, maxresiderr, mass, x0, v0, a0, u0):
    n = len(x0);
    nion = n / 3
    eout = 0.0;
    keout = 0.0;
    uout = 0.0
    x1 = [0] * n;
    v1 = [0] * n;
    a1 = [0] * n
    it = 0

    ## calculate kinetic energy ##
    ke0 = 0.0
    for ii in range(nion):
        ke0 += 0.5 * mass[ii] * (
            v0[3 * ii] * v0[3 * ii] + v0[3 * ii + 1] * v0[3 * ii + 1] + v0[3 * ii + 2] * v0[3 * ii + 2])
    e0 = u0 + ke0

    errorx = [1.0] * (M + 1);
    errorv = [1.0] * (M + 1);
    erroroldx = [1.0] * (M + 1);
    erroroldv = [1.0] * (M + 1)
    keall = [ke0] * (M + 1);
    uall = [u0] * (M + 1)
    xall = [];
    vall = [];
    aall = []
    fxall = [];
    fvall = []
    Jfxall = [];
    Jfvall = []
    Jfxall2 = [];
    Jfvall2 = []
    Jfxall3 = [];
    Jfvall3 = []
    for m1 in range(M + 1):
        xall.append([0] * n);
        vall.append([0] * n);
        aall.append([0] * n)
        fxall.append([0] * n);
        fvall.append([0] * n)
        Jfxall.append([0] * n);
        Jfvall.append([0] * n)
        Jfxall2.append([0] * n);
        Jfvall2.append([0] * n)
        Jfxall3.append([0] * n);
        Jfvall3.append([0] * n)

    for m1 in range(M + 1):
        for i in range(n):
            xall[m1][i] = x0[i]
            vall[m1][i] = v0[i]
            aall[m1][i] = a0[i]

    ### initial quasihessians ###
    qall = []
    for m1 in range(M + 1):
        qall.append(QHess2(n, 15))

    mstart = 1
    residx0 = 1.0e9;
    residv0 = 1.0e9;
    residx = 0.0;
    residv = 0.0;
    residmax = 0.0

    done = False
    while (not done):
        residx = 0.0;
        residv = 0.0;
        residmax = 0.0
        dmstart = True
        mstart0 = mstart
        for m1 in range(mstart, M + 1):
            qall[m1].AddHistory(xall[m1], aall[m1])

            xres = 0.0;
            vres = 0.0
            for i in range(n):
                dx = xall[m1][i] - (
                    xall[m1 - 1][i] + vall[m1 - 1][i] * timestep + 0.5 * aall[m1 - 1][i] * timestep * timestep)
                dv = vall[m1][i] - (vall[m1 - 1][i] + 0.5 * (aall[m1 - 1][i] + aall[m1][i]) * timestep)
                fxall[m1][i] = dx;
                fvall[m1][i] = dv
                xres += abs(dx);
                vres += abs(dv)
            residx += xres;
            residv += vres

            if ((xres + vres) > residmax):
                residmax = xres + vres

            if (((vres + xres) < (1.0e-4 * maxresiderr / (M + 1))) and dmstart):
                mstart0 = m1
            else:
                dmstart = False

        done = ( (it > M + 2) or ((residx + residv) < maxresiderr)  )
        if (not done):
            normax = 1.0e8
            ### compute B^{-1} * F ###
            for m1 in range(mstart, M + 1):
                Bf2 = qall[m1].Broyden(fxall[m1])
                for i in range(n):
                    Jfxall2[m1][i] = fxall[m1][i];
                    Jfvall2[m1][i] = fvall[m1][i] + 0.5 * timestep * Bf2[i]
                    Jfxall[m1][i] = Jfxall2[m1][i];
                    Jfvall[m1][i] = Jfvall2[m1][i]
            L = 15
            for k in range(L):
                ### compute Q * JF2 ; Q = 1- B^{-1}*DF ###

                ### JF3=C*JF2 ###
                for m1 in range(mstart, M + 1):
                    Bf2 = qall[m1 - 1].Broyden(Jfxall2[m1 - 1])
                    for i in range(n):
                        Jfxall3[m1][i] = Jfxall2[m1 - 1][i] + timestep * Jfvall2[m1 - 1][
                            i] + 0.5 * timestep * timestep * Bf2[i]
                        Jfvall3[m1][i] = Jfvall2[m1 - 1][i] + 0.5 * timestep * Bf2[i]

                ### JF2=B^{-1}*JF3 ###
                for m1 in range(mstart, M + 1):
                    Bf2 = qall[m1].Broyden(Jfxall3[m1])
                    for i in range(n):
                        Jfxall2[m1][i] = Jfxall3[m1][i]
                        Jfvall2[m1][i] = Jfvall3[m1][i] + 0.5 * timestep * Bf2[i]

                for m1 in range(mstart, M + 1):
                    for i in range(n):
                        Jfxall[m1][i] += Jfxall2[m1][i]
                        Jfvall[m1][i] += Jfvall2[m1][i]

            for m1 in range(mstart, M + 1):
                for i in range(n):
                    xall[m1][i] -= Jfxall[m1][i]
                    vall[m1][i] -= Jfvall[m1][i]

            #norm1 = 0.0;
            #norm2 = 0.0;
            #norm3 = 0.0;
            #for m1 in range(mstart,M+1):
            #   for i in range(n):
            #      norm1 += fxall[m1][i]*fxall[m1][i]    + fvall[m1][i]*fvall[m1][i]
            #      norm2 += Jfxall[m1][i]*Jfxall[m1][i]  + Jfvall[m1][i]*Jfvall[m1][i]
            #      norm3 += Jfxall[m1][i]*fxall[m1][i]   + Jfvall[m1][i]*fvall[m1][i]
            #print "norm1,norm2,norm3,residx,residv,mstart=",norm1,norm2,norm3,residx,residv,mstart

            result = runcalcs(hostsports, job, xall, mstart)
            it += 1

            for m1 in range(mstart, M + 1):
                uall[m1] = result[m1 - mstart][0]
                for ii in range(nion):
                    aall[m1][3 * ii] = result[m1 - mstart][1][3 * ii] / mass[ii]
                    aall[m1][3 * ii + 1] = result[m1 - mstart][1][3 * ii + 1] / mass[ii]
                    aall[m1][3 * ii + 2] = result[m1 - mstart][1][3 * ii + 2] / mass[ii]
            for m1 in range(mstart, M + 1):
                keall[m1] = 0.0
                for ii in range(nion):
                    keall[m1] += 0.5 * mass[ii] * (
                        vall[m1][3 * ii] * vall[m1][3 * ii] + vall[m1][3 * ii + 1] * vall[m1][3 * ii + 1] + vall[m1][
                            3 * ii + 2] * vall[m1][3 * ii + 2])

            if (abs(residx0 + residv0 - residx - residv) > maxresiderr):
                mstart = mstart0;
                diff = (M + 1 - mstart + 1) - 2
                if (diff < 0): mstart += diff;
            else:
                mstart = 1;
            residx0 = residx;
            residv0 = residv

            #done  = ( (it>M+2) or ((residx+residv) < maxresiderr)  )

    iterations = it
    for i in range(n):
        x1[i] = xall[M][i];
        v1[i] = vall[M][i];
        a1[i] = aall[M][i]
    ke1 = keall[M];
    u1 = uall[M];
    e1 = u1 + ke1

    return (x1, v1, a1, e1, ke1, u1, iterations)


##################################################
#                                                #
#               mdverlet_parallel2               #
#                                                #
##################################################
#
#  Preconditioned quasi-Newton parallel in time velocity verlet integration.
#
#  Entry - hostports: list of hostport strings
#          job: python dictionary for the fine grain energy-force calculation
#          jobcoarse: python dictionary for the coarse grain energy-force calculation.
#          M: length of parallel in time segment
#          timestep: velocity verlet time step 
#          Perror: not used
#          maxresiderr: maximum path residual error
#          mass: array of atomic masses in au
#          x0,v0,a0: arrays of initial positions, velocities,and fine grain accelerations in au
#          u0: initial energy, not used
#          acoarse0: initial coarse grain accelerations in au
#  Exit - returns (x1,v1,a1,e1,ke1,u1,iterations) tuple where
#          x1,v1,a1: arrays of final positions, velocities,and accelerations in au
#          e1,ke1,u1: final total energy, kinetic energy, and potential energy
#          iterations: number of quasi-Newton steps taken
#
def mdverlet_parallel2(hostsports, job, jobcoarse, M, timestep, Perror, maxresiderr, mass, x0, v0, a0, u0, acoarse0):
    n = len(x0);
    nion = n / 3
    it = 0
    eout = 0.0;
    keout = 0.0;
    uout = 0.0
    x1 = [0] * n;
    v1 = [0] * n;
    a1 = [0] * n;
    acoarse1 = [0] * n

    ### calculate energy force and kinetic energy ###
    errorx = [1.0] * (M + 1);
    errorv = [1.0] * (M + 1);
    erroroldx = [1.0] * (M + 1);
    erroroldv = [1.0] * (M + 1)
    keall = [0] * (M + 1);
    uall = [0] * (M + 1)
    xall = [];
    vall = [];
    aall = []
    xall2 = [];
    vall2 = [];
    aall2 = []
    wall = [];
    zall = []
    fxall = [];
    fvall = []
    Jfxall = [];
    Jfvall = []
    Jfxall1 = [];
    Jfvall1 = []
    Jfxall2 = [];
    Jfvall2 = []
    Jfxall3 = [];
    Jfvall3 = []
    Jfxall4 = [];
    Jfvall4 = []
    for m1 in range(M + 1):
        xall.append([0] * n);
        vall.append([0] * n);
        aall.append([0] * n)
        xall2.append([0] * n);
        vall2.append([0] * n);
        aall2.append([0] * n)
        wall.append([0] * n);
        zall.append([0] * n)
        fxall.append([0] * n);
        fvall.append([0] * n)
        Jfxall.append([0] * n);
        Jfvall.append([0] * n)
        Jfxall1.append([0] * n);
        Jfvall1.append([0] * n)
        Jfxall2.append([0] * n);
        Jfvall2.append([0] * n)
        Jfxall3.append([0] * n);
        Jfvall3.append([0] * n)
        Jfxall4.append([0] * n);
        Jfvall4.append([0] * n)

    for m1 in range(M + 1):
        for i in range(n):
            xall[m1][i] = x0[i]
            vall[m1][i] = v0[i]

    ### initial quasihessians ###
    qall = []
    for m1 in range(M + 1):
        qall.append(QHess2(n, 15))

    qcoarse = []
    for m1 in range(M + 1):
        qcoarse.append(QHess2(n, 3))

    ### compute initial coarse grain path ###
    for i in range(n):
        aall2[0][i] = acoarse0[i]
    for m1 in range(M):
        for i in range(n):
            xall[m1 + 1][i] = xall[m1][i] + vall[m1][i] * timestep + 0.5 * aall2[m1][i] * timestep * timestep
        result = runcalcs(hostsports, jobcoarse, [xall[m1 + 1]])  #result：from server cal
        f0 = result[0][1]
        for ii in range(nion):
            aall2[m1 + 1][3 * ii] = f0[3 * ii] / mass[ii]
            aall2[m1 + 1][3 * ii + 1] = f0[3 * ii + 1] / mass[ii]
            aall2[m1 + 1][3 * ii + 2] = f0[3 * ii + 2] / mass[ii]
        for i in range(n):
            vall[m1 + 1][i] = vall[m1][i] + 0.5 * (aall2[m1][i] + aall2[m1 + 1][i]) * timestep

    ### compute F(x) ###
    for i in range(n):
        aall[0][i] = a0[i]
    uall[0] = u0
    mstart = 1
    result = runcalcs(hostsports, job, xall, mstart)  #result：from server cal
    it += 1
    for m1 in range(mstart, M + 1):
        uall[m1] = result[m1 - mstart][0]
        for ii in range(nion):
            aall[m1][3 * ii] = result[m1 - mstart][1][3 * ii] / mass[ii]
            aall[m1][3 * ii + 1] = result[m1 - mstart][1][3 * ii + 1] / mass[ii]
            aall[m1][3 * ii + 2] = result[m1 - mstart][1][3 * ii + 2] / mass[ii]
    for m1 in range(M + 1):
        keall[m1] = 0.0
        for ii in range(nion):
            keall[m1] += 0.5 * mass[ii] * (
                vall[m1][3 * ii] * vall[m1][3 * ii] + vall[m1][3 * ii + 1] * vall[m1][3 * ii + 1] + vall[m1][
                    3 * ii + 2] *
                vall[m1][3 * ii + 2])

    mstart = 1
    residx0 = 1.0e9;
    residv0 = 1.0e9;
    residx = 0.0;
    residv = 0.0;
    residmax = 0.0

    done = False
    while (not done):
        residx = 0.0;
        residv = 0.0;
        residmax = 0.0
        dmstart = True
        mstart0 = mstart
        for m1 in range(mstart, M + 1):
            qall[m1].AddHistory(xall[m1], aall[m1])
            qcoarse[m1].AddHistory(xall[m1], aall2[m1])

            xres = 0.0;
            vres = 0.0
            for i in range(n):
                dx = xall[m1][i] - (
                    xall[m1 - 1][i] + vall[m1 - 1][i] * timestep + 0.5 * aall[m1 - 1][i] * timestep * timestep)
                dv = vall[m1][i] - (vall[m1 - 1][i] + 0.5 * (aall[m1 - 1][i] + aall[m1][i]) * timestep)
                fxall[m1][i] = dx;
                fvall[m1][i] = dv
                xres += abs(dx);
                vres += abs(dv)
            residx += xres;
            residv += vres

            if ((xres + vres) > residmax):
                residmax = xres + vres

            if (((vres + xres) < (maxresiderr / (M))) and dmstart):
                mstart0 = m1
            else:
                dmstart = False

        #print "residx,residv,mstart0=",residx,residv,mstart0
        done = ( (it > M + 2) or ((residx + residv) < maxresiderr)  )
        if (not done):

            normax = 1.0e8

            ### compute JF=B^{-1} * F ###
            for m1 in range(mstart, M + 1):
                Bf2 = qall[m1].Broyden(fxall[m1])
                #Bf2 = [0]*n
                for i in range(n):
                    Jfxall[m1][i] = fxall[m1][i];
                    Jfvall[m1][i] = fvall[m1][i] + 0.5 * timestep * Bf2[i]
                    Jfxall1[m1][i] = Jfxall[m1][i];
                    Jfvall1[m1][i] = Jfvall[m1][i]
                    Jfxall2[m1][i] = Jfxall[m1][i];
                    Jfvall2[m1][i] = Jfvall[m1][i]

            L = 15
            for k in range(L):
                ### compute JF2=Q * JF2 ; Q = (1- B^{-1})*DF ###

                ### JF3=C*JF2 ###
                for m1 in range(mstart, M + 1):
                    Bf2 = qall[m1 - 1].Broyden(Jfxall2[m1 - 1])
                    #Bf2 = [0]*n
                    for i in range(n):
                        Jfxall3[m1][i] = Jfxall2[m1 - 1][i] + timestep * Jfvall2[m1 - 1][
                            i] + 0.5 * timestep * timestep * Bf2[i]
                        Jfvall3[m1][i] = Jfvall2[m1 - 1][i] + 0.5 * timestep * Bf2[i]

                ### JF2=B^{-1}*JF3 ###
                for m1 in range(mstart, M + 1):
                    Bf2 = qall[m1].Broyden(Jfxall3[m1])
                    #Bf2 = [0]*n
                    for i in range(n):
                        Jfxall2[m1][i] = Jfxall3[m1][i]
                        Jfvall2[m1][i] = Jfvall3[m1][i] + 0.5 * timestep * Bf2[i]

                ### JF1 += JF2
                for m1 in range(mstart, M + 1):
                    for i in range(n):
                        Jfxall1[m1][i] += Jfxall2[m1][i]
                        Jfvall1[m1][i] += Jfvall2[m1][i]

            ### compute Q*JF1 ###
            for m1 in range(mstart, M + 1):
                Bf2 = qall[m1 - 1].Broyden(Jfxall1[m1 - 1])
                for i in range(n):
                    Jfxall3[m1][i] = Jfxall1[m1 - 1][i] + timestep * Jfvall1[m1 - 1][i] + 0.5 * timestep * timestep * \
                                                                                          Bf2[i]
                    Jfvall3[m1][i] = Jfvall1[m1 - 1][i] + 0.5 * timestep * Bf2[i]
            for m1 in range(mstart, M + 1):
                Bf2 = qall[m1].Broyden(Jfxall3[m1])
                for i in range(n):
                    Jfxall2[m1][i] = Jfxall3[m1][i]
                    Jfvall2[m1][i] = Jfvall3[m1][i] + 0.5 * timestep * Bf2[i]

            ### compute P*JF1 ###
            for m1 in range(mstart, M + 1):
                Bf2 = qcoarse[m1 - 1].Broyden(Jfxall1[m1 - 1])
                for i in range(n):
                    Jfxall4[m1][i] = Jfxall1[m1 - 1][i] + timestep * Jfvall1[m1 - 1][i] + 0.5 * timestep * timestep * \
                                                                                          Bf2[i]
                    Jfvall4[m1][i] = Jfvall1[m1 - 1][i] + 0.5 * timestep * Bf2[i]
            for m1 in range(mstart, M + 1):
                Bf2 = qcoarse[m1].Broyden(Jfxall4[m1])
                for i in range(n):
                    Jfxall3[m1][i] = Jfxall4[m1][i]
                    Jfvall3[m1][i] = Jfvall4[m1][i] + 0.5 * timestep * Bf2[i]
            for m1 in range(mstart, M + 1):
                for i in range(n):
                    Jfxall[m1][i] += Jfxall2[m1][i] - Jfxall3[m1][i]
                    Jfvall[m1][i] += Jfvall2[m1][i] - Jfvall3[m1][i]

            ### compute A*JF ###
            for m1 in range(mstart, M + 1):
                Bf2 = qcoarse[m1].Broyden(Jfxall[m1])
                for i in range(n):
                    Jfvall[m1][i] -= 0.5 * timestep * Bf2[i]

            ### update wall and zall ###
            for m1 in range(mstart, M + 1):
                for i in range(n):
                    wall[m1][i] -= Jfxall[m1][i]
                    zall[m1][i] -= Jfvall[m1][i]

            #norm1 = 0.0;
            #norm2 = 0.0;
            #norm3 = 0.0;
            #for m1 in range(mstart,M+1):
            #   for i in range(n):
            #      norm1 += fxall[m1][i]*fxall[m1][i]    + fvall[m1][i]*fvall[m1][i]
            #      norm2 += Jfxall[m1][i]*Jfxall[m1][i]  + Jfvall[m1][i]*Jfvall[m1][i]
            #      norm3 += Jfxall[m1][i]*fxall[m1][i]   + Jfvall[m1][i]*fvall[m1][i]
            #print "norm1,norm2,norm3,residx,residv,mstart0=",norm1,norm2,norm3,residx,residv,mstart0

            ### compute x=invG(w) - compute initial coarse grain path ###
            for m1 in range(M):
                for i in range(n):
                    xall[m1 + 1][i] = wall[m1 + 1][i] + xall[m1][i] + vall[m1][i] * timestep + 0.5 * aall2[m1][
                        i] * timestep * timestep
                result = runcalcs(hostsports, jobcoarse, [xall[m1 + 1]])  #result：from server cal
                for ii in range(nion):
                    aall2[m1 + 1][3 * ii] = result[0][1][3 * ii] / mass[ii]
                    aall2[m1 + 1][3 * ii + 1] = result[0][1][3 * ii + 1] / mass[ii]
                    aall2[m1 + 1][3 * ii + 2] = result[0][1][3 * ii + 2] / mass[ii]
                for i in range(n):
                    vall[m1 + 1][i] = zall[m1 + 1][i] + vall[m1][i] + 0.5 * (aall2[m1][i] + aall2[m1 + 1][i]) * timestep


            ### compute F(x) ###
            result = runcalcs(hostsports, job, xall, mstart)  #result：from server cal
            it += 1
            for m1 in range(mstart, M + 1):
                uall[m1] = result[m1 - mstart][0]
                for ii in range(nion):
                    aall[m1][3 * ii] = result[m1 - mstart][1][3 * ii] / mass[ii]
                    aall[m1][3 * ii + 1] = result[m1 - mstart][1][3 * ii + 1] / mass[ii]
                    aall[m1][3 * ii + 2] = result[m1 - mstart][1][3 * ii + 2] / mass[ii]
            for m1 in range(mstart, M + 1):
                keall[m1] = 0.0
                for ii in range(nion):
                    keall[m1] += 0.5 * mass[ii] * (
                        vall[m1][3 * ii] * vall[m1][3 * ii] + vall[m1][3 * ii + 1] * vall[m1][3 * ii + 1] + vall[m1][
                            3 * ii + 2] * vall[m1][3 * ii + 2])

            if (abs(residx0 + residv0 - residx - residv) > maxresiderr):
                mstart = mstart0;
                diff = (M + 1 - mstart + 1) - 2
                if (diff < 0): mstart += diff;
            else:
                mstart = 1;
            residx0 = residx
            residv0 = residv

            #done  = ( (it>M+2) or ((residx+residv) < maxresiderr)  )

    iterations = it
    for i in range(n):
        x1[i] = xall[M][i]
        v1[i] = vall[M][i]
        a1[i] = aall[M][i]
        acoarse1[i] = aall2[M][i]
    ke1 = keall[M];
    u1 = uall[M];
    e1 = u1 + ke1

    return (x1, v1, a1, acoarse1, e1, ke1, u1, iterations)


###########################################################################
########################### main program ##################################
###########################################################################

hostsports = eval(raw_input("Enter hostsports: "))  #['localhost:50001','localhost:50002']
xyzfilename0 = raw_input("Enter initial xyz filename: ")  #HCl_4water.00.xyz
xyzfilename1 = raw_input("Enter final xyz filename: ")  #HCl_4water.01.xyz
xyzfilename2 = raw_input("Enter trajectory xyz filename: ")  #HCl_4water.traj.xyz
theory = raw_input("Enter theory: ")  #scf/mp2
if (theory == 'lj'):
    epsilon = eval(raw_input("Enter epsilon: "))
    rmin = eval(raw_input("Enter rmin: "))
    basis = 'md potential'
elif (theory == 'spring'):
    epsilon = eval(raw_input("Enter epsilon: "))
    rmin = eval(raw_input("Enter rmin: "))
    basis = 'md potential'
else:
    basis = raw_input("Enter basis: ")  #3-21G
    mult = eval(raw_input("Enter mult: "))  #1
    charge = eval(raw_input("Enter charge: "))  #0

islinear = eval(raw_input("Is this a linear molecule? "))  #False
print
timestep = eval(raw_input("Enter timestep in au: "))  #5.0
m = eval(raw_input("Enter length of parallel time segment (M): "))  #10
Ntimesteps = eval(raw_input("Enter total number of timesteps (Ntimesteps): "))  #20
residerr = eval(raw_input("Enter maximum residual error (residerr): "))  #1.0e-4
ParallelTimeJob = eval(raw_input("Is this a parallel time job? "))  #True/False
print
precondition = eval(raw_input("Is this a preconditioned parallel time job? "))  #False
if precondition:
    pretheory = raw_input("Enter preconditioning theory: ")
    if (pretheory == 'lj'):
        preepsilon = eval(raw_input("Enter preconditioning epsilon: "))
        prermin = eval(raw_input("Enter preconditioning rmin: "))
        prebasis = 'md potential'
    elif (pretheory == 'spring'):
        #preKspring = eval(raw_input("Enter preconditioning Kspring: "))
        preepsilon = eval(raw_input("Enter preconditioning epsilon: "))
        prermin = eval(raw_input("Enter preconditioning rmin: "))
        prebasis = 'md potential'
    else:
        prebasis = raw_input("Enter preconditioning basis: ")


## read in xyzfile ##
xyzfile = open(xyzfilename0, 'r')
symbols = [];
x0 = [];
v0 = [];
mass = []
nion = eval(xyzfile.readline())
xyzfile.readline()
for ii in range(nion):
    line = xyzfile.readline().split()
    sym = line[0].lower().capitalize();
    symbols.append(sym)
    i = def_symbols.index(sym);
    mass.append(def_masses[i] * 1822.89)
    x0.append(eval(line[1]) / 0.529177);
    x0.append(eval(line[2]) / 0.529177);
    x0.append(eval(line[3]) / 0.529177)
    if (len(line) == 7):
        v0.append(eval(line[4]) / 0.529177);
        v0.append(eval(line[5]) / 0.529177);
        v0.append(eval(line[6]) / 0.529177)
    else:
        v0.append(0.0);
        v0.append(0.0);
        v0.append(0.0)
xyzfile.close()

#ParallelTimeJob = True
#Ntimesteps = 400
#m    = 20
#timestep   = 5.0
#residerr   = 1.0e-4
Perror = 1.0
tempset = 300.0

cpu1 = time.time()

nwjob = {}
nwjob['theory'] = theory
nwjob['symbols'] = symbols
nwjob['nion'] = nion
if (theory == 'lj'):
    sigma = rmin / 2 ** (1.0 / 6.0)
    nwjob['epsilon'] = []
    nwjob['sigma'] = []
    for ii in range(nion):
        nwjob['epsilon'].append(epsilon)
        nwjob['sigma'].append(sigma)
elif (theory == 'spring'):
    Kspring = 72.0 * epsilon / (rmin * rmin)
    nwjob['Kspring'] = Kspring
    nwjob['r0spring'] = rmin
else:
    nwjob['basis'] = "basis \"ao basis\" cartesian print\n" + "  * library " + basis + "\n" + "end\n"
    nwjob['mult'] = mult
    nwjob['charge'] = charge

#precondition ：false
if precondition:
    coarsenwjob = {}
    coarsenwjob['theory'] = pretheory
    coarsenwjob['symbols'] = symbols
    coarsenwjob['nion'] = nion
    if (pretheory == 'lj'):
        presigma = prermin / 2 ** (1.0 / 6.0)
        coarsenwjob['epsilon'] = []
        coarsenwjob['sigma'] = []
        for ii in range(nion):
            coarsenwjob['epsilon'].append(preepsilon)
            coarsenwjob['sigma'].append(presigma)
    elif (pretheory == 'spring'):
        preKspring = 72.0 * preepsilon / (rmin * rmin)
        coarsenwjob['Kspring'] = preKspring
        coarsenwjob['r0spring'] = prermin
    else:
        coarsenwjob['basis'] = "basis \"ao basis\" cartesian print\n" + "  * library " + prebasis + "\n" + "end\n"
        coarsenwjob['mult'] = mult
        coarsenwjob['charge'] = charge

print "          ****************************************************"
print "          *                                                  *"
print "          *      Parallel in Time Molecular Dynamics         *"
print "          *                                                  *"
print "          *     [      Velocity Verlet Integration    ]      *"
print "          *     [     Master/Slave Parallelization    ]      *"
print "          *     [     Python Socket Implementation    ]      *"
print "          *                                                  *"
print "          *            version #1.00   02/03/12              *"
print "          *                                                  *"
print "          ****************************************************"
print "          >>> job started at       ", util_date(), " <<<"
print
print "number of HOSTSPORTS =", len(hostsports)
print "HOSTSPORTS =", hostsports
print
print "Parallel in Time Job             =", ParallelTimeJob
print
if precondition:
    print "theory/basis=", theory + "/" + basis + "   (preconditioner " + pretheory + "/" + prebasis + ")"
else:
    print "theory/basis=", theory + "/" + basis
print
print "number of atoms                  =", nion
print "length of parallel time segments =", m
print "timestep                         =", timestep
print "Total number of timesteps        =", Ntimesteps
print "relative energy threshold        =", Perror
print "maximum residual error           =", residerr
print "linear molecule                  =", islinear
#print "Temperature scale                =", tempset
print
print "initial xyz filename    = ", xyzfilename0
print "final   xyz filename    = ", xyzfilename1
print "trajectory xyz filename = ", xyzfilename2
print
print "intitial geometry (a.u.):"
for ii in range(nion):
    print '%d  %s  %8.4f %8.4f %8.4f  %8.4f %8.4f %8.4f - mass = %8.2f' % (
        ii + 1, nwjob['symbols'][ii], x0[3 * ii], x0[3 * ii + 1], x0[3 * ii + 2], v0[3 * ii], v0[3 * ii + 1],
        v0[3 * ii + 2], mass[ii] / 1822.89)
print
print
print "          ================ iteration ========================="
print
print "          >>> iteration started at ", util_date(), " <<<"
print
print "   iter.             Energy    Potential      Kinetic   Temperature  Iterations"
print "   ----------------------------------------------------------------------------"
sys.stdout.flush()

xyzfile = open(xyzfilename2, 'w')
xyzfile.write('%d\n\n' % nion)
for ii in range(nion):
    xyzfile.write('%s  %e %e %e\n' % (
        nwjob['symbols'][ii], x0[3 * ii] * 0.529177, x0[3 * ii + 1] * 0.529177, x0[3 * ii + 2] * 0.529177))
xyzfile.flush()

itave = 0
Nit = Ntimesteps / m
kesum = 0.0
cpu2 = time.time()

## calculate initial acceleration ##
a0 = [0] * 3 * nion
result = runcalcs(hostsports, nwjob, [x0])  #client->server with data nwjob
u0 = result[0][0]
for ii in range(nion):
    a0[3 * ii] = result[0][1][3 * ii] / mass[ii]
    a0[3 * ii + 1] = result[0][1][3 * ii + 1] / mass[ii]
    a0[3 * ii + 2] = result[0][1][3 * ii + 2] / mass[ii]
cpu2a0 = time.time()

## calculate initial coarse acceleration ##
acoarse0 = [0] * 3 * nion
#preconditon：False
if precondition:
    result = runcalcs(hostsports, coarsenwjob, [x0])
    for ii in range(nion):
        acoarse0[3 * ii] = result[0][1][3 * ii] / mass[ii]
        acoarse0[3 * ii + 1] = result[0][1][3 * ii + 1] / mass[ii]
        acoarse0[3 * ii + 2] = result[0][1][3 * ii + 2] / mass[ii]

for it in range(Nit):
    if ParallelTimeJob:
        if precondition:
            #execute this
            x1, v1, a1, acoarse1, e1, ke1, u1, iterations = mdverlet_parallel2(hostsports, nwjob, coarsenwjob, m,
                                                                               timestep, Perror, residerr, mass, x0, v0,
                                                                               a0, u0, acoarse0)
        else:
            x1, v1, a1, e1, ke1, u1, iterations = mdverlet_parallel(hostsports, nwjob, m, timestep, Perror, residerr,
                                                                    mass, x0, v0, a0, u0)
    else:
        x1, v1, a1, e1, ke1, u1, iterations = mdverlet_serial(hostsports, nwjob, m, timestep, Perror, residerr, mass,
                                                              x0, v0, a0, u0)
    itave += iterations
    kesum += ke1;
    keave = kesum / float(it + 1)
    if islinear:
        temp = 2.0 * keave / (3.0 * nion - 5.0);
    else:
        temp = 2.0 * keave / (3.0 * nion - 6.0);
    temp /= 3.16679e-6;
    #if (it>20):
    #   stemp = 1.0
    #else:
    #   stemp = math.sqrt(tempset/temp)
    stemp = 1.0

    tup = ((it + 1) * m, e1, u1, ke1, temp, iterations)
    print '%8d        %11.6f  %11.6f  %11.6f       %7.1f %11d' % tup
    sys.stdout.flush()

    for i in range(3 * nion):
        x0[i] = x1[i]
        v0[i] = v1[i] * stemp
        a0[i] = a1[i]
    if precondition:
        for i in range(3 * nion):
            acoarse0[i] = acoarse1[i]

    xyzfile.write('%d\n\n' % nion)
    for ii in range(nion):
        xyzfile.write('%s  %e %e %e\n' % (
            nwjob['symbols'][ii], x0[3 * ii] * 0.529177, x0[3 * ii + 1] * 0.529177, x0[3 * ii + 2] * 0.529177))
    xyzfile.flush()

xyzfile.close()
cpu3 = time.time()
print "          >>> iteration ended at   ", util_date(), " <<<"
print
print
print "final geometry (a.u.):"
for ii in range(nion):
    print '%d  %s  %8.4f %8.4f %8.4f  %8.4f %8.4f %8.4f - mass = %8.2f' % (
        ii + 1, nwjob['symbols'][ii], x0[3 * ii], x0[3 * ii + 1], x0[3 * ii + 2], v0[3 * ii], v0[3 * ii + 1],
        v0[3 * ii + 2], mass[ii] / 1822.89)
print

## print out final geometry ##
xyzfile = open(xyzfilename1, 'w')
xyzfile.write('%d\n\n' % nion)
for ii in range(nion):
    xyzfile.write('%s  %18.9e %18.9e %18.9e %18.9e %18.9e %18.9e\n' % (
        nwjob['symbols'][ii], x0[3 * ii] * 0.529177, x0[3 * ii + 1] * 0.529177, x0[3 * ii + 2] * 0.529177,
        v0[3 * ii] * 0.529177, v0[3 * ii + 1] * 0.529177, v0[3 * ii + 2] * 0.529177))
xyzfile.close()

cpu4 = time.time()
t1 = cpu2 - cpu1
t2 = cpu3 - cpu2
t3 = cpu4 - cpu3
t4 = cpu4 - cpu1
av = t2 / (m * Nit)

maxspeedup = m / (itave / float(Nit))
speedup = (cpu2a0 - cpu2) / av

print "average iterations     =", itave / float(Nit)
print "maximum speedup        =", maxspeedup
print "approximate speedup    =", speedup
print "approximate efficiency =", speedup / maxspeedup
print
print "-----------------"
print "cputime in seconds"
print "prologue       : ", t1
print "main loop      : ", t2
print "epilogue       : ", t3
print "total          : ", t4
print "cputime/step   : ", av
print "cpuserial/step : ", cpu2a0 - cpu2
print

