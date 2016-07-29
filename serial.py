__author__ = 'ShuD'

import sys,time,pickle,math
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

def util_date():
    tt = time.localtime()
    dd = "%d-%d-%d-%d:%d" % (tt[0], tt[1], tt[2], tt[3], tt[4])
    return dd


def run_simulation(x):
    if (x['theory'] == 'lj'):
        data = run_ljgradient(x)
    if (x['theory'] == 'spring'):
        data = run_spring(x)
    return data


def run_jobs(jobs):
    egrads = []
    njobs = len(jobs)
    if (njobs == 1):
        egrads += [run_simulation(jobs[0])]
    else:
        ## route to ss sockets, or route to ss sockets and calculate locally ##
        ok = (not onlyroute) or (len(ss) > 0)
        if (ok):
            if (onlyroute):
                ns = len(ss)
                jobcount = [0] * ns
                for i in range(ns):                  jobcount[i] = int(njobs * sscpus[i] / float(ncpus))
                for i in range(njobs - sum(jobcount)): jobcount[i % ns] += 1
                #jobcount = [njobs/ns]*ns
                #for i in range(njobs%ns): jobcount[i] += 1

                ## send jobs to ss ##
                i = -1
                istart = 0
                for s in ss:
                    iend = istart + jobcount[i + 1]
                    s.sendall(pickle.dumps(jobs[istart:iend]) + 'ericdone')
                    istart = iend
                    i += 1
            else:
                ns = len(ss) + 1
                sscpus = [1] + sscpus
                ncpus += 1
                jobcount = [0] * ns
                for i in range(ns):                  jobcount[i] = int(njobs * sscpus[i] / float(ncpus))
                for i in range(njobs - sum(jobcount)): jobcount[i % ns] += 1

                ## send jobs to ss ##
                i = 0
                istart = jobcount[i]
                for s in ss:
                    iend = istart + jobcount[i + 1]
                    s.sendall(pickle.dumps(jobs[istart:iend]) + 'ericdone')
                    istart = iend
                    i += 1

                ## cacluate jobs locally ##
                for i in range(jobcount[0]): egrads += [run_simulation(PORT, jobs[i])]

            ## read from ss and close. Also keep track of badjobs and goodports ##
            if (onlyroute):
                i = -1
                istart = 0
            else:
                i = 0
                istart = jobcount[i]
            j = 0
            badjobs = []
            newhostsports = []
            for s in ss:
                cpu1 = time.time()
                iend = istart + jobcount[i + 1]
                msg = ''
                done = False
                while not done:
                    data = s.recv(1024)
                    msg = msg + data
                    done = msg[len(msg) - 8:len(msg)] == 'ericdone'
                msg = msg[:len(msg) - 8]
                data = pickle.loads(msg)
                egrads += data
                if (data == [None]):
                    badjobs.append(jobs[istart:iend])
                else:
                    newhostsports.append(hostsports[j])
                istart = iend
                i += 1
                j += 1
                s.close()
                cpu2 = time.time()
                sstimes[j - 1] = cpu2 - cpu1

            print 'sstimes=', sstimes[:ns]
            ## rerun badjobs ##
            for j in range(len(badjobs)):
                data = run_jobs(onlyroute, PORT, newhostsports, badjobs[j])
                i = egrads.index(None)
                egrads = egrads[:i] + data + egrads[i + 1:]

        ## failed because onlyroute==True and ss==[] ##
        else:
            egrads = [None]

    return egrads

#############################################
#                                           #
#             run_spring                    #
#                                           #
#############################################
#
def run_spring(job):
    x = job['xyz']
    K = job['Kspring']
    r0 = job['r0spring']
    n = len(x)
    nion = n / 3
    f = [0] * n
    e = 0
    for jj in range(nion - 1):
        ii = jj + 1
        dx = x[3 * ii] - x[3 * jj]
        dy = x[3 * ii + 1] - x[3 * jj + 1]
        dz = x[3 * ii + 2] - x[3 * jj + 2]
        r = math.sqrt(dx * dx + dy * dy + dz * dz)
        e += 0.5 * K * (r - r0) * (r - r0)
        dvdr = K * (r - r0)
        f[3 * ii] -= dx * (dvdr / r)
        f[3 * ii + 1] -= dy * (dvdr / r)
        f[3 * ii + 2] -= dz * (dvdr / r)
        f[3 * jj] += dx * (dvdr / r)
        f[3 * jj + 1] += dy * (dvdr / r)
        f[3 * jj + 2] += dz * (dvdr / r)
    return (e, f)


#############################################
#                                           #
#             run_ljgradient                #
#                                           #
#############################################
#
def run_ljgradient(job):
    x = job['xyz']
    eps = job['epsilon']
    sig = job['sigma']
    n = len(x)
    nion = n / 3
    f = [0] * n
    e = 0
    for jj in range(nion):
        for ii in range(jj + 1, nion):
            epsilon = math.sqrt(eps[ii] * eps[jj])
            sigma = 0.5 * (sig[ii] + sig[jj])
            dx = x[3 * ii] - x[3 * jj]
            dy = x[3 * ii + 1] - x[3 * jj + 1]
            dz = x[3 * ii + 2] - x[3 * jj + 2]
            r = math.sqrt(dx * dx + dy * dy + dz * dz)
            u = (sigma / r)
            u6 = u * u * u * u * u * u
            u12 = u6 * u6
            e += 4.0 * epsilon * (u12 - u6)
            dvdr = -4.0 * (epsilon / r) * (12 * u12 - 6 * u6)
            f[3 * ii] -= dx * (dvdr / r)
            f[3 * ii + 1] -= dy * (dvdr / r)
            f[3 * ii + 2] -= dz * (dvdr / r)
            f[3 * jj] += dx * (dvdr / r)
            f[3 * jj + 1] += dy * (dvdr / r)
            f[3 * jj + 2] += dz * (dvdr / r)
    return (e, f)

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
        for i in range(ns):
            jobcount[i] = int(njobs * sscpus[i] / float(ncpus))
        for i in range(njobs - sum(jobcount)):
            jobcount[i % ns] += 1
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



###########################################################################
########################### main program ##################################
###########################################################################
xyzfilename0 = raw_input("Enter initial xyz filename: ")  #HCl_4water.00.xyz
xyzfilename1 = raw_input("Enter final xyz filename: ")  #HCl_4water.01.xyz
xyzfilename2 = raw_input("Enter trajectory xyz filename: ")  #HCl_4water.traj.xyz
theory = raw_input("Enter theory: ")  #scf/mp2/lj/spring
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

Perror = 1.0
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

print "          >>> job started at       ", util_date(), " <<<"
print "theory/basis=", theory + "/" + basis
print "number of atoms                  =", nion
print "Total number of timesteps        =", Ntimesteps
print "relative energy threshold        =", Perror
print "maximum residual error           =", residerr
print "linear molecule                  =", islinear
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
xyzfile.flus

Nit = Ntimesteps / m
for it in range(Nit):
     x1, v1, a1, e1, ke1, u1, iterations = mdverlet_serial(hostsports, nwjob, m, timestep, Perror, residerr, mass,
                                                              x0, v0, a0, u0)