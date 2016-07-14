#!/usr/bin/python

##### other choices
#### #!/usr/local/bin/python
#### #!/usr/bin/python
#### #!/Library/Frameworks/Python.framework/Versions/2.6/bin/python


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

import sys, os, socket, pickle, math, time

# nwchem_binary = '/home/bylaska/bin/nwchem '
#nwchem_binary = '/Users/bylaska/bin/nwchem-6.1 '
#nwchem_binary = '/usr/bin/nwchem '
#nwchem_binary = '/home/yic017/nwchem-releases/nwchem/bin/LINUX64/nwchem '

nwchem_binary = '/Users/bylaska/bin/nwchem '


#############################################
#                                           #
#             my_sockconnect                #
#                                           #
#############################################
#
# Returns connected socket to hostport or None if failed
# where hostport = "ipaddress:port[:ncpus]" (e.g. hostport="localhost:50001:4")
#
def my_sockconnect(hostport):
    l = hostport.index(':')
    r = hostport.rindex(':')
    if (r == l):
        hp = (hostport[:l], eval(hostport[l + 1:]))
    else:
        hp = (hostport[:l], eval(hostport[l + 1:r]))
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.settimeout(1.0)
    try:
        s.connect(hp)
        s.settimeout(None)
        return s
    except:
        s.close()
        s = None
    return s


#############################################
#                                           #
#             run_grgradient                #
#                                           #
#############################################
#
def wbw_splint(xa, ya, y2a, nx, x):
    khi = nx + 1
    klo = nx
    while ((xa[klo] > x) or (xa[khi] < x)):
        if (xa[klo] > x):
            klo = klo - 1
            khi = khi - 1
        if (xa[khi] < x):
            klo = klo + 1
            khi = khi + 1
    h = xa[khi] - xa[klo]
    a = (xa[khi] - x) / h
    b = (x - xa[klo]) / h
    y = a * ya[klo] + b * ya[khi] + ((a ** 3 - a) * y2a[klo] + (b ** 3 - b) * y2a[khi]) * h ** 2 / 6.0
    return y


def wbw_dsplint(xa, ya, y2a, nx, x):
    khi = nx + 1
    klo = nx
    while ((xa[klo] > x) or (xa[khi] < x)):
        if (xa[klo] > x):
            klo = klo - 1
            khi = khi - 1
        if (xa[khi] < x):
            klo = klo + 1
            khi = khi + 1
    h = xa[khi] - xa[klo]
    a = (xa[khi] - x) / h
    b = (x - xa[klo]) / h
    da = -1.0 / h
    db = 1.0 / h
    dy = da * ya[klo] + db * ya[khi] + ((da * 3 * a ** 2 - da) * y2a[klo] + (db * 3 * b ** 2 - db) * y2a[
        khi]) * h ** 2 / 6.0
    return dy


def run_grgradient(job):
    x = job['xyz']
    xgrid = job['xgrid']
    ygrid = job['ygrid']
    ygrid2 = job['ygrid2']
    cutoff = job['cutoff']
    l = job['l']
    n = len(x)
    nion = n / 3
    f = [0] * n
    e = 0.0

    ### compute nlist and nn ###
    cutoff2 = cutoff * cutoff
    nlist = [0] * nion
    nn = [-1] * nion
    for jj in range(nion - 1):
        for ii in range(jj + 1, nion):
            dx = x[3 * ii] - x[3 * jj]
            dy = x[3 * ii + 1] - x[3 * jj + 1]
            dz = x[3 * ii + 2] - x[3 * jj + 2]
            dx = dx - l * round(dx / l)
            dy = dy - l * round(dy / l)
            dz = dz - l * round(dz / l)
            dij = dx * dx + dy * dy + dz * dz
            if (dij < cutoff2):
                nlist[ii] += 1
                nlist[jj] += 1
                if (nn[ii] == -1):
                    nn[ii] = [jj]
                else:
                    nn[ii].append(jj)
                if (nn[jj] == -1):
                    nn[jj] = [ii]
                else:
                    nn[jj].append(ii)
    for ii in range(nion):
        if (nlist[ii] > 0):
            for jj in nn[ii]:
                dx = x[3 * ii] - x[3 * jj]
                dy = x[3 * ii + 1] - x[3 * jj + 1]
                dz = x[3 * ii + 2] - x[3 * jj + 2]
                dx = dx - l * round(dx / l)
                dy = dy - l * round(dy / l)
                dz = dz - l * round(dz / l)
                dij = math.sqrt(dx * dx + dy * dy + dz * dz)
                nr = int(dij / (xgrid[1] - xgrid[0]))

                tmp = wbw_splint(xgrid, ygrid, ygrid2, nr, dij)
                e += tmp

                tmp2 = wbw_dsplint(xgrid, ygrid, ygrid2, nr, dij) / dij
                f[3 * ii] -= tmp2 * dx
                f[3 * ii + 1] -= tmp2 * dy
                f[3 * ii + 2] -= tmp2 * dz
    e *= 0.5
    return (e, f)


#############################################
#                                           #
#             run_swgradient                #
#                                           #
#############################################
#
def run_swgradient(job):
    x = job['xyz']
    a0 = job['a0']
    b0 = job['b0']
    eps = job['eps']
    sig = job['sig']
    lamb = job['lamb']
    sa0 = job['sa0']
    gam = job['gam']
    p = job['p']
    l = job['l']

    n = len(x)
    nion = n / 3
    f = [0] * n
    e = 0.0

    ### compute nlist and nn ###
    cutoff = sig * sa0
    cutoff2 = cutoff * cutoff
    nlist = [0] * nion
    nn = [-1] * nion
    for ii in range(nion - 1):
        for jj in range(ii + 1, nion):
            dx = x[3 * ii] - x[3 * jj]
            dy = x[3 * ii + 1] - x[3 * jj + 1]
            dz = x[3 * ii + 2] - x[3 * jj + 2]
            dx = dx - l * round(dx / l)
            dy = dy - l * round(dy / l)
            dz = dz - l * round(dz / l)
            dij = dx * dx + dy * dy + dz * dz
            if (dij < cutoff2):
                nlist[ii] += 1
                nlist[jj] += 1
                if (nn[ii] == -1):
                    nn[ii] = [jj]
                else:
                    nn[ii].append(jj)
                if (nn[jj] == -1):
                    nn[jj] = [ii]
                else:
                    nn[jj].append(ii)


    ### compute the two-body potential ###
    const = eps * a0
    for ii in range(nion):
        if (nlist[ii] > 0):
            for j in range(nlist[ii]):
                jj = nn[ii][j]
                dx = x[3 * ii] - x[3 * jj]
                dy = x[3 * ii + 1] - x[3 * jj + 1]
                dz = x[3 * ii + 2] - x[3 * jj + 2]
                dx = dx - l * round(dx / l)
                dy = dy - l * round(dy / l)
                dz = dz - l * round(dz / l)
                dij = math.sqrt(dx * dx + dy * dy + dz * dz)

                tmp = dij / sig
                tmp2 = b0 * (tmp ** (-p)) - 1.0
                tmp3 = math.exp(1.0 / (tmp - sa0))
                e = e + const * tmp2 * tmp3

                tmp2 = -tmp3 * p * b0 * (tmp ** (-p - 1.0)) - tmp2 * tmp3 / ((tmp - sa0) ** 2)
                tmp2 = const * tmp2 / (sig * dij)
                f[3 * ii] += -tmp2 * dx
                f[3 * ii + 1] += -tmp2 * dy
                f[3 * ii + 2] += -tmp2 * dz
    e *= 0.5

    ### compute the three-body potential ###
    const = eps * lamb
    for ii in range(nion):
        if (nlist[ii] > 1):
            for j in range(nlist[ii] - 1):
                jj = nn[ii][j]
                dxij = x[3 * ii] - x[3 * jj]
                dyij = x[3 * ii + 1] - x[3 * jj + 1]
                dzij = x[3 * ii + 2] - x[3 * jj + 2]
                dxij = dxij - l * round(dxij / l)
                dyij = dyij - l * round(dyij / l)
                dzij = dzij - l * round(dzij / l)
                dij = math.sqrt(dxij * dxij + dyij * dyij + dzij * dzij)
                tmpij = dij / sig

                for k in range(j + 1, nlist[ii]):
                    kk = nn[ii][k]
                    dxik = x[3 * ii] - x[3 * kk]
                    dyik = x[3 * ii + 1] - x[3 * kk + 1]
                    dzik = x[3 * ii + 2] - x[3 * kk + 2]
                    dxik = dxik - l * round(dxik / l)
                    dyik = dyik - l * round(dyik / l)
                    dzik = dzik - l * round(dzik / l)
                    dik = math.sqrt(dxik * dxik + dyik * dyik + dzik * dzik)
                    tmpik = dik / sig

                    costheta = (dxij * dxik + dyij * dyik + dzij * dzik) / (dij * dik)

                    tmp2 = costheta + 1.0 / 3.0
                    tmp3 = gam / (tmpij - sa0) + gam / (tmpik - sa0)
                    e0 = const * tmp2 * tmp2 * math.exp(tmp3)
                    e += e0

                    tmp4 = -gam * e0 / ((tmpij - sa0) * (tmpij - sa0) * sig * dij)
                    tmp5 = -gam * e0 / ((tmpik - sa0) * (tmpik - sa0) * sig * dik)

                    f[3 * ii] -= (tmp4 * dxij + tmp5 * dxik)
                    f[3 * ii + 1] -= (tmp4 * dyij + tmp5 * dyik)
                    f[3 * ii + 2] -= (tmp4 * dzij + tmp5 * dzik)

                    f[3 * jj] += tmp4 * dxij
                    f[3 * jj + 1] += tmp4 * dyij
                    f[3 * jj + 2] += tmp4 * dzij

                    f[3 * kk] += tmp5 * dxik
                    f[3 * kk + 1] += tmp5 * dyik
                    f[3 * kk + 2] += tmp5 * dzik

                    tmp4 = 2.0 * const * tmp2 * math.exp(tmp3)
                    tmp5 = costheta

                    f[3 * ii] += -tmp4 * (
                    (dxij + dxik) / (dij * dik) - tmp5 * (dxij / (dij * dij) + dxik / (dik * dik)))
                    f[3 * ii + 1] += -tmp4 * (
                    (dyij + dyik) / (dij * dik) - tmp5 * (dyij / (dij * dij) + dyik / (dik * dik)))
                    f[3 * ii + 2] += -tmp4 * (
                    (dzij + dzik) / (dij * dik) - tmp5 * (dzij / (dij * dij) + dzik / (dik * dik)))

                    f[3 * jj] += tmp4 * (dxik / (dij * dik) - tmp5 * dxij / (dij * dij))
                    f[3 * jj + 1] += tmp4 * (dyik / (dij * dik) - tmp5 * dyij / (dij * dij))
                    f[3 * jj + 2] += tmp4 * (dzik / (dij * dik) - tmp5 * dzij / (dij * dij))

                    f[3 * kk] += tmp4 * (dxij / (dij * dik) - tmp5 * dxik / (dik * dik))
                    f[3 * kk + 1] += tmp4 * (dyij / (dij * dik) - tmp5 * dyik / (dik * dik))
                    f[3 * kk + 2] += tmp4 * (dzij / (dij * dik) - tmp5 * dzik / (dik * dik))

    return (e, f)


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


#############################################
#                                           #
#             make_nwinput                  #
#                                           #
#############################################
#
# returns an nwchem input deck for pspw, dft, or mp2 gradient calculation
#
def make_nwinput(PORT, job):
    ss = ""
    if (job['theory'] == 'pspw'):
        ss += "title \"automatic pspw nwchemjob\"\n\n" + "start nwchempspw%d\n\n" % PORT
    elif (job['theory'] == 'dft'):
        ss += "title \"automatic dft  nwchemjob\"\n\n" + "start nwchemdft%d\n\n" % PORT
    elif (job['theory'] == 'mp2'):
        ss += "title \"automatic mp2  nwchemjob\"\n\n" + "start nwchemmp2%d\n\n" % PORT
    elif (job['theory'] == 'scf'):
        ss += "title \"automatic scf  nwchemjob\"\n\n" + "start nwchemscf%d\n\n" % PORT
    else:
        ss += "title \"automatic nwchemjob\"\n\n" + "start nwchemother%d\n\n" % PORT
    ss += "memory 1500 mb\n\n"
    if (job.has_key('charge')): ss += "charge %d\n" % (job['charge'])
    ss += "geometry units bohrs noautosym noautoz nocenter\n"
    for ii in range(len(job['symbols'])):
        ss += "%s  %le %le %le\n" % (
        job['symbols'][ii], job['xyz'][3 * ii], job['xyz'][3 * ii + 1], job['xyz'][3 * ii + 2])
    ss += "end\n"
    if (job['theory'] == 'pspw'):
        ss += "nwpw\n" + "  simulation_cell\n" + "    fcc 38.0\n" + "  end\n"
        if (job.has_key('cutoff')): ss += "  cutoff %f\n" % job['cutoff']
        if (job.has_key('mult')):   ss += "  mult %d\n" % job['mult']
        if (job.has_key('xc')):
            if (job['xc'] == 'lda'):  ss += "  xc vosko\n"
            if (job['xc'] == 'pbe'):  ss += "  xc pbe96\n"
            if (job['xc'] == 'pbe0'): ss += "  xc pbe0\n"
        ss += "  lmbfgs\n" + "end\n"
        ss += "task pspw gradient\n"
    if (job['theory'] == 'dft'):
        if (job.has_key('basis')):
            ss += job['basis'] + "\n\n"
        else:
            ss += "basis \"ao basis\" cartesian print\n" + "  * library \"6-311++G(2d,2p)\"\n" + "end\n\n"
        ss += "dft\n"
        if (job.has_key('mult')): ss += "direct\n" + "  mult %d\n" % job['mult']
        if (job.has_key('xc')):
            if (job['xc'] == 'pbe'):   ss += "  xc xpbe96 cpbe96\n"
            if (job['xc'] == 'pbe0'):  ss += "  xc pbe0\n"
            if (job['xc'] == 'b3lyp'): ss += "  xc b3lyp\n"
        ss += "  iterations 5001\n" + "end\n"
        ss += "task dft gradient\n"
    if (job['theory'] == 'mp2'):
        if (job.has_key('basis')):
            ss += job['basis'] + "\n\n"
        else:
            ss += "basis \"ao basis\" cartesian print\n" + "  * library \"6-311++G(2d,2p)\"\n" + "end\n"
        if (job.has_key('mult')):
            if (job['mult'] == 2): ss += "scf\n" + "uhf\n" + "doublet\n" + "end\n"
            if (job['mult'] == 3): ss += "scf\n" + "uhf\n" + "triplet\n" + "end\n"
        ss += "mp2\n" + "end\n" + "set cphf:maxiter 1500\n"
        ss += "task mp2 gradient\n"
    if (job['theory'] == 'scf'):
        if (job.has_key('basis')):
            ss += job['basis'] + "\n\n"
        else:
            ss += "basis \"ao basis\" cartesian print\n" + "  * library \"6-311++G(2d,2p)\"\n" + "end\n"
        ss += "scf\n"
        if (job.has_key('mult')):
            if (job['mult'] == 2): ss += "scf\n" + "uhf\n" + "doublet\n" + "end\n"
            if (job['mult'] == 3): ss += "scf\n" + "uhf\n" + "triplet\n" + "end\n"
        ss += "end\n"
        ss += "task scf gradient\n"
    return ss


#############################################
#                                           #
#              parse_nwout                  #
#                                           #
#############################################
#
# The function returns (e, force) from parsing an nwchem output deck 
# of a pspw, dft, or mp2 gradient calculation
#
def parse_nwout(jobtheory, nwofilename):
    e = 0.0
    force = []
    if (jobtheory == 'pspw'):
        done = False
        ofound = False
        count = -1
        xyzdat = []
        edat = []
        ofile = open(nwofilename, 'r')
        for line in ofile:
            if (line.find("Total PSPW energy") != -1):
                edat.append(line)
            if (not done):
                if (count > 0):
                    if (line.find("C.O.M.") != -1):
                        done = True
                    else:
                        xyzdat.append(line)
                if (not done):
                    if (count >= 0):
                        count += 1
                    if (line.find("===  Ion Gradients ===") != -1):
                        ofound = True
                        count = 0
        ofile.close()
        split = edat[0].split()
        e = eval(split[4])
        for xyz in xyzdat:
            split = xyz.split()
            force.append(eval(split[3]))
            force.append(eval(split[4]))
            force.append(eval(split[5]))
    if (jobtheory == 'dft'):
        done = False
        ofound = False
        count = -1
        xyzdat = []
        edat = []
        ofile = open(nwofilename, 'r')
        for line in ofile:
            if (line.find("Total DFT energy") != -1):
                edat.append(line)
            if (not done):
                if (count > 2):
                    if (len(line) <= 5):
                        done = True
                    else:
                        xyzdat.append(line)
                if (not done):
                    if (count >= 0):
                        count += 1
                    if (line.find("DFT ENERGY GRADIENTS") != -1):
                        ofound = True
                        count = 0
        ofile.close()
        split = edat[0].split()
        e = eval(split[4])
        for xyz in xyzdat:
            split = xyz.split()
            force.append(-eval(split[5]))
            force.append(-eval(split[6]))
            force.append(-eval(split[7]))
    if (jobtheory == 'mp2'):
        done = False
        ofound = False
        count = -1
        xyzdat = []
        edat = []
        ofile = open(nwofilename, 'r')
        for line in ofile:
            if (line.find("Total MP2 energy") != -1):
                edat.append(line)
            if (not done):
                if (count > 2):
                    if (len(line) <= 5):
                        done = True
                    else:
                        xyzdat.append(line)
                if (not done):
                    if (count >= 0):
                        count += 1
                    if (line.find("mp2 ENERGY GRADIENTS") != -1):
                        ofound = True
                        count = 0
        ofile.close()
        split = edat[0].split()
        e = eval(split[3])
        for xyz in xyzdat:
            split = xyz.split()
            force.append(-eval(split[5]))
            force.append(-eval(split[6]))
            force.append(-eval(split[7]))
    if (jobtheory == 'scf'):
        done = False
        ofound = False
        count = -1
        xyzdat = []
        edat = []
        ofile = open(nwofilename, 'r')
        for line in ofile:
            if (line.find("Total SCF energy") != -1):
                edat.append(line)
            if (not done):
                if (count > 2):
                    if (len(line) <= 5):
                        done = True
                    else:
                        xyzdat.append(line)
                if (not done):
                    if (count >= 0):
                        count += 1
                    if (line.find("ENERGY GRADIENTS") != -1):
                        ofound = True
                        count = 0
        ofile.close()
        split = edat[0].split()
        e = eval(split[4])
        for xyz in xyzdat:
            split = xyz.split()
            force.append(-eval(split[5]))
            force.append(-eval(split[6]))
            force.append(-eval(split[7]))
    return (e, force)


#############################################
#                                           #
#              run_simulation               #
#                                           #
#############################################
#  This routine runs the specified simulation.  This routine should be modified
#  to the simulation your interested in running.  This example computes 
#    1) the lj  energy and force otoms if
#             x['theory']='lj'
#             x['xyz']=[x1,y1,z1, x2,y2,z2, ....],
#    2) the nwchem energy and force otoms at pspw, dft, scf, and mp2 if 
#             x['theory']='pspw','dft','scf','mp2'
#             x['xyz']=[x1,y1,z1, x2,y2,z2, ....],
#
#  The input and output be modified to be any python data types appropiate to the simulation
#
#  Entry - x : python data type,      e.g. dictionary
#  Exit - returns a python data type - e.g. tuple (e,forces) 
#
def run_simulation(PORT, x):
    if (x['theory'] == 'lj'):
        data = run_ljgradient(x)
    elif (x['theory'] == 'sw'):
        data = run_swgradient(x)
    elif (x['theory'] == 'gr'):
        data = run_grgradient(x)
    elif (x['theory'] == 'spring'):
        data = run_spring(x)
    else:
        ss = make_nwinput(PORT, x)
        curdir = os.getcwd()
        infile = curdir + '/wbwnwchem%d.nw' % PORT
        outfile = curdir + '/wbwnwchem%d.nwout' % PORT
        nwfile = open(infile, 'w')
        nwfile.write(ss)
        nwfile.close()
        cmd = nwchem_binary + infile + " > " + outfile
        os.system(cmd)
        data = parse_nwout(x['theory'], outfile)
    return data

#############################################
#                                           #
#              run_jobs                     #
#                                           #
#############################################
#  This routine runs a list of jobs. where each of the jobs 
#  is a dictionary specifying the computation to be done.
#
sstimes = [0] * 50


def run_jobs(onlyroute, PORT, hostsports, jobs):
    egrads = []
    njobs = len(jobs)
    if ((njobs == 1) and (not onlyroute)):
        egrads += [run_simulation(PORT, jobs[0])]
    else:
        nhp = len(hostsports)

        ## set up connections ###
        ss = []
        goodhostsports = []
        sscpus = []
        ncpus = 0
        for hp in hostsports[:njobs]:
            s = my_sockconnect(hp)
            if (s != None):
                ss.append(s)
                goodhostsports.append(hp)
                l = hp.index(':')
                r = hp.rindex(':')
                if (l == r):
                    cpus = 1
                else:
                    cpus = eval(hp[r + 1:])
                sscpus.append(cpus)
                ncpus += cpus

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
#              run_wbwserver                #
#                                           #
#############################################
#  This server routine listens for a socket requests from clients until a
#  client tells it to stop.  Once a connection is accepted 
#   1) it receives x
#   2) calls result = run_simulation(x) or stops if x='stop'
#   3) sends result back to the client
#   4) closes the connection
#
#  Entry - PORT, port number to listen to
#  Exit - none
#
def run_wbwserver(onlyroute, PORT, hostsports):
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    HOST = ''
    s.bind((HOST, PORT))
    while 1:
        s.listen(5)
        conn, addr = s.accept()
        print 'Connected to ', addr, ' on port ', PORT

        cpu1 = time.time()
        msg = ''
        done = False
        while not done:
            data = conn.recv(1024)
            msg = msg + data
            done = msg[len(msg) - 8:len(msg)] == 'ericdone'

        # unpickle the msg
        msg = msg[:len(msg) - 8]
        x = pickle.loads(msg)
        if (x == 'stop'): break

        # calculate result, e.g. energy and gradient
        cpu2 = time.time()
        result = run_jobs(onlyroute, PORT, hostsports, x)
        cpu3 = time.time()

        # send back result
        data = pickle.dumps(result)
        conn.sendall(data + 'ericdone')
        conn.close()
        cpu4 = time.time()
        print "recv,calc,send times=", cpu2 - cpu1, cpu3 - cpu2, cpu4 - cpu3, ' on port ', PORT

    conn.close()
    s.close()


#############################################
#                                           #
#              kill_wbwserver               #
#                                           #
#############################################
#  This routine acts as a client and sends "stop" to the PORT.

def kill_wbwserver(PORT):
    HOST = 'localhost'
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect((HOST, PORT))
    s.sendall(pickle.dumps('stop') + 'ericdone')
    s.close()


#############################################################################################
##################################### main program ##########################################
#############################################################################################
usage = \
    """
      The bwwserver_nwchem program runs and routes pspw, dft, scf, and mp2
      nwchem gradient jobs. The program is run as follows:

      bwwserver_nwchem start|route|stop [portnumber [hostport0 hostport1 hostport2 ....]]

      where hostport = 'ipaddr':port:ncpu:w

      e.g.
      bwwserver_nwchem start
      bwwserver_nwchem stop
      bwwserver_nwchem start 50001
      bwwserver_nwchem stop 50001
      bwwserver_nwchem start 50001 localhost:50002 localhost:50003:2 localhost:50004:3
      bwwserver_nwchem route 50001 localhost:50002 localhost:50003

    """
tt = time.localtime()
dd = "%d-%d-%d-%d:%d" % (tt[0], tt[1], tt[2], tt[3], tt[4])

PORT = 50001  #default port if none specified
hostsports = []
onlyroute = False

if (len(sys.argv) >= 2):

    if (len(sys.argv) == 2):
        runtype = sys.argv[1]
    elif (len(sys.argv) == 3):
        runtype = sys.argv[1]
        PORT = eval(sys.argv[2])
    else:
        runtype = sys.argv[1]
        PORT = eval(sys.argv[2])
        for i in range(3, len(sys.argv)):
            hostsports.append(sys.argv[i])

    if (runtype.find('stop') != -1):
        print 'wbwserver program stopped at ', dd, ' on port ', PORT
        kill_wbwserver(PORT)
    else:
        if (runtype.find('route') != -1): onlyroute = True
        print 'wbwserver program started at ', dd, ' on port ', PORT
        if onlyroute:
            print '              - port ', PORT, ' routing to ', hostsports
        else:
            print '              - port ', PORT, ' using ', hostsports
        run_wbwserver(onlyroute, PORT, hostsports)
        tt = time.localtime()
        dd = "%d-%d-%d-%d:%d" % (tt[0], tt[1], tt[2], tt[3], tt[4])
        print 'wbwserver program ended at ', dd, ' on port ', PORT
else:
    print usage

