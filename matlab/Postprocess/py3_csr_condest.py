#!/usr/bin/python3
import numpy as np
import h5py
from scipy.sparse import csr_matrix
from pyamg.util.linalg import condest, cond
import time




"""
Script to compute the condition number estimate of the matrix.

Please use before the matlab routine (HDF5readCSRmat) with the option save=1 to compute the mat_ijk

In:
nproc       : (int) number of parallel processes
pathName    : (string) path to the Mat_ijk
niter       : (int) number of iteration for condest()
sym         : (bool) is the matrix symmetric?
tol         : (float) tolerance for while loop
err         : (float) set > tol for while loop (very long runtime)

n: number of (locally own) columns in the matrix
nnz: number of (locally own) non-zeros elements
cols: column of indices for each elements
rowptr: indices such that columns indices for row i are cols[rowptr[i]:rowptr[i+1]] with corresponding values in val[rowptr[i]:rowptr[i+1]]
loc2glob: global column number of the local column in cols
vals: (local) values

"""
nproc = 16
#pathName = '/home/bluce/MHDG_sim/West1090_P4/2021_07_01_NGT2D/'
#pathName = '/marconi_scratch/userexternal/bluce000/MHDG/West1090_P4/2021_07_01_NGT2D/'
pathName = '/scratch/bluce/MHDG/West1090_P4/2021_07_01_NGT2D/'
sfile = '2D.txt'
niter = 25
sym = False
tol = 0.5
err = tol*10


def display_time(seconds, granularity=5):
    """
    Make readable time

    Example(s)
    -------

    Parameter(s)
    ----------
    seconds: float
    time in seconds

    Parameters (optional)
    ----------
    granularity: int
    smallest unit.
    1 is for day and 5 is for second

    Return
    -------
    return: type(return)
    A short description

    Author(s)
    -------
    Benjamin LUCE (benjamin.luce[at]univ-amu.fr)
    """
    result = []
    time_intervals = (
        ('weeks', 604800),  # 60 * 60 * 24 * 7
        ('days', 86400),    # 60 * 60 * 24
        ('hours', 3600),    # 60 * 60
        ('minutes', 60),
        ('seconds', 1))

    for name, count in time_intervals:
        value = seconds // count
        if value:
            seconds -= value * count
            if value == 1:
                name = name.rstrip('s')
            result.append("{} {}".format(value, name))
    return ', '.join(result[:granularity])


I,J,K = np.array([]),np.array([]),np.array([])
for i in range(nproc):
    if nproc == 1:
        fileName ='Mat_ijk.h5'
    else:
        fileName ='Mat_ijk_'+str(i+1)+'_'+str(nproc)+'.h5'
        try:
            with h5py.File(pathName+fileName, 'r') as fname:
                Iloc = fname['I'][()]
                Jloc = fname['J'][()]
                Kloc = fname['K'][()]
                I = np.append(I,Iloc)
                J = np.append(J,Jloc)
                K = np.append(K,Kloc)
        except OSError or FileNotFoundError:
            raise('File not found in {0}'.format(pathName+fileName))
        # This commented part is too long in python due to the loop
        #            with h5py.File(pathName+fileName) as fname:
        #                n        = fname['n'][0]        
        #                nnz      = fname['nnz'][0]      
        #                cols     = fname['cols'][()]    
        #                rowptr   = fname['rowptr'][()]  
        #                loc2glob = fname['loc2glob'][()]
        #                vals     = fname['vals'][()]
        #                Iloc, Jloc, Kloc = np.zeros(nnz),np.zeros(nnz),np.zeros(nnz)
        #                for i in range(np.size(rowptr)-1):
        #                    Iloc[rowptr[i]-1:rowptr[i+1]] = loc2glob[i]
        #                    Jloc[rowptr[i]-1:rowptr[i+1]] = cols[rowptr[i]-1:rowptr[i+1]]
        #                    Kloc[rowptr[i]-1:rowptr[i+1]] = vals[rowptr[i]-1:rowptr[i+1]]
        #                I = np.append(I,Iloc)
        #                J = np.append(J,Jloc)
        #                K = np.append(K,Kloc)
        I = I.astype(int)
        J = J.astype(int)
        matrix = csr_matrix((K,(I,J)))

while err > tol:
    t1 = time.time()
    a = condest(matrix,maxiter=niter,symmetric=sym)
    b = condest(matrix,maxiter=2*niter,symmetric=sym)
    err = abs((b-a)/a)
    niter = 2*niter
    t2 = time.time()
    print('Condition number estimate loop (time={0}, err={1}, niter={2}, sym={3}): {4:e}'.format(display_time(t2-t1),err,niter,sym,b))
    f = open(sfile, 'a')
    f.write('Condition number estimate loop (time={0}, err={1}, niter={2}, sym={3}): {4:e}\n'.format(display_time(t2-t1),err,niter,sym,b))
    f.close()
    if niter > 25000:
        err = tol/10
else:
    t1 = time.time()
    b = condest(matrix,maxiter=niter,symmetric=sym)
    t2 = time.time()
    print('Condition number estimate no loop (time={0}, err={1}, niter={2}, sym={3}): {4:e}'.format(display_time(t2-t1),err,niter,sym,b))
    f = open(sfile, 'a')
    f.write('Condition number estimate loop (time={0}, err={1}, niter={2}, sym={3}): {4:e}\n'.format(display_time(t2-t1),err,niter,sym,b))
    f.close()

# cond(matrix) ONLY for small matrix
# print('Condition number: {0}'.format(cond(matrix)))
