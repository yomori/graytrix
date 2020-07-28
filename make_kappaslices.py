#This file generates kappa files for each redshift
# essentially Equation C2 of https://arxiv.org/pdf/1504.05672.pdf

#from cosmotools import utils
#from cosmotools import theory
#from cosmotools import cosmo
import numpy as np
import healpy as hp
import camb
from camb import model, initialpower
import sys

h  = 0.6777

dir_in   = sys.argv[1]
dir_out  = sys.argv[2]
cambpar  = sys.argv[3]#'mdpl2_params.ini'

pars     = camb.read_ini(cambpar)
results  = camb.get_results(pars)
 
#---------------------------------------------------------------------------
Nsnap    = 146 #highest slicenumber-4+1
dchish   = 25

chish   = (np.arange(Nsnap+1)+0.5)*dchish
chishuu = (np.arange(Nsnap+1)+1.0)*dchish
chishll = (np.arange(Nsnap+1)+0.0)*dchish

chis    = chish/h       # mid-z in Mpc
chisuu  = chishuu/h     # upper edge
chisll  = chishll/h     # lower edges
dchis   = dchish/h      # tmpchis[1:]-tmpchis[:-1]

midzs   = results.redshift_at_comoving_radial_distance(chis)
uuzs    = results.redshift_at_comoving_radial_distance(chisuu)
llzs    = results.redshift_at_comoving_radial_distance(chisll)

for i in range(4,Nsnap):
    print(i,llzs[i],midzs[i],uuzs[i])


tmp=np.zeros(hp.nside2npix(8192))
ll=np.arange(6001)



for i in range(0,4):
    print("adding slice %d zmid=%.3f zuedge=%.3f"%(i,midzs[i],uuzs[i]))
    print(i)
    print(1./(1+midzs[i]))
    d=hp.read_map(dir_in+'mdpl2_density_%d_%d_%d.fits'%(i,i*25,(i+1)*25))
    avg = np.mean(d)
    o   = (d-avg)

    omm      = 0.307115 # Omega_matter
    boxlen   = 1000.    # box length of sim in [Mpc/h]
    npart    = 3840.    # number of particles**(1/3)
    chubble0 = 2.9979e3 # 
    nside    = 8192 

    # note: comoving radial distance is the SAME as comoving angular diameter distance in a flat universe
    # Don't confuse with the "standard" angular diamter distance, which is dA=chi/(1+z)
    # A good reference http://www.tapir.caltech.edu/~chirata/ph217/lec02.pdf 
    print (chisuu[i],chisll[i])
    dlplane  = 0.75*((chisuu[i]*h)**4-(chisll[i]*h)**4)/((chisuu[i]*h)**3-(chisll[i]*h)**3) # comoving radial distance NOT comoving angular diameter distance
    zlens    = results.redshift_at_comoving_radial_distance(dlplane/h)
    alens    = 1./(zlens+1.0)
    
    norm0    = 1.5*omm*boxlen**3*hp.nside2npix(nside)/(alens*(dlplane)*chubble0**2*npart**3*4.0*np.pi)
    print("alens %f"%alens)
    print("chiuu %f"%(chisuu[i]*h))
    print("chill %f"%(chisll[i]*h)) 
    print("dlplane %f"%dlplane)
    print(norm0)
    hp.write_map(dir_out+'mdpl2_graytrix_sdens_%d.fits'%(i),norm0*o,overwrite=True)

