#!/usr/bin/env python
import os,re,glob,pyfits,gc
import numpy as np

zmax=1.0
zmin=0.0
nbins=100
zbins=[i*(zmax-zmin)/nbins+zmin for i in range(0,nbins+1)]
match_string='*.fit'
dir='./'
outfile_name='z_distr_redm'
files=glob.glob(os.path.join(dir,match_string))

zn=[0]*nbins
zave=np.zeros(nbins)
for i,file in enumerate(files):
    print 'Reading '+file
    try:
        pyfile=pyfits.open(file)
        data=pyfile[1].data['ZREDMAGIC']
    except:
        print "Can't read file "+file+" skipping"
        continue

    # count the entries in each bin
    nbin=np.histogram(data,bins=zbins)[0]
    # sum the number of entries in each bin
    ave=np.histogram(data,bins=zbins,weights=data)[0]
    
    ave=[a/b if b>0 else 0 for a,b in zip(ave,nbin)]

    # average the total and the current histograms together
    for j,(cbin,cave,tbin,tave) in enumerate(zip(nbin,ave,zn,zave)):
        tot=cbin+tbin
        zave[j]=cave*cbin/tot+tave*tbin/tot

    zn+=nbin

    pyfile.close()

    # This seems to be needed for older versions of pyfits
    gc.collect()
norm = sum(zn*1.0)*(zbins[1]-zbins[0])
outfile=open(outfile_name,'w')
for i in range(0,nbins):
    print '%e %e\n%e %e'%(zbins[i],zn[i],zbins[i+1],zn[i])
    outfile.write('%e %e %e %e\n'%(zbins[i],(zbins[i]+zbins[i+1])/2.,zbins[i+1],zn[i]*1.0/norm))
outfile.write("\n")
outfile.close()
