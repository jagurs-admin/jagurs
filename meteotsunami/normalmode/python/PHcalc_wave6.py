JAGURS=1
NONAN_NOFILT=0

import sys
import numpy as np
import pandas
from scipy import signal
from scipy.special import hankel2
from scipy.io import FortranFile
import json

def sub2ind(shape, rows, cols):
    return (cols*shape[1] + rows).astype(np.int64)

f = open('normalmode.json','r')
json_dict = json.load(f)
srclon = json_dict['srclon']
srclat = json_dict['srclat']
dist_s = json_dict['dist_s']
dist_e = json_dict['dist_e']
tmax   = json_dict['tmax']
dt     = json_dict['dt']
if 'sind' in json_dict:
    sind = json_dict['sind']
else:
    sind = 10
if 'SDt' in json_dict:
    SDt = json_dict['SDt']
else:
    SDt = 2000 # source duration (s)
if 'DCflag' in json_dict:
    DCflag = json_dict['DCflag']
else:
    DCflag = 1     # 1: Calculate dispersion curve
                   # 0: Use a priori dispersion curve save as DCfile
if 'SYNflag' in json_dict:
    SYNflag = json_dict['SYNflag']
else:
    SYNflag = 1    # 1: Calculate synthetic waveform
                   # 0: Not calculate
if 'OH' in json_dict:
    OH = json_dict['OH']
else:
    OH = 5 # depth (km)
if 'atmfile' in json_dict:
    atmfile = json_dict['atmfile']
else:
    atmfile = 'AtmModel.dat'

pi = 3.14159265

if JAGURS == 0:
    maxmode = 13 # maximum mode number for dispersion curve and synthetic wave calculation
    minmode = 8  # minimum mode for synthetic wave calculation
else:
    if 'maxmode' in json_dict:
        maxmode = json_dict['maxmode']
    else:
        maxmode = 13
    if 'minmode' in json_dict:
        minmode = json_dict['minmode']
    else:
        minmode = 10

with_Ocean = 1 # 1: Calculate with ocean layer
               # 0: Without ocean
Only_Ocean = 0 # 1: Cacluate only ocean layer
               # 0: Calculate with atmospheric model
DCsaveflag = 1 # 1: Save dispersion curve in DCfile
               # 0: Not save
DCplotflag = 0 # 1: Plot dispersion curve (phase velocity)
               # 0: Not plot

DCfile = 'DC_PH5_' + str(minmode) + '-' + str(maxmode) + '.npz' # dispersion curve
DCtext = 'DC_PH5_' + str(minmode) + '-' + str(maxmode) + '.txt' # dispersion curve

distlist = np.arange(dist_s,dist_e+1) # [km]

outfile = 'atm_NM_US_' + str(dist_s) + '-' + str(dist_e) + 'km_' + str(tmax) + 'sec_PH5_GR'

eps = np.finfo(np.float64).eps
t = np.arange(0,tmax+eps*tmax,dt) # time of synthetic wave (s)

synwave1D = np.zeros((t.shape[0],distlist.shape[0]))
synwave1D2 = np.zeros((t.shape[0],distlist.shape[0]))
synwave1D3 = np.zeros((t.shape[0],distlist.shape[0]))
synwave1D4 = np.zeros((t.shape[0],distlist.shape[0]))
synwave1D5 = np.zeros((t.shape[0],distlist.shape[0]))

# -- read atmospheric model --
AtmModel = pandas.read_csv(atmfile,header=None)
zl = AtmModel[0].values
H = AtmModel[1].values
alpha = AtmModel[2].values
rho = AtmModel[3].values

# -- parameters of water layer --
Oalpha = 1.5 # speed of sound (km/s)
Orho = 1030 # density (kg/m3)

# -- renew model with ocean --
g0 = 9.8e-3 # gravitational acceleration (km/s)
if with_Ocean == 1:
    alpha = np.append(Oalpha, alpha)
    rho = np.append(Orho, rho)
    H = np.append(OH, H)
    zl = np.append(0, OH+zl)

nz = alpha.shape[0] # number of layers
im = 1.j # imaginary unit

if Only_Ocean == 1:
    maxmode = 1
    minmode = 1
    nz = 1

# -- phase velocity (km/s)--
vp = np.arange(0.1,0.4+eps*0.4/10,0.00001)

# -- angular frequency (rad/s) --
Df = 2*pi/(t.shape[0]*dt)
omega = np.arange(Df,0.05+eps*0.05,Df)

fs = 1/dt
fc = 1/5000

b_filt,a_filt = signal.butter(2, fc/(fs/2), 'high')

# calculation parameters
Re = 6371.0 # earth radius (km)
g = g0*(Re/(Re+zl+H/2))**2
ganma = np.ones(alpha.shape[0])
# -- sin source --
if 'Samp' in json_dict:
    Samp = json_dict['Samp']
else:
    Samp = 1000 # source amplitude

if with_Ocean == 1:
    ganma[1:] = 1.4 # specific heat ratio (1 at ocean, 1.4 at atm)
    BVf = g*np.sqrt(ganma-1)/alpha # Brund-Vaisala
    BVf[0] = 0.003 # Brunt-Vaisala frequency at ocean
else:
    ganma[:] = 1.4
    BVf = g*np.sqrt(ganma-1)/alpha
    sind = 1

if Only_Ocean == 1:
    sind = 1

lamda = ganma*g/(2*alpha**2)
beta = 2*alpha**2*BVf/(g*ganma)

# Calculate dispersion curve
if DCflag == 1:
    print('Calculate dispersion curve')
    
    # -- Set wavenumber --
    VPmat,OMGmat = np.meshgrid(vp,omega)
    Kmat = OMGmat/VPmat
    nk = Kmat.size
    Kvec = Kmat.reshape(-1)
    OMGvec = OMGmat.reshape(-1)
    VPvec = VPmat.reshape(-1)
    
    Kvec2 = Kvec**2
    OMGvec4 = OMGvec**4
    VPvec2 = VPvec**2
    kvp3 = Kvec**3*VPvec**3
    
    # -- Calculate matrix am (PH62-Eq.12) --
    # Press and Harkrider(1962)
    am = np.zeros((nk,nz,2,2),dtype=np.complex128)
    
    for l in range(nz):
        delta = g[l]**2*Kvec2 - OMGvec4
        
        tmp1 = VPvec2/alpha[l]**2
        tmp2 = VPvec2/beta[l]**2
        kr = -im*np.sqrt(Kvec2*(1 - tmp1) - BVf[l]**2/VPvec2*(1 - tmp2) + 0.j)
        
        Pm = kr*H[l]
        
        tmp3 = g[l]/alpha[l]**2
        tmp4 = alpha[l]**2/VPvec2 - ganma[l]/2
        tmp5 = np.sin(Pm)/kr
        am[:,l,0,0] = np.exp(lamda[l]*H[l])*(np.cos(Pm) + tmp3*tmp4*tmp5)
        am[:,l,1,1] = np.exp(-lamda[l]*H[l])*(np.cos(Pm) - tmp3*tmp4*tmp5)
        
        tmp6 = g[l]**2*tmp4**2
        tmp7 = alpha[l]**4*kr**2
        tmp8 = rho[l]*alpha[l]**4*delta*kr
        am[:,l,0,1] = im*kvp3*(tmp6 + tmp7)*np.sin(Pm)/tmp8
        am[:,l,1,0] = im*rho[l]*delta*np.sin(Pm)/(kvp3*kr)
    
    # -- Calculate normal mode dispersion function A22 (PH62-Eq.16) --
    A = am[:,0,:,:]
    
    for l in range(1,nz):
        A = am[:,l,:,:]@A
    
    if NONAN_NOFILT == 0:
        for l in range(nk):
            Acheck = np.linalg.det(A[l,:,:])
            if np.abs(1 - Acheck) > 0.1:
                A[l,:,:] = np.nan
    
    A22 = A[:,1,1].real
    
    # -- Dispersion curve estimation --
    row = np.zeros(omega.shape[0])
    col = np.arange(omega.shape[0])
    sz = Kmat.shape
    
    DC = np.zeros((3,maxmode,vp.shape[0]))
    DC[:,:,:] = np.nan
    # DC(:,:,1) -> phase speed for each mode
    # DC(:,:,2) -> wavenumber for each mode
    # DC(:,:,3) -> angular frequency for each mode
    DCind = np.zeros(maxmode, dtype=np.int64)
    
    nmode = 0 # mode number
    for l in range(vp.shape[0]):
        row[:] = l
        
        kind = sub2ind(sz,row,col) # search k for each vp
        ktmp = Kvec[kind]
        omgtmp = OMGvec[kind]
        ksort = np.sort(ktmp) # sort k in ascending order
        ksi = np.argsort(ktmp) # sort k in ascending order
        
        A22tmp = np.sign(A22[kind])
        A22sort = A22tmp[ksi]
        omgsort = omgtmp[ksi]
        
        if l == 0:
            init_sign = np.sign(A22sort[0]).astype(np.int64) # set initial sign of each vp
            old_sign = init_sign
        else:
            init_sign = np.sign(A22sort[0]).astype(np.int64)
            if init_sign != old_sign: # update initial mode
                nmode = nmode + 1
                old_sign = init_sign
        
        zcind = np.where(A22sort[:-1]*A22sort[1:] < 0)[0] # zero-cross index;
        
        for ll in range(zcind.shape[0]):
            if nmode+ll+1 > maxmode:
                break
            DC[0,nmode+ll,DCind[nmode+ll]] = vp[l]
            DC[1,nmode+ll,DCind[nmode+ll]] = ksort[zcind[ll]]
            DC[2,nmode+ll,DCind[nmode+ll]] = omgsort[zcind[ll]]
            DCind[nmode+ll] = DCind[nmode+ll] + 1
    
    # search for each omega (<0.01)
    omg_searchmax = 0.01
    row = np.arange(vp.shape[0])
    col = np.zeros(vp.shape[0])
    omglist = np.where(omega <= omg_searchmax)[0]
    DComg = np.zeros((3,maxmode,omglist.shape[0]))
    DComg[:,:,:] = np.nan
    DCind = np.zeros(maxmode, dtype=np.int64)
    
    nmode = 0
    for l in range(omglist.shape[0]):
        col[:] = l
        kind = sub2ind(sz,row,col)
        vptmp = VPvec[kind]
        ksort = np.sort(Kvec[kind])
        ksi = np.argsort(Kvec[kind])
        A22tmp = (np.sign(A22[kind]).real).astype(np.int64)
        A22sort = A22tmp[ksi]
        vpsort = vptmp[ksi]
        if l == 0:
            init_sign = np.sign(A22sort[0]).astype(np.int64) # set initial sign of each vp
            old_sign = init_sign
        else:
            init_sign = np.sign(A22sort[0]).astype(np.int64)
            if init_sign != old_sign: # update initial mode
                nmode = nmode + 1
                old_sign = init_sign
        
        zcind = np.where(A22sort[:-1]*A22sort[1:] < 0)[0] # zero-cross index;
        for ll in range(zcind.shape[0]):
            if nmode+ll+1 > maxmode:
                break
            DComg[0,nmode+ll,DCind[nmode+ll]] = vpsort[zcind[ll]]
            DComg[1,nmode+ll,DCind[nmode+ll]] = ksort[zcind[ll]]
            DComg[2,nmode+ll,DCind[nmode+ll]] = omega[l]
            DCind[nmode+ll] = DCind[nmode+ll] + 1
    
    DCall = np.zeros(DC.shape)
    DCall[:,:,:] = np.nan
    DComg = DComg[:,::-1,:]
    for l in range(minmode-1,maxmode):
        DCmin = np.nanmin(DC[2,l,:])
        DCallind1 = np.nanargmin(DC[2,l,:])
        DCallind2 = np.where(DComg[2,l,:] == DCmin)[0]
        if DCallind2.shape[0] == 0:
            DCallind2 = 0
        else:
            DCallind2 = DCallind2[0]
        if DCallind2 == 0:
            DCall[:,l,:DCallind1+1] = DC[:,l,:DCallind1+1]
        else:
            DCall[:,l,:DCallind1+1] = DC[:,l,:DCallind1+1]
            DCall[:,l,DCallind1+1:DCallind1+DCallind2+2] = DComg[:,l,DCallind2::-1]
    
    # make & save dispersion curve
    nkr = np.zeros(maxmode-minmode+1,dtype=np.int64) # number of k root in each mode
    ll = 0
    for l in range(minmode-1,maxmode):
        nkr[ll] = np.count_nonzero(~np.isnan(DCall[0,l,:]))
        ll = ll + 1;
    
    kroot = np.zeros(np.sum(nkr)) # k root
    vpkr = np.zeros(kroot.shape[0]) # phase velocity at k root
    omgkr = np.zeros(kroot.shape[0]) # angular freq. at k root
    tmp = DCall[1,minmode-1,:]
    kroot[:nkr[0]] = tmp[~np.isnan(tmp)]
    tmp = DCall[0,minmode-1,:]
    vpkr[:nkr[0]] = tmp[~np.isnan(tmp)]
    tmp = DCall[2,minmode-1,:]
    omgkr[:nkr[0]] = tmp[~np.isnan(tmp)]
    ll = 1
    for l in range(minmode,maxmode):
        ind1 = np.sum(nkr[:ll])
        ind2 = np.sum(nkr[:ll])+nkr[ll]
        tmp = DCall[1,l,:]
        kroot[ind1:ind2] = tmp[~np.isnan(tmp)]
        tmp = DCall[0,l,:]
        vpkr[ind1:ind2] = tmp[~np.isnan(tmp)]
        tmp = DCall[2,l,:]
        omgkr[ind1:ind2] = tmp[~np.isnan(tmp)]
        ll = ll + 1
     
    if DCsaveflag == 1:
        print('Save dispersion curve')
        np.savez(DCfile, nkr=nkr, kroot=kroot, vpkr=vpkr, omgkr=omgkr) # save dispersion curve and matrix A
        
        f = open(DCtext, 'w')
        f.write('----------------\n')
        f.write('mode    nkr     \n')
        f.write('----------------\n')
        for l in range(minmode,maxmode+1):
            f.write('%8d%8d\n'%(l,nkr[l-minmode]))
        
        f.write('-------------------------------------------------------------\n')
        f.write('mode    nkr     kroot          vpkr           omgkr          \n')
        f.write('-------------------------------------------------------------\n')
        ll = 0
        for l in range(minmode,maxmode+1):
            for i in range(nkr[l-minmode]):
                f.write('%8d%8d%15.8f%15.8f%15.8f\n'%(l,i,kroot[ll],vpkr[ll],omgkr[ll]))
                ll = ll + 1
        f.close()
    
else:
    print('Load dispersion curve')
    npz = np.load(DCfile)
    nkr = npz['nkr']
    kroot = npz['kroot']
    vpkr = npz['vpkr']
    omgkr = npz['omgkr']

# Calculate synthetic waveforms
if SYNflag == 1:
    print('Calcualte synthetic waveforms')
    
    flag_zero = False
    jj = 0
    for dist in distlist:
        if dist == 0:
            flag_zero = True
            jj = jj + 1
            continue
        
        baro = np.zeros((maxmode-minmode+1,t.shape[0])) # synthetic barogram
        
        # -- parameters for fourier transform --
        Nf = 2*pi/(2*dt)
        FFTf = np.arange(0,Nf*2-Df+eps*(Nf*2-Df),Df)
        
        # -----------------------
        #     figure % Figure for check
        #     hold on
        # -----------------------
        
        for ll in range(nkr.shape[0]):
            if nkr[ll] == 0:
                break
            
            FFTbaro = np.zeros(FFTf.shape[0],dtype=np.complex128)
            
            if ll == 0:
                omgkrUni,indUni = np.unique(omgkr[0:nkr[0]],return_index=True) # omega for first mode
            else:
                ind1 = np.sum(nkr[:ll])
                ind2 = np.sum(nkr[:ll])+nkr[ll]
                omgkrUni,indUni = np.unique(omgkr[ind1:ind2],return_index=True)
                indUni = ind1+indUni
            
            indtmp1 = np.where(FFTf-omega[0] >= 0)[0]
            Find1 = indtmp1[0]
            indtmp2 = np.where(0.0124-FFTf >= 0)[0] # longer than 506 sec
            Find2 = indtmp2[-1]
            
            # -- interpolate each variable for FFT
            kint = np.interp(FFTf[Find1:Find2+1], omgkrUni, kroot[indUni])
            vpint = np.interp(FFTf[Find1:Find2+1], omgkrUni, vpkr[indUni])
            omgint = FFTf[Find1:Find2+1]
            
            # -- Calculate d(am)/dk (H64-Eq.69)--
            dam = np.zeros((kint.shape[0],nz,2,2),dtype=np.complex128) # d(am)/d(k)
            damo = np.zeros((kint.shape[0],nz,2,2),dtype=np.complex128) # d(am)/d(omega)
            amkr = np.zeros(dam.shape,dtype=np.complex128) # am for k root
            
            kroot2 = kint**2
            kroot3 = kint**3
            vpkr2 = vpint**2
            vpkr3 = vpint**3
            omgkr2 = omgint**2
            omgkr3 = omgint**3
            omgkr4 = omgint**4
            kvpkr3 = kroot3*vpkr3
            
            for l in range(nz):
                delta = g[l]**2*kroot2 - omgkr4
                ddelk = 2*g[l]**2*kint # d(delta)/dk
                ddelo = -4*omgkr3 # d(delta)/d(omega)
                
                dtmp1 = kroot2 - omgkr2/alpha[l]**2 - BVf[l]**2*kroot2/omgkr2 + BVf[l]**2/beta[l]**2
                kr = -im*np.sqrt(dtmp1 + 0.j)
                dkrk = -im*(1 - BVf[l]**2/omgkr2)*kint/np.sqrt(dtmp1 + 0.j) # d(kr)/dk
                dkro = -im*(BVf[l]**2*kroot2/omgkr3 - omgint/alpha[l]**2)/np.sqrt(dtmp1 + 0.j) # d(kr)/d(omgea)
                
                Pm = kr * H[l]
                dPmk = dkrk*H[l] # d(Pm)/dk
                dPmo = dkro*H[l] # d(Pm)/d(omega)
                
                dPkrk = (dPmk*np.cos(Pm)*kr - dkrk*np.sin(Pm))/kr**2 # d(sin(Pm)/kr)/dk
                dPkro = (dPmo*np.cos(Pm)*kr - dkro*np.sin(Pm))/kr**2 # d(sin(Pm)/kr)/d(omega)
                
                # -- calculate each element of d(am)/dk --
                dtmp2 = 2*g[l]*np.sin(Pm)/kr
                dtmp3 = g[l]/alpha[l]**2*(alpha[l]**2*kroot2/omgkr2 - ganma[l]/2)
                dam[:,l,0,0]  = np.exp( lamda[l]*H[l])*(-dPmk*np.sin(Pm) + dtmp2*kint/omgkr2 + dtmp3*dPkrk)
                dam[:,l,1,1]  = np.exp(-lamda[l]*H[l])*(-dPmk*np.sin(Pm) - dtmp2*kint/omgkr2 - dtmp3*dPkrk)
                damo[:,l,0,0] = np.exp( lamda[l]*H[l])*(-dPmo*np.sin(Pm) - dtmp2*kroot2/omgkr3 + dtmp3*dPkro)
                damo[:,l,1,1] = np.exp(-lamda[l]*H[l])*(-dPmo*np.sin(Pm) + dtmp2*kroot2/omgkr3 - dtmp3*dPkro)
                
                dam[:,l,1,0] = im*rho[l]/omgkr3*(ddelk*np.sin(Pm)/kr + delta*dPkrk)
                
                dtmp4o = (ddelo*omgkr3 - 3*delta*omgkr2)*np.sin(Pm)/(omgint**6*kr)
                damo[:,l,1,0] = im*rho[l]*(dtmp4o + delta*dPkro/omgkr3)
                
                dtmp5 = alpha[l]**2*kroot2/omgkr2 - ganma[l]/2
                dtmp6 = (g[l]**2*dtmp5**2 + alpha[l]**4*kr**2)/delta
                dtmp7k = 4*alpha[l]**2*g[l]**2*dtmp5*kint/omgkr2 + 2*alpha[l]**4*kr*dkrk
                dtmp8k = (-ddelk*dtmp6 + dtmp7k)/delta
                dam[:,l,0,1] = im*omgkr3/(rho[l]*alpha[l]**4)*(dtmp8k*np.sin(Pm)/kr + dtmp6*dPkrk)
                
                dtmp7o = -4*alpha[l]**2*g[l]**2*dtmp5*kroot2/omgkr3 + 2*alpha[l]**4*kr*dkro
                dtmp8o = (-ddelo*dtmp6 + dtmp7o)/delta
                dtmp9o = dtmp8o*np.sin(Pm)/kr + dtmp6*dPkro
                damo[:,l,0,1] = im/(rho[l]*alpha[l]**4)*(3*omgkr2*dtmp6*np.sin(Pm)/kr + omgkr3*dtmp9o)
                
                # -- calculate each element of am for k root --
                tmp3 = g[l]/alpha[l]**2
                tmp4 = alpha[l]**2/vpkr2 - ganma[l]/2
                tmp5 = np.sin(Pm)/kr
                amkr[:,l,0,0] = np.exp(lamda[l]*H[l])*(np.cos(Pm) + tmp3*tmp4*tmp5)
                amkr[:,l,1,1] = np.exp(-lamda[l]*H[l])*(np.cos(Pm) - tmp3*tmp4*tmp5)
                
                tmp6 = g[l]**2*tmp4**2
                tmp7 = alpha[l]**4*kr**2
                tmp8 = rho[l]*alpha[l]**4*delta*kr
                amkr[:,l,0,1] = im*kvpkr3*(tmp6+tmp7)*np.sin(Pm)/tmp8
                amkr[:,l,1,0] = im*rho[l]*delta*np.sin(Pm)/(kvpkr3*kr)
            
            # -- Calculate dA/dk (H64-Eq.70) --
            dA = dam[:,0,:,:] + amkr[:,0,:,:] # dA/dk
            dAo = damo[:,0,:,:] + amkr[:,0,:,:] # dA/d(omega)
            Akr = amkr[:,0,:,:] # A for k root
            for l in range(1,nz):
                dA = dam[:,l,:,:]@Akr + amkr[:,l,:,:]@dA
                dAo = damo[:,l,:,:]@Akr + amkr[:,l,:,:]@dAo
                Akr = amkr[:,l,:,:]@Akr
            
            # -- Calculate group velocity U (H64-Eq.71)--
            U = -dA[:,1,1]/dAo[:,1,1]
            
            # -- Calculate synthetic wave  --
            Lt = Samp*np.sin(2*pi*t/SDt) # sinusoidal source in time domain
            Lt[int(np.round(SDt/dt)):] = 0
            Lf = np.fft.fft(Lt) # sinusoidal source in frequency domain
            L = Lf[Find1:Find2+1]
            
            As1 = amkr[:,0,:,:]
            if sind > 2:
                for l in range(1,sind-1):
                    As1 = amkr[:,l,:,:]@As1
            syntmp1 = (-Akr[:,1,0]).conjugate()* \
                ((As1[:,1,1]).conjugate() - im*g0*rho[sind-1]*(As1[:,0,1]).conjugate()/omgint) \
                /(dA[:,1,1]).conjugate() # propagation term (H64-Eq.51)
            
            theta = dist/Re # station distance in angle [rad]
            if dist != 0:
                syntmp2 = np.sqrt(dist/(Re*np.sin(theta)))*hankel2(0,dist*kint)
            else:
                syntmp2 = np.sqrt(1/Re)*hankel2(0,dist*kint)
            
            P2W = -im*amkr[:,0,0,1] # vetical velocity (H64-Eq.56)
            P2Ppm = amkr[:,0,1,1] #  pressure at upper interface
            P2Pm = P2Ppm + rho[0]*g0*P2W/omgint # upper interface + rho*g*eta
            
            # -- sin wave source (HP67) --
            FFTbaro[Find1:Find2+1] = P2Ppm.conjugate()*L*syntmp1*syntmp2 # (HP67-Eq.14)
            
            ind = int((FFTf.shape[0]-1)/2)+1
            FFTbaro[ind:] = np.conjugate(FFTbaro[ind-1:0:-1])
            
            # -- apply group velocity taper --
            tind1 = np.argmin(np.abs(t-dist/(np.max(U)*1.2)))
            tind2 = np.argmin(np.abs(t-dist/(np.min(U)*0.8)))
            taper = np.zeros(t.shape[0])
            if tind2 < t.shape[0]-1:
                tptmp = signal.windows.tukey(t[tind1:tind2+1].shape[0],0.1)
                taper[tind1:tind2+1] = tptmp
            else:
                tptmp = signal.windows.tukey(t[tind1:].shape[0],0.1)
                taper[tind1:] = tptmp
            
            if NONAN_NOFILT == 0:
                #baro[ll,:]= signal.filtfilt(b_filt,a_filt,taper*np.fft.ifft(FFTbaro).real)
                baro[ll,:]= signal.filtfilt(b_filt,a_filt,np.fft.ifft(FFTbaro).real)
            else:
                #baro[ll,:] = taper*np.fft.ifft(FFTbaro).real
                baro[ll,:] = np.fft.ifft(FFTbaro).real
        
        if JAGURS == 0:
            synwave1D[:,jj] = np.sum(baro[2:-1,:],axis=0)
            synwave1D2[:,jj] = baro[2,:]
            synwave1D3[:,jj] = baro[3,:]
            synwave1D4[:,jj] = baro[4,:]
            synwave1D5[:,jj] = baro[5,:]
        else:
            synwave1D[:,jj] = np.sum(baro[:,:],axis=0)
        jj = jj + 1
    
    if flag_zero == True:
        if JAGURS == 0:
            synwave1D[:,0] = synwave1D[:,1]
            synwave1D2[:,0] = synwave1D2[:,1]
            synwave1D3[:,0] = synwave1D3[:,1]
            synwave1D4[:,0] = synwave1D4[:,1]
            synwave1D5[:,0] = synwave1D5[:,1]
        else:
            synwave1D[:,0] = synwave1D[:,1]
    
    if JAGURS == 0:
        np.save(outfile+'all.npy', 100*synwave1D)
        np.save(outfile+'0.npy', 100*synwave1D5)
        np.save(outfile+'1.npy', 100*synwave1D4)
        np.save(outfile+'2.npy', 100*synwave1D3)
        np.save(outfile+'3.npy', 100*synwave1D2)
    else:
        f = FortranFile('normalmode.dat', 'w')
        f.write_record(100*synwave1D)
        f.close()
        
        f = open('normalmode.namelist', 'w')
        f.write('&normalmode\n')
        f.write('srclon=%fd0\n'%srclon)
        f.write('srclat=%fd0\n'%srclat)
        f.write('dist_s=%d\n'%dist_s)
        f.write('dist_e=%d\n'%dist_e)
        f.write('tmax=%fd0\n'%tmax)
        f.write('dt=%fd0\n'%dt)
        f.write('/\n')
        f.close()
