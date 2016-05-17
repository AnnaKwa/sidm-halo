import numpy as np
import py3_dir.atomic_plots_AK as ap
import sys
sys.path.append('/home/anna/Codes/Cosmology')       #or change to wherever you keep manoj's cosmology code
from cosmo.halos import haloModel 
from scipy.integrate import quad
from scipy import interpolate
from scipy.optimize import brentq
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib as mpl
import astropy.cosmology as cosmology
from astropy import units as u
from decimal import *

mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14



# needs to be run in python 3 (halos.py)



## fractional energy losses for upscattering ~ Ehf / scattering


## to be specified
mH=100.    # GeV
alpha=0.01
Ehf=1e-5

DM_params=[mH, alpha, Ehf]
M = 1e8   # M200 at identification redshift zi
z = 10	

print ('\n Please be patient, CAMB takes a while here.')

cosmo_model = cosmology.FlatLambdaCDM(name="WMAP9", H0=69.3 * u.km / (u.Mpc * u.s), 
                      Om0=0.286, Tcmb0=2.725 * u.K, Neff=3.04, 
                      m_nu=[0.02, 0.02, 0.02] * u.eV, Ob0=0.0463)
power_params = {'sigma8' : 0.77, 'ns' : 0.9603, 'baryonic_effects' : True}
halo_params = {'cM_model': 'Bolshoi'}
hm = haloModel(cosmo_model, power_params, halo_params)

print ('\n Done.')

## constants [cgs] ##

g_per_GeV = 1.78e-24   # [g/GeV/c^2]
g_per_solarM = 2.e33   # [g/solar mass]
cm_per_kpc = 3.1e21    # [cm/kpc]
km_per_Mpc = 3.086e19  # [km/Mpc]
kB = 1.38e-16          # [erg/K]
s_per_year = 365*24*60*60.

A_mw = (0.3 * g_per_GeV) * (8./20)**2 * (270.e5)**(-3)    # phase space density normalization using Milky Way [g s^3/cm^6]


## functions ##


def scale_r(Mi,zi):
	# returns in kpc
	return hm.rs( z=zi, m200=Mi )

def rho0(Mi,zi):
	# returns in g/cm^3 
	r_s = scale_r(Mi,zi)
	r200 = hm.r200(z=zi,m200=Mi)
	int_tot = quad(lambda r: 4.*np.pi* (r**2)* ((r/r_s)*(1+(r/r_s))**2)**(-1), 1e-4*r200, r200)[0]
	return (Mi*g_per_solarM) / (int_tot * cm_per_kpc**3 )

def rho(x, Mi, zi):
	# returns in g/cm^3 
	rho_0 = rho0(Mi,zi)
	return rho_0/( x*(1+x)**2) 

def T(x, mH, A=A_mw):
	mDM = mH * g_per_GeV
	return (mDM/(3.*kB)) * (rho(x,mi,zi) * x**2 / A)**(2./3.)

def v(x, Mi, zi, A=A_mw):
	# returns in cm/s
	return (rho(x, Mi, zi) * x**2 / A)**(1./3.)

def R_scatter(x, Mi, zi, sigma_m, A=A_mw):
	# returns s^-1
	# assumes v input in cm/s, sigma_m in cm^2/g
	return rho(x,Mi,zi) * np.sqrt(32. * (v(x, Mi, zi))**2 / (3.*np.pi)) * sigma_m 
#	return np.sqrt(16./(3.*np.pi)) * A**(-1./3.) * sigma_m * rho(x,Mi,zi)**(7./3) * x**(2./3.)

def deltaT(zstart,zstop):
	return cosmo_model.age(zstop)-cosmo_model.age(zstart)

def frac_KE_loss_per_scatter(x,Mi,zi,mH, alpha, Ehf, A=A_mw):
	KE_init = 0.5 * mH * ( v(x,Mi,zi)*1e-5*ap.km_per_sec )**2 
	ge, gp = 2., 2.  
	fR = fR_from_Ehf(Ehf,alpha,ge,gp)
	E0 = alpha**2 * mH/fR
	Ehf_lost = Ehf * E0
	return Ehf_lost/KE_init

def tscatter_vs_r(dataRec, Mi,zi, mH, alpha, Ehf10, sig_type='total', A=A_mw):
	# t_scatter in yrs
	#mH,alpha,Ehf10 = DM_params[0], DM_params[1], DM_params[2]
	#sigma_m = 0.1      # = ap.sigma_from_m_v( ... )    need to write this into that module
	conc = hm.bolshoiConc(m200=Mi,z=zi)
	x_arr=np.logspace(-3, np.log10(conc),70)
	tscatter_arr=[]

	for i in range(len(x_arr)):
		x = x_arr[i]
		vel = v(x, Mi, zi)*1e-5*km_per_sec
		if sig_type=='total':    #returns total viscosity cross section
			sigma_m = sigM_from_m(dataRec,mH,alpha,vel,Ehf10, sig_type='total')[1] 
		if sig_type=='inelastic':   #returns total inelastic upscattering (not viscosity averaged) cross section
			sigma_m = sigM_from_m(dataRec,mH,alpha,vel,Ehf10, sig_type='inelastic')[0] 
		tscatter_arr.append( R_scatter(x_arr[i], Mi, zi, sigma_m, A)**(-1.) / s_per_year  )

	return ([x_arr, tscatter_arr])

def ionization_constraint(v):
	# returns minimum Ehf to not have KE>binding E in halo of v_rms=v
	v = v*km_per_sec
	ge, gp=2.,2.
	return ( 2.*v**2/(3*ge*gp) )

def t_decay(mH, R, Ehf10):
	#sigma_over_m = sig_m
	#mH = sidm_mass_single(sig_m, alpha,v,Ehf10, Ehf10_m, Etot_m, R_m, sig_interp6, sig_interp5, sig_interp4, sig_interp3, fast_decay=1) #[GeV]
	Ehf = pow(10.,Ehf10)
	ge, gp = 2., 2.
	fR = fR_from_R(R)
	alpha = np.sqrt( fR*Ehf / (2.*ge*gp/3.) )
	mp = (R + 1.) * mH / (2. + R + (1./R) ) 
	return ( (6.58e-25) * 3. * mp**2 / ( alpha * Ehf**3 * R**2 ) )      # [seconds]


def x1(dataRec, Mi, zi, mH, alpha, Ehf10):
	# returns [x1 = r1/rs, rho(x1), v(x1)] 
	# outside  r1, avg particle does not scatter within a Hubble time, i.e. approximately non-interacting
	
	hubbleT = (cosmo_model.H(zi).value / km_per_Mpc)**-1 / s_per_year      #since tscatter_vs_r returns in units of [yr]
	tscatter_output =  tscatter_vs_r(dataRec, Mi,zi, mH, alpha, Ehf10)		#using total viscosity (default in tscatter_vs_r)) to determine scatterings that lead to core
	x_arr, tscatter_arr = np.ndarray.flatten(np.array(tscatter_output[0])), np.ndarray.flatten(np.array(tscatter_output[1]))    
	tscatter_interp = interpolate.interp1d(x_arr, tscatter_arr)

	def f(x):
		return ( hubbleT - tscatter_interp(x) )
	x1_root = brentq(f, np.min(x_arr), np.max(x_arr)) 

	return (x1_root, rho(x1_root,Mi,zi), v(x1_root, Mi, zi))



def plot_tscatter_vs_r(dataRec, Mi,zi, mH, alpha, Ehf10, A=A_mw):
	#sigM = sigma / m
	fig= plt.figure(figsize=(5,5))
	ax = fig.add_subplot(111)
	ax.loglog(nonposy='clip')


	x_arr, tscatter_visc_arr = tscatter_vs_r(dataRec, Mi,zi, mH, alpha, Ehf10, sig_type='total')
	x_arr, tscatter_inel_arr = tscatter_vs_r(dataRec, Mi,zi, mH, alpha, Ehf10, sig_type='inelastic')

	x_1 = x1(dataRec, Mi, zi, mH, alpha, Ehf10)[0]

	for i in range(len(tscatter_inel_arr)):
		tscatter_visc_arr[i] = tscatter_visc_arr[i] * 1e-9
		tscatter_inel_arr[i] = tscatter_inel_arr[i] * 1e-9

	ax.plot([x_1, x_1],[1e-5,1e10], ':', color='grey', linewidth=1.1, label='$r_1$')
	ax.plot(x_arr, tscatter_visc_arr,'-', color='Red', label='Viscosity-weighted total scattering')
	ax.plot(x_arr, tscatter_inel_arr,'-', color='Blue', label='Upscattering')
	ax.set_xlabel(r'r/r$_\mathrm{s}$', fontsize=16)
	ax.set_ylabel(r'$\tau_{\mathrm{scatter}}$ [Gyr]', fontsize=16)
	ax.set_xlim([np.min(x_arr), np.max(x_arr)])
	ax.set_ylim([0.2*np.min(tscatter_inel_arr), 1.2*np.max(tscatter_inel_arr)])

	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	ax.legend(prop={'size':14},
				borderpad=.5,
				numpoints=1,
				markerscale=1,
				ncol=1,
				columnspacing=.25,
				handletextpad=.15,
				#bbox_to_anchor=(2.5, 0.75) ,)
				loc='upper left')
	

def plot_core_scattering_times(dataRec, zi, mH, alpha, Ehf10, Mmin=1e8, Mmax=1e14, A=A_mw):
	# warning: this takes forever to run
	# for given mH, alpha, Ehf10, show avg time for 1 scattering in halo cores of different mass/v
	#assume that the core density is roughly rho(r1), in reality isothermal profile goes a bit higher than this

	M_arr=np.logspace(np.log10(Mmin), np.log10(Mmax), 10)
	tscale_arr=[]     # v(r1) estimate of v in core (probably ~2x low)

	for i in range(len(M_arr)):
		Mi = M_arr[i]
		x_1,rho1,v1 = x_1(dataRec, Mi, zi, mH, alpha, Ehf10)
		v1 = v1 * 1e-5 * ap.km_per_sec            # x1 func returns v1 in units of [cm/s] but ap.sigM_from_m needs input in natural units
		sigma_m = ap.sigM_from_m(dataRec,mH,alpha,v1,Ehf10, sig_type='total')[1] 
		core_timescale = R_scatter(x_1, Mi, zi, sigma_m, A=A_mw)**-1 * 1e9 / s_per_year
		tscale_arr.append(core_timescale)
			
	fig= plt.figure(figsize=(5,5))
	ax = fig.add_subplot(111)
	ax.loglog(nonposy='clip')
	ax.plot(M_arr, tscale_arr, 'o')
	ax.set_xlabel('M$_{200}$ [M$_\odot$]')
	ax.set_ylabel('scattering time [Gyr]')

def plot_sigmaM(dataRec_m, mH, alpha, Ehf10):
	"""
	input mH in GeV
	plots total viscosity and elastic sigma/m vs velocity for a given mH, alpha, Ehf
	similar to Kim's original plot_sigma, but axes units translated for astro interpretations
	"""
	Ehf = pow(10.,Ehf10)
	cond = (dataRec_m.Ehf10)==Ehf10
	data = dataRec_m[cond]
	# determine range of alpha we can plot (based on how much data we have)
	ge,gp = (2.,2.)
	Etot = pow(10.,data.E10)
	Emin = np.min(Etot)
	Emax = np.max(Etot)

	vmin = np.sqrt( 4.*Emin*Ehf / (2.*ge*gp/3.) ) /ap.km_per_sec
	vmax = np.sqrt( 4.*(Emax-2.*Ehf)*Ehf / (2.*ge*gp/3.) ) /ap.km_per_sec

	v_arr = np.logspace(np.log10(vmin), np.log10(vmax), 100)

	sigM_tot_arr=[]
	sigM_Vtot_arr=[]
	sigM_inel_arr=[]

	for i in range(len(v_arr)):
		v = v_arr[i] * ap.km_per_sec
		sigM_tot_arr.append( ap.sigM_from_m(dataRec_m,mH,alpha,v,Ehf10, sig_type='total')[0] )     #total cross section, relevant for ionization
		sigM_Vtot_arr.append( ap.sigM_from_m(dataRec_m,mH,alpha,v,Ehf10, sig_type='total')[1] )     #second element returns viscosity sigma/m
		sigM_inel_arr.append( ap.sigM_from_m(dataRec_m,mH,alpha,v,Ehf10, sig_type='inelastic')[0] ) #first element returns total sigma/m

	fig= plt.figure(figsize=(5,5))
	ax=fig.add_subplot(111)
	#ax.set_title(r'$E_\mathrm{{hf}} = 1e{:.0f}$'.format(Ehf10))
	ax.annotate(r'$\mathrm{m}_\mathrm{{H}}=$'+repr(int(mH))+' GeV' ,
				xy=(0.05,0.17),
				xycoords='axes fraction',horizontalalignment='left',size=14)
	ax.annotate(r'$E_\mathrm{{hf}}=10^{{- {:.0f} }}$'.format(abs(Ehf10))+' E$_\mathrm{0}$',
				xy=(0.05,0.11),
				xycoords='axes fraction',horizontalalignment='left',size=14)
	ax.annotate(r'$\alpha={:.2f}$'.format(alpha),xy=(0.05,0.05),
				xycoords='axes fraction',horizontalalignment='left',size=14)
	ax.set_xscale('log')
	ax.set_yscale('log')

	ax.plot(v_arr, sigM_tot_arr, '-', color='limegreen', label='Total')
	ax.plot(v_arr, sigM_Vtot_arr, '-', color='red', label='Total viscosity')
	ax.plot(v_arr, sigM_inel_arr, '-', color='blue', label='Inelastic upscattering')
	ax.set_xlabel('v [km/s]')
	ax.set_ylabel( '$\sigma / m$ [cm$^2$ s$^{-1}$]' )
	ax.set_xlim([10., vmax])
	ax.set_ylim([1e-3, 2.])

	ax.legend()


def plot_timescales(dataRec_m, Ehf10, sigM_arr, v, rho, sig_type='total'):
	"""
	plots a single scattering timescale plot for chosen Ehf10, halo mass + redshift --> core velocity
	"""
	fig= plt.figure(figsize=(5,5))
	if sig_type=='total':
		st = fig.suptitle('Viscosity-weighted scattering',
					  fontsize="x-large")	
	if sig_type=='inelastic':
		st = fig.suptitle('Upscattering',
					  fontsize="x-large")	
	Ehf = pow(10.,Ehf10)
	cond = (dataRec_m.Ehf10)==Ehf10
	data = dataRec_m[cond]

	# determine range of alpha we can plot (based on how much data we have)
	ge,gp = (2.,2.)
	Etot = pow(10.,data.E10)
	Emin = np.min(Etot)
	Emax = np.max(Etot)

	vmin = np.sqrt( 4.*Emin*Ehf / (2.*ge*gp/3.) )
	vmax = np.sqrt( 4.*(Emax-2.*Ehf)*Ehf / (2.*ge*gp/3.) )
	if v < vmin or v > vmax:
		print ('Only have data for v between {:e} and {:e} km/s'.format(vmin /ap.km_per_sec, vmax / ap.km_per_sec))

        
	fR_min = np.min(ap.fR_from_R(data.R))
	fR_max = np.max(ap.fR_from_R(data.R))
	alpha_min = np.sqrt( fR_min*Ehf / (2.*ge*gp/3.) )
	alpha_max = np.sqrt( fR_max*Ehf / (2.*ge*gp/3.) )
	alpha_list = np.linspace(alpha_min,alpha_max,150)
	mH_arr = []
	scattering_time_arr = list( np.zeros( len(sigM_arr) ) )

	for k in range(len(sigM_arr)):
		mH_arr.append( np.array([]) )
		scattering_time_arr[k] = ap.coll_time(sigM_arr[k], rho, v)*1e-9
	
	#print np.shape(mH_arr)
	alpha_list_plot=[]
	for a in range(len(alpha_list)):
		alpha=alpha_list[a]
		if sig_type=='total':
			mass = ap.sidm_mass_arr(sigM_arr, dataRec_m,alpha,v,Ehf10, sig_type, sig_min=0.01, sig_max=100., Nsig=100)[1]    #want viscosity cross section if interested in puffing type scattering
		if sig_type=='inelastic':
			mass = ap.sidm_mass_arr(sigM_arr, dataRec_m,alpha,v,Ehf10, sig_type, sig_min=0.01, sig_max=100., Nsig=100)[0]   #want total cross section if interested in upscattering-->decays

		if mass.size>1.:
			alpha_list_plot.append(alpha)
			for k in range(len(mH_arr)):
				mH_arr[k] = np.append( mH_arr[k], mass[k] )

	ax=fig.add_subplot(111)
	#ax.set_title(r'$E_\mathrm{{hf}} = 1e{:.0f}$'.format(Ehf10))
	ax.set_xlabel(r'$m_H$ [GeV]')
	ax.set_ylabel(r'$\alpha$')

	for k in range(len(sigM_arr)):
		label='%.1E' % Decimal(scattering_time_arr[k])
		x=np.log10(scattering_time_arr[k] / np.min(scattering_time_arr) ) /np.log10(np.max(scattering_time_arr) / np.min(scattering_time_arr))
		sigM_label='%s' % float('%.1g' % sigM_arr[k])
		ax.plot(mH_arr[k], alpha_list_plot, 
			'-',
			color=cm.jet(x), 
#			'-',
			linewidth=2,
			#label=r'$\sigma$/m='+sigM_label)
			label=r'$\sigma$/m='+sigM_label +' cm$^2$/g, t= '+label+' Gyr' )

	ax.legend(prop={'size':10},
				borderpad=.5,numpoints=1,
				markerscale=1,ncol=1,columnspacing=.25,handletextpad=.15,
				bbox_to_anchor=(2.5, 0.75) )
	ax.set_xlim([0, np.max(mH_arr)*1.05])
	ax.set_ylim([alpha_min, alpha_max])
    
	scattering_time_arr.reverse()
	normalize = mcolors.LogNorm(vmin=np.max(scattering_time_arr), vmax= np.min(scattering_time_arr) )
	scalarmappable = cm.ScalarMappable(norm=normalize, cmap=cm.jet)
	scalarmappable.set_array(scattering_time_arr)
    
	cbar =fig.colorbar(scalarmappable)
	cbar.ax.set_ylabel('Scattering timescale [Gyr]')
	ax.annotate(r'v={0:d} km/s'.format(int(v/ap.km_per_sec)) ,xy=(0.95,0.10), xycoords='axes fraction', ha='right', size=14)

	ax.annotate(r'E$_{\mathrm{hf}}=1$e'+repr(Ehf10),xy=(0.95,0.03), xycoords='axes fraction', ha='right', size=14)    


##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

"""
def plot_timescales_old(dataRec_m, sigM_arr, v, rho, sig_type='total'):
	"""
	old version: produces multiple plots for all Ehf10s
	Plots lines of constant scattering timescales as defined in
	sidm_mass() in the plane of alpha vs mH. These lines have endpoints only
	due to having finite data.
	"""
	Ehf10_arr = np.unique(dataRec_m.Ehf10)
	#Ehf10_arr = [-6,-5]

	color_range=cm.jet(np.linspace(0, 1.0, len(sigM_arr)+1) )

	xdim=4.*len(Ehf10_arr)
	fig= plt.figure(figsize=(xdim,10))
	st = fig.suptitle(r'v={0:d} km/s'.format(int(v/ap.km_per_sec)),
					  fontsize="x-large")
	fig.subplots_adjust(hspace=0.5)
	fig.subplots_adjust(wspace=0.5)

	for i in range(len(Ehf10_arr)):
		# pick out data with input Ehf10
		Ehf10 = Ehf10_arr[i]
		Ehf = pow(10.,Ehf10)
		cond = (dataRec_m.Ehf10)==Ehf10
		data = dataRec_m[cond]

		# determine range of alpha we can plot (based on how much data we have)
		ge,gp = (2.,2.)
		Etot = pow(10.,data.E10)
		Emin = np.min(Etot)
		Emax = np.max(Etot)

		vmin = np.sqrt( 4.*Emin*Ehf / (2.*ge*gp/3.) )
		vmax = np.sqrt( 4.*(Emax-2.*Ehf)*Ehf / (2.*ge*gp/3.) )
		if v < vmin or v > vmax:
			print ('Only have data for v between {:e} and {:e} km/s'.format(vmin/ap.km_per_sec,vmax/ap.km_per_sec))
			ax=fig.add_subplot(2,len(Ehf10_arr),i+1)
			ax.set_title(r'$E_\mathrm{{hf}} = 1e{:.0f}$'.format(Ehf10))
			ax=fig.add_subplot(2,len(Ehf10_arr),i+1+len(Ehf10_arr))
			ax.set_title(r'$E_\mathrm{{hf}} = 1e{:.0f}$'.format(Ehf10))
			continue
		fR_min = np.min(ap.fR_from_R(data.R))
		fR_max = np.max(ap.fR_from_R(data.R))
		alpha_min = np.sqrt( fR_min*Ehf / (2.*ge*gp/3.) )
		alpha_max = np.sqrt( fR_max*Ehf / (2.*ge*gp/3.) )
		alpha_list = np.linspace(alpha_min,alpha_max,100)
		# for each alpha, find masses that match sig/mH = 0.1,1,10 cm^2/g
		# do so for total sigma and sigmaV

		mH_arr, mH_V_arr = [], []
		scattering_time_arr = np.zeros( len(sigM_arr) )

		for k in range(len(sigM_arr)):
			mH_arr.append( np.array([]) )

			mH_V_arr.append( np.array([]) )
			scattering_time_arr[k] = ap.coll_time(sigM_arr[k], rho, v)
		
		#print np.shape(mH_arr)

		for a in range(len(alpha_list)):
			alpha=alpha_list[a]
			mass,massV = ap.sidm_mass_arr(sigM_arr, dataRec_m,alpha,v,Ehf10, sig_type, sig_min=0.01, sig_max=100., Nsig=100)
			if mass.size>1.:
				for k in range(len(mH_arr)):
					mH_arr[k] = np.append( mH_arr[k], mass[k] )
					mH_V_arr[k] = np.append( mH_V_arr[k], massV[k] )
	
		ax=fig.add_subplot(2,len(Ehf10_arr),i+1)
		ax.set_title(r'$E_\mathrm{{hf}} = 1e{:.0f}$'.format(Ehf10))
		ax.set_xlabel(r'$m_H$ [GeV]')
		ax.set_ylabel(r'$\alpha$')

		for k in range(len(sigM_arr)):
			label='%.2E' % Decimal(scattering_time_arr[k])
			sigM_label='%s' % float('%.1g' % sigM_arr[k])
			print label, sigM_label
			ax.plot(mH_arr[k], alpha_list, 
				color=color_range[k], 
				marker='o',
				#label=r'$\sigma$/m='+sigM_label)
				label=r'$\sigma$/m='+sigM_label +' cm$^2$/g, t= '+label+' yr' )
		

		ax=fig.add_subplot(2,len(Ehf10_arr),i+1+len(Ehf10_arr))
		ax.set_title(r'$E_\mathrm{{hf}} = 1e{:.0f}$'.format(Ehf10))
		ax.set_xlabel(r'$m_H$ [GeV]')
		ax.set_ylabel(r'$\alpha$')
	
		for k in range(len(sigM_arr)):
			label=repr( scattering_time_arr[k] ) 
			ax.plot(mH_V_arr[k], alpha_list, 
				color=color_range[k], 
				#color=s_m.to_rgba(crange_params[k]),
				marker='o',
				label=r'Scattering time = '+label+' yr' )

	ax=fig.add_subplot(2,len(Ehf10_arr),len(Ehf10_arr))
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	ax.legend(prop={'size':10},
				borderpad=.5,numpoints=1,
				markerscale=1,ncol=1,columnspacing=.25,handletextpad=.15,
				bbox_to_anchor=(2.5, 0.75) )
				#loc='lower right',)
	st.set_y(0.95)
	#fig.subplots_adjust(top=0.83)
	#fig.tight_layout(
"""