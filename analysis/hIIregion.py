import numpy as np
import pylab as pl
from scipy import optimize
from mpfit import mpfit
from astropy import constants
from astropy import units as u

kb = constants.k_B
c = constants.c
mu = 1.4
mh = constants.m_p
msun = u.M_sun
hplanck = constants.h
emu = u.cm**-6*u.pc

def val_with_unit(var, unit):
    if hasattr(var, 'unit'):
        return var.to(unit).value
    else:
        return var

def with_unit(var, unit):
    if hasattr(var, 'unit'):
        return var.to(unit)
    else:
        return var*unit

def tnu(Te,nu,EM):
    """
    Te - excitation temperature
    nu - frequency in GHz
    EM - Emission Measure

    Calculates optical depth as a function of temperature, frequency, and
    emission measure from Rohlfs and Wilson 2000's eqns 9.33 and 9.34.

    """
    Te = val_with_unit(Te, u.K)
    nu = val_with_unit(nu, u.GHz)
    EM = val_with_unit(EM, emu)

    nu0 = Te**1.5 / 1000
    answer_highnu = (nu > nu0) * 3.014e-2 * Te**-1.5 * nu**-2 * EM  
    gff_lownu = ( np.log(4.955e-2 * nu**-1) + 1.5 * np.log(Te) )  # <gff> Gaunt factor for free-free
    answer_lownu  = (nu < nu0) * 3.014e-2 * Te**-1.5 * nu**-2 * EM * gff_lownu
    tau = answer_lownu+answer_highnu
    # altenhoff version
    tau = 8.235e-2 * Te**-1.35 * nu**-2.1 * EM
    return tau

def Inu(nu,tau,Te,I0=0):
    """
    Calculates flux for a given optical depth, frequency, and temperature
    assuming Rayleigh-Jeans

    Parameters
    ----------
    I0 : float
        Scale factor.  If zero, will be determined.  If not, will still be
        determined.
    """
    Te = with_unit(Te, u.K)
    nu = with_unit(nu, u.GHz)

    if I0==0 and not np.isscalar(nu):
        whtau1 = np.argmin(np.abs(tau-1))
        nutau1 = nu[whtau1]
        taufactor = 1
    else:
        nutau1 = nu
        taufactor = tau
        """ assumes I0 is set"""
    #I0 = 2 * kb * Te * nutau1**2 / c**2 * taufactor
    #thin = (tau < 1) * (exp(1-tau)) * I0

    I0 = 2 * kb * Te * nutau1**2 / c**2 
    thin = (tau < 5) * (1-np.exp(-taufactor)) * I0
    thick = 2 * kb * Te * (nu * (tau > 5))**2 / c**2
    return (thin+thick).to(u.Jy)

def inorm(em,nu,nu0,intens0,Te=8500*u.K):
    nu = with_unit(nu, u.GHz)
    nu0 = with_unit(nu0, u.GHz)
    em = with_unit(em, emu)
    I0 = 2 * kb * Te * nu0**2 / c**2
    model_intensity0 = Inu(nu0,tnu(Te,nu0,em),Te,I0=I0)
    model_intensity = Inu(nu,tnu(Te,nu,em),Te,I0=I0)
    model_norm = intens0/model_intensity0 * model_intensity
    return model_norm

def inufit(nu,em,normfac,Te=8500*u.K):
    nu = with_unit(nu, u.GHz)
    em = with_unit(em, emu)
    Te = with_unit(Te, u.K)
    I0 = 2 * kb * Te * nu[0]**2 / c**2
    model_intensity = Inu(nu,tnu(Te,nu,em),Te,I0=I0) 
    model_norm = normfac * model_intensity
    return model_norm

def inufit_dust(nu,em,normfac,alpha,normfac2,Te=8500*u.K):
    nu = with_unit(nu, u.GHz)
    em = with_unit(em, emu)
    Te = with_unit(Te, u.K)
    I0 = 2 * kb * Te * nu[0]**2 / c**2
    model_intensity = Inu(nu,tnu(Te,nu,em),Te,I0=I0) 
    model_norm = normfac * model_intensity + normfac2*nu**alpha
    return model_norm

def inufit_dustT(nu,em,normfac,beta,normfac2,dustT,Te=8500*u.K):
    Te = with_unit(Te, u.K)
    dustT = with_unit(dustT, u.K)
    nu = with_unit(nu, u.GHz)
    em = with_unit(em, emu)

    I0 = 2 * kb * Te * nu[0]**2 / c**2
    model_intensity = Inu(nu,tnu(Te,nu,em),Te,I0=I0) 
    dustem = 2*hplanck*(nu)**(3+beta) / c**2 * (np.exp(hplanck*nu/(kb*(dustT))) - 1)**-1
    model_norm = normfac * model_intensity + normfac2/(dustT)*dustem
    return model_norm

def mpfitfun(x,y,err=None,dust=False,dustT=False):
    if dust:
        if err == None:
            def f(p,fjac=None): return [0,(y-inufit_dust(x,*p))]
        else:
            def f(p,fjac=None): return [0,(y-inufit_dust(x,*p))/err]
        return f
    elif dustT:
        if err == None:
            def f(p,fjac=None): return [0,(y-inufit_dustT(x,*p))]
        else:
            def f(p,fjac=None): return [0,(y-inufit_dustT(x,*p))/err]
        return f
    else:
        if err == None:
            def f(p,fjac=None): return [0,(y-inufit(x,*p))]
        else:
            def f(p,fjac=None): return [0,(y-inufit(x,*p))/err]
        return f
