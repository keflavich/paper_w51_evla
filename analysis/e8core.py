"""
Simple analysis of the e8 core.  Results:

Mass assuming n=300000.0 1 / cm3: 1.1395251902046273 solMass
Mass assuming n=1000000.0 1 / cm3: 3.798417300682091 solMass
Q(20) = 83.78643448498852, Q(150) = 1178.4217244185054
Q(20) = 49.49153310105068, Q(150) = 1019.97053159
Column density assuming LTE and T=Tex=25.0 K: 6.298933e+15 1 / cm2
Column density assuming LTE and T=Tex=50.0 K: 1.075059e+16 1 / cm2
Column density assuming LTE and T=20.0 K and Tex=20.0 K: 5.748722e+15 1 / cm2
Column density assuming LTE and T=100.0 K and Tex=100.0 K: 2.353634e+16 1 / cm2
Column density assuming LTE and T=200.0 K and Tex=200.0 K: 5.861634e+16 1 / cm2
Mass given column, assuming T=200.0 K X=1e-09: 1893.813870236937 solMass
Mass given column, assuming T=200.0 K X=1e-07: 18.938138702369372 solMass
Mass given column, assuming T=200.0 K X=3e-07: 6.3127129007897915 solMass
Mass given column, assuming T=200.0 K X=5e-08: 37.876277404738744 solMass

"""
import numpy as np
from astropy import units as u
from astropy import constants
import pyradex

integrated_intensity = 18.4*u.mJy # from centroid_e8
line_integral = 57.247*u.K*u.km/u.s # from centroid_e8; 2 beams
radius = 0.0152*u.pc

# assume it is a spherical gaussian, so sqrt(2 pi)^3
density = 3e5 * u.cm**-3
mass_n3e5 = density * ((2*np.pi)**0.5*radius)**3 * (2.8*u.Da)
print("Mass assuming n={0}: {1}".format(density, mass_n3e5.to(u.M_sun)))
density = 1e6 * u.cm**-3
mass_n1e6 = density * ((2*np.pi)**0.5*radius)**3 * (2.8*u.Da)
print("Mass assuming n={0}: {1}".format(density, mass_n1e6.to(u.M_sun)))

# column density computation based on Mangum & Shirley 2015
dipole_moment = 2.331*u.Debye*(1e-18*u.esu/u.Debye)
dipole_moment = 2.311e-18*u.esu*u.cm
# 1 Debye = 1e-18 esu

def Jnu(nu, T):
    return constants.h*nu/constants.k_B / (np.exp(constants.h*nu/(constants.k_B*T))-1)

R = pyradex.Radex(species='oh2co-h2', density=1e3, temperature=20, column=1e12)

from h2co_params import Qrot
print("Q(20) = {0}, Q(150) = {1}".format(Qrot(20), Qrot(150)))

from h2co_modeling import lte_model
Qrot = lambda x: np.interp(x.to(u.K).value if hasattr(x,'to') else x,
                           lte_model.tem, lte_model.Q)
print("Q(20) = {0}, Q(150) = {1}".format(Qrot(20), Qrot(150)))


def Nthintot(nu, line_integral, degeneracy_u, dipole_moment, temperature, Eu, Tex, Tbg=2.7315*u.K, fillingfactor=1.0):
    """
    Eqn 80, where S mu^2 R_i = line_strength^1/2
    because R_i = relative strength of hyperfines = 1
    Copied from ColumnDensity.ipynb
    """
    term1 = (3*constants.h/(8*np.pi**3 * dipole_moment**2))
    term2 = Qrot(temperature)/degeneracy_u
    term3 = np.exp(Eu/(constants.k_B * Tex)) / (np.exp(constants.h*nu/(constants.k_B * Tex)) - 1)
    term4 = 1./(Jnu(nu,Tex)-Jnu(nu,Tbg))
    term5 = line_integral.to(u.K*u.km/u.s) / fillingfactor
    return term1 * term2 * term3 * term4 * term5

# def Nthintot(nu, line_integral, degeneracy_u, dipole_moment, temperature, Eu,
#              Tex, Ju=2, K=1, Tbg=2.7315*u.K):
#     Ri = 1
#     term1 = (3*constants.h/(8*np.pi**3 * dipole_moment**2 * Ri))
#     term1a = (Ju*(Ju+1))/K**2
#     gJ = 2*Ju+1
#     gK = 2
#     gI = 3/4.
#     term2 = Qrot(temperature) / (gJ*gK*gI)
#     term3 = np.exp(Eu/(constants.k_B * Tex)) / (np.exp(constants.h*nu/(constants.k_B * Tex)) - 1)
#     term4 = 1./(Jnu(nu,Tex)-Jnu(nu,Tbg))
#     term5 = line_integral.to(u.K*u.km/u.s)
#     return term1 * term1a * term2 * term3 * term4 * term5

temperature = Tex = 25*u.K

J = np.array([int(x.split(b"_")[0]) for x in R.quantum_number[R.upperlevelindex]])
gI = 3/4. # see h2co_modeling.lte_model
degeneracy_u = (2*J+1)*gI

column = Nthintot(nu=R.frequency[2], line_integral=line_integral,
                  degeneracy_u=degeneracy_u[2],
                  dipole_moment=dipole_moment,
                  temperature=temperature,
                  Eu=R.upperstateenergy[2]*u.K*constants.k_B,
                  Tex=Tex,).to(u.cm**-2)

print("Column density assuming LTE and T=Tex={0}: {1:e}".format(temperature, column))

temperature = Tex = 50*u.K
column = Nthintot(nu=R.frequency[2], line_integral=line_integral,
                  degeneracy_u=degeneracy_u[2],
                  dipole_moment=dipole_moment, temperature=temperature,
                  Eu=R.upperstateenergy[2]*u.K*constants.k_B,
                  Tex=Tex,).to(u.cm**-2)
print("Column density assuming LTE and T=Tex={0}: {1:e}".format(temperature, column))

temperature = 20*u.K
Tex = 20*u.K

column = Nthintot(nu=R.frequency[2], line_integral=line_integral,
                  degeneracy_u=degeneracy_u[2],
                  dipole_moment=dipole_moment, temperature=temperature,
                  Eu=R.upperstateenergy[2]*u.K*constants.k_B,
                  Tex=Tex,).to(u.cm**-2)
print("Column density assuming LTE and T={2} and Tex={0}: {1:e}".format(temperature, column, temperature))

temperature = 100*u.K
Tex = 100*u.K

column = Nthintot(nu=R.frequency[2], line_integral=line_integral,
                  degeneracy_u=degeneracy_u[2],
                  dipole_moment=dipole_moment, temperature=temperature,
                  Eu=R.upperstateenergy[2]*u.K*constants.k_B,
                  Tex=Tex,).to(u.cm**-2)
print("Column density assuming LTE and T={2} and Tex={0}: {1:e}".format(temperature, column, temperature))

temperature = 200*u.K
Tex = 200*u.K

column = Nthintot(nu=R.frequency[2], line_integral=line_integral,
                  degeneracy_u=degeneracy_u[2],
                  dipole_moment=dipole_moment, temperature=temperature,
                  Eu=R.upperstateenergy[2]*u.K*constants.k_B,
                  Tex=Tex,).to(u.cm**-2)
print("Column density assuming LTE and T={2} and Tex={0}: {1:e}".format(temperature, column, temperature))

Xh2co = 1e-9
mass_col = column * ((2*np.pi)**0.5 * radius)**2 * 2.8*u.Da / Xh2co
print("Mass given column, assuming T={tex} X={Xh2co}: {0}".format(mass_col.to(u.M_sun), Xh2co=Xh2co, tex=Tex))
Xh2co = 1e-7
mass_col = column * ((2*np.pi)**0.5 * radius)**2 * 2.8*u.Da / Xh2co
print("Mass given column, assuming T={tex} X={Xh2co}: {0}".format(mass_col.to(u.M_sun), Xh2co=Xh2co, tex=Tex))
Xh2co = 3e-7
mass_col = column * ((2*np.pi)**0.5 * radius)**2 * 2.8*u.Da / Xh2co
print("Mass given column, assuming T={tex} X={Xh2co}: {0}".format(mass_col.to(u.M_sun), Xh2co=Xh2co, tex=Tex))
Xh2co = 5e-8
mass_col = column * ((2*np.pi)**0.5 * radius)**2 * 2.8*u.Da / Xh2co
print("Mass given column, assuming T={tex} X={Xh2co}: {0}".format(mass_col.to(u.M_sun), Xh2co=Xh2co, tex=Tex))
