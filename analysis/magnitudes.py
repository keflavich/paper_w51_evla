
distance_modulus 5*np.log10(5400)-5
# 13.661968799114842

# http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
mkabso9v = - 3.20
mko9v = distance_modulus + mkabso9v
# 10.461968799114842

# goldader 1994
ak = 2.6
mko9vext = mko9v + ak
# 13.061968799114842
