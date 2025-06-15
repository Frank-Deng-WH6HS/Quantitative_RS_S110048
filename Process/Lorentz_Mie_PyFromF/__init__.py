#Lorentz-Mie散射代码改编包. 原始代码由Mishchenko M I以Fortran写成. 

#参考文献: 
#[1] Mishchenko MI, Travis LD, Lacis AA. Scattering, Absorption, and  
#    Emission of Light by Small Particles[M]. Cambridge: Cambridge 
#    University Press, 2002. 
#[2] Mishchenko MI. Double-precision Lorenz-Mie program for the case 
#    of a nonabsorbing host medium[EB/OL]. (2005-08-06) 
#    https://www.giss.nasa.gov/staff/mmishchenko/ftpcode/spher.f
#[3] Mishchenko MI. Double-precision Lorenz-Mie program for the case 
#    of an absorbing host medium[EB/OL]. (2018-04-13) 
#    https://www.giss.nasa.gov/staff/mmishchenko/ftpcode/spher_abs_host_2.f
#[4] Schmunk RB, Schmidt GA. Electromagnetic Scattering by Particles 
#    and Surfaces[EB/OL]. (2023-08-08) 
#    https://www.giss.nasa.gov/staff/mmishchenko/Lorenz-Mie.html.

import numpy as np; 
import pyfunc_type; 
from . import math_func; 
from .scattering_io import LorentzMieScatteringConfig; 
from .scattering_io import LorentzMieScatteringResult; 
