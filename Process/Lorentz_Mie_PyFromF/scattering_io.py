#Lorentz-Mie散射输入输出结果对象模型. 

from . import np; 
from numpy import uint16, uint32; 
from numpy import float64 as fp64, complex128 as cmplx128; 

class LorentzMieScatteringConfig(): 
    
    def __init__(self): 
        self.distr_mode = 3; #NDISTR, use power-law distr by default. 
        self.radius_min = np.nan; #R1
        self.radius_max = np.nan; #R2
        self.wavelength = np.nan; #LAM
        self.refractive_idx_medium = np.nan * cmplx128(); #MRR, MRI or CM1
        self.refractive_idx_particle = np.nan * cmplx128(); #CM2
        self.integ_intvl = uint32(1); #N
        self.integ_intvl_power_law = uint32(1); #NP, use n_p instead. 
        self.integ_intvl_gauss_divis = uint32(1); #NK, use n_k instead. 
        self.scattering_angle_divis = uint32(1); #NPNA, use n_pna instead. 
        self.scattering_matrix_accracy = np.nan; #DDELT
        

#定义输出结果对象
class LorentzMieScatteringResult(): 
    
    def __init__(self, n_mie, n_pl, n_drdi): 
        self.radius_min = np.nan; #R1
        self.radius_max = np.nan; #R2
        self.radius_effective = np.nan; #REFF
        self.variance_effective = np.nan; #VEFF
        self.cross_sectn_extinction = np.nan; #CEXT
        self.cross_sectn_scattering = np.nan; #CSCA
        self.asymmetry_param = np.nan; #<COS>
        self.area_projected_avg = np.nan; #<G>
        self.volume_avg = np.nan; #<V>
        self.radius_avg_volume_weighted = np.nan; #Rvw
        self.radius_avg = np.nan; #<R>
        self.matr_elem \
            = np.full((4, n_pl), np.nan, dtype=fp64); #F11, F33, F12, F34
        self.coeff_alpha \
            = np.full((4, n_pl), np.nan, dtype=fp64); #ALPHA1, ..., ALPHA4
        self.coeff_beta \
            = np.full((2, n_pl), np.nan, dtype=fp64); #BETA1, BETA2
        