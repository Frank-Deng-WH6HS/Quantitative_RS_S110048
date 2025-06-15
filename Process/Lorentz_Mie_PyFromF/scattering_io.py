#Lorentz-Mie散射输入输出结果对象模型. 原始代码由Mishchenko M I以Fortran写成. 

from . import np; 
from numpy import float64 as fp64, complex128 as cmplx128; 
from . import pyfunc_type; 
from .compatibility import do, span, real; 

#定义输入参数对象
class LorentzMieScatteringConfig(): 
    
    def __init__(self): 
        self.distr_mode = 3; #NDISTR, use power-law distr by default. 
        self.distr_AA = np.full(3, np.nan, dtype=fp64) #AA, AA1, AA2
        self.distr_BB = np.full(3, np.nan, dtype=fp64) #BB, BB1, BB2
        self.distr_gamma = np.nan; #GAM
        self.radius_min = np.nan; #R1
        self.radius_max = np.nan; #R2
        self.wavelength = np.nan; #LAM
        self.refractive_idx_medium = np.nan * cmplx128(); #MRR, MRI or CM1
        self.refractive_idx_particle = np.nan * cmplx128(); #CM2
        self.integ_intvl = 1; #N
        self.integ_intvl_gauss_divis = 1; #NK, use n_k instead. 
        self.integ_intvl_power_law = 1; #NP, use n_p instead. 
        self.scattering_angle_divis = 1; #NPNA, use n_pna instead. 
        self.scattering_matrix_accracy = np.nan; #DDELT
        

#定义输出结果对象
class LorentzMieScatteringResult(): 
    
    def __init__(self, n_mie: int): 
        self.radius_min = np.nan; #R1
        self.radius_max = np.nan; #R2
        self.radius_effective = np.nan; #REFF
        self.variance_effective = np.nan; #VEFF
        self.cross_sectn_extinction = np.nan; #CEXT
        self.cross_sectn_scattering = np.nan; #CSCA
        self.asymmetry_param = np.nan; #<COS>, also known as Q1
        self.area_projected_avg = np.nan; #<G>, also known as AREA
        self.volume_avg = np.nan; #<V>, also known as VOL
        self.radius_avg_volume_weighted = np.nan; #Rvw
        self.radius_avg = np.nan; #<R>, also known as RMEAN
        self.n_mie = n_mie; #NMIE; 
        n_pl = 2 * n_mie; 
        self.n_pl = n_pl; 
        self.matr_elem = real(4, n_pl); #F11, F33, F12, F34
        self.coeff_alpha = real(4, n_pl); #ALPHA1, ..., ALPHA4
        self.coeff_beta = real(2, n_pl); #BETA1, BETA2
        
#计算粒子半径的分布模式
#由于python传参方式的特殊性, 其用法改为distrb(nnk, yy, wy, config, result)
#强制要求对yy, wy执行类型检查, 要求其须为可变(mutable)对象np.ndarray; 
#强制要求对config执行类型检查, 要求其须为可变对象LorentzMieScatteringConfig; 
#强制要求对result执行类型检查, 要求其须为可变对象LorentzMieScatteringResult
#方可通过传引用方式修改. 
#返回的结果为None
def distrb(nnk: int, yy: np.ndarray, wy: np.ndarray, 
    config: LorentzMieScatteringConfig, 
    result: LorentzMieScatteringResult): 
    '''
    Usage: distrb(nnk, yy, wy, config, result)
        yy and wy are both ndarrays; 
        config is an object of LorentzMieScatteringConfig; 
        result is an object of LorentzMieScatteringResult; 
        wy and result will be modified. 
    '''
    pyfunc_type.type_check(); 
    #将传入的config中的属性值传送至函数内的局部变量
    n_distr = config.distr_mode; 
    aa = config.distr_AA.copy(); 
    bb = config.distr_BB.copy(); 
    gam = config.distr_gamma; 
    #将传入的result中的属性值传送至函数内的局部变量
    r_1 = result.radius_min; 
    r_2 = result.radius_max; 
    #开始计算
    if n_distr == 1: #Modified gamma distribution
        a2 = aa[0] / gam; 
        db = fp64(1e0) / bb[0]; 
        x = yy.copy(); 
        y = x ** aa[0]; 
        x *= db; 
        y *= np.exp(-a2 * (x ** gam)); 
        wy *= y; 
    elif n_distr == 2: #Log normal distribution
        #Label 100: 
        da = fp64(1e0) / aa[0]; 
        x = yy.copy(); 
        y = np.log(x * da); 
        y = np.exp(-(y ** 2) * fp64(0.5e0) / bb[0]) / x; 
        wy *= y; 
    elif n_distr == 3: #Power law distribution
        #Label 200: 
        x = yy.copy(); 
        wy /= (x ** 3); 
    elif n_distr == 4: #Gamma distribution
        #Label 300: 
        b2 = (fp64(1e0) - fp64(3e0) * bb[0]) / bb[0]; 
        dab = fp64(1e0) / (aa[0] * bb[0]); 
        x = yy.copy(); 
        x = (x ** b2) * np.exp(-x * dab); 
        wy *= x; 
    elif n_distr == 5: #Modified power law distribution
        #Label 360: 
        x = yy.copy(); 
        idx_x_gt_r1 = np.where(x > r_1); 
        wy[idx_x_gt_r1] *= (x[idx_x_gt_r1] / r_1) ** bb[0]; 
    elif n_distr == 6: #Bimodal volume log normal distribution
        #Label 380: 
        da1 = fp64(1e0) / aa[1]; 
        da2 = fp64(1e0) / aa[2]; 
        x = yy.copy(); 
        y1 = np.log(x * da1); 
        y2 = np.log(x * da2); 
        y1 = np.exp(-(y1 ** 2) * fp64(0.5e0) / bb[1]); 
        y2 = np.exp(-(y2 ** 2) * fp64(0.5e0) / bb[2]); 
        wy *= (y1 + gam * y2) / (x ** 4); 
    else: #Invalid n_distr
        raise ValueError("No. of distribution mode is invalid. "); 
    #Label 400: 
    sum_wy = np.nansum(wy); 
    #Label 450: 
    wy /= sum_wy; 
    #Label 500: 
    x = yy.copy(); 
    g = np.nansum((x ** 2) * wy); 
    #Label 550: 
    #x = yy.copy(); #重复计算
    r_eff = np.nansum((x ** 3) * wy); 
    #Label 600: 
    r_eff /= g; 
    #x = yy.copy(); #重复计算
    xi = x - r_eff; 
    v_eff = np.nansum((xi ** 2) * (x ** 2) * wy); 
    vol = np.nansum((x ** 3) * wy); 
    r_vw = np.nansum((x ** 4) * wy); 
    r_mean = np.nansum(x * wy); 
    #Label 650: 
    v_eff /= (g * r_eff ** 2); 
    area = g * np.pi; 
    r_vw /= vol; 
    vol *= fp64(4e0) * np.pi / fp64(3e0); 
    #用局部变量修改result中的特定属性值
    result.radius_min = r_1; 
    result.radius_max = r_2; 
    result.radius_effective = r_eff; 
    result.variance_effective = v_eff; 
    result.area_projected_avg = area; 
    result.volume_avg = vol; 
    result.radius_avg_volume_weighted = r_vw; 
    result.radius_avg = r_mean; 