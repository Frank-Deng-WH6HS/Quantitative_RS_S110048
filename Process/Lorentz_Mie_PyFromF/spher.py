#Lorentz-Mie代码中spher.f专有函数. 原始代码由Mishchenko M I以Fortran写成. 

from . import np; 
from numpy import float64 as fp64; 
from . import pyfunc_type; 
from .math_func import gener, gauss, power; 
from .scattering_io import LorentzMieScatteringConfig; 
from .scattering_io import LorentzMieScatteringResult; 
from .scattering_io import distrb; 
from .compatibility import do, span, real; 

#angl函数

#spher函数, 用于计算α1, ..., α4, β1, β2
#由于python传参方式的特殊性, 其用法改为spher(config, result)
#强制要求对config执行类型检查, 要求其须为可变对象LorentzMieScatteringConfig; 
#强制要求对result执行类型检查, 要求其须为可变对象LorentzMieScatteringResult
#方可通过传引用方式修改. 
#返回的结果为None
def spher(
    config: LorentzMieScatteringConfig, 
    result: LorentzMieScatteringResult): 
    '''
    Usage: spher(config, result)
        config is an object of LorentzMieScatteringConfig; 
        result is an object of LorentzMieScatteringResult; 
        result will be modified. 
    '''
    pyfunc_type.type_check(); 
    #将传入的config中的属性值传送至函数内的局部变量
    n_distr = config.distr_mode; 
    aa = config.distr_AA; 
    bb = config.distr_BB; 
    gam = config.distr_gamma; 
    r_1 = config.radius_min; 
    r_2 = config.radius_max; 
    lam = config.wavelength; 
    mr_r = config.refractive_idx_medium.real; 
    mr_ii = config.refractive_idx_medium.imag; 
    n = config.integ_intvl; 
    n_k = config.integ_intvl_gauss_divis; 
    n_p = config.integ_intvl_power_law; 
    n_pna = config.scattering_angle_divis; 
    d_delt = config.scattering_matrix_accracy; 
    #将传入的result中的属性值传送至函数内的局部变量
    n_mie = result.n_mie; 
    n_grad = 100000; #暂时与源代码中保持一致, 后续可能扩展
    n_pl = result.n_pl; 
    n_drdi = 3 * n_mie; 
    #声明数组
    psi = real(n_mie); hi = real(n_drdi); rpsi = real(n_drdi); 
    p = real(4, n_pl); 
    d_r = real(n_drdi); d_i = real(n_drdi); 
    a_r = real(n_mie); a_i = real(n_mie); 
    b_r = real(n_mie); b_i = real(n_mie); 
    x = real(n_pl); w = real(n_pl); 
    xx = real(n_pl); yy = real(n_grad); wy = real(n_grad); 
    coeff = real(3, n_mie); 
    coef = real(8, n_pl);  
    pin = real(n_mie); taun = real(n_mie); 
    #开始计算
    mr_i = -mr_ii; 
    if n_distr == 3: 
        (r_1, r_2) = power(aa[0], bb[0], r_1, r_2); 
    result.radius_min = r_1; 
    result.radius_max = r_2; 
    #计算高斯区间划分点和权重. 
    if n_k > n_pl: 
        raise ValueError("n_k > n_pl. Evaluation terminated. "); 
    breakpoint(); 
    gauss(n_k, 0, 0, x, w); 
    n_nk = n * n_k; 
    if n_distr != 5: 
        n_npk = 0; 
    else: 
        n_npk = n_p * n_k; 
    n_nk += n_npk; 
    if n_nk > n_grad: 
        raise ValueError("n_nk > n_grad. Evaluation terminated. "); 
    wn = fp64(2e0) * np.pi / lam; 
    rx = r_2 * wn; 
    m = int(rx + fp64(4.05e0) * np.cbrt(rx) + fp64(8e0)); 
    if m > n_mie: 
        raise ValueError("Too many Mie coefficients, increase n_mie. "); 
    dd = np.array(list(do(1, m)), dtype=fp64); 
    dd1 = dd + fp64(1e0); 
    coeff[1, span(1, m)] = dd1 / dd; 
    dd2 = fp64(2e0) * dd + fp64(1e0); 
    coeff[2, span(1, m)] = dd2.copy(); 
    coeff[3, span(1, m)] = fp64(0.5e0) * dd2 / (dd * dd1); 
    #Label 27: 
    n_g = 2 * m - 1; 
    l1_max = 2 * m; 
    cm = fp64(1e0) / (mr_r ** 2 + mr_i ** 2); 
    rm_r = mr_r * cm; 
    rm_i = -mr_i * cm; 
    if not n_distr != 5: 
        z_1 = r_1 / fp64(n_p); 
        z_2 = z_1 * fp64(0.5e0); 
        xx[span(1, n_k)] = z_2 * x[span(1, n_k)] + z_2; 
        #Label 28: 
        for j in do(1, n_p): 
            j_1 = j - 1; 
            z_j = z_1 * fp64(j_1); 
            i_j = j_1 * n_k; 
            for i in do(1, n_k): 
                i_1 = i + i_j; 
                yy[i_1] = xx[i] + z_j; 
                wy[i_1] = w[i] * z_2; 
    #Label 30: 
    z_1 = (r_2 - r_1) / fp64(n); 
    z_2 = z_1 * fp64(0.5e0); 
    xx[span(1, n_k)] = z_2 * x[span(1, n_k)] + z_2; 
    #Label 32: 
    for j in do(1, n): 
        j_1 = j - 1; 
        z_j = z_1 * fp64(j_1) + r_1; 
        i_j = j_1 * n_k + n_npk; 
        for i in do(1, n_k): 
            i_1 = i + i_j; 
            yy[i_1] = xx[i] + z_j; 
            wy[i_1] = w[i] * z_2; 
    breakpoint(); 
    #Label 34: 
    distrb(n_nk, yy, wy, config, result); 