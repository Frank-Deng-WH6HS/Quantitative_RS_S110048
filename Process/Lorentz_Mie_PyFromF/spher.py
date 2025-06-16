#Lorentz-Mie代码中spher.f专有函数. 原始代码由Mishchenko M I以Fortran写成. 

import sys; 

from . import np; 
from numpy import float64 as fp64; 
from . import pyfunc_type; 
from .math_func import gener, gauss, power; 
from .scattering_io import LorentzMieScatteringConfig; 
from .scattering_io import LorentzMieScatteringResult; 
from .scattering_io import distrb; 
from .compatibility import do, span, real; 

#angl函数
#由于python传参方式的特殊性, 其用法改为angl(n_max, u, pin, taun, coeff)
def angl(n_max: int, u: np.float64, 
    pin: np.ndarray, taun: np.ndarray, coeff: np.ndarray):
    '''
    Usage: angl(n_max, u, pin, taun, coeff)
        pin, taun and coeff are ndarrays, 
        pin and taun will be modified. 
    '''
    #由于本函数将在spher执行期间的循环中大量调用, 因此类型检查功能在调试无误后移除
    p_1 = fp64(0e0); p_2 = fp64(1e0); 
    for n in do(1, n_max): 
    #Label 5: 
        s = u * p_2; t = s - p_1; 
        taun[n] = fp64(n) * t - p_1; pin[n] = p_2; 
        p_1 = p_2; p_2 = s + coeff[1, n] * t;   

#spher函数, 用于计算α1, ..., α4, β1, β2
#由于python传参方式的特殊性, 其用法改为spher(config, result)
#强制要求对config执行类型检查, 要求其须为可变对象LorentzMieScatteringConfig; 
#强制要求对result执行类型检查, 要求其须为可变对象LorentzMieScatteringResult
#方可通过传引用方式修改. 
#返回的结果为int
def spher(
    config: LorentzMieScatteringConfig, 
    result: LorentzMieScatteringResult) -> int: 
    '''
    Usage: l1_max = spher(config, result)
        config is an object of LorentzMieScatteringConfig; 
        result is an object of LorentzMieScatteringResult; 
        result will be modified. 
    '''
    pyfunc_type.type_check(); 
    sys.stderr.write(
        "Performing initial variables...\n"
    ); sys.stderr.flush(); 
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
    f_ij = result.matr_elem; 
    alpha = result.coeff_alpha; 
    beta = result.coeff_beta; 
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
    sys.stderr.write(
        "Calculating Gaussian division points and weights...\n"
    ); sys.stderr.flush(); 
    if n_k > n_pl: 
        raise ValueError("n_k > n_pl. Evaluation terminated. "); 
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
    
    #计算粒子半径分布
    sys.stderr.write(
        "Calculating distribution of particle radii...\n"
    ); sys.stderr.flush(); 
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
    #Label 34: 
    distrb(n_nk, yy, wy, config, result); 
    c_ext = fp64(0e0); c_sca = fp64(0e0); 
    gauss(n_g, 0, 0, x, w); 
    f_ij[span(1, 4), span(1, n_g)] = fp64(0e0); 
    
    #Label 40: 
    #根据粒径计算均值
    sys.stderr.write(
        "Averaging over sizes...\n"
    ); sys.stderr.flush(); 
    y = yy.copy(); zw = wy.copy(); 
    rx = y * wn; rx_r = mr_r * rx; rx_i = mr_i * rx; 
    cdd = rx_r ** 2 + rx_i ** 2; 
    cd = np.sqrt(cdd); 
    dc = np.cos(rx); ds = np.sin(rx); 
    cx = 1 / rx; cx_r = rx_r / cdd; cx_i = -rx_i / cdd; 
    
    for i in do(1, n_nk): 
        #在控制台显示进度
        sys.stderr.write("Cycle: %5d of %5d" % (i, n_nk)); 
        sys.stderr.write("\x0d"); 
        #正式执行循环
        m_1 = int(rx[i] + fp64(4.05e0) * np.cbrt(rx[i]) + fp64(8e0)); 
        m_2 = m_1 + 2 + int(fp64(1.2e0) * np.sqrt(rx[i])) + 5; 
        if m_2 > n_drdi: 
            raise ValueError("m_2 > n_drdi. Evaluation terminated. "); 
        m_3 = m_2 - 1; 
        q_max = np.max([fp64(m_1), cd[i]]); 
        m_4 = int(fp64(6.4e0) * np.cbrt(q_max) + q_max) + 8; 
        if m_4 > n_drdi: 
            raise ValueError("m_4 > n_drdi. Evaluation terminated. "); 
        d_4 = fp64(m_4 + 1); 
        d_r[m_4] = d_4 * cx_r[i]; d_i[m_4] = d_4 * cx_i[i]; 
        hi[1] = ds[i] + dc[i] * cx[i]; 
        hi[2] = fp64(3e0) * hi[1] * cx[i] - dc[i]; 
        psi[1] = cx[i] * ds[i] - dc[i]; 
        rpsi[m_2] = rx[i] / fp64(2 * m_2 + 1); 
        for j in do(2, m_3): 
            j_1 = m_2 - j + 1; 
            rpsi[j_1] = fp64(1e0) / fp64(
                (2 * j_1 + i) * cx[i] - rpsi[j_1 + 1]
            ); 
        for j in do(2, m_4): 
            j_1 = m_4 - j + 2; j_2 = j_1 - 1; 
            d_j = fp64(j_1); 
            f_r = d_j * cx_r[i]; f_i = d_j * cx_i[i]; 
            o_r = d_r[j_1] + f_r; o_i = d_i[j_1] + f_i; 
            o_ri = fp64(1e0) / (o_r ** 2 + o_i ** 2); 
            d_r[j_2] = f_r - o_r * o_ri; 
            d_i[j_2] = f_i + o_i * o_ri; 
        m_2 = m_1 - 1; 
        for j in do(2, m_2): 
            j_1 = j - 1; j_2 = j + 1; 
            hi[j_2] = fp64(2 * j + 1) * hi[j] * cx[i] - hi[j_1]; 
        for j in do(2, m_1): 
            psi[j] = rpsi[j] * psi[j - 1]; 
        psi_1 = psi[1]; 
        d_r_1 = d_r[1]; 
        d_i_1 = d_i[1]; 
        hi_1 = hi[1]; 
        o_r = d_r_1 * rm_r - d_i_1 * rm_i + cx[i];
        o_r_1 = o_r * psi_1 - ds[i]; 
        o_i = d_r_1 * rm_i + d_i_1 * rm_r; 
        o_i_1 = o_i * psi_1; 
        o_r_2 = o_r * psi_1 - o_i * hi_1 - ds[i]; 
        o_i_2 = o_r * hi_1 + o_i * psi_1 - dc[i]; 
        #以下部分是复数除法, 但原作者甚至未将除法公式封装成函数...
        #复数除法在这个循环体当中出现了四次...
        o_ab = fp64(1e0) / (o_r_2 ** 2 + o_i_2 ** 2); 
        a_r[1] = (o_r_1 * o_r_2 + o_i_1 * o_i_2) * o_ab; 
        a_i[1] = (o_r_2 * o_i_1 - o_r_1 * o_i_2) * o_ab; 
        o_r = mr_r * d_r_1 - mr_i * d_i_1 + cx[i]; 
        o_i = mr_r * d_i_1 + mr_i * d_r_1; 
        o_r_1 = o_r * psi_1 - ds[i]; 
        o_r_2 = o_r * psi_1 - o_i * hi_1 - ds[i]; 
        o_i_1 = o_i * psi_1; 
        o_i_2 = o_r * hi_1 + o_i * psi_1 - dc[i]; 
        o_ab = fp64(1e0) / (o_r_2 ** 2 + o_i_2 ** 2); 
        b_r[1] = (o_r_1 * o_r_2 + o_i_1 * o_i_2) * o_ab; 
        b_i[1] = (o_r_2 * o_i_1 - o_r_1 * o_i_2) * o_ab; 
        for j in do(2, m_1): 
            j_1 = j - 1; 
            d_j = fp64(j) * cx[i]; 
            psi_1 = psi[j]; psi_2 = psi[j_1]; 
            hi_1 = hi[j]; hi_2 = hi[j_1]; 
            d_r_1 = d_r[j]; d_i_1 = d_i[j]; 
            o_r = d_r_1 * rm_r - d_i_1 * rm_i + d_j; 
            o_i = d_r_1 * rm_i + d_i_1 * rm_r; 
            o_r_1 = o_r * psi_1 - psi_2; 
            o_i_1 = o_i * psi_1; 
            o_r_2 = o_r * psi_1 - o_i * hi_1 - psi_2; 
            o_i_2 = o_r * hi_1 + o_i * psi_1 - hi_2; 
            o_ab = fp64(1e0) / (o_r_2 ** 2 + o_i_2 ** 2); 
            ya_r = (o_r_1 * o_r_2 + o_i_1 * o_i_2) * o_ab; 
            ya_i = (o_r_2 * o_i_1 - o_r_1 * o_i_2) * o_ab; 
            a_r[j] = ya_r; a_i[j] = ya_i; 
            o_r = mr_r * d_r_1 - mr_i * d_i_1 + d_j; 
            o_i = mr_r * d_i_1 + mr_i * d_r_1; 
            o_r_1 = o_r * psi_1 - psi_2; 
            o_r_2 = o_r * psi_1 - o_i * hi_1 - psi_2; 
            o_i_1 = o_i * psi_1; 
            o_i_2 = o_r * hi_1 + o_i * psi_1 - hi_2; 
            o_ab = fp64(1e0) / (o_r_2 ** 2 + o_i_2 ** 2); 
            yb_r = (o_r_1 * o_r_2 + o_i_1 * o_i_2) * o_ab; 
            yb_i = (o_r_2 * o_i_1 - o_r_1 * o_i_2) * o_ab; 
            b_r[j] = yb_r; b_i[j] = yb_i;                           
            ya_r = ya_r ** 2 + ya_i ** 2 + yb_r ** 2 + yb_i ** 2; 
        #Label 50: 
        ce = fp64(0e0); cs = fp64(0e0); 
        for j in do(1, m_1): 
            cj = coeff[2, j]; 
            a_rj = a_r[j]; a_ij = a_i[j]; 
            b_rj = b_r[j]; b_ij = b_i[j]; 
            cda = a_rj ** 2 + a_ij ** 2; 
            cdb = b_rj ** 2 + b_ij ** 2; 
            ce += cj * (a_rj + b_rj); 
            cs += cj * (cda + cdb); 
            cj = coeff[3, j]; 
            a_r[j] = cj * (a_rj + b_rj); 
            a_i[j] = cj * (a_ij + b_ij); 
            b_r[j] = cj * (a_rj - b_rj); 
            b_i[j] = cj * (a_ij - b_ij); 
        #Label 60: 
        c_ext += zw[i] * ce; c_sca += zw[i] * cs; 
        for k in do(1, n_g): 
            angl(m_1, x[k], pin, taun, coeff); 
            sp_r = sp_i = sm_r = sm_i = fp64(0e0); 
            for j in do(1, m_1): 
                pj = pin[j]; tj = taun[j]; 
                pp = tj + pj; pm = tj - pj; 
                sp_r += a_r[j] * pp; 
                sp_i += a_i[j] * pp; 
                sm_r += b_r[j] * pp; 
                sm_i += b_i[j] * pp; 
            #Label 65: 
            d_1 = (sp_r ** 2 + sp_i ** 2) * zw[i]; 
            d_2 = (sm_r ** 2 + sm_i ** 2) * zw[i]; 
            f_ij[1, k] += d_1 + d_2; 
            f_ij[2, k] += d_1 - d_2; 
            f_ij[3, k] += (sp_r * sm_r + sp_i * sm_i) * zw[i] * fp64(2e0); 
            f_ij[4, k] += (sp_r * sm_i - sp_i * sm_r) * zw[i] * fp64(2e0); 
        #Label 70: 
    #Label 100: 
    dd = fp64(2e0) / c_sca; 
    f_ij[span(1, 4), span(1, n_g)] *= dd; 
    #Label 120: 
    #计算截面和单次散射反照率
    vol_1 = fp64(2e0) * np.pi / (wn ** 2); 
    c_ext *= vol_1; c_sca *= vol_1; 
    alb = c_sca / c_ext; 
    #计算展开系数
    sys.stderr.write(
        "Calculating expansion coefficients...\n"
    ); sys.stderr.flush(); 
    for l1 in do(3, l1_max): 
        l = l1 - 1;                                                                 
        coef[span(1, 8), l1] = np.array([
            fp64(1e0) / fp64(l + 1), 
            fp64(2 * l + 1), 
            fp64(1e0) / np.sqrt(fp64((l + 1) ** 2 - 4)), 
            np.sqrt(fp64(l ** 2 - 4)), 
            fp64(1e0) / (fp64(l) * fp64((l + 1) ** 2 - 4)), 
            fp64(2* l + 1) * fp64(l * (l + 1)), 
            fp64((2 * l + 1) * 4), 
            fp64(l + 1) * fp64(l ** 2 - 4) 
        ], dtype=fp64); 
    #Label 150: 
    alpha[span(1, 4), span(1, l1_max)] = fp64(0e0); 
    beta[span(1, 2), span(1, l1_max)] = fp64(0e0); 
    d6 = fp64(0.25e0) * np.sqrt(fp64(6e0)); 
    for i in do(1, n_g): 
        gener(x[i], l1_max, p, coef, d6); 
        wi = w[i]; 
        (ff11, ff33, ff12, ff34) = f_ij[span(1, 4), i] * wi; 
        fp = ff11 + ff33; fm = ff11 - ff33; 
        p1l1 = p[1, span(1, l1_max)]; 
        p2l1 = p[2, span(1, l1_max)]; 
        p3l1 = p[3, span(1, l1_max)]; 
        p4l1 = p[4, span(1, l1_max)]; 
        alpha[1, span(1, l1_max)] += ff11 * p1l1
        alpha[2, span(1, l1_max)] += fp * p2l1; 
        alpha[3, span(1, l1_max)] += fm * p3l1; 
        alpha[4, span(1, l1_max)] += ff33 * p1l1; 
        beta[1, span(1, l1_max)] += ff12 * p4l1; 
        beta[2, span(1, l1_max)] += ff34 * p4l1; 
    #Label 300: 
    for l1 in do(1, l1_max): 
        cl = fp64(l1 - 1) + fp64(0.5e0); 
        l = l1; 
        alpha[1, l1] *= cl; 
        a_2 = alpha[2, l1] * cl * fp64(0.5e0); 
        a_3 = alpha[3, l1] * cl * fp64(0.5e0); 
        alpha[2, l1] = a_2 + a_3; 
        alpha[3, l1] = a_2 - a_3;   
        alpha[4, l1] *= cl; 
        beta[1, l1] *= cl; 
        beta[2, l1] *= -cl; 
        if np.abs(alpha[1, l1] < d_delt): 
            break; 
    #Label 380: 
    l1_max = l; 
    #用局部变量修改result中的特定属性值
    result.cross_sectn_extinction = c_ext; 
    result.cross_sectn_scattering = c_sca; 
    result.single_scattering_albedo = alb; 
    result.asymmetry_param = alpha[1, 2] / fp64(3); 
    sys.stderr.write(
        "Calculation finished. \n"
    ); sys.stderr.flush(); 
    return l1_max; 