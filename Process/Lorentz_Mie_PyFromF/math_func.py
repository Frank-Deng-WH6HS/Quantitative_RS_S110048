#Lorentz-Mie散射代码主要数学计算函数. 

#参考文献: 
#[1]上海计算技术研究所. 电子计算机算法手册[M]. 上海: 上海教育出版社, 1982. 
#[2]喻文健. 数值分析与算法[M]. 3版. 北京: 清华大学出版社, 2020. 

from . import np; 
from numpy import float64 as fp64; 

#获取当前系统浮点数的机器精度
def machine_precision(): 
    eps = fp64(1e0); tol1 = np.inf; 
    while tol1 > fp64(1e0): 
        eps *= fp64(0.5e0); 
        tol1 = fp64(1e0) + eps; 
    return eps; 

#求一元函数零点的ZeroIn算法(结合了二分法, 线性插值和逆二次插值[1-2])
def zeroin(a_x: np.float64, b_x: np.float64, 
    f: callable, tol: np.float64) -> np.float64: 
    #Label 10: 获取浮点数的机器精度
    eps = machine_precision(); 
    #Label 15: 初始化
    a = fp64(a_x); b = fp64(b_x); 
    f_a = f(a); f_b = f(b); 
    tail_loop_lbl20 = True; 
    while True: 
        #Label 20: 调整区间端点
        if tail_loop_lbl20: 
            c = a; f_c = f_a; 
            d = b - a; e = d; 
        #Label 30: 交换a, b, c的值, 使得abs(f(b))始终较小
        if not (abs(f_a) >= abs(f_b)): 
            a = b; b = c; c = a; 
            f_a = f_b; f_b = f_c; f_c = f_a; 
        #Label 40: 计算过程中实际使用的精度
        tol1 = fp64(2e0) * eps * abs(b) + fp64(0.5e0) * tol; 
        xm = fp64(0.5e0) * (c - b); 
        #Label 44: 跳出循环的条件
        if abs(xm) <= tol1 or f_b == 0e0: 
            break; 
        #Label 45, 46: 
        if abs(e) < tol1 or abs(f_a) <= abs(f_b): 
            #Label 70: 
            d = xm; e = d; 
        else: 
            #Label 47: 
            if a != c: 
                #Label 50: 
                q = f_a / f_c; r = f_b / f_c; 
                s = f_b / f_a; 
                p = s * (fp64(2e0) * xm * q * (q - r) \
                    - (b - a) * (r - fp64(1e0))); 
                q = (q - fp64(1e0)) * (r - 1e0) * (s - 1e0); 
            else: 
                #Label 48: 
                s = f_b / f_a; 
                p = fp64(2e0) * xm * s; 
                q = fp64(1e0) - s; 
            #Label 60: 
            if p > fp64(0e0): 
                q = -q; 
            p = abs(p); 
            #Label 64: 
            if (fp64(2e0) * p >= fp64(3e0) * xm * q - abs(tol1 * q)) \
                or (p >= abs(fp64(0.5e0) * e * q)): 
                #Label 70: 
                d = xm; e = d; 
            else: 
                #Label 65: 
                e = d; d = p / q; 
        #label 80: 
        a = b; f_a = f_b; 
        if abs(d) > tol1: 
            b += d; 
        else: 
            b += np.copysign(tol1, xm); 
        f_b = f(b); 
        #Label 85: 
        tail_loop_lbl20 = bool(f_b * (f_c / abs(f_c)) > fp64(0e0)); 
    return b; 

#F函数
#在源代码中是被POWER(A, B, R1, R2)调用的一元函数(与ZEROIN合用)
#由于python传参方式的特殊性, 其用法改为r = f_trivar(a, b, r1)
#在后续的power函数中, f_trivar将通过lambda表达式重新构造一元函数f
def f_trivar(a: np.float64, b: np.float64, 
    r_1: np.float64) -> np.float64: 
    r_2 = (fp64(1e0) + b) * fp64(2e0) * a - r_1; 
    f = (r_2 - r_1) / np.log(r_2 / r_1) - a; 
    return f; 

#Power函数
#由于python传参方式的特殊性, 其用法改为(r_1, r_2) = power(a, b, r_1, r_2)
#返回的结果是包含两个np.float64的tuple
def power(a: np.float64, b: np.float64, 
    r_1: np.float64, r_2: np.float64) -> tuple: 
    '''Usage: (r_1, r_2) = power(a, b, r_1, r_2)'''
    f = lambda r_1: f_trivar(a, b, r_1); 
    a_x = fp64(1e-12); b_x = a; 
    r_min = zeroin(a_x, b_x, f, fp64(0e0)); 
    r_max = (fp64(1e0) + b) * 2 * a - r_min; 
    return (r_min, r_max); 

