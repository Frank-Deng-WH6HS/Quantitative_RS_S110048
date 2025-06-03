#Lorentz-Mie散射代码主要数学计算函数. 

#参考文献: 
#[1]上海计算技术研究所. 电子计算机算法手册[M]. 上海: 上海教育出版社, 1982. 
#[2]喻文健. 数值分析与算法[M]. 3版. 北京: 清华大学出版社, 2020. 

import math; 

#获取当前系统浮点数的机器精度
def machine_precision(): 
    eps = 1e0; tol1 = math.inf; 
    while tol1 > 1e0: 
        eps *= 0.5e0; 
        tol1 = 1e0 + eps; 
    return eps; 

#求一元函数零点的ZeroIn算法(结合了二分法, 线性插值和逆二次插值[1-2])
def zeroin(a_x: float, b_x: float, f: callable, tol: float) -> float: 
    #Label 10: 获取浮点数的机器精度
    eps = machine_precision(); 
    #Label 15: 初始化
    a = a_x; b = b_x; 
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
        tol1 = 2e0 * eps * abs(b) + 0.5e0 * tol; 
        xm = 0.5e0 * (c - b); 
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
                p = s * (2e0 * xm * q * (q - r) - (b - a) * (r - 1e0)); 
                q = (q - 1e0) * (r - 1e0) * (s - 1e0); 
            else: 
                #Label 48: 
                s = f_b / f_a; 
                p = 2e0 * xm * s; 
                q = 1e0 - s; 
            #Label 60: 
            if p > 0: 
                q = -q; 
            p = abs(p); 
            #Label 64: 
            if (2e0 * p >= 3e0 * xm * q - abs(tol1 * q)) or \
                (p >= abs(0.5e0 * e * q)): 
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
            b += math.copysign(tol1, xm); 
        f_b = f(b); 
        #Label 85: 
        tail_loop_lbl20 = bool(f_b * (f_c / abs(f_c)) > 0e0); 
    return b; 