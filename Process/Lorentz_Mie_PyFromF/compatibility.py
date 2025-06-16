#Lorentz-Mie散射代码中兼容Fortran行文的部分调整. 
#不要轻易在本项目模块以外的开发工作中使用! 

from . import np; 
from numpy import float64 as fp64; 

#Fortran中的DO循环生成器
#DO I=START,END
#I从START开始每次循环加1, 直到END. END处仍循环最后一次(两端都包含)
#对比: Python中range(start, end)只包含start, 不包含end
def do(start: int, end: int): 
    idx = start; 
    while idx <= end: 
        yield idx; 
        idx += 1; 
        
#与DO类似的数组切片对象, 同时包含start和end
def span(start: int, end: int): 
    return slice(start, end + 1, 1); 
        
#Fortran中声明实数数组
#REAL*8 ARR(DIM) 或者 REAL*8 ARR(DIM1, DIM2, ...)
#由于Fortran中数组下标从1开始, 因此使用np.array模拟时, 需要对每个维度多加一行(列)
def real(*shape): 
    shape_arr = np.array(shape) + 1; 
    uninit_arr = np.full(shape_arr, np.nan, dtype=fp64); 
    return uninit_arr; 