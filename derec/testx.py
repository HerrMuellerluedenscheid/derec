from pyrocko.trace import *
import numpy as num
import ctypes
import time


fun=ctypes.cdll.LoadLibrary('./liblxnorm.so')
size = 1500
norm = 1
u = num.random.uniform(-1e21, 1e21,size)
v = num.random.uniform(-1e21, 1e21,size)

loops = 10000

t1 = time.time()

for i in range(loops):
    size_c = ctypes.c_int(size)
    norm_c = ctypes.c_double(norm)
    u_c = u.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    v_c = v.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    fun.lxnorm_m.restype = ctypes.c_double
    fun.lxnorm_n.restype = ctypes.c_double
    m = fun.lxnorm_m(v_c, u_c, norm_c, size_c )
    n = fun.lxnorm_n(v_c, norm_c, size_c )

t2 = time.time()
print t2-t1
print m, n

t1 = time.time()
for i in range(loops):
    if norm ==1:
        m,n = L1_norm(u,v)

t2 = time.time()
print t2-t1
print m, n

for i in range(loops):
    m,n = Lx_norm(u,v, 1)

t2 = time.time()
print t2-t1
print m, n
