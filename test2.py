import numpy as np
from Parameters import *
from Initial_conditions import *

def temperature_root(p):
    t_array = np.roots(p)
    print(t_array)
    t = T_amb
    for i in range(len(t_array)):
        if isinstance(t_array[i], complex):
            continue
        elif t_array[i] <= 0:
            continue
        else:
            t = t_array[0]
    return t


w = -0.02
T = (temperature_root([r / (1 / w - b), -P_0i, -a / ((1 / w) * (b + 1 / w))])) ** 2
print (T)


