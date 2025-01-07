import sys
from math import *

success = True
max_error = [1.e-4, 2.e-3]
radius = (1. + sqrt(2)*1.e-3)


delta = []
method = {
    1: {"name": "smoothed_VOF", "L2": [], "L_inf": [], "max_L2": 1.e-3, "max_L_inf": .2, "min_order_L2": 1., "min_order_L_inf": 0.},
    2: {"name": "extended_heights", "L2": [], "L_inf": [], "max_L2": 3.e-4, "max_L_inf": 1.e-2, "min_order_L2": 2., "min_order_L_inf": 1.},
    3: {"name": "extact_LS", "L2": [], "L_inf": [], "max_L2": 2.e-7, "max_L_inf": 2.e-5, "min_order_L2": 3., "min_order_L_inf": 3.}
}
file = open("out",'r')
for line in file.readlines():
    try:
        key = int(line[0])
    except:
        key = -1
    if key == 0:
        delta.append(2*radius/float(line.split(" ")[1]))
    elif key in [1,2,3]:
        data = line.strip('\n').split(" ")
        method[key]["L2"].append(float(data[2]))
        method[key]["L_inf"].append(float(data[3]))
file.close()

delta_range = log(delta[0]) - log(delta[-1])
for k in range(1,4):
    order = {"L2": [], "L_inf": []}
    for norm in ["L2", "L_inf"]:
        order[norm] = (log(method[k][norm][0]) - log(method[k][norm][-1]))/delta_range
        if order[norm] < method[k]["min_order_"+norm]:
            success = False
        if method[k][norm][-1] > method[k]["max_"+norm]:
            success = False
    print(method[k]["name"]+" : order_L2="+str(round(order["L2"],1))+", order_L_inf="+str(round(order["L_inf"],1)))

file = open("log",'r')
for line in file.readlines():
    if "WARNING: maximum iterations reached in normal_scalar_extension" in line:
        success = False
file.close()

if (success):
    print(0)
else:
    print(1)
