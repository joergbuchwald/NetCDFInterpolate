import numpy as np
import NetCDFInterpolate
import matplotlib.pyplot as plt


f = NetCDFInterpolate.NetCDFInterpolate("pointheatsource_linear-mesh.h5",dim=2)

print(f"point field names: {f.get_point_field_names()}")
print(f"data arrays: {f.data_arrays}")

pts = {"pt0": (0.2,0.1,0.0)}

# extract a time series at given point
temperature = f.read_time_series("temperature",pts=pts)
plt.plot(f.times, temperature["pt0"])
plt.ylabel("$T$ / K")
plt.xlabel("$t$ / s")

plt.show()

# define an array
x = np.linspace(start=0.0, stop=10.0, num=100)
r =  [(i,0,0) for i in x]

# extract data along given point set array and time step
p = f.get_set_data("pressure", 1, pointsetarray=r)

# extract data along given point set array and time
p2 = f.read_set_data("pressure", 45000, pointsetarray=r)

plt.plot(x,p)
plt.plot(x,p2)
plt.show()
