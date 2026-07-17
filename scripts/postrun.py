import sarracen as sar
import matplotlib.pyplot as plt

print("Hello World!")

sdf = sar.read_csv("hydro32.csv")
sdf.params['hfact'] = 1.4
sdf.params['mass'] = 3.0517578125e-05
sdf.calc_density()
sdf.render(target="rho")

plt.show()
