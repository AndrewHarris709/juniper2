import sarracen as sar
import matplotlib.pyplot as plt
import pandas as pd

acc = pd.read_csv("./output/kh_00000_accenergy.jrep", header=None, names=['i', 'u', 'ax', 'ay', 'az'])

print("Hello World!")

sdf = sar.read_csv("kh_00000.csv")
sdf.params['hfact'] = 1.4
sdf.params['mass'] = 4.040675383339015e-06
sdf.calc_density()

sdf[['ax', 'ay', 'az']] = acc[['ax', 'ay', 'az']]

fig, ax = plt.subplots()

sdf.render(target="rho", ax=ax)
sdf.arrowplot(target=['ax', 'ay', 'az'], ax=ax, x_arrows=30, color="blue")

print(sdf[['x', 'y', 'z']].describe())

plt.show()
