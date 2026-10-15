from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import sarracen as sar

mass = 4.040675383339015e-06

def to_sarracen(file_location: str):
    sdf: sar.SarracenDataFrame = sar.read_csv(file_location)
    sdf.params["mass"] = mass
    sdf.params["hfact"] = 1.2
    sdf.calc_density()
    return sdf

config = "../kh_00000.toml"
snapshot_folder = Path("../output").resolve()
snapshots = [to_sarracen("../kh_00000.csv")] + [to_sarracen(path) for path in sorted(snapshot_folder.glob(r"kh_?????.csv"))]
densities = [sar.interpolate_3d_proj(snapshot, "rho", x_pixels=500, y_pixels=500) for snapshot in snapshots]

print(snapshot_folder)

fig, ax = plt.subplots(figsize=(9, 9))

density_image = ax.imshow(densities[0], origin="lower", extent=[0, 1, 0, 1])
fig.colorbar(density_image, ax=ax, label=r"$\rho$")
ax.set(title="Density", xlabel="x", ylabel="y")

title = fig.suptitle("")

def draw(frame):
    density_image.set_data(densities[frame])
    title.set_text(f"Kelvin-Helmholtz: step {frame} of {len(snapshots) - 1}")
    return density_image, title

animation = FuncAnimation(fig, draw, frames=len(snapshots), interval=100)
animation.save("../output/output_animation.gif", writer=PillowWriter(fps=10))
print(f"Animation Saved!")

plt.show()
