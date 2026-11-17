import argparse
import re
import sys
import tomllib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import sarracen as sar
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Button, Slider

HFACT = 1.2
PIXELS = 500


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate an animated figure from the frames of a Juniper simulation.")
    parser.add_argument("config", type=Path,
                        help="simulation configuration file (.toml)")
    parser.add_argument("frames_dir", type=Path,
                        help="folder containing the simulation frames (<prefix>NNNNN.csv)")
    parser.add_argument("output", type=Path,
                        help="where to save the animation (.gif uses Pillow, other formats use ffmpeg)")
    parser.add_argument("-n", "--num-frames", type=int, default=None,
                        help="number of frames to draw (default: all simulation frames)")
    parser.add_argument("-t", "--target", default="rho",
                        help="particle quantity to render (default: rho)")
    parser.add_argument("-z", "--z-slice", type=float, default=None,
                        help="render a cross-section at this z value via interpolate_3d_cross "
                             "(default: column integration via interpolate_3d_proj)")
    parser.add_argument("--fps", type=int, default=10,
                        help="playback rate of the animation (default: 10)")
    parser.add_argument("--no-show", action="store_true",
                        help="save the animation without opening the interactive viewer")
    args = parser.parse_args(argv)

    if args.config.suffix != ".toml":
        parser.error(f"config must be a .toml file, got '{args.config}'")
    if not args.config.is_file():
        parser.error(f"config file '{args.config}' does not exist")
    if not args.frames_dir.is_dir():
        parser.error(f"frames folder '{args.frames_dir}' does not exist")
    if args.num_frames is not None and args.num_frames < 1:
        parser.error("--num-frames must be at least 1")
    if args.fps < 1:
        parser.error("--fps must be at least 1")
    return args


def load_config(path: Path) -> dict:
    with open(path, "rb") as file:
        return tomllib.load(file)["simconfig"]


def find_frames(config_path: Path, name: str, frames_dir: Path) -> list[Path]:
    """Initial conditions (``<name>.csv`` beside the config) followed by every numbered frame."""
    prefix = re.sub(r"\d+$", "", name)
    pattern = re.compile(re.escape(prefix) + r"(\d+)\.csv")

    numbered = {}
    for path in frames_dir.glob(f"{prefix}*.csv"):
        match = pattern.fullmatch(path.name)
        if match:
            numbered[int(match.group(1))] = path
    frames = [numbered[step] for step in sorted(numbered)]

    initial = config_path.parent / f"{name}.csv"
    if initial.is_file() and initial.name not in {frame.name for frame in frames}:
        frames.insert(0, initial)
    return frames


def render_frame(path: Path, args: argparse.Namespace, mass: float, xlim, ylim) -> np.ndarray:
    sdf: sar.SarracenDataFrame = sar.read_csv(path)
    sdf.params["mass"] = mass
    sdf.params["hfact"] = HFACT
    sdf.calc_density()

    if args.target not in sdf.columns:
        sys.exit(f"\nerror: target '{args.target}' not found in {path.name}; "
                 f"available: {', '.join(sdf.columns)}")

    common = dict(x_pixels=PIXELS, y_pixels=PIXELS, xlim=xlim, ylim=ylim)
    if args.z_slice is None:
        image = sar.interpolate_3d_proj(sdf, args.target, **common)
    else:
        image = sar.interpolate_3d_cross(sdf, args.target, z_slice=args.z_slice, **common)
    return image.astype(np.float32)


def build_figure(images, frames, args, name, extent, figsize):
    """Figure with the image, colourbar and title; returns it with a ``draw(index)`` callback."""
    fig, ax = plt.subplots(figsize=figsize)

    vmin = min(float(np.nanmin(image)) for image in images)
    vmax = max(float(np.nanmax(image)) for image in images)
    artist = ax.imshow(images[0], origin="lower", extent=extent, vmin=vmin, vmax=vmax)

    target_label = rf"$\rho$" if args.target == "rho" else args.target
    if args.z_slice is None:
        mode = "column integrated"
        colorbar_label = rf"$\int$ {target_label} dz"
    else:
        mode = f"cross-section at z = {args.z_slice:g}"
        colorbar_label = target_label
    fig.colorbar(artist, ax=ax, label=colorbar_label)
    ax.set(title=f"{args.target} ({mode})", xlabel="x", ylabel="y")
    title = fig.suptitle("")

    def draw(index):
        artist.set_data(images[index])
        title.set_text(f"{name}: frame {index} of {len(images) - 1} ({frames[index].name})")
        return artist, title

    draw(0)
    return fig, draw


class Player:
    """Interactive playback: step back, play, pause and step forward buttons plus a frame slider."""

    def __init__(self, fig, draw, frame_count, fps):
        self.fig = fig
        self.draw = draw
        self.frame_count = frame_count
        self.index = 0

        fig.subplots_adjust(bottom=0.18)
        slider_ax = fig.add_axes((0.15, 0.09, 0.7, 0.03))
        self.slider = Slider(slider_ax, "Frame", 0, max(frame_count - 1, 1),
                             valinit=0, valstep=1, valfmt="%d")
        self.slider.on_changed(self._on_slider)

        self.buttons = {}
        for i, (label, callback) in enumerate([("◀", self.step_back), ("Play", self.play),
                                               ("Pause", self.pause), ("▶", self.step_forward)]):
            button = Button(fig.add_axes((0.3 + i * 0.105, 0.02, 0.095, 0.05)), label)
            button.on_clicked(callback)
            self.buttons[label] = button

        self.animation = FuncAnimation(fig, self._tick, init_func=lambda: self.draw(self.index),
                                       interval=1000 / fps, cache_frame_data=False)

    def _tick(self, _):
        self.show((self.index + 1) % self.frame_count)

    def _on_slider(self, value):
        self.index = min(int(value), self.frame_count - 1)
        self.draw(self.index)
        self.fig.canvas.draw_idle()

    def show(self, index):
        self.slider.set_val(index)

    def play(self, _=None):
        self.animation.resume()

    def pause(self, _=None):
        self.animation.pause()

    def step_back(self, _=None):
        self.pause()
        self.show((self.index - 1) % self.frame_count)

    def step_forward(self, _=None):
        self.pause()
        self.show((self.index + 1) % self.frame_count)


def main(argv=None):
    args = parse_args(argv)
    config = load_config(args.config)
    name = config["name"]
    (x1, x2), (y1, y2), (z1, z2) = config["limits"]

    if args.z_slice is not None and not z1 <= args.z_slice <= z2:
        print(f"warning: z-slice {args.z_slice} is outside the simulation limits [{z1}, {z2}]",
              file=sys.stderr)

    frames = find_frames(args.config, name, args.frames_dir)
    if not frames:
        sys.exit(f"error: no frames matching '{name}' found in {args.frames_dir}")
    if args.num_frames is not None:
        frames = frames[:args.num_frames]

    images = []
    for i, path in enumerate(frames):
        print(f"\rInterpolating frame {i + 1}/{len(frames)} ({path.name})", end="", flush=True)
        images.append(render_frame(path, args, config["mass"], (x1, x2), (y1, y2)))
    print()

    extent = (x1, x2, y1, y2)

    fig, draw = build_figure(images, frames, args, name, extent, figsize=(9, 9))
    animation = FuncAnimation(fig, draw, frames=len(images), interval=1000 / args.fps)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    animation.save(args.output, writer="pillow" if args.output.suffix.lower() == ".gif" else None,
                   fps=args.fps)
    plt.close(fig)
    print(f"Animation saved to {args.output.resolve()}")

    if not args.no_show:
        fig, draw = build_figure(images, frames, args, name, extent, figsize=(9, 10))
        player = Player(fig, draw, len(images), args.fps)  # noqa: F841 (keeps widgets alive)
        plt.show()


if __name__ == "__main__":
    main()
