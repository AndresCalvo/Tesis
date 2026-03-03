"""
Animated 3-D surface plot of |A(t, r)|² read from field.dat.

The file is written by Fortran as blocks of (t, r, |A|²) separated by
blank lines – one block per time-step.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 – needed for 3D projection


# ── 1. Load data ─────────────────────────────────────────────────────────────
# numpy's loadtxt with unpack handles Fortran blank-line separators poorly,
# so we parse manually for robustness.

print("Loading field.dat …")
# Fast load: genfromtxt skips blank lines and returns (N, 3)
raw = np.loadtxt("field.dat")          # shape (Nt*Nr, 3)
t_col, r_col, f_col = raw[:, 0], raw[:, 1], raw[:, 2]

# Detect Nr: count consecutive rows that share the first t value
Nr = 1
while Nr < len(t_col) and t_col[Nr] == t_col[0]:
    Nr += 1
Nt = len(t_col) // Nr

t_vals = t_col[::Nr][:Nt]
r_vals = r_col[:Nr]
field  = f_col[: Nt * Nr].reshape(Nt, Nr)

print(f"  Nt = {Nt},  Nr = {Nr}")
print(f"  t ∈ [{t_vals[0]:.4g}, {t_vals[-1]:.4g}]")
print(f"  r ∈ [{r_vals[0]:.4g}, {r_vals[-1]:.4g}]")
print(f"  |A|² ∈ [{field.min():.4g}, {field.max():.4g}]")

# ── 2. Build meshgrid ────────────────────────────────────────────────────────
R, T = np.meshgrid(r_vals, t_vals)

# ── 3. Decide how many frames to show ────────────────────────────────────────
# If there are too many time‑steps the animation would be very slow to render;
# subsample to at most ~200 frames.
MAX_FRAMES = 200
if Nt > MAX_FRAMES:
    frame_indices = np.linspace(0, Nt - 1, MAX_FRAMES, dtype=int)
else:
    frame_indices = np.arange(Nt)
n_frames = len(frame_indices)

# ── 4. Set up figure ─────────────────────────────────────────────────────────
fig = plt.figure(figsize=(10, 7))
ax  = fig.add_subplot(111, projection="3d")

# Global z‑limits so the camera doesn't jump each frame
z_max = field.max()
z_min = field.min()

# Initial surface (will be replaced each frame)
surf = [ax.plot_surface(R[0:1, :], T[0:1, :], field[0:1, :],
                         cmap=cm.inferno, alpha=0.0)]

ax.set_xlabel("r")
ax.set_ylabel("t")
ax.set_zlabel("|A|²")
ax.set_zlim(z_min, z_max * 1.05 if z_max > 0 else 1)
title = ax.set_title("")


def update(frame_num):
    """Draw the surface up to the current time index."""
    idx = frame_indices[frame_num]

    # Remove the old surface
    surf[0].remove()

    # Plot the surface as a 2-D slice: all t up to current time vs r
    # For a "growing" surface use t_vals[:idx+1] × r_vals
    # Or for a single-slice "wave front" animation, show only the current
    # profile as a wireframe and the accumulated surface behind it.
    end = idx + 1
    surf[0] = ax.plot_surface(
        R[:end, :], T[:end, :], field[:end, :],
        cmap=cm.inferno,
        rstride=max(1, end // 60),   # keep triangle count reasonable
        cstride=max(1, Nr // 60),
        antialiased=False,
    )
    title.set_text(f"t = {t_vals[idx]:.4g}  (frame {frame_num+1}/{n_frames})")
    return (surf[0],)


print(f"Animating {n_frames} frames …")
ani = FuncAnimation(fig, update, frames=n_frames,
                    interval=80, blit=False, repeat=True)

plt.tight_layout()
plt.show()
