import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

parser = argparse.ArgumentParser()
parser.add_argument('--save')
parser.add_argument('--file')
args = parser.parse_args()

if not sys.stdin.isatty():
    data = np.loadtxt(sys.stdin, delimiter=',')
elif args.file:
    data = np.loadtxt(args.file, delimiter=',')

t, b1x, b1y, b2x, b2y, b3x, b3y = data.T

fig, ax = plt.subplots(figsize=(8, 8))
ax.set_aspect('equal')
ax.grid(True, linestyle='--', alpha=0.6)

all_coords = data[:, 1:]
limit = np.max(np.abs(all_coords)) * 1.1
ax.set_xlim(-limit, limit)
ax.set_ylim(-limit, limit)

b1, = ax.plot([], [], 'ro', ms=8)
b2, = ax.plot([], [], 'bo', ms=8)
b3, = ax.plot([], [], 'go', ms=8)

ax.plot(b1x, b1y, color='lightgrey', alpha=0.6)
ax.plot(b2x, b2y, color='lightgrey', alpha=0.6)
ax.plot(b3x, b3y, color='lightgrey', alpha=0.6)

def update(i):
    b1.set_data([b1x[i]], [b1y[i]])
    b2.set_data([b2x[i]], [b2y[i]])
    b3.set_data([b3x[i]], [b3y[i]])

    return b1, b2, b3

ani = FuncAnimation(fig, update, frames=range(0, len(t), 24), interval=20, blit=True)

if args.save:
    ani.save(args.save, writer='pillow')
else:
    plt.show()