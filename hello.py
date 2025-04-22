import numpy as np
import polars as pl
import matplotlib.pyplot as plt


def euler(f, t0, y0, tf, dt):
    ts = np.arange(t0, tf, dt)
    ys = np.empty((ts.size, y0.size))
    dys = np.empty((ts.size, y0.size))
    ys[0] = y0
    for i in range(1, ts.size):
        dys[i - 1] = f(ts[i - 1], ys[i - 1])
        ys[i] = ys[i - 1] + dt * dys[i - 1]

    dys[-1] = f(ts[-1], ys[-1])
    return ts, ys, dys


# 3) dy/dt = y^2-4t
def f(t, y):
    return y**2 - 4 * t


t0 = 0
y0 = np.array([0.5])
tf = 2.1
dt = 0.25

ts, ys, dys = euler(f, t0, y0, tf, dt)
df = pl.DataFrame(
    {
        "t": ts,
        "y": ys[:, 0].tolist(),
        "dy/dt": dys[:, 0].tolist(),
    }
)

df.write_csv("problem_3_table.csv")

plt.plot(ts, ys[:,0], marker="o", linestyle="-")
plt.title("Euler’s Method Approximation\nfor $dy/dt = y^2 - 4t$")
plt.xlabel("Time $t$")
plt.ylabel("Solution $y$")
plt.grid(True)
plt.tight_layout()
plt.show()