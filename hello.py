import numpy as np
import polars as pl
import matplotlib.pyplot as plt
import inspect
import re 

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

def f2(t,y):
    return y**2 - y**3 

t0 = 0
y0 = np.array([.2])
tf = 10.01
dt = .1

ts, ys, dys = euler(f2, t0, y0, tf, dt)
df2 = pl.DataFrame(
    {
        "t": ts,
        "y": ys[:, 0].tolist(),
        "dy/dt": dys[:, 0].tolist(),
    }
)
df2.write_csv("data/problem_9_table.csv")

import matplotlib.pyplot as plt
import inspect
import re

def plot_euler_solution(ts, ys, f=None, ode_str=None):
    """
    Plot y vs. t and auto‑generate a title.
    
    ts      – array of time‑points
    ys      – array of solution values (shape [len(ts), state_dim])
    f       – (optional) the dy/dt function you passed to your solver
    ode_str – (optional) a string to describe the ODE (e.g. "y^2 - 4*t")
    """
    # plot the first component
    plt.plot(ts, ys[:,0], marker="o", linestyle="-")
    
    # build title
    if ode_str:
        expr = ode_str
    elif f:
        # try to extract the “return …” expression from the source
        try:
            src = inspect.getsource(f)
            m   = re.search(r"return\s+(.+)", src)
            expr = m.group(1).strip()
        except (OSError, AttributeError):
            # fallback to the function’s name
            expr = f.__name__
    else:
        expr = "dy/dt = ?"
    
    plt.title(f"Euler’s Method Approximation\nfor $d y/dt = {expr}$")
    plt.xlabel("Time $t$")
    plt.ylabel("Solution $y$")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

plot_euler_solution(ts, ys, ode_str="y**2 - 4*t")
