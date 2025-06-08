import sympy as sp
import numpy as np

# Time and dynamic variables
t = sp.Symbol('t')
theta = sp.Function('theta')(t)
phi = sp.Function('phi')(t)

theta_dot = theta.diff(t)
phi_dot = phi.diff(t)
theta_ddot = theta.diff(t, t)
phi_ddot = phi.diff(t, t)

# Physical parameters
m1, m2 = sp.symbols('m1 m2')
L1, L2 = sp.symbols('L1 L2')
I1, I2 = sp.symbols('I1 I2')
g = sp.Symbol('g')
r_shoulder, r_wrist = sp.symbols('r_shoulder r_wrist')

# Define non-inertial force terms
N1 = (-m2 * L1 * L2 * phi_dot**2 * sp.sin(phi)
      + m1 * g * (L1 / 2) * sp.sin(theta)
      + m2 * g * (L1 * sp.sin(theta) + (L2 / 2) * sp.sin(theta + phi)))

N2 = (m2 * L1 * L2 * theta_dot**2 * sp.sin(phi)
      + m2 * g * (L2 / 2) * sp.sin(theta + phi))

# Define full shoulder muscle force equation
F_shoulder = ((I1 + m2 * L1**2 + I2 + m2 * L1 * L2 * sp.cos(phi)) * theta_ddot
              + (I2 + m2 * L1 * L2 * sp.cos(phi)) * phi_ddot
              + N1) / r_shoulder

# Define full wrist muscle force equation
F_wrist = ((I2 + m2 * L1 * L2 * sp.cos(phi)) * theta_ddot
           + I2 * phi_ddot
           + N2) / r_wrist

# Create Python-callable functions
F_shoulder_func = sp.lambdify((theta, phi, theta_dot, phi_dot, theta_ddot, phi_ddot,
                               m1, m2, L1, L2, I1, I2, g, r_shoulder), F_shoulder, 'numpy')

F_wrist_func = sp.lambdify((theta, phi, theta_dot, phi_dot, theta_ddot, phi_ddot,
                            m1, m2, L1, L2, I1, I2, g, r_wrist), F_wrist, 'numpy')

"Shoulder and wrist force functions created for simulation or symbolic use."

# --- 1. Numeric parameters, using *_val so we don’t shadow the Sympy symbols ---
m1_val        = 4.0       # kg, upper-arm + forearm + hand
m2_val        = 0.35      # kg, driver
L1_val        = 0.70      # m, shoulder→wrist
L2_val        = 1.10      # m, club length
I1_val        = (1.0/3.0) * m1_val * L1_val**2
I2_val        = (1.0/3.0) * m2_val * L2_val**2
g_val         = 9.81      # m/s²
r_shoulder_val = 0.05     # m, moment arm at shoulder
r_wrist_val    = 0.03     # m, moment arm at wrist

# --- 2. Initial‐condition values for the very start of the downswing ---
theta0       = 1.20      # rad
phi0         = -1.57     # rad
theta_dot0   = 0.0       # rad/s
phi_dot0     = 0.0       # rad/s
theta_ddot0  = 25.0      # rad/s²
phi_ddot0    = 0.0       # rad/s²

# --- 3. Call your lambdified functions with these new names ---
# (Assuming F_shoulder_func and F_wrist_func have already been defined as in your previous code.)

F_shoulder_val = F_shoulder_func(
    theta0,      # θ
    phi0,        # φ
    theta_dot0,  # θ̇
    phi_dot0,    # φ̇
    theta_ddot0, # θ̈
    phi_ddot0,   # φ̈
    m1_val,
    m2_val,
    L1_val,
    L2_val,
    I1_val,
    I2_val,
    g_val,
    r_shoulder_val
)

F_wrist_val = F_wrist_func(
    theta0,
    phi0,
    theta_dot0,
    phi_dot0,
    theta_ddot0,
    phi_ddot0,
    m1_val,
    m2_val,
    L1_val,
    L2_val,
    I1_val,
    I2_val,
    g_val,
    r_wrist_val
)

print(f"Shoulder force (N) at downswing start: {F_shoulder_val:.1f}")
print(f"Wrist   force (N) at downswing start: {F_wrist_val:.1f}")

F_shoulder_val*r_shoulder_val, F_wrist_val*r_wrist_val  # Return the torques at shoulder and wrist