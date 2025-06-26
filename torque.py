import numpy as np
import pandas as pd

# 1) load and window between t=3.667 and t=4.1
df = pd.read_csv('data/coordinate_speeds_clubgolfswing.csv')
df = df[(df['time'] >= 3.667) & (df['time'] <= 4.1)]
k = "arm_rot_l"
# 2) convert absolute left-shoulder angle (deg → rad)
theta = np.deg2rad(df['arm_flex_l'].values)
theta = np.deg2rad(df['arm_flex_r'].values)
theta = np.deg2rad(df[k].values)
t = df['time'].values

# 3) compute derivatives
omega = np.gradient(theta, t)  # rad/s
# alpha = np.gradient(omega, t)  # rad/s²

# 4) pick inertia values
m_arm, L_arm = 2.7, 0.34     # kg, m
I_arm = (1/3) * m_arm * L_arm**2

m_club, R = 0.2, 1.0          # kg, m
I_club = m_club * R**2

I_total = I_arm + I_club

# 5) compute torque and peak
tau = I_total * omega
peak_tau = np.max(np.abs(tau))

print(f"Estimated peak shoulder torque (t ≤ 4.1 s): {k} : {peak_tau:.1f} N·m")
