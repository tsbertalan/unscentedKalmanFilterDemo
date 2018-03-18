import matplotlib.pyplot as plt
import numpy as np
from pprint import pprint 

datafolder = '/home/tsbertalan/Documents/Projects/Extended Kalman Filter/data/'
datafilename = 'obj_pose-laser-radar-synthetic-input.txt'

with open(datafolder + datafilename, 'r') as f:
    d = f.readlines()

events = []
for l in d:
    l = l.split()
    if l[0] == 'L': # LIDAR
        keys = [k.strip() for k in '''x_measured, y_measured,
        timestamp, x_groundtruth, y_groundtruth, 
        vx_groundtruth, vy_groundtruth, yaw_groundtruth, yawrate_groundtruth'''.replace(',', '').split()
        ]
        
    else: # RADAR
        keys = [k.strip() for k in '''rho_measured, phi_measured, rhodot_measured,
        timestamp, x_groundtruth, y_groundtruth, 
        vx_groundtruth, vy_groundtruth, yaw_groundtruth, yawrate_groundtruth'''.replace(',', '').split()
        ]

    ev = {
        k: float(x)
        for k, x in zip(keys, l[1:len(keys)+1])
    }

    ev['sensor_type'] = l[0]

    events.append(ev)

pprint(events[-2:])


def gv(key, sensor='any'):
    if sensor in 'LR':
        cond = lambda ev: ev['sensor_type'] == sensor
    else:
        cond = lambda ev: True
    return np.array([
        ev[key] for ev in events
        if cond(ev)
    ])


# Plot trajectory.
fig, ax = plt.subplots()
T = gv('timestamp')
T -= T.min()
T /=  1e6
X = gv('x_groundtruth')
Y = gv('y_groundtruth')
ax.plot(X, Y, label='truth')

lasX = gv('x_measured', sensor='L')
lasY = gv('y_measured', sensor='L')
ax.scatter(lasX, lasY, label='LIDAR', color='black')

radPhi = gv('phi_measured', sensor='R')
radRho = gv('rho_measured', sensor='R')
radX = np.cos(radPhi) * radRho
radY = np.sin(radPhi) * radRho
ax.scatter(radX, radY, label='RADAR', color='red')

ax.legend()
fig.savefig('y_vs_x.png')

# Plot accelerations.
fig, ax = plt.subplots()
a2 = ax.twinx()

DT = np.diff(T)
V = np.stack((gv('vx_groundtruth'), gv('vy_groundtruth')))
Vmag = np.linalg.norm(V, axis=0)
Amag = np.diff(Vmag) / DT
a2.plot(
    T[:-1], Amag, 
    label=r'$\Delta  v / \Delta t$, $\sigma_a = %.2g$' % Amag.std()
)

Pddmag = np.diff(gv('yawrate_groundtruth')) / DT
a2.plot(
    T[:-1], Pddmag, 
    label=r'$\Delta  \dot\phi / \Delta t$, $\sigma_{\ddot\phi} = %.2g$' % Pddmag.std()
    )

a2.set_ylabel('accelerations')
a2.legend()

ax.plot(T, X, label='$x$', color='red')
ax.plot(T, Y, label='$y$', color='blue')
ax.set_xlabel('$t$')
ax.set_ylabel('$x$ (red) and $y$ (blue)')

fig.tight_layout()
fig.savefig('data_vs_time.png')

# Show plots.
plt.show()