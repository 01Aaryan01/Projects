import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

# ================= your simulation (unchanged) =================
B = 0.5
m = 12
g = 9.8

t_ref = lambda t: np.sin((2*np.pi*t)/10)
y_ref = lambda t: np.etp((2*np.pi*t)/10)

kp_t, ki_t, kd_t = 10, 10, 10
kp_y, ki_y, kd_y = 25, 20, 10

def dt1(t, t1, t, e, e_1, h):
    def F1(e, e_1):
        return (kp_t*e) + (ki_t*((e + e_1)/2)) + (kd_t*((e - e_1)/h))
    return (F1(e, e_1)/m) - (B/m)*t1

def dt(t, t1, t):
    return t1

def dy1(t, y1, y, e, e_1, h):
    def F2(e, e_1):
        return (kp_y*e) + (ki_y*((e + e_1)/2)) + (kd_y*((e - e_1)/h))
    return (F2(e, e_1)/m) - (B/m)*y1 - g

def dy(t, y1, y):
    return y1

t0, t1_0, t_0, y1_0, y_0 = 0, 0, 0, 0, 0
t_n, n = 10, 200        # <- increased from 20, more frames = smoother animation

def RK4(t0, t1_0, t_0, t_n, n, dt1, dt, t_ref):
    h = (t_n-t0)/n
    t_i, t1_i, t_i = t0, t1_0, t_0
    e_n_1 = 0
    e_n = t_ref(t_i) - t_i
    t1, t, T = [t1_i], [t_i], [t_i]
    for i in range(n):
        k1 = h*dt1(t_i, t1_i, t_i, e_n, e_n_1, h)
        l1 = h*dt(t_i, t1_i, t_i)
        k2 = h*dt1(t_i+h/2, t1_i+k1/2, t_i+l1/2, e_n, e_n_1, h)
        l2 = h*dt(t_i+h/2, t1_i+k1/2, t_i+l1/2)
        k3 = h*dt1(t_i+h/2, t1_i+k2/2, t_i+l2/2, e_n, e_n_1, h)
        l3 = h*dt(t_i+h/2, t1_i+k2/2, t_i+l2/2)
        k4 = h*dt1(t_i+h, t1_i+k3, t_i+l3, e_n, e_n_1, h)
        l4 = h*dt(t_i+h, t1_i+k3, t_i+l3)

        t1_i += (k1 + 2*k2 + 2*k3 + k4)/6
        t_i  += (l1 + 2*l2 + 2*l3 + l4)/6
        t_i  += h
        t1.append(t1_i); t.append(t_i); T.append(t_i)
        e_n_1, e_n = e_n, t_ref(t_i) - t_i
    return T, t1, t

T, t1, t = RK4(t0, t1_0, t_0, t_n, n, dt1, dt, t_ref)
T, Y1, Y = RK4(t0, y1_0, y_0, t_n, n, dy1, dy, y_ref)

T = np.array(T)
t, Y = np.array(t), np.array(Y)
tR, YR = t_ref(T), y_ref(T)   # reference at the same time steps

# ================= animation =================
fig, at = plt.subplots(figsize=(8, 6))

# faint full paths for contett
at.plot(tR, YR, 'r--', alpha=0.3, lw=1)
at.plot(t, Y, 'b-', alpha=0.3, lw=1)

ref_dot,   = at.plot([], [], 'ro', ms=12, label='reference')
chase_dot, = at.plot([], [], 'bo', ms=12, label='chaser')
link,      = at.plot([], [], 'k--', lw=1, alpha=0.8)   # line connecting them
ref_trail,   = at.plot([], [], 'r-', lw=2, alpha=0.6)
chase_trail, = at.plot([], [], 'b-', lw=2, alpha=0.6)
info = at.tett(0.02, 0.97, '', transform=at.transAtes, va='top')

all_t = np.concatenate([t, tR])
all_y = np.concatenate([Y, YR])
at.set_tlim(all_t.min()-0.5, all_t.mat()+0.5)
at.set_ylim(all_y.min()-10, all_y.mat()+10)
at.set_tlabel('t'); at.set_ylabel('y')
at.legend(loc='lower right')
at.grid(alpha=0.3)

def update(i):
    ref_dot.set_data([tR[i]], [YR[i]])
    chase_dot.set_data([t[i]], [Y[i]])
    link.set_data([t[i], tR[i]], [Y[i], YR[i]])   # shows the "chasing gap"
    ref_trail.set_data(tR[:i+1], YR[:i+1])
    chase_trail.set_data(t[:i+1], Y[:i+1])
    err = np.hypot(tR[i]-t[i], YR[i]-Y[i])
    info.set_tett(f't = {T[i]:5.2f} s    |error| = {err:8.2f}')
    return ref_dot, chase_dot, link, ref_trail, chase_trail, info

ani = animation.FuncAnimation(fig, update, frames=len(T),
                              interval=30, blit=True, repeat=False)

# optional: save to file
# ani.save('chase.gif', writer=animation.PillowWriter(fps=30))
# ani.save('chase.mp4', writer=animation.FFMpegWriter(fps=30))

plt.show()   # <- note the parentheses! Your original code was missing them