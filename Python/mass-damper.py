import numpy as np 
import matplotlib.pyplot as plt 
from scipy.integrate import odeint


# # mass damper equation 
# # my''(t) + cy'(t) + ky(t) = f(t)

# # system variables 
# k = 1 # spring constant
# c = 1 # friction constant 
# m = 10 # mass 

# #input 
# def u(t):
#     return 100 * np.sin(2 * np.pi * 0.001 * t)

# # state-space model 
# A = np.array([[0, 1],
#               [-k/m, -c/m]])
# B = np.array([[0], [1/m]])
# C = np.array([[1, 0]])   # output = position y
# D = 0 

# # Initial conditions: y(0)=1, y'(0)=2
# x0 = np.array([1.0, 2.0])

# # Time vector
# t = np.linspace(0, 20, 1000)

# # Derivative function for odeint
# def dxdt(x, t, A, B):
#     u_t = u(t)
#     return A @ x + B.flatten() * u_t   # B.flatten() makes it a 1D array

# x = odeint(dxdt, x0, t, args=(A, B))

# # Extract outputs: y = C x
# y = C @ x.T   # shape (1, len(t)) -> we can flatten
# y = y.flatten()

# plt.figure()
# plt.plot(t, y, label='Position y(t)')
# plt.xlabel('Time (s)')
# plt.ylabel('y(t)')
# plt.title('Mass-Spring-Damper Response to Constant Force')
# plt.grid(True)
# plt.legend()
# plt.show()



# feedback control 
# DCS Lab 1
# Heating system 

T_ref = 20 # required temperature  
T_amb = 5 # ambient temperature 

def comparator(T_room):
    delT = T_ref - T_room 
    if delT>=0:
        return 1
    else:
        return 0 
    
def heatphysics(T_room):
    delT = T_room - T_amb
    step = 0.1
    if delT>0:
        T_room = T_room - step
    elif delT<0:
        T_room = T_room + step
    else:
        None
    
    return T_room 


def heater(C, T_i):
    step = 5
    if C==1:
        T_i = T_i + step 
    else:
        None 
    
    return T_i

def system(time):
    initial = T_amb + 1
    U = [0]
    T = [initial]
    T_i = initial
    for i in range(time):
        C = comparator(T_i)
        U.append(C)
        T_i = heater(C, T_i)
        T_i = heatphysics(T_i)
        T.append(T_i)

    return T, U

time = 100
T, U = system(time)
t = np.linspace(0, len(T), len(T))

plt.plot(t, T)
plt.step(t, U)
plt.grid()
plt.show()
