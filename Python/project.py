# import signal/control system packages
import signal
import sys
# import time package
import time

# define a DC motor dynamics 
class DCMotor:
    def __init__(self, J, b, Kt, Ke, R, L):
        self.J = J  # moment of inertia
        self.b = b  # damping coefficient
        self.Kt = Kt  # torque constant
        self.Ke = Ke  # back EMF constant
        self.R = R  # resistance
        self.L = L  # inductance

    def dynamics(self, state, voltage):
        theta_dot = state[1]
        i_dot = (voltage - self.Ke * theta_dot - self.R * state[2]) / self.L
        theta_ddot = (self.Kt * state[2] - self.b * theta_dot) / self.J
        return [theta_dot, theta_ddot, i_dot], i_dot    
    
# define a control system for the DC motor
class DCMotorController:
    def __init__(self, motor, Kp, Ki, Kd):
        self.motor = motor
        self.Kp = Kp  # proportional gain
        self.Ki = Ki  # integral gain
        self.Kd = Kd  # derivative gain
        self.integral_error = 0
        self.previous_error = 0

    def control(self, desired_theta, current_state):
        error = desired_theta - current_state[0]
        self.integral_error += error
        derivative_error = error - self.previous_error
        voltage = (self.Kp * error + 
                   self.Ki * self.integral_error + 
                   self.Kd * derivative_error)
        self.previous_error = error
        return voltage
    
# signal handler for graceful shutdown
def signal_handler(sig, frame):
    print('Shutting down gracefully...')
    sys.exit(0) 

# register signal handler
signal.signal(signal.SIGINT, signal_handler)
signal.signal(signal.SIGTERM, signal_handler)

# main function to run the control loop
def main():
    # create a DC motor instance
    motor = DCMotor(J=0.01, b=0.1, Kt=0.01, Ke=0.01, R=1, L=0.5)
    # create a controller instance
    controller = DCMotorController(motor, Kp=2.0, Ki=0.5, Kd=0.1)
    
    # initial state [theta, theta_dot, current]
    state = [0, 0, 0]
    desired_theta = 1.0  # desired position in radians
    
    while True:
        voltage = controller.control(desired_theta, state)
        state_dot = motor.dynamics(state, voltage)
        # simple Euler integration for state update
        state = [state[i] + state_dot[i] * 0.01 for i in range(3)]
        print(f"Current State: {state}, Voltage: {voltage}")
        time.sleep(0.01)  # control loop at 100 Hz

if __name__ == "__main__":    
    main()    
