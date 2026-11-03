// demo_simulation.cpp - Complete DC Motor Control Simulation
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using namespace std;

// Motor Parameters (First-order model)
struct Motor {
    double K = 0.05;      // Gain (rad/s per volt)
    double tau = 0.1;     // Time constant (seconds)
    double speed = 0.0;   // Current speed (rad/s)
    double max_voltage = 12.0; // Maximum input voltage
};

// PI Controller with Anti-Windup
struct PIController {
    double Kp = 0.8;      // Proportional gain
    double Ki = 2.0;      // Integral gain
    double integral = 0.0; // Integral term
    double max_output = 12.0; // Maximum output (volts)
    double min_output = 0.0;  // Minimum output (volts)
    
    double calculate(double error, double dt) {
        // Calculate integral with anti-windup
        double new_integral = integral + error * dt;
        
        // Calculate output without limits
        double output = (Kp * error) + (Ki * new_integral);
        
        // Apply anti-windup
        if (output > max_output) {
            output = max_output;
        } else if (output < min_output) {
            output = min_output;
        } else {
            // Only update integral if output is within limits
            integral = new_integral;
        }
        
        return output;
    }
    
    void reset() {
        integral = 0.0;
    }
};

// Simulate motor response (First-order system)
void simulateMotor(Motor &motor, double voltage, double dt) {
    // Motor dynamics: τ * ds/dt + s = K * V
    // Discrete approximation: s[k+1] = s[k] + (K*V - s[k]) * dt/τ
    
    // Limit input voltage
    if (voltage > motor.max_voltage) voltage = motor.max_voltage;
    if (voltage < 0) voltage = 0;
    
    // Update motor speed
    motor.speed = motor.speed + (motor.K * voltage - motor.speed) * (dt / motor.tau);
    
    // Add some realistic noise (optional)
    motor.speed += ((rand() % 100) / 100.0 - 0.5) * 0.1;
}

// Convert RPM to rad/s
double rpmToRadPerSec(double rpm) {
    return rpm * (2 * M_PI / 60);
}

// Convert rad/s to RPM
double radPerSecToRpm(double rad_s) {
    return rad_s * (60 / (2 * M_PI));
}

int main() {
    cout << "DC Motor Speed Control Simulation\n";
    cout << "================================\n\n";
    
    // Initialize motor and controller
    Motor motor;
    PIController controller;
    
    // Simulation parameters
    double dt = 0.01;           // Control loop period (10 ms)
    double simulation_time = 5.0; // Total simulation time (seconds)
    int steps = static_cast<int>(simulation_time / dt);
    
    // Setpoints (RPM) - create a speed profile
    vector<double> setpoints(steps);
    for (int i = 0; i < steps; i++) {
        double t = i * dt;
        if (t < 1.0) {
            setpoints[i] = 0;           // Start at 0 RPM
        } else if (t < 2.0) {
            setpoints[i] = 500;         // Step to 500 RPM
        } else if (t < 3.0) {
            setpoints[i] = 1000;        // Step to 1000 RPM
        } else if (t < 4.0) {
            setpoints[i] = 700;         // Step down to 700 RPM
        } else {
            setpoints[i] = 0;           // Back to 0 RPM
        }
    }
    
    // Data logging
    vector<double> time(steps);
    vector<double> actual_speed(steps);
    vector<double> target_speed(steps);
    vector<double> control_signal(steps);
    vector<double> error_signal(steps);
    
    // Reset motor and controller
    motor.speed = 0.0;
    controller.reset();
    
    // Run simulation
    for (int i = 0; i < steps; i++) {
        time[i] = i * dt;
        target_speed[i] = setpoints[i];
        
        // Convert setpoint from RPM to rad/s for control
        double setpoint_rad_s = rpmToRadPerSec(setpoints[i]);
        
        // Calculate error
        double error = setpoint_rad_s - motor.speed;
        error_signal[i] = radPerSecToRpm(error);
        
        // Calculate control signal (voltage)
        double voltage = controller.calculate(error, dt);
        control_signal[i] = voltage;
        
        // Apply to motor
        simulateMotor(motor, voltage, dt);
        actual_speed[i] = radPerSecToRpm(motor.speed);
        
        // Display progress every 100 steps
        if (i % 100 == 0) {
            cout << "Time: " << time[i] << "s | ";
            cout << "Target: " << target_speed[i] << " RPM | ";
            cout << "Actual: " << actual_speed[i] << " RPM | ";
            cout << "Error: " << error_signal[i] << " RPM | ";
            cout << "Voltage: " << voltage << "V\n";
        }
    }
    
    // Save data to CSV file for plotting
    ofstream data_file("motor_simulation.csv");
    data_file << "Time,Target_Speed,Actual_Speed,Control_Signal,Error\n";
    for (int i = 0; i < steps; i++) {
        data_file << time[i] << ","
                  << target_speed[i] << ","
                  << actual_speed[i] << ","
                  << control_signal[i] << ","
                  << error_signal[i] << "\n";
    }
    data_file.close();
    
    cout << "\nSimulation complete! Data saved to 'motor_simulation.csv'\n";
    cout << "\nTry plotting the data:\n";
    cout << "1. Open the CSV in Excel or Google Sheets\n";
    cout << "2. Create line charts for Target_Speed vs Actual_Speed\n";
    cout << "3. Adjust controller gains (Kp, Ki) and rerun\n";
    
    return 0;
}