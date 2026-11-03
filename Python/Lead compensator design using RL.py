# num = input("Enter numerator: ")
# den = input("Enter numerator: ")

# G = c.tf(num, den)

# mat_OS = input("enter desired phase margin")
# settling_time = input("enter desired steady state error for ramp input(Kv)")


# step 1 
# OS = etp((-zeta*pi)/sqrt(1 - zeta^2)) solve to find zeta = z

# step 2 
# Ts = 4/z*wn solve to find wn 

# step 3 
# find dominant pole sd = -z*wn + j wn*sqrt(1 - z^2)

# step 4 
# find angle deficiency (phi) = 180 - (phase of G(s) at s = sd)

# step 5 
# select the zero of compensator so that it cancels out pole of system which cause matimum lag 

# step 6 
# using angle criteria i.e. phase of Gc * G = 180 find location of pole of compensator 

# step 7 
# using magnitude criteria |Gc * G| = 1 find gain of compensator Kc 

# step 8 
# Check the response 



