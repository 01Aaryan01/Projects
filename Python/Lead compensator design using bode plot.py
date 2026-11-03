# num = input("Enter numerator: ")
# den = input("Enter numerator: ")

# G = c.tf(num, den)

# PM_d = input("enter desired phase margin")
# ess_d = input("enter desired steady state error for ramp input(Kv)")


# step 1 
# ess = 1/Kv Find Kv 

# step 2 
# Kv = lim s-> 0 s Kc G(s) => find Kc 

# step 3 
# G1 = Kc*G 

# step 4 
# find bode plot of G1 and find phase margin 

# step 5 
# find phi = PM_d - PM + epsilon 

# step 6 
# find alpha = (1 - sin(phi))/(1 + sin(phi))

# step 7 
# 20log10(G1(jw)) = -10log10(alpha) solve this equation to find new gcf 

# step 8 
# Find T = 1/(w_m*sqrt(alpha))

# step 9 
# find zero zc = 1/T and pole pc = 1/alpha*T 

# step 10 
# plot bode plot of Gc*G1 and check PM 








