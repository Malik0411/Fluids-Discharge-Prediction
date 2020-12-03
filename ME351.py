import sympy as sym

Vout, hf, f = sym.symbols('Vout hf f')

# Starting level (maximum)
z = 0
t = 0
tinc = 0.1

# First pipe
L = 0.2

# Known Values
d = 7.24/1000
e = 0.0025/1000
g = 9.81
u = 0.001002
l = 32/100
w = 26/100

# Array of Vout solutions
A = []
B = []
C = []
i = 1

while z >= -0.08:
    # Initial guess
    f0 = 1
    f1 = 0.03
    
    # Iterating to find Vout
    while abs(f0 - f1) > 0.00001:
        f0 = f1

        # Solving for hf (no summation needed)
        hf = (f0*L*Vout^2)/(2*d*g)

        eqn = Vout == sym.sqrt(19.6*z+1.96+19.6*L/150-19.6*hf)
        Uavg = float(sym.solveSet(eqn, Vout))
        Re = (998.19*Uavg*d)/u

        # Assuming Turbulent flow
        eqn = 1/sym.sqrt(f) == -2*sym.log(e/(d*3.7)+2.51/(Re*sym.sqrt(f)))
        f1 = float(sym.solveSet(eqn, f))
    
    C.append(z)
    z = z - (Uavg*tinc*sym.pi*(d/2)^2)/(l*w)

    B.append(t)
    t = t + tinc

    A.append(Uavg)
    i = i + 1

print(len(B)*tinc)

# figure(1); % opens a figure window
# plot(B, A, '-r'); % plots acceleration versus time
# title('Output velocity vs. Water Level'); % creates a title for the plot
# ylabel('Output Velocity, Vout [m/s]');
# xlabel('Water Level, z [m]') % labels the x-axis

