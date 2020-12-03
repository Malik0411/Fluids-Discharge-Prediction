syms Vout hf f

% Starting level (maximum)
z = 0;
t = 0;
tinc = 1;

% First pipe
L = 0.2;

% Known Values
d = 7.24/1000;
e = 0.0025/1000;
g = 9.81;
u = 0.001002;
l = 32/100;
w = 26/100;
K = 0.5;
A1 = l*w;
A2 = pi*(d/2)^2;

% Array of Vout solutions
A = zeros(1, 2000);
B = zeros(1, 2000);
C = zeros(1, 2000);
i = 1;

while z >= -0.08
    % Initial guess
    f0 = 1;
    f1 = 0.03;
    
    % Iterating to find Vout
    while abs(f0 - f1) > 0.001
        f0 = f1;

        % Solving for hf (no summation needed)
        hf = (f0*L*Vout^2)/(2*d*g);
        hm = (K*Vout^2)/(2*g);

        eqn = Vout == sqrt(19.6*z+1.96+19.6*L/150-19.6*hf-19.6*hm);
        Uavg = double(solve(eqn, Vout));
        Re = (998.19*Uavg*d)/u;

        % Assuming Turbulent flow
        eqn = 1/sqrt(f) == -2*log(e/(d*3.7)+2.51/(Re*sqrt(f)));
        f1 = double(solve(eqn, f));
    end
    
    % Position array
    C(i) = z;
    z = z - Uavg*tinc*A2/A1
    
    % Time array
    B(i) = t;
    t = t + tinc;
    
    % Velocity array
    A(i) = Uavg;
    i = i + 1;
end

figure(1); % opens a figure window
plot(B, A, '-r'); % plots acceleration versus time
title('Output velocity vs. Water Level'); % creates a title for the plot
ylabel('Output Velocity, Vout [m/s]');
xlabel('Water Level, z [m]') % labels the x-axis

