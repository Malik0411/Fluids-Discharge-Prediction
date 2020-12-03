syms Vout hf f

% Starting level (maximum)
z = 0.1;
t = 0;
tinc = 1;

% First pipe
L = 0.4;

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
A = zeros(1, 500);
B = zeros(1, 500);
C = zeros(1, 500);
i = 1;

while z >= 0.02
    % Initial guess
    f0 = 1;
    f1 = 0.03;
    
    % Iterating to find Vout
    while abs(f0 - f1) > 0.001
        f0 = f1;

        % Solving for hf (no summation needed)
        hf = (f0*L*Vout^2)/(2*d*g);
        hm = (K*Vout^2)/(2*g);

        eqn = Vout == sqrt(19.62*z+19.62*L/150-19.62*hf-19.62*hm);
        Uavg = double(solve(eqn, Vout));
        Re = (998.19*Uavg*d)/u;

        % Assuming Turbulent flow
        eqn = f == 64/Re;
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
title('Output Velocity vs. Water Level'); % creates a title for the plot
ylabel('Output Velocity, Vout [m/s]');
xlabel('Water Level, z [m]') % labels the x-axis
