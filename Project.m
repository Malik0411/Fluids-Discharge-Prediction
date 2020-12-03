syms Vout hf hm f

%% Variables - All units in  standard SI
% Starting level (maximum)
z = 0.1;
t = 0;
tinc = 1;

% First pipe
L = 0.6;

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

% Arrays of solutions
Vel = zeros(1, 500); % Velocity
Tim = zeros(1, 500); % Time
Pos = zeros(1, 500); % Position
Rey = zeros(1, 500); % Reynold's number
Hf = zeros(1, 500); % hf

% For iterating through solution arrays
i = 1;

%% Iterative solution
while z >= 0.02
    % Initial guess
    f0 = 0.1;
    f1 = 0.03;
    
    % Iterating to find Vout
    while abs(f0 - f1) > 0.001
        % Current f to previous f
        f0 = f1;

        % Solving for hf (no summation needed)
        hf = (f0*L*Vout^2)/(2*d*g);
        hm = (K*Vout^2)/(2*g);

        % Defining implicit equation for Vout = Uavg
        eqn = Vout == sqrt(19.62*z+19.62*L/150-19.62*hf-9.81*hm);
        Uavg = double(solve(eqn, Vout));
        Re = 998.19*Uavg*d/u;

        % Assuming Turbulent flow
        eqn = 1/sqrt(f) == -2*log(e/(d*3.7)+2.51/(Re*sqrt(f)));
        
        % % Assuming Laminar flow
        % eqn = f == 64/Re;
        
        f1 = double(solve(eqn, f));
    end
    
    % Position array
    Pos(i) = z;
    z = z - Uavg*tinc*A2/A1
    
    % Time array
    Tim(i) = t;
    t = t + tinc;
    
    % Velocity array
    Vel(i) = Uavg;
    
    % Reynold's Number array
    Rey(i) = Re;
    
    % Hf array
    hf_val = (f1*L*Uavg^2)/(2*d*g);
    Hf(i) = hf_val
    
    % Iterating the index
    i = i + 1;
    
end

%% Strip arrays
% Finding where the zeroes in B start
idx = 1;
val = -1;
while val == -1
    % Start from index 2 since 1 is zero
    idx = idx + 1;
    if (Tim(idx) == 0)
        val = 0;
    end
end

%% Plots of Velocity, Position (z), Reynold's number, hf with time, t
figure(1); % opens a figure window
% Vel vs Time
subplot(4, 1, 1)
plot(Tim(1:idx-1), Vel(1:idx-1), '-r')
ylabel('Velocity, [m/s]');
% Pos vs Time
subplot(4, 1, 2)
plot(Tim(1:idx-1), Pos(1:idx-1), '-b')
ylabel('Position, [m]');
% Rey vs Time
subplot(4, 1, 3)
plot(Tim(1:idx-1), Rey(1:idx-1), '-g')
ylabel('Re');
% Hf vs Time
subplot(4, 1, 4)
plot(Tim(1:idx-1), Hf(1:idx-1), '-c')
ylabel('hf');
xlabel('Time, [s]');

title('Parameters vs time'); % creates a title for the plot
hold off

figure(2)
plot(Pos(1:idx), Vel(1:idx), '-r'); % plots velocity on y vs position on x
ylabel('Velocity, [m/s]');
xlabel('Position, [m]');
title('Velocity vs Position'); % creates a title for the plot
