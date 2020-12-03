syms Vout Re hf f

%% Variables - All units in  standard SI
% Starting level (maximum)
z = 0.1;
t = 0;
tinc = 5;

% Pipe length
L = 0.1;

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

% for iterating through solution arrays
i = 1;

%% Euler method solution
% While the water level is still above 2cm (end height), starting at 0.1m
while z >= 0.02
    % Initial guess
    f0 = 1;
    f1 = 0.03;
    
    % Iterating to find Vout
    while abs(f0 - f1) > 0.001
        f0 = f1;

        % Solving for hf and hm (no summation needed)
        hf = (f0*L*Vout^2)/(2*d*g);
        hm = (K*Vout^2)/(2*g);

        % Defining implicit equation for Vout = Uavg
        eqn = Vout == sqrt(19.62*z+19.62*L/150-19.62*hf-19.62*hm);
        Uavg = double(solve(eqn, Vout));
        Re = (998.19*Uavg*d)/u;
        
        hf_val = (f0*L*Uavg^2)/(2*d*g);
        
        % Assuming Turbulent flow
        %eqn = (1/sqrt(f)) == -2*log(((e/d)/3.7) + (2.51/(Re*sqrt(f))));
        % Assuming Laminar flow
        eqn = f == 64/Re;
        f1 = double(solve(eqn, f));
    end
    
    %% Filling arrays
    % Position array
    Pos(i) = z;
    z = z - Uavg*tinc*A2/A1
    
    % Time array
    Tim(i) = t;
    t = t + tinc;
    
    % Velocity array
    Vel(i) = Uavg;
    
    % Reynold's number array
    Rey(i) = Re;
    
    % Hf array
    Hf(i) = hf_val;
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
figure(6); % opens a figure window
% Vel vs Time
plot(Tim(1:idx-1), Vel(1:idx-1), '-r')
ylabel('Velocity, [m/s]');
xlabel('Time, [s]');
title('Velocity vs time (Laminar)'); % creates a title for the plot

figure(7)
% Pos vs Time
plot(Tim(1:idx-1), Pos(1:idx-1), '-b')
ylabel('Position, [m]');
xlabel('Time, [s]');
title('Position vs time (Laminar)'); % creates a title for the plot

figure(8)
% Rey vs Time
plot(Tim(1:idx-1), Rey(1:idx-1), '-g')
ylabel('Re');
xlabel('Time, [s]');
title('Re vs time (Laminar)'); % creates a title for the plot

figure(9)
% Hf vs Time
plot(Tim(1:idx-1), Hf(1:idx-1), '-c')
ylabel('hf');
xlabel('Time, [s]');
title('hf vs time (Laminar)'); % creates a title for the plot

figure(10)
% Flip the arrays
Pos_snip = Pos(1:idx-1);
Pos_flip = flip(Pos_snip);
Vel_snip = Vel(1:idx-1);
Vel_flip = flip(Vel_snip);

plot(Pos_flip, Vel_flip, '-r'); % plots velocity on y vs position on x
set ( gca, 'xdir', 'reverse' )
ylabel('Velocity, [m/s]');
xlabel('Position, [m]');
title('Velocity vs Position (Laminar)'); % creates a title for the plot