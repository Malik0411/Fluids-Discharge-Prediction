syms Vout hf hm f

%% Variables - All units in  standard SI
% Starting level (maximum)
z = 0.08;
t = 0;
tinc = 1;

% Length of the tube
L = 50/100;

% Known Values
d = 7.94/1000;
e = 0.0025/1000;
g = 9.81;
rho = 998.19;
u = 0.001002;
l = 32/100;
w = 26/100;
K = 0.8;
A1 = l*w;
A2 = pi*(d/2)^2;

% Declaring the arrays for each plotted parameter
Vel = zeros(1, 300); % Velocity
Tim = zeros(1, 300); % Time
Pos = zeros(1, 300); % Position
Rey = zeros(1, 300); % Reynold's number
F = zeros(1, 300); % Friction Factor

% For iterating through the above arrays
i = 1;

%% Iterating through z from 0.08 -> 0
while z >= 0
    % Initial guess for converging friction factors
    f0 = 0.03;
    f1 = 0.05;
    
    % Iterating until f converges based on recalculated Vout = Uavg
    while abs(f0 - f1) > 0.0001
        % Setting the current f to the previous f
        f0 = f1;

        % Defining implicit equation for Vout = Uavg
        eqn = Vout == sqrt((2*g*z+(g*L)/75+0.04*g-(L*f0*Vout^2)/d-K*Vout^2)/(1-(A2/A1)^2));
        Uavg = double(solve(eqn, Vout));
        Re = rho*Uavg*d/u;

        % Piecewise assumption based on the calculated Reynold's number
        if Re >= 4000
            eqn = 1/sqrt(f) == -2*log(e/(d*3.7)+2.51/(Re*sqrt(f)));
        elseif Re < 2300
            eqn = f == 64/Re;
        else
            eqn = f == 0.045;
        end
            
        % Calculating the new friction factor to compare against the previous one
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
    F(i) = f1;
    
    % Iterating the index
    i = i + 1;
    
end

%% Strip arrays of zeros
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

%% Plots of Velocity, Position (z), Reynold's number, f with time, t
figure(1); % opens a figure window
% Vel vs Time
plot(Tim(1:idx-1), Vel(1:idx-1), '-b'); % plots time on x vs output velocity on y
ylabel('Velocity, [m/s]');
xlabel('Time, [s]');
title('Output Velocity vs. Time'); % creates a title for the plot

figure(2); % opens a figure window
% Rey vs Time
plot(Tim(1:idx-1), Rey(1:idx-1), '-g'); % plots time on x vs reynold's number on y
ylabel('Re');
xlabel('Time, [s]');
title("Reynold's Number vs. Time");

figure(3); % opens a figure window
plot(Pos(1:idx-1), Vel(1:idx-1), '-r'); % plots velocity on y vs position on x
ylabel('Velocity, [m/s]');
xlabel('Position, [m]');
title('Output Velocity vs. Position');

figure(4); % opens a figure window
plot(Tim(1:idx-1), F(1:idx-1), '-m'); % plots time on x vs friction factor on y
ylabel('Friction Factor');
xlabel('Time, [s]');
title('Friction Factor vs. Time');