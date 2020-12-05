syms Vx hf hm f

%% Variables - All units in  standard SI
% Starting level (maximum)
z = 0.08;
t = 0;
tinc = 1;

% Length of the tube
L_tube = 20/100;
L_tee = 4/100;

% Known Values
d_tube = 7.94/1000;
d_tee = 11.1125/1000;
e = 0.0025/1000;
g = 9.81;
rho = 998.19;
u = 0.001002;
l = 32/100;
w = 26/100;
K_tube = 0.8;
K_tee = 0.962772;
A1 = l*w;
A2 = pi*(d_tube/2)^2;

% Declaring the arrays for each plotted parameter
Vel = zeros(1, 800); % Output Velocity (Vy)
Tim = zeros(1, 800); % Time
Pos = zeros(1, 800); % Position
Reyp = zeros(1, 800); % Reynold's number pipe
Reyt = zeros(1, 800); % Reynold's number tee
fp = zeros(1, 800); % Friction coefficient for pipe
ft = zeros(1, 800); % Friction coefficient for tee

% For iterating through solution arrays
i = 1;

%% Iterating through z from 0.08 -> 0
while z >= 0
    % Initial guess for converging friction factors
    f0_tube = 0.03;
    f1_tube = 0.05;
    f0_tee = 0.03;
    f1_tee = 0.05;
    
    % Iterating until either f converges based on recalculated Vin = Vx
    while abs(f0_tube - f1_tube) > 0.0001 || abs(f0_tee - f1_tee) > 0.0001
        % Iterating our next friction coefficient
        f0_tube = f1_tube;
        f0_tee = f1_tee;
        
        % Defining direct equation for Vx
        eqn = Vx == sqrt((z+L_tube/150+0.02)/(0.1422377+L_tube*f0_tube/(d_tube*2*g)+f0_tee*L_tee*0.25526276^2/(d_tee*g)));
        Vx0 = double(solve(eqn, Vx));
        
        % Computing Vy based on derived ratio
        Vy = 0.25526276*Vx0;
        
        % Calculating the Reynold's number for each tube
        Re_tube = rho*Vx0*d_tube/u;
        Re_tee = rho*Vy*d_tee/u;
        
        % Piecewise assumption based on the calculated Reynold's number
        if Re_tube >= 4000
            eqn_tube = 1/sqrt(f) == -2*log(e/(d_tube*3.7)+2.51/(Re_tube*sqrt(f)));
        elseif Re_tube < 2300
            eqn_tube = f == 64/Re_tube;
        else
            eqn_tube = f == 0.045;
        end
        
        % Piecewise assumption based on the calculated Reynold's number
        if Re_tee >= 4000
            eqn_tee = 1/sqrt(f) == -2*log(e/(d_tee*3.7)+2.51/(Re_tee*sqrt(f)));
        elseif Re_tee < 2300
            eqn_tee = f == 64/Re_tee;
        else
            eqn_tee = f == 0.045;
        end
        
        % Calculating the new friction factor to compare against the previous one
        f1_tube = double(solve(eqn_tube, f));
        f1_tee = double(solve(eqn_tee, f));
    end
    
    % Position array
    Pos(i) = z;
    z = z - Vx0*tinc*A2/A1
    
    % Time array
    Tim(i) = t;
    t = t + tinc;
    
    % Velocity array
    Vel(i) = Vy;
    
    % Reynold's array
    Reyp(i) = Re_tube;
    Reyt(i) = Re_tee;
    
    % Friction coefficient array
    fp(i) = f1_tube;
    ft(i) = f1_tee;
    
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
subplot(2, 1, 1);
plot(Tim(1:idx-1), Reyp(1:idx-1), '-g'); % plots time on x vs reynold's number on y
ylabel('Re');
xlabel('Time, [s]');
title("Tube Reynold's Number vs. Time");
subplot(2, 1, 2);
plot(Tim(1:idx-1), Reyt(1:idx-1), '-g'); % plots time on x vs reynold's number on y
ylabel('Re');
xlabel('Time, [s]');
title("Tee Reynold's Number vs. Time");

figure(3);
plot(Pos(1:idx-1), Vel(1:idx-1), '-r'); % plots output velocity on y vs position on x
ylabel('Velocity, [m/s]');
xlabel('Position, [m]');
title('Output Velocity vs. Position');

figure(4);
subplot(2, 1, 1);
plot(Tim(1:idx-1), fp(1:idx-1), '-m'); % plots time on x vs friction factor on y
ylabel('Friction Factor');
xlabel('Time, [s]');
title('Tube Friction Factor vs. Time');
subplot(2, 1, 2)
plot(Tim(1:idx-1), ft(1:idx-1), '-m'); % plots time on x vs friction factor on y
ylabel('Friction Factor');
xlabel('Time, [s]');
title("Tee Friction Factor vs. Time");
