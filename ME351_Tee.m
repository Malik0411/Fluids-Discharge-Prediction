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

% Arrays of solutions
Vel = zeros(1, 800); % Velocity
Tim = zeros(1, 800); % Time
Pos = zeros(1, 800); % Position
Reyp = zeros(1, 800); % Reynold's number pipe
Reyt = zeros(1, 800); % Reynold's number tee
fp = zeros(1, 800); % Friction coefficient for pipe
ft = zeros(1, 800); % Friction coefficient for tee

% For iterating through solution arrays
i = 1;

% Iterative solution
while z >= 0
    
    % Initial guess
    f0_tube = 0.03;
    f1_tube = 0.05;
    f0_tee = 0.03;
    f1_tee = 0.05;
    
    % Iterating to find f_tube and f_tee
    while abs(f0_tube - f1_tube) > 0.001 || abs(f0_tee - f1_tee) > 0.001
        
        % Iterating our next friction coefficient
        f0_tube = f1_tube;
        f0_tee = f1_tee;
        
        % Defining implicit equation for Vx
        eqn = Vx == sqrt((z+L_tube/150+0.02)/(0.1422377+L_tube*f0_tube/(d_tube*2*g)+f0_tee*L_tee*0.25526276^2/(d_tee*g)));
        Vx0 = double(solve(eqn, Vx));
        
        % compute Vy
        Vy = 0.25526276*Vx0;
        
        % Friction coefficients
        Re_tube = rho*Vx0*d_tube/u;
        Re_tee = rho*Vy*d_tee/u;
        
        % Depending on the type of flow for tube
        if Re_tube >= 4000
            eqn_tube = 1/sqrt(f) == -2*log(e/(d_tube*3.7)+2.51/(Re_tube*sqrt(f)));
        elseif Re_tube < 2300
            eqn_tube = f == 64/Re_tube;
        else
            eqn_tube = f == 0.045;
        end
        
        % Depending on the type of flow for tee
        if Re_tee >= 4000
            eqn_tee = 1/sqrt(f) == -2*log(e/(d_tee*3.7)+2.51/(Re_tee*sqrt(f)));
        elseif Re_tee < 2300
            eqn_tee = f == 64/Re_tee;
        else
            eqn_tee = f == 0.045;
        end
        
        % Solving for the new friction coefficients
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
    Vel(i) = Vx0;
    
    % Reynold's array
    Reyp(i) = Re_tube;
    Reyt(i) = Re_tee;
    
    % Friction coefficient array
    fp(i) = f1_tee;
    ft(i) = f1_tube;
    
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
plot(Tim(1:idx-1), Reyp(1:idx-1), '-g')
ylabel('Pipe Reynolds');
% Rey vs Time
subplot(4, 1, 4)
plot(Tim(1:idx-1), Reyt(1:idx-1), '-m')
ylabel('Tee Reynolds');
hold off

figure(2)
subplot(3, 1, 1)
plot(Pos(1:idx), Vel(1:idx), '-r'); % plots velocity on y vs position on x
ylabel('Velocity, [m/s]');
xlabel('Position, [m]');
title('Velocity vs Position'); % creates a title for the plot
% Friction coefficient vs Time
subplot(3, 1, 2)
plot(Tim(1:idx-1), fp(1:idx-1), '-b')
ylabel('Pipe f');
xlabel('Time, [s]');
% Friction coefficient vs Time
subplot(3, 1, 3)
plot(Tim(1:idx-1), ft(1:idx-1), '-m')
ylabel('Tee f');
xlabel('Time, [s]');
hold off
