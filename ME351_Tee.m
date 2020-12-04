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
Hfp = zeros(1, 800); % hf for pipe
Hft = zeros(1, 800); % hf for tee

% For iterating through solution arrays
i = 1;

%% Iterative solution
while z >= 0
    % Initial guess
    Vx0 = 0.5;
    Vx1 = 1;
    
    % Iterating to find Vx
    while abs(Vx0 - Vx1) > 0.001
        % Current Vx to previous Vx
        Vx0 = Vx1;
        Vy = 0.25526276*Vx0;
        
        % Friction coefficients
        Re_tube = rho*Vx0*d_tube/u;
        Re_tee = rho*Vx0*d_tee/u;
        
        % Depending on the type of flow for tube
        if Re_tube >= 4000
            eqn_tube = 1/sqrt(f) == -2*log(e/(d_tube*3.7)+2.51/(Re_tube*sqrt(f)));
        elseif Re_tube < 2300
            eqn_tube = f == 64/Re;
        else
            eqn_tube = f == 0.045;
        end
        
        % Depending on the type of flow for tee
        if Re_tee >= 4000
            eqn_tee = 1/sqrt(f) == -2*log(e/(d_tee*3.7)+2.51/(Re_tee*sqrt(f)));
        elseif Re_tee < 2300
            eqn_tee = f == 64/Re;
        else
            eqn_tee = f == 0.045;
        end
        
        % Solving for the friction coefficients
        f_tube = double(solve(eqn_tube, f));
        f_tee = double(solve(eqn_tee, f));
        
        % Friction equations
        hf_tube = (L_tube*f_tube*Vx0^2)/(d_tube*2*g);
        hf_tee = (L_tee*f_tee*Vy^2)/(d_tee*g);
        hm_entrance = K_tube*Vx^2/(2*g);
        hm_tee = K_tee*Vy^2/g;
        
        % Defining implicit equation for Vout = Uavg
        eqn = Vx == sqrt(301.110876*(z+L_tube/150+0.02-hf_tube-hf_tee-hm_entrance-hm_tee));
        Vx1 = double(solve(eqn, Vx));
        
    end
    
    % Position array
    Pos(i) = z;
    z = z - Vx1*tinc*A2/A1
    
    % Time array
    Tim(i) = t;
    t = t + tinc;
    
    % Velocity array
    Vel(i) = Vx1;
    
    % Reynold's array
    Reyp(i) = Re_tube;
    Reyt(i) = Re_tee;
    
    % Hf array
    Hfp(i) = hf_tube;
    Hft(i) = hf_tee;
    
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
subplot(3, 1, 1)
plot(Tim(1:idx-1), Vel(1:idx-1), '-r')
ylabel('Velocity, [m/s]');
% Pos vs Time
subplot(3, 1, 2)
plot(Tim(1:idx-1), Pos(1:idx-1), '-b')
ylabel('Position, [m]');
hold off

figure(2)
% Rey vs Time
subplot(4, 1, 1)
plot(Tim(1:idx-1), Reyp(1:idx-1), '-r')
ylabel('Re Tube');
% Rey vs Time
subplot(4, 1, 2)
plot(Tim(1:idx-1), Reyt(1:idx-1), '-b')
ylabel('Re Tee');
% Hf vs Time
subplot(4, 1, 3)
plot(Tim(1:idx-1), Hfp(1:idx-1), '-c')
ylabel('hf Pipe');
xlabel('Time, [s]');
% Hf vs Time
subplot(4, 1, 4)
plot(Tim(1:idx-1), Hft(1:idx-1), '-m')
ylabel('hf Tee');
xlabel('Time, [s]');