%% Declaring the Datasets for Length, Experimental and Calculated Drain Times
Len = [20, 30, 40, 60];
Experimental = [199, 214, 266, 288];
Calculated = [197, 210, 224, 256];

%% Creating the Experimental vs. Simulated Plot for Simple Model Comparison
figure(1);
plot(Len, Experimental, '-b', Len, Calculated, '-m');
title('Experimental & Simulated Run Times vs. Tube Length');
xlabel('Tube Length, [cm]');
ylabel('Run Time, [s]');
legend('Experimental', 'Simulated');
