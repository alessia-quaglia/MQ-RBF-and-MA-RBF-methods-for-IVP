% The following script generates the absolute stability regions 
% for the analyzed MA-RBF methods.

addpath('Stability region')

% Figure 3.2(a)
x = linspace(-4, 2, 400);
y = linspace(-3, 3, 400);
[X,Y] = meshgrid(x,y);
Z = X + 1i*Y;
R = 1 + Z;
figure
plot_stability(x, y, R, 'k');
hold on 
R3 = 1 + Z + Z.^2/2 - Z.^3/2;
plot_stability(x, y, R3, 'g');
R4 = 1 + Z + (Z.^2 .* (3 + Z)) ./ (6 .* (1 - Z)) - (Z.^3 .*...
(3 + Z)) ./(6 .* (1 - Z));
plot_stability(x, y, R4, 'r');
legend('Euler', 'MA-RBF Euler 2', 'MA-RBF Euler 3', 'Location', 'best');
title('Stability regions');
grid on
axis([-4 3 -3 3])
hold off;

% Figure 3.2(b)
x = linspace(-4, 2, 400);
y = linspace(-3, 3, 400);
[X,Y] = meshgrid(x,y);
Z = X + 1i*Y;
R = 1 + Z;
figure
plot_stability(x, y, R, 'k');
hold on 
R1 = 1 + Z + Z.^2/2 + Z.^3/2; % MQ-Euler 2
plot_stability(x, y, R1, 'b');
R2 = 1 + Z + (Z.^2 .* (3 + Z)) ./ (6 .* (1 + Z)) + (Z.^3 .*...
(3 + Z)) ./(6 .* (1 + Z)); % MQ-Euler 3
plot_stability(x, y, R2, 'y');
R3 = 1 + Z + Z.^2/2 - Z.^3/2;
plot_stability(x, y, R3, 'g');
R4 = 1 + Z + (Z.^2 .* (3 + Z)) ./ (6 .* (1 - Z)) - (Z.^3 .*...
(3 + Z)) ./ (6 .* (1 - Z));
plot_stability(x, y, R4, 'r');
legend('Euler', 'MQ-RBF Euler 2', 'MQ-RBF Euler 3', 'MA-RBF Euler 2', ...
'MA-RBF Euler 3', 'Location', 'best');
title('Stability regions');
grid on
axis([-4 3 -3 3])
hold off;

% Figure 3.3
x = linspace(-0.5, 0, 400);  
y = linspace(-3, 3, 1000);   
coeffs = @(z)[ 1, -(2*z+z^3/3), -1];
figure
RBF_multistep(x, y, coeffs, 'b')
title('Stability region');
legend('MA-RBF Midpoint', 'Location', 'best')
axis ([-3 3 -3 3]);

% Figure 3.4(a)
figure
AB2_classic
hold on
x = linspace(-3, 3, 800);  
y = linspace(-3, 3, 800);   
coeffsMA = @(z)[ 1, -(1+3/2*z+41/24*z^3), 1/2*z+31/24*z^3];
RBF_multistep(x, y, coeffsMA, 'g')
title('Stability regions')
legend('AB2', 'MA-RBF AB2', 'Location', 'best')
axis([-2 1 -1.5 1.5])
hold off;

% Figure 3.4(b)
figure
AB2_classic
hold on
x = linspace(-3, 3, 800);  
y = linspace(-3, 3, 800);   
coeffsMQ = @(z)[ 1, -(1+3/2*z-7/24*z^3), 1/2*z-17/24*z^3];
RBF_multistep(x,y,coeffsMQ,'r')
coeffsMA = @(z)[ 1, -(1+3/2*z+41/24*z^3), 1/2*z+31/24*z^3];
RBF_multistep(x,y,coeffsMA,'g')
title('Stability regions')
legend('AB2', 'MQ-RBF AB2', 'MA-RBF AB2', 'Location', 'best')
axis([-2 1 -1.5 1.5])
hold off;

% Figure 3.5
figure
RBF_AM1
title('Stability region')

legend('MA-RBF AM1')
