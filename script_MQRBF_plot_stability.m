% The following script generates the absolute stability regions 
% for the analyzed MQ-RBF methods.

addpath('Stability region')

% Figure 2.1(a)
x = linspace(-4, 2, 400);
y = linspace(-3, 3, 400);
[X,Y] = meshgrid(x,y);
Z = X + 1i*Y;
R = 1 + Z;
figure
plot_stability(x, y ,R, 'k');
hold on 
R1 = 1 + Z + Z.^2/2 + Z.^3/2;
plot_stability(x, y, R1, 'r');
R2 = 1 + Z + (Z.^2 .* (3 + Z)) ./ (6 .* (1 + Z)) + ...
(Z.^3 .* (3 + Z)) ./(6 .* (1 + Z));
plot_stability(x, y, R2, 'b');
legend('Euler', 'MQ-RBF Euler 2', 'MQ-RBF Euler 3', 'Location', 'best');
title('Stability regions');
grid on
axis([-4 3 -3 3])
hold off;

% Figure 2.1(b)
x = linspace(-4, 2, 400);
y = linspace(-3, 3, 400);
[X,Y] = meshgrid(x,y);
Z = X + 1i*Y;
R1 = 1 + Z + Z.^2/2 + Z.^3/2;
figure
plot_stability(x, y, R1, 'r');
hold on 
AB2_classic
title('Stability regions')
legend('MQ-RBF Euler 2', 'AB2','Location', 'best')
axis([-3 3 -3 3])
hold off;

% Figure 2.2
x = linspace(-0.5, 0, 400);  
y = linspace(-3, 3, 1000);   
coeffs = @(z)[ 1, -(2*z+z^3/3), -1];
figure
RBF_multistep(x, y, coeffs, 'b')
title('Stability region');
legend('MQ-RBF Midpoint', 'Location', 'best')
axis ([-3 3 -3 3]);

% Figure 2.3(a)
figure
AB2_classic
hold on
x = linspace(-3, 3, 400);  
y = linspace(-3, 3, 400);   
coeffs = @(z)[ 1, -(1+3/2*z-7/24*z^3), 1/2*z-17/24*z^3];
RBF_multistep(x, y, coeffs, 'r')
title('Stability regions')
legend('AB2', 'MQ-RBF AB2', 'Location', 'best')
axis([-2 2 -2 2])
hold off; 

% Figure 2.3(b)
figure
RBF_AM1
title('Stability region')

legend('MQ-RBF AM1')

