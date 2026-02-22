% The following script computes and shows the global errors 
% versus various N and local convergence orders for the midpoint method,
% the MQ/MA-RBF and the GA-RBF midpoint method applied to the IVP: 
% u'= (2t^2-u)/(t^2u-t), u(1) = 2

addpath('Classic methods IVP')
addpath('RBF methods IVP')

f = @(t,u) (2.*t.^2-u)/(t.^2.*u-t);
u_esatta = @(t) 1./t + sqrt(1./(t.^2)+4.*t-4);
a = 1;
b = 3;
u0 = 2;
N = [10, 20, 40, 80, 160, 320];

err_midpoint = zeros(size(N));
err_MQmidpoint = zeros(size(N));  
err_GAmidpoint = zeros(size(N));

for i = 1:length(N)
    n = N(i);
    err_midpoint(i) = midpoint(f, u_esatta, a, b, u0, n);
    err_MQmidpoint(i) = RBF_midpoint(f, u_esatta, a, b, u0, n, 1);
    err_GAmidpoint(i) = RBF_midpoint(f, u_esatta, a, b, u0, n, 2);
end

ord_midpoint = [NaN, log2(err_midpoint(1:end-1)./err_midpoint(2:end))];
ord_MQmidpoint = [NaN, log2(err_MQmidpoint(1:end-1)./err_MQmidpoint(2:end))];
ord_GAmidpoint = [NaN, log2(err_GAmidpoint(1:end-1)./err_GAmidpoint(2:end))];

%% Table generation
met_name = {'Midpoint', 'MA-RBF Midpoint', 'GA-RBF Midpoint'};
met_err = {err_midpoint, err_MQmidpoint, err_GAmidpoint};
met_ord  = {ord_midpoint, ord_MQmidpoint, ord_GAmidpoint};
table_latex(met_name, N, met_err, met_ord);

%% Figure generation
figure
hold on;  
slope_N = [10, 320];
slope = 1e-2 * (slope_N/slope_N(1)).^-3; 
g_s = loglog(slope_N, slope, 'k--', 'LineWidth', 1.2);
g_midpoint   = loglog(N, err_midpoint, 'b-o', 'LineWidth', 1, 'MarkerSize', 7);
g_MAmidpoint  = loglog(N, err_MQmidpoint,   'k-s', 'LineWidth', 1, 'MarkerSize', 7);
g_GAmidpoint = loglog(N, err_GAmidpoint,   'r-*', 'LineWidth', 1, 'MarkerSize', 7);
legend([g_midpoint, g_MAmidpoint, g_GAmidpoint, g_s], ...
       {'Midpoint', 'MA-RBF Midpoint', 'GA-RBF Midpoint', 'Slope -3'}, ...
       'Location', 'southwest', 'FontSize', 9);
xlabel('N'); ylabel('Global error');
set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'XTick', N, 'XTickLabel', string(N));
axis([10 400 10^(-7) 10^(1)])
title('Global errors for $u'' = \frac{2t^2-u}{t^2u-t}$', 'Interpreter', 'latex');
axis square; grid on; box on;
hold off;
