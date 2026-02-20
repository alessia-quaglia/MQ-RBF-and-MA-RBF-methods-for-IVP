% The following script computes and shows the global errors 
% versus various N and local convergence orders for the midpoint method,
% the MQ/MA-RBF and the GA-RBF midpoint method applied to the IVP: 
% u'(t) = -u^2, u(0) = 1

addpath('Classic methods IVP')
addpath('RBF methods IVP')

f = @(t,u) -u.^2;
u_esatta = @(t) 1./(1 + t); 
a = 0;
b = 1;
u0 = 1;
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
slope1 = 1e-3 * (slope_N/slope_N(1)).^-2; 
slope2 = 0.3e-3 * (slope_N/slope_N(1)).^-3; 
g_s1 = loglog(slope_N, slope1, 'b--', 'LineWidth', 1.2);
g_s2 = loglog(slope_N, slope2, 'k--', 'LineWidth', 1.2);
g_midpoint   = loglog(N, err_midpoint, 'b-o', 'LineWidth', 1, 'MarkerSize', 7);
g_MAmidpoint  = loglog(N, err_MQmidpoint,   'k-s', 'LineWidth', 1, 'MarkerSize', 7);
g_GAmidpoint = loglog(N, err_GAmidpoint,   'r-*', 'LineWidth', 1, 'MarkerSize', 7);
legend([g_midpoint, g_MAmidpoint, g_GAmidpoint, g_s1, g_s2], ...
       {'Midpoint', 'MA-RBF Midpoint', 'GA-RBF Midpoint', 'Slope -2', 'Slope -3'}, ...
       'Location', 'northeast', 'FontSize', 9);
xlabel('N'); ylabel('Global error');
set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'XTick', N, 'XTickLabel', string(N));
axis([10 400 10^(-10) 10^(-2)])
title('Global errors for $u'' = -u^2$', 'Interpreter', 'latex');
axis square; grid on; box on;
hold off;


