% The following script computes and shows the global errors 
% versus various N and local convergence orders for Euler's method,
% MQ-RBF and MA-RBF Euler's method applied to the IVP: 
% u'(t) = -u^2, u(0) = 1

addpath('Classic methods IVP')
addpath('RBF methods IVP')

f = @(t,u) -u.^2;
u_esatta = @(t) 1./(1 + t); 
a = 0;
b = 1;
u0 = 1;
N = [10, 20, 40, 80, 160, 320];

err_euler = zeros(size(N));
err_MQeuler = zeros(size(N));
err_MAeuler = zeros(size(N));

for i = 1:length(N)
    n = N(i); 
    err_euler(i) = euler(f, u_esatta, a, b, u0, n);
    err_MQeuler(i) = RBF_euler(f, u_esatta, a, b, u0, n, 1);
    err_MAeuler(i) = RBF_euler(f, u_esatta, a, b, u0, n, 2);
end

ord_euler = [NaN, log2(err_euler(1:end-1)./err_euler(2:end))];
ord_MQeuler = [NaN, log2(err_MQeuler(1:end-1)./err_MQeuler(2:end))];
ord_MAeuler = [NaN, log2(err_MAeuler(1:end-1)./err_MAeuler(2:end))];

%% Table generation
met_name = {'Euler', 'MQ-RBF Euler 2', 'MA-RBF Euler 2'};
met_err = {err_euler, err_MQeuler, err_MAeuler};
met_ord  = {ord_euler, ord_MQeuler, ord_MAeuler};
table_latex(met_name, N, met_err, met_ord);

%% Figure generation
figure
hold on; 
slope_N = [10, 320];
slope1 = 1e-2 * (slope_N/slope_N(1)).^-1; 
slope2 = 1e-3 * (slope_N/slope_N(1)).^-2; 
g_s1 = loglog(slope_N, slope1, 'b--', 'LineWidth', 1.2);
g_s2 = loglog(slope_N, slope2, 'k--', 'LineWidth', 1.2);
g_euler   = loglog(N, err_euler, 'b-o', 'LineWidth', 1, 'MarkerSize', 7);
g_MQeuler  = loglog(N, err_MQeuler,   'k-s', 'LineWidth', 1, 'MarkerSize', 7);
g_MAeuler = loglog(N, err_MAeuler,   'r-*', 'LineWidth', 1, 'MarkerSize', 7);
legend([g_euler, g_MQeuler, g_MAeuler, g_s1, g_s2], ...
       {'Euler', 'MQ-RBF Euler 2', 'MA-RBF Euler 2', 'Slope -1', 'Slope -2'}, ...
       'Location', 'southwest', 'FontSize', 9);
xlabel('N'); ylabel('Global error');
set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'XTick', N, 'XTickLabel', string(N));
axis([10 400 10^(-7) 10^(-1)])
title('Global errors for $u'' = -u^2$', 'Interpreter', 'latex');
axis square; grid on; box on;
hold off;


