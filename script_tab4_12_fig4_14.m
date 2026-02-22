% The following script computes and shows the global errors 
% versus various N and local convergence orders for Euler's method, 
% MQ-RBF Euler's methods without the condition C2 and with C2 for 
% p = 1/2 and L=h^(-1/2), p=1 and L=0, p=1 and L=h^(-1) applied to the IVP:
% u' = u + 2, u(0) = -1

addpath('Classic methods IVP')
addpath('RBF methods IVP')

f = @(t,u) u + 2;
u_esatta = @(t) exp(t) - 2; 
a = 0;
b = 1;
u0 = -1;
N = [10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10000];

err_euler = zeros(size(N));
err_MQeuler = zeros(size(N));
err_MQeulerC2_1 = zeros(size(N));
err_MQeulerC2_2 = zeros(size(N));
err_MQeulerC2_3 = zeros(size(N));

for i = 1:length(N)
    n = N(i);
    h = (b-a)/n;
    err_euler(i) = euler(f, u_esatta, a, b, u0, n);
    err_MQeuler(i) = RBF_euler(f, u_esatta, a, b, u0, n, 1);
    err_MQeulerC2_1(i) = RBF_eulerC2(f, u_esatta, a, b, u0, n, 1, 0, 1);
    err_MQeulerC2_2(i) = RBF_eulerC2(f, u_esatta, a, b, u0, n, 1, h^(-1), 1);
    err_MQeulerC2_3(i) = RBF_eulerC2(f, u_esatta, a, b, u0, n, 0.5, h^(-0.5), 1);
end

ord_euler = [NaN, log2(err_euler(1:end-1)./err_euler(2:end))];
ord_MQeuler = [NaN, log2(err_MQeuler(1:end-1)./err_MQeuler(2:end))];
ord_MQeulerC2_1 = [NaN, log2(err_MQeulerC2_1(1:end-1)./err_MQeulerC2_1(2:end))];
ord_MQeulerC2_2 = [NaN, log2(err_MQeulerC2_2(1:end-1)./err_MQeulerC2_2(2:end))];
ord_MQeulerC2_3 = [NaN, log2(err_MQeulerC2_3(1:end-1)./err_MQeulerC2_3(2:end))];

%% Table generation
met_name = {'Euler', 'MQ-RBF Euler 2', 'MQ-RBF Euler 2 $p=1$, $L=0$', 'MQ-RBF Euler 2 $p=1$, $L=h^(-1)$',...
           'MQ-RBF Euler 2 $p=\frac{1}{2}$, $L=h^{\frac{1}{2}}$'};
met_err = {err_euler, err_MQeuler, err_MQeulerC2_1, err_MQeulerC2_2, err_MQeulerC2_3};
met_ord = {ord_euler, ord_MQeuler, ord_MQeulerC2_1, ord_MQeulerC2_2, ord_MQeulerC2_3};
table_latex(met_name, N, met_err, met_ord);

%% Figure generation
figure
hold on; 
slope_N = [10, 10000];
slope1 = 0.8e-1 * (slope_N/10).^-1; 
slope2 = 0.5e-2 * (slope_N/10).^-2; 
g_s1 = loglog(slope_N, slope1, 'b--', 'LineWidth', 1.2);
g_s2 = loglog(slope_N, slope2, 'k--', 'LineWidth', 1.2);
g_euler = loglog(N, err_euler, 'b-x', 'LineWidth', 1, 'MarkerSize', 7);
g_MQeuler = loglog(N, err_MQeuler, 'm-s', 'LineWidth', 1, 'MarkerSize', 7);
g_MQeulerC2_1 = loglog(N, err_MQeulerC2_1, 'k-*', 'LineWidth', 1, 'MarkerSize', 7);
g_MQeulerC2_2 = loglog(N, err_MQeulerC2_2, 'g-d', 'LineWidth', 1, 'MarkerSize', 7);
g_MQeulerC2_3 = loglog(N, err_MQeulerC2_3, 'r-o', 'LineWidth', 1, 'MarkerSize', 7);
legend([g_euler, g_MQeuler, g_MQeulerC2_1, g_MQeulerC2_2, g_MQeulerC2_3, g_s1, g_s2], ...
       {'Euler', 'MQ-RBF Euler 2', 'MQ-RBF Euler 2 $p=1$, $L=0$', 'MQ-RBF Euler 2 $p=1$, $L=h^{-1}$',...
       'MQ-RBF Euler 2 $p=\frac{1}{2}$, $L=h^{-\frac{1}{2}}$', 'Slope -1', 'Slope -2'}, ...
       'Interpreter','latex','Location', 'southwest', 'FontSize', 9);
xlabel('N'); ylabel('Global error');
set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'XTick', N, 'XTickLabel', string(N));
axis([10 10000 10^(-10) 10^(0)])
title('Global errors for $u'' = u+2$', 'Interpreter', 'latex');
axis square; grid on; box on;
hold off;

