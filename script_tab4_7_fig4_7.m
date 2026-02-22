% The following script computes and shows the global errors 
% versus various N and local convergence orders for the two-step 
% Adams-Bashforth method, the MQ-RBF and MA-RBF two-step Adams-Bashforth methods, 
% the Adams-Bashforth-Moulton predictor-corrector method and its 
% corresponding MQ-RBF and MA-RBF modifications applied to the IVP: 
% u'= (2t^2-u)/(t^2u-t), u(1) = 2

addpath('Classic methods IVP')
addpath('RBF methods IVP')

f = @(t,u) (2.*t.^2-u)/(t.^2.*u-t);
u_esatta = @(t) 1./t + sqrt(1./(t.^2)+4.*t-4);
a = 1;
b = 3;
u0 = 2;
N = [10, 20, 40, 80, 160, 320];

err_ab2 = zeros(size(N));
err_MQab2 = zeros(size(N));
err_MAab2 = zeros(size(N));
err_ab2am1 = zeros(size(N));
err_MQab2am1 = zeros(size(N));
err_MAab2am1 = zeros(size(N));

for i = 1:length(N)
    n = N(i); 
    err_ab2(i) = ab2(f, u_esatta, a, b, u0, n);
    err_MQab2(i) = RBF_adams(f, u_esatta, a, b, u0, n, 1);
    err_MAab2(i) = RBF_adams(f, u_esatta, a, b, u0, n, 3);
    err_ab2am1(i) = ab2am1(f, u_esatta, a, b, u0, n);
    err_MQab2am1(i) = RBF_adams(f, u_esatta, a, b, u0, n, 2);
    err_MAab2am1(i) = RBF_adams(f, u_esatta, a, b, u0, n, 4);
end

ord_ab2 = [NaN, log2(err_ab2(1:end-1)./err_ab2(2:end))];
ord_MQab2 = [NaN, log2(err_MQab2(1:end-1)./err_MQab2(2:end))];
ord_MAab2 = [NaN, log2(err_MAab2(1:end-1)./err_MAab2(2:end))];
ord_ab2am1 = [NaN, log2(err_ab2am1(1:end-1)./err_ab2am1(2:end))];
ord_MQab2am1 = [NaN, log2(err_MQab2am1(1:end-1)./err_MQab2am1(2:end))];
ord_MAab2am1 = [NaN, log2(err_MAab2am1(1:end-1)./err_MAab2am1(2:end))];

%% Table generation
met_name = {'AB2', 'MQ-RBF AB2', 'MA-RBF AB2', 'AB2-AM1', 'MQ-RBF AB2-AM1', 'MA-RBF AB2-AM1'};
met_err = {err_ab2, err_MQab2, err_MAab2, err_ab2am1, err_MQab2am1, err_MAab2am1};
met_ord  = {ord_ab2, ord_MQab2, ord_MAab2, ord_ab2am1, ord_MQab2am1, ord_MAab2am1};
table_latex(met_name, N, met_err, met_ord);

%% Figure generation
figure
hold on; 
slope_N = [10, 320];
slope1 = 0.15e-1 * (slope_N/slope_N(1)).^-2; 
slope2 = 0.3e-2 * (slope_N/slope_N(1)).^-3; 
g_s1 = loglog(slope_N, slope1, 'b--', 'LineWidth', 1.2);
g_s2 = loglog(slope_N, slope2, 'k--', 'LineWidth', 1.2);
g_ab2 = loglog(N, err_ab2, 'b-o', 'LineWidth', 1, 'MarkerSize', 7);
g_MQab2  = loglog(N, err_MQab2, 'k-s', 'LineWidth', 1, 'MarkerSize', 7);
g_MAab2 = loglog(N, err_MAab2, 'r-*', 'LineWidth', 1, 'MarkerSize', 7);
g_ab2am1 = loglog(N, err_ab2am1, 'g-o', 'LineWidth', 1, 'MarkerSize', 7);
g_MQab2am1  = loglog(N, err_MQab2am1, 'm-s', 'LineWidth', 1, 'MarkerSize', 7);
g_MAab2am1 = loglog(N, err_MAab2am1, 'c-*', 'LineWidth', 1, 'MarkerSize', 7);
legend([g_ab2, g_MQab2, g_MAab2, g_ab2am1, g_MQab2am1, g_MAab2am1, g_s1, g_s2], ...
       {'AB2', 'MQ-RBF AB2', 'MA-RBF AB2', 'AB2-AM1', 'MQ-RBF AB2-AM1', ...
       'MA-RBF AB2-AM1', 'Slope -2', 'Slope -3'}, 'Location', 'southwest', 'FontSize', 9);
xlabel('N'); ylabel('Global error');
set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'XTick', N, 'XTickLabel', string(N));
axis([10 400 10^(-8) 10^(-1)])
title('Global errors for $u'' = \frac{2t^2-u}{t^2u-t}$', 'Interpreter', 'latex');
axis square; grid on; box on;
hold off;
