% Exact solution u(t) = exp(t) - 2 of the IVP 
% u' = u + 2 for 0 <= t <= 1

t = linspace(0,1,1000);
u_exact = exp(t) - 2;
plot(t, u_exact, 'r', 'LineWidth', 1.2)
hold on;
yline(0, 'k', 'LineWidth', 1.2); 
t_int = log(2); 
line([t_int t_int], [-1 0], 'Color', 'k', 'LineStyle', '--');
xlabel('t'); ylabel('u(t)');
axis([0 1 -1 1])
val_axis = 0:0.1:1;
set(gca, 'XTick', val_axis, 'XTickLabel', string(val_axis));
title('Exact solution $u(t)$ of $u''= u + 2$', 'Interpreter', 'latex');
axis square; grid on; box on;
hold off;