% Exact solution u(t) = 1/(t^4+1) of the IVP 
% u' = -4t^3u^2 for -10 <= t <= 10

t = linspace(-10,10,1000);
u_exact = 1./(t.^4+1);
plot(t, u_exact, 'r', 'LineWidth', 1.2)
xlabel('t'); ylabel('u(t)');
axis([-10 10 0 1.2])
val_axis = [-10 -8 -6 -4 -2 0 2 4 6 8 10];
set(gca, 'XTick', val_axis, 'XTickLabel', string(val_axis));
title('Exact solution $u(t)$ of $u''= -4t^3u^2$', 'Interpreter', 'latex');
axis square; grid on; box on;
