function RBF_AM1()
%
% Usage:   RBF_AM1
% Purpose: the function plots the stability region of the 
%          MQ-RBF (or equivalently MA-RBF) Adams-Moulton 
%          one-step method in the complex plane, based on  
%          its characteristic polynomial
% Input:   none
% Output:  this function does not return outputs, it 
%          produces a filled contour plot showing the 
%          absolute stability region where all roots of
%          the characteristic polynomial satisfy |xi| <= 1
%
x = linspace(-3, 3, 400);   
y = linspace(-3, 3, 400);   
[X,Y] = meshgrid(x,y);
Z = X + 1i*Y;
Stable = zeros(size(Z));
for i = 1:numel(Z)
    z = Z(i);
    coeffs = [ 0, (1-z/2+z^3/24), -(1+z/2-z^3/24)]; 
    xi = roots(coeffs);
    if all(abs(xi) <= 1)
        Stable(i) = 1;
    end
end
contourf(X, Y, Stable, [0.5 1.5], 'LineColor','none');
colormap([1 1 1; 0.6 0.8 1]);
hold on
contour(X, Y, Stable, [0.5 0.5], 'LineWidth', 2, 'Color', 'b'); 
xlabel('\Re(h\lambda)');
ylabel('\Im(h\lambda)');
axis equal;
grid on;

