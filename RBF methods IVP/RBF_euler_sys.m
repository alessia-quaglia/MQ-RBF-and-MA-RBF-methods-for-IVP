function [err1, err2] = RBF_euler_sys(f, u_es1, u_es2, a, b, u0, N, Meth)
%
% Usage:     [err1, err2] = RBF_euler_sys(f, u_es1, u_es2, a, b, u0, N, Meth)
% Purpose:   it solves a system of 2 first order differential equations 
%            u'(t) = f(u(t),t) with an initial condition u(a) = u0 
%            and computes the global error for both components
%            using a RBF Euler's method of the second order
% Input:     f = given function f(t,u) of the problem
%        u_es1 = exact solution of the first component
%        u_es2 = exact solution of the second component
%            a = initial point
%            b = end point 
%           u0 = row vector containing the 2 initial values
%            N = controls the step size h
%         Meth = RBF Euler's method: '1' MQ-RBF, '2' MA-RBF
% Output: err1 = global error of the first component at the final point t = b
%         err2 = global error of the second component at the final point t = b
%
h = (b-a) / N;
t = linspace(a, b, N+1)';
u = u0' * ones(1, N+1); 
u(1, 2) = u_es1(t(2));
u(2, 2) = u_es2(t(2));
for n = 2:N
    fn = f(t(n), u(:,n)); 
    fn_1 = f(t(n-1), u(:,n-1));
    if Meth == 1
        eps2 = (fn - fn_1) ./ (h * u(:,n));
        u(:,n+1) = (1 + (eps2.*h^2)/2) .* (u(:,n) + h*fn);
    elseif Meth == 2
        eps2 = -(3*fn - 3*fn_1) ./ (h * u(:,n));
        u(:,n+1) = (1 - (eps2.*h^2)/6) .* u(:,n) + (h + (eps2.*h^3)/6) .* fn;
    else
        warning('Method not available')
    end
end
u = u';
err1 = abs(u(end,1) - u_es1(b));
err2 = abs(u(end,2) - u_es2(b));