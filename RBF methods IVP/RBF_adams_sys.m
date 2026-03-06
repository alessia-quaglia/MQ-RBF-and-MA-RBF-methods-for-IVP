function [err1, err2] = RBF_adams_sys(f, u_es1, u_es2, a, b, u0, N, Meth)
%
% Usage:     [err1, err2] = RBF_adams_sys(f, u_es1, u_es2, a, b, u0, N, Meth)
% Purpose:   it solves a system of 2 first order differential equations
%            u'(t) = f(u(t),t) with an initial condition u(a) = u0 and 
%            computes the global error for both components using a 
%            RBF Adams method of the third order
% Input:     f = given function f(t,u) of the problem
%        u_es1 = exact solution of the first component
%        u_es2 = exact solution of the second component
%            a = initial point
%            b = end point 
%           u0 = row vector containing the 2 initial values
%            N = controls the step size h
%         Meth = RBF Adams method: '1' MQ-RBF two-step Adams-Bashforth, 
%                '2' MQ-RBF two step Adams-Bashforth one-step Adams-Moulton
%                predictor-corrector, '3' MA-RBF two-step Adams-Bashforth, 
%                '4' MA-RBF two step Adams-Bashforth one-step Adams-Moulton
%                predictor-corrector
% Output: err1 = global error of the first component at the final point t = b
%         err2 = global error of the second component at the final point t = b
%
h = (b-a) / N;
t = linspace(a, b, N+1)';
u = u0' * ones(1, N+1); 
u(1, 2) = u_es1(t(2));
u(2, 2) = u_es2(t(2));
u(1, 3) = u_es1(t(3));
u(2, 3) = u_es2(t(3));
t3 = f(t(1), u(:,1)); 
t2 = f(t(2), u(:,2)); 
for n = 3:N
    t1 = f(t(n), u(:,n));
    if Meth == 1 
        eps2 = (t1 - 2.*t2 + t3) ./ (h^2 .* t2);
        u(:,n+1) = u(:,n) + (3/2*h - 7/24.*eps2*h^3).*t1 + (-1/2*h + 17/24.*eps2*h^3).*t2;
    elseif Meth == 2 
        eps2 = (t1 - 2.*t2 + t3) ./ (h^2 .* t2);
        u_p = u(:,n) + (3/2*h - 7/24.*eps2*h^3).*t1 + (-1/2*h + 17/24.*eps2*h^3).*t2;
        t0 = f(t(n+1), u_p);
        u(:,n+1) = u(:,n) + (h/2 - eps2./24*h^3) .* (t0 + t1); 
    elseif Meth == 3 
        eps2 = -(3*t1 - 6.*t2 + 3*t3) ./ (h^2 .* t2);
        u(:,n+1) = u(:,n) + (3/2*h - 41/72.*eps2*h^3).*t1 + (-1/2*h + 31/72.*eps2*h^3).*t2;
    elseif Meth == 4 
        eps2 = -(3*t1 - 6.*t2 + 3*t3) ./ (h^2 .* t2);
        u_p = u(:,n) + (3/2*h - 41/72.*eps2*h^3).*t1 + (-1/2*h + 31/72.*eps2*h^3).*t2;
        t0 = f(t(n+1), u_p);
        u(:,n+1) = u(:,n) + (h/2 + eps2./72*h^3) .* (t0 + t1); 
    else 
        warning('Method not available')
    end
    if n < N
       t3 = t2; 
       t2 = t1; 
    else
       break
    end
end
u = u';
err1 = abs(u(end,1) - u_es1(b));
err2 = abs(u(end,2) - u_es2(b));


