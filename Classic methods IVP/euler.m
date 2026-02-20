function err = euler(f, u_es, a, b, u0, N)
%
% Usage:    err = euler(f, u_es, a, b, u0, N)
% Purpose:  it solves the differential equation u'(t) = f(u(t),t) with
%           an initial condition u(a) = u0 and computes the global error
%           using the Euler's method
% Input:    f = given function f(t,u) of the problem
%        u_es = exact solution
%           a = initial point
%           b = end point 
%          u0 = initial value
%           N = controls the step size h
% Output: err = global error at the final point t = b
%
h = (b-a) / N;
t = linspace(a, b, N+1)';
u = zeros(N+1, 1); 
u(1) = u0;
for n = 1:N
    u(n+1) = u(n) + h*f(t(n),u(n));
end
err = abs(u(end) - u_es(b));
