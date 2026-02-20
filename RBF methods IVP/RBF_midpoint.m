function err = RBF_midpoint(f, u_es, a, b, u0, N, Meth)
%
% Usage:    err = RBF_midpoint(f, u_es, a, b, u0, N, Meth)
% Purpose:  it solves the differential equation u'(t) = f(u(t),t) with
%           an initial condition u(a) = u0 and computes the global error
%           using a RBF Midpoint method of the third order
% Input:    f = given function f(t,u) of the problem
%        u_es = exact solution
%           a = initial point
%           b = end point 
%          u0 = initial value
%           N = controls the step size h
%        Meth = RBF Midpoint method: '1' MQ-RBF (or equivalently
%               MA-RBF), '2' GA-RBF 
% Output: err = global error at the final point t = b
%
h = (b-a)/N;
t = linspace(a,b,N+1)';
u = zeros(N+1, 1); 
u(1) = u0;
u(2) = u_es(t(2));
u(3) = u_es(t(3));
% u(2) = rk4_step(f, u(1), h);
% u(3) = rk4_step(f, u(2), h);
for n = 2:N-1
    fn = f(t(n),u(n)); 
    fn_m1 = f(t(n-1),u(n-1));
    fn_p1 = f(t(n+1),u(n+1));
    if Meth == 1
        eps2 = -(fn_p1 - 2*fn + fn_m1) / (3 * h^2 * fn);
        u(n+2) = u(n) + (2*h - eps2*h^3) * fn_p1;
    elseif Meth == 2
        eps2 = -(fn_p1 - 2*fn + fn_m1) / (6 * h^2 * fn);
        u(n+2) = u(n) + 2*h*fn_p1*exp(-eps2*h^2);
    else 
        warning('Method not available')
    end
end
err = abs(u(end) - u_es(b));
