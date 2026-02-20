function err = RBF_eulerC2(f, u_es, a, b, u0, N, p, L, Meth)
%
% Usage:    err = RBF_eulerC2(f, u_es, a, b, u0, N, p, L, Meth)
% Purpose:  it solves the differential equation u'(t) = f(u(t),t) with
%           an initial condition u(a) = u0 and computes the global error
%           using a RBF Euler's method with the condition C2
% Input:    f = given function f(t,u) of the problem
%        u_es = exact solution
%           a = initial point
%           b = end point 
%          u0 = initial value
%           N = controls the step size h
%        p, L = nonnegative constants in the condition C2
%        Meth = RBF Euler's method with C2: '1' MQ-RBF, '2' MA-RBF
% Output: err = global error at the final point t = b
%
h = (b-a)/N;
t = linspace(a,b,N+1)';
u = zeros(N+1, 1); 
u(1) = u0;
u(2) = u_es(t(2));
for n = 2:N
    fn = f(t(n),u(n)); 
    fn_1 = f(t(n-1),u(n-1));
    if Meth == 1
        if abs(u(n)) <= h^p
            eps2 = sign((fn - fn_1) / (h * u(n))) * L;
        else
            eps2 = (fn - fn_1) / (h * u(n));
        end
        u(n+1) = (1 + (eps2*h^2)/2) * (u(n) + h*fn);
    elseif Meth == 2
        if abs(u(n)) <= h^p
            eps2 = sign(-(3*fn - 3*fn_1) / (h * u(n))) * L;
        else
            eps2 = -(3*fn - 3*fn_1) / (h * u(n));
        end
        u(n+1) = (1 - (eps2*h^2)/6) * u(n) + (h + (eps2*h^3)/6) * fn;
    else 
        warning('Method not available')
    end
end
err = abs(u(end) - u_es(b));
