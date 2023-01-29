function unlap_vals = invlaplace_dehoog2(lap_fun,t,tmax_acc,Err,varargin)

%dehoog2: Function for calculating the inverse laplace transform given a
%solution in the laplace domain. Original source is de Hoog et al. (1982)
%"An improved method for numerical inversion of Laplace transforms."
%
%Usage:
%unlap_vals = invlaplace_dehoog2(lap_fun,t,tmax,Err)
%where:
%   -unlap_vals (scalar) is the calculated time-domain solution at the
%   times specified in t
%   -Lap_fun (anonymous function) is a function that takes as input a
%   single value, the laplace parameter, and gives as output the value of
%   the laplace transform solution
%   -t (scalar) are the times at which to compute the solution
%   -tmax (scalar) is a "late-time" defining parameter used to estimate the
%   error in the solution
%
%Solution by de Hoog et al., and originally used in code by Malama et al.
%for slug test solutions.
%Generalization of deHoog algorithm by M. Cardiff

nreqin = 4;
if nargin > nreqin
    Nk = varargin{1};
else
    Nk = 20;
end

%These parameters control the error associated with the solution. tmax
%should be on the order of the final / latest time for which a solution is
%desired. Nk determines the number of terms that will be included in the
%series summation - larger Nk will generally lead to better converged
%solutions.
%Values originally specified in code:
%Nk = 20;
%tmax = 20.0;
% Err = 1e-6;
%
%TODO: Discuss with Bwalya reasoning for choosing different Err values.


T = 0.8*tmax_acc;
beta = -log(Err)/(2*T);
alpha = 1i*pi*t/T;
z = exp(alpha);
lapf = zeros(Nk+1,1);
for k=1:Nk+1
    p = beta + 1i*(k-1)*pi/T;
    lapf(k) = lap_fun(p);
end
M = Nk/2;
eps1 = zeros(Nk+1,1) + 1i*zeros(Nk+1,1);
a = lapf; a(1) = lapf(1)/2.0;

q1 = zeros(Nk,1);
q1(1:1:Nk) = a(2:1:(Nk+1))./a(1:1:Nk);

d(1) = a(1);
d(2) = -q1(1);
eps2 = zeros(Nk-3,1);
q2 = zeros(Nk-4,1);
for r=2:M+1
    for n=1:(Nk-2*r+1)
        eps2(n)=q1(n+1)-q1(n)+eps1(n+1);
    end
    for n=1:(Nk-2*r)
        q2(n)=q1(n+1)*eps2(n+1)/eps2(n);
    end
    d(2*r-1) = -eps2(1);
    if(r<M+1)
        d(2*r) = -q2(1);
    end
    eps1(1:1:(Nk-2*r+1)) = eps2(1:1:(Nk-2*r+1));
    q1(1:1:(Nk-2*r)) = q2(1:1:(Nk-2*r));
end
A=zeros(Nk+2,1) + 1i*zeros(Nk+2,1);
A(2)=d(1);
B=ones(Nk+2,1) + 1i*zeros(Nk+2,1);
for n=3:(Nk+2)
    A(n)=A(n-1)+d(n-1)*z*A(n-2);
    B(n)=B(n-1)+d(n-1)*z*B(n-2);
end

unlap_vals = exp(beta*t)*(real(A(Nk+2)/B(Nk+2)))/T;