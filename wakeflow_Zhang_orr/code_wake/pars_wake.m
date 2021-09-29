function [method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,k,Fr2,Re,eps,c0,h,f] = pars_wake(varargin)
%% Set Solver & Algorithm
diff_meth = ["Schimd", "Trefethen"];
method = diff_meth(1);
solveGEPmeth = ["qr", "qz", "eig", "invB"];
alg = solveGEPmeth(1);
baseflowlist = ["exponential", "error function"];
bflow = baseflowlist(1);
de_singularize = 'y';
do_balancing = 'y';
eig_spectrum = 'all';
N = 600;
k = 1.5;
Re = inf;
Fr2 = 1.5^2;
% zL = 0.74708299;
h = 3*2*pi/real(k);
ddm_number = 44;
c0 = 1./sqrt(k*Fr2);
f = wZhang_ddm.ddmtype(ddm_number);
eps = 0.1;
if ~isempty(varargin)
    switch(varargin{1})
        case(2)
            eig_spectrum = 'max';
        case(3)
            eig_spectrum = 'max';
            k = linspace(0.01,4,100);
            h = 3*2*pi./real(k);
    end
end
end