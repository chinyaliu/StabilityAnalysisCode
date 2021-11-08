function [method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,H,k,Fr2,Re,eps,c0,h,f] = pars_wake(varargin)
%% Set Solver & Algorithm
diff_meth = ["Schimd", "Trefethen"];
method = diff_meth(1);
solveGEPmeth = ["qr", "qz", "eig", "invB"];
alg = solveGEPmeth(1);
baseflowlist = ["cosh"];
bflow = baseflowlist(1);
de_singularize = 'y';
do_balancing = 'y';
eig_spectrum = 'all';
N = 800;
k = 0.04;
Re = inf;
Fr2 = 1.5^2;
% zL = 0.74708299;
H = 5;
h = 3*2*pi/real(k)+H;
% h = 12;
ddm_number = 3;
c0 = min(0.99,1./sqrt(k*Fr2));
f = wSubmerged.ddmtype(ddm_number);
eps = 0.1;
if ~isempty(varargin)
    switch(varargin{1})
        case(2)
            eig_spectrum = 'max';
            k = [0.1 0.4 0.8 1.5];
            h = 3*2*pi./real(k)+H;
        case(3)
            eig_spectrum = 'max';
            k = linspace(0.04,2.5,400);
            h = 3*2*pi./real(k)+H;
    end
end
end