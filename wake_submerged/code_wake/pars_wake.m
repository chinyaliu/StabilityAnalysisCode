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
k = 1;
Re = inf;
% Re = 1000;
Fr2 = 3^2;
H = 0;
ddm_number = 44;
c0 = min(0.99, 1./sqrt(k*Fr2));
f = wSubmerged.ddmtype(ddm_number);
eps = 0.1;
if ~isempty(varargin)
    switch(varargin{1})
        case(2)
            k = [0.05 0.25 0.8];
        case(3)
            eig_spectrum = 'max';
            k = [linspace(0.01,2.6,100) linspace(2.8, 4, 10)];
%             k = linspace(0.01, 4, 100);
    end
end
h = 3*2*pi./real(k)+H;
end