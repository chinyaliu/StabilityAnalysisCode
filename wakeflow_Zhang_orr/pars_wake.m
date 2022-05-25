function [method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,k,Fr2,Re,eps,c0,h,f] = pars_wake(varargin)
%% Set Solver & Algorithm
diff_meth = ["Schimd", "Trefethen"];
method = diff_meth(1);
solveGEPmeth = ["qr", "qz", "eig", "invB"];
alg = solveGEPmeth(1);
baseflowlist = ["cosh"];
bflow = baseflowlist(1);
de_singularize = 'n';
do_balancing = 'n';
eig_spectrum = 'max';
N = 600;
k = 1.5;
% Re = 1000;
Re = inf;
Fr2 = 1.5^2;
% zL = 0.74708299;
h = 2*2*pi/real(k);
% ddm_number = 1;
ddm_number = 44;
c0 = 1./sqrt(k*Fr2);
f = wZhang_ddm.ddmtype(ddm_number);
eps = 0.1;
if ~isempty(varargin)
    switch(varargin{1})
        case(2)
            k = [0.3 1.1 1.5 3];
            eig_spectrum = 'max';
        case(3)
%             eig_spectrum = 'max';
%             k = linspace(0.01,4,150);
            k = [linspace(0.01,1,40) linspace(1.02,3,50) linspace(3.05,4,20)];
            h = 3*2*pi./real(k);
    end
end
end