function [method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,H,k,Fr2,Re,eps,c0,h,f] = pars_wake(varargin)
%% Set Solver & Algorithm
diff_meth = ["Schimd", "Trefethen"];
method = diff_meth(1);
solveGEPmeth = ["qr", "qz", "eig", "invB"];
alg = solveGEPmeth(1);
baseflowlist = ["cosh"];
bflow = baseflowlist(1);
de_singularize = 'n';
do_balancing = 'y';
eig_spectrum = 'max';
N = 600;
k = 3;
Re = 1000;
% Re = inf;
Fr2 = 1.5^2;
H = 0;
ddm_number = 1; % Re = 1000
% ddm_number = 45; % Re = inf
f = wSubmerged.ddmtype(ddm_number);
eps = 0.1;
if ~isempty(varargin)
    switch(varargin{1})
        case(2)
%             k = [0.3 1.1 1.5 2]; % Re = 1000
            k = [0.3 1.1 1.5 3]; % Re = inf
        case(3)
            k = linspace(0.01, 4, 200);
    end
    switch(lower(varargin{2}))
        case 'inv'
            Re = inf;
            f = wSubmerged.ddmtype(45);
            alg = solveGEPmeth(1);
            do_balancing = 'y';
            de_singularize = 'n';
            N = 600;
        case 'vis'
            Re = 1000;
            f = wSubmerged.ddmtype(1);
            alg = solveGEPmeth(1);
            do_balancing = 'y';
            de_singularize = 'n';
            N = 600;
    end
end
c0 = min(0.99, 1./sqrt(real(k)*Fr2));
h = 3*2*pi./real(k)+H;
end