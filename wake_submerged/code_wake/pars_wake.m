function [method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,H,k,Fr2,Re,eps,c0,h,f] = pars_wake(varargin)
%% Set Solver & Algorithm
diff_meth = ["Schimd", "Trefethen"];
method = diff_meth(1);
solveGEPmeth = ["qr", "qz", "eig", "invB"];
alg = solveGEPmeth(1);
baseflowlist = ["cosh", "pwlinear", "tanh"];
bflow = baseflowlist(1);
de_singularize = 'y';
do_balancing = 'y';
eig_spectrum = 'all';
N = 600;
k = 1.5;
% Re = 1000;
Re = inf;
Fr2 = 1.5^2;
H = 0;
% ddm_number = 1; % Re = 1000
ddm_number = 45; % Re = inf
f = wSubmerged.ddmtype(ddm_number);
eps = 0.1;
if ~isempty(varargin)
    switch(varargin{1})
        case(2)
%             k = [0.3 1.1 1.5 2]; % Re = 1000
            k = [0.3 1.1 1.5 3]; % Re = inf
%             k = [0.1 0.5 1.1 1.8];
        case(3)
%             k = [linspace(0.01,0.2,20) linspace(0.21,1,250) linspace(1.01,4,100)];
%             k = linspace(0.01, 4, 100);
            k = [linspace(0.01, 2, 150) linspace(2, 4, 100)];
    end
    if length(varargin)>1
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
end
c0 = min(0.99, 1./sqrt(real(k)*Fr2));
lambda = 2*pi./real(k);
h = 2*lambda+H;
end