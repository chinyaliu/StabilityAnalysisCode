classdef wZhang_complex < handle
    properties
        k = 1; h = -6; Re = inf; Fr2 = 0; del;
    end
    properties(SetObservable)
        N = 400; method;
    end
    properties (SetAccess = private)
        z; zc; phi;
    end
    properties (Access = private)
        zeta; Din; ord;
        g = @(x) (5000*acosh((-2497/(625*(x - 1)))^(1/2)/2))/4407;
    end
    methods
        function obj = wZhang_complex(N,k,h,Re,Fr2,meth,delt)
            if (nargin == 7)
                obj.N = N; obj.k = k; obj.h = h; obj.Re = Re;
                obj.Fr2 = Fr2; obj.method = lower(meth); obj.del = delt;
            end
            obj.diffmat();
            addlistener(obj,'N','PostSet',@obj.chgDM);
            addlistener(obj,'method','PostSet',@obj.chgDM);
        end
        [o, an, cA, errGEP, dob] = solver(obj, alg, bal);
    end
    methods (Access = private)
        function chgDM(obj, varargin)
            obj.diffmat();
        end
        diffmat(obj);
        U = baseflow(obj);
        [A, B] = matAB(obj, D, U);
    end
end