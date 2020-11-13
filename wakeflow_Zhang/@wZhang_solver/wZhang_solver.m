classdef wZhang_solver < handle
    properties
        k = 1; h = -6; Re = inf; Fr2 = 0; 
    end
    properties( SetObservable )
        N = 400; method = ['','',''];
    end
    properties (SetAccess = private)
        z; zc; phi; zL = 1;
    end
    properties (Access = private)
        zeta; Din; BC; ord;
    end
    methods
        function obj = wZhang_solver(N,k,h,Re,Fr2,meth)
            if (nargin == 6)
                obj.N = N; obj.k = k; obj.h = h; obj.Re = Re;
                obj.Fr2 = Fr2; obj.method = meth; 
            end
            obj.diffmat();
            addlistener(obj,'N','PostSet',@obj.chgDM);
            addlistener(obj,'method','PostSet',@obj.chgDM);
        end
        [o, an, cA, errGEP, dob] = solver(obj, zL1, iter, alg, bal);
    end
    methods (Access = private)
        function chgDM(obj, varargin)
            obj.diffmat();
        end
        diffmat(obj);
        U = baseflow_zhang2(obj);
        [A, B] = matAB_zhang2(obj, D, U, w1, w2);
    end
end