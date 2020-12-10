classdef wZhang_block < handle
    properties
        k = 1; h = -6; Re = inf; Fr2 = 0; eps;
    end
    properties(SetObservable)
        N = 400; method;
    end
    properties (SetAccess = private)
        z; zc; phi; zL; dm; ord;
        g = @(x) (5000*acosh((-2497/(625*(x - 1)))^(1/2)/2))/4407;
    end
    methods
        function obj = wZhang_block(N,k,h,Re,Fr2,meth)
            if (nargin == 6)
                obj.N = N; obj.k = k; obj.h = h; obj.Re = Re;
                obj.Fr2 = Fr2; obj.method = lower(meth); 
            end
            obj.chgM();
            addlistener(obj,'N','PostSet',@obj.chgN);
            addlistener(obj,'method','PostSet',@obj.chgM);
        end
        [o, an, cA, errGEP, dob] = solver(obj, zL1, iter, alg, bal,eps);
    end
    methods (Access = private)
        function chgM(obj, varargin)
            switch obj.method(1)
                case 'ray'
                    obj.ord = 2;
                otherwise
                    error('Invalid method(1) name');
            end
            obj.diffmat();
        end
        function chgN(obj, varargin)
            obj.diffmat();
        end
        diffmat(obj);
        [A, B] = matAB(obj, D, U);
        [N, arr] = setN4sub(obj);
    end
end