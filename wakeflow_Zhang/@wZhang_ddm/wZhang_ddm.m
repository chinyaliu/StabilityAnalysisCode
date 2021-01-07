classdef wZhang_ddm < handle
    properties
        k = 1; h = -6; Re = inf; Fr2 = 0;
    end
    properties(SetObservable)
        N = 400; method;
    end
    properties (SetAccess = private)
        zc;
    end
    properties (Constant)
        g = @(x) (5000*acosh((-2497./(625*(x - 1))).^(1/2)/2))./4407;
    end
    properties (Access = private)
        ord; dm; subD;
    end
    methods
        function obj = wZhang_ddm(N,k,h,Re,Fr2,meth)
            if (nargin == 6)
                obj.N = N; obj.k = k; obj.h = h; obj.Re = Re;
                obj.Fr2 = Fr2; obj.method = lower(meth); obj.chgM();
            end
            addlistener(obj,'method','PostSet',@obj.chgM);
        end
        [o, an, cA, errGEP, dob] = solver(obj, alg, bal, funcN, zL1, eps);
        function [z, phi] = findmodeshape(obj, an)
            num = 0; z = []; phi = [];
            for i = 1:length(obj.subD)
                z = [z;obj.subD(i).z];
                phi = [phi;obj.subD(i).modeshape(an(num+1:num+obj.subD(i).N+1))];
                num = num + obj.subD(i).N + 1;
            end
        end
        function out = getarr(obj)
            for i = 1:length(obj.subD)
                out(i) = obj.subD(i).ztop;
            end
            out(i+1) = obj.subD(end).zbot;
        end
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
        diffmat(obj);
    end
    methods (Static)
        function out = ddmtype(name)
            switch(name)
                case 1
                    out = @wZhang_ddm.setN1sub;
                case 2
                    out = @wZhang_ddm.setN2sub;
                case 4
                    out = @wZhang_ddm.setN4sub;
                otherwise
                    error('Function for domain number not specified.');
            end
        end
        [N, arr] = setN4sub(obj,init,o,eps,varargin);
        function [N, arr] = setN2sub(obj,varargin)
            arr = [0 obj.zc obj.h];
            N = 0.5*obj.N*ones(1,2);
        end
        function [N, arr] = setN1sub(obj,varargin)
            arr = [0 obj.h];
            N = obj.N;
        end
    end
end