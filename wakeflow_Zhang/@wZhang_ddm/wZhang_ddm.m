classdef wZhang_ddm < handle
    properties
        k = 1; h = -6; Fr2 = 0; N = 400;
    end
    properties (SetAccess = private)
        zc; Re = inf; 
        method; ord; dm; subD;
    end
    properties (Constant)
        g = @(x) (5000*acosh((-2497./(625*(x - 1))).^(1/2)/2))./4407;
    end
    methods
        function obj = wZhang_ddm(N,k,h,Re,Fr2,meth)
            if (nargin == 6)
                obj.N = N; obj.k = k; obj.h = h; obj.Re = Re;
                obj.Fr2 = Fr2; obj.chgM(meth);
            end
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
        function out = getcut(obj)
            for i = 1:length(obj.subD)
                out(i) = obj.subD(i).ztop;
            end
            out(i+1) = obj.subD(end).zbot;
        end
        chgM(obj, meth);
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