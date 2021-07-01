classdef wZhang_ddm < handle
    properties
        k = 1; h = -6; Fr2 = 2.25; N = 400;
    end
    properties (SetAccess = private)
        zc; Re = inf; 
        method; ord; dm; subD;
    end
    methods
        function obj = wZhang_ddm(N,k,h,Re,Fr2)
            if (nargin >= 5)
                obj.N = N; obj.k = k; obj.h = h; obj.Re = Re;
                obj.Fr2 = Fr2; 
            end
        end
        [o, an, cA, errGEP, dob] = solver(obj, alg, bal, eigspec, funcN, addvar);
        numMeth(obj,meth);
        function chgRe(obj,Re)
            obj.Re = Re;
            if isinf(obj.Re)
                obj.method(1) = 'Ray';
                obj.ord = 2;
            else
                if strcmpi(obj.method(2),'trefethen')
                    error('Trefethen''s differential method can''t be used for M = N-2\n');
                end
                obj.method(1) = 'd4';
                obj.ord = 4;
            end
        end
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
    end
    methods (Static)
        out = ddmtype(name);
        function out = criticalH(x)
            out = -(5000*acosh((-2497./(625*(x - 1))).^(1/2)/2))./4407;
        end
    end
end