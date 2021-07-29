classdef wMorland < handle
    properties
        k = 1; h = -6; N = 400;
    end
    properties (SetAccess = private)
        zc; ud; delta; 
        method; ord; dm; 
        subD; subDclass;
        criticalH; 
    end
    methods
        function obj = wMorland(N, h, ud, delta, lambda, meth, bf)
            if (nargin >= 7)
                obj.N = N; obj.k = 2*pi/lambda; obj.h = h;
                obj.ud = ud; obj.delta = delta; 
                obj.criticalH = @(c) 0.5*delta*log(c/ud);
                obj.numMeth(meth);
                obj.setbaseflow(bf);
            end
        end
        function setbaseflow(obj,bftype)
            switch(lower(bftype))
                case 'exponential'
                    obj.subDclass = @subRay;
                case 'error function'
                    obj.subDclass = @subErf;
                otherwise
                    error('Undifined base flow name.');
            end
        end
        [o, an, cA, errGEP, dob] = solver(obj, alg, bal, eigspec, funcN, addvar);
        numMeth(obj,meth);
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
    end
end