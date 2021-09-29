classdef wMorland < handle
    properties
        ud = 2; delta = 1; k = 1; % Non-dimensional flow parameters
        h = 6;                    % Truncation height
        N = 400;                  % Number of collocation points
    end
    properties (SetAccess = private)
        criticalH;                % Inverse function of base flow (function handle)
        zc;                       % Critical height
        subD; subDclass;          % Subdomains and its corresponding class function 
        method;                   % Numerical methods, set by @numMeth
        ord;                      % Order of GE, 2 for Rayleigh, 4 for Orr-Sommerfeld
        dm;                       % Differential matrix constructing method (function handle)
    end
    methods
        % Class constructor, set flow properties
        function obj = wMorland(N, h, ud, delta, lambda, meth, bf)
            if (nargin >= 7)
                obj.N = N; obj.k = 2*pi/lambda; obj.h = h;
                obj.ud = ud; obj.delta = delta; 
                obj.numMeth(meth);
                obj.setbaseflow(bf);
            end
        end
        % Set velocity profile used in this case (subDclass, criticalH)
        function setbaseflow(obj,bftype)
            switch(lower(bftype))
                case 'exponential'
                    obj.subDclass = @subRay;
                    obj.criticalH = @(c) 0.5*obj.delta*log(c/obj.ud);
                case 'error function'
                    obj.subDclass = @subErf;
                    obj.criticalH = @(c) -erfcinv(c/obj.ud)*sqrt(pi)*obj.delta/2;
                otherwise
                    error('Undefined base flow name.');
            end
        end
        % Solve the stability equations of the flow
        [c, an, cA, errGEP, dob] = solver(obj, alg, bal, eigspec, funcN, addvar);
        [c, an, cA, errGEP, dob] = solvers(obj, alg, des, bal, eigspec, funcN, addvar);
        % Select numerical methods and governing equations
        numMeth(obj,meth);
        % Construct and modify subdomains with specified DD methods
        function setsubd(obj,init,funcN,addvar)
            [Na, arr] = funcN(obj, init, addvar);
            % Create subdomains according to (N, arr)
            if strcmpi(init,'init')
                obj.subD = repmat(obj.subDclass(),length(Na),1); % initialize
                for j = 1:length(Na)
                    obj.subD(j) = obj.subDclass(Na(j),arr(j),arr(j+1),obj);
                    obj.subD(j).makeAB();
                end
            else
                for j = 1:length(Na)
                    if j > length(obj.subD)
                        obj.subD(j) = obj.subDclass(Na(j),arr(j),arr(j+1),obj);
                        obj.subD(j).makeAB();
                    else
                        obj.subD(j).chgL(Na(j),arr(j),arr(j+1));
                    end
                end
                % Remove redundant subdomains
                if length(Na) < length(obj.subD)
                    obj.subD = obj.subD(1:length(Na));
                end
            end
        end
        % Combine the eigenvectors of the subdomains
        function [z, phi] = findmodeshape(obj, an)
            num = 0; z = []; phi = [];
            for i = 1:length(obj.subD)
                z = [z;obj.subD(i).z];
                phi = [phi;obj.subD(i).modeshape(an(num+1:num+obj.subD(i).N+1))];
                num = num + obj.subD(i).N + 1;
            end
        end
        % Return the domain decompose positions
        function out = getcut(obj)
            for i = 1:length(obj.subD)
                out(i) = obj.subD(i).ztop;
            end
            out(i+1) = obj.subD(end).zbot;
        end
        % Return the condition number of inverse for A, B
        [cA,cB,cAb,cBb] = calcond(obj, des, bal, funcN, addvar);
    end
    methods (Static)
        % Select the domain decompose strategy (meth, ord, dm)
        out = ddmtype(name);
    end
    methods (Access = private)
        % Construct matrix A, B for the GEP
        function [A, B] = makeAB(obj)
            % Governing equation
            [Age, Bge] = obj.subD(1).match(obj.subD);
            % BC (truncated, free surface)
            [Abc1, Bbc1] = obj.subD(1).BC0(size(Age,2)-1);
%             % BC (free slip)
%             [Abc2, Bbc2] = obj.subD(end).BCh(size(Age,2)-1);
%             A = [Age; Abc1; Abc2];
%             B = [Bge; Bbc1; Bbc2];
            % BC (truncated, exponential decay)
            [Abc2, Bbc2] = obj.subD(end).BCh3(size(Age,2)-1);
            A = [Age; Abc1; Abc2];
            B = [Bge; Bbc1; Bbc2];    
        end
    end
end