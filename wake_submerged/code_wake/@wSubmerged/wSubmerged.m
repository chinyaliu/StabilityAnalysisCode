classdef wSubmerged < handle
    properties
        H = 5;                    % Submerge depth
        k = 1;                    % Wavenumber
        Fr2 = 2.25;               % Squared Froude number
        h = 6;                    % Truncation height
        N = 400;                  % Number of collocation points
    end
    properties (SetAccess = private)
        Re = inf;                 % Reynolds number
        zc;                       % Critical height
        subD; subDclass;          % Subdomains and its corresponding class function 
        method;                   % Numerical methods, set by @numMeth
        ord;                      % Order of GE, 2 for Rayleigh, 4 for Orr-Sommerfeld
        dm;                       % Differential matrix constructing method (function handle)
    end
    methods
        function obj = wSubmerged(N,H,k,h,Re,Fr2)
            if (nargin >= 6)
                obj.N = N; obj.k = k; obj.h = h; obj.Re = Re;
                obj.Fr2 = Fr2; obj.H = H; 
            end
            obj.method = strings(1,2);
        end
        % Solve the stability equations of the flow
        [c, an, cA, errGEP, dob] = solver(obj, alg, des, bal, eigspec, funcN, addvar);
        [c, an, cA, errGEP, dob] = solver2(obj, alg, des, bal, eigspec, funcN, addvar);
        % Solve the stability equations of the flow
        [c, an, cA, errGEP, dob] = solverGPU(obj, alg, des, bal, eigspec, funcN, addvar);
        % Select numerical methods and governing equations
        numMeth(obj,meth);
        % Critical height of given phase velocity
        function out = criticalH(obj,x)
            if x <= 0.0012
                out = 0;
            else
                out = obj.H+abs((5000*acosh((-2497./(625*(x-1))).^(1/2)/2))./4407);
            end
        end
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
                        obj.subD(j).chgL(Na(j),arr(j),arr(j+1),obj.H);
                    end
                end
                % Remove redundant subdomains
                if length(Na) < length(obj.subD)
                    obj.subD = obj.subD(1:length(Na));
                end
            end
        end
        % Set the Reynolds number of the flow
        function chgRe(obj,Re)
            obj.method(1) = "d4";
            obj.Re = Re;
            if isinf(obj.Re)
                obj.method(1) = "Ray";
                obj.subDclass = @subRay;
                obj.ord = 2;
            else
                obj.subDclass = @subOrr;
                obj.ord = 4;
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
        % Return the baseflow
        function [z, U] = getU(obj)
            num = 0; z = []; U = [];
            for i = 1:length(obj.subD)
                z = [z;obj.subD(i).z];
                U = [U;obj.subD(i).U];
            end
        end
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
            [Abc1, Bbc1] = obj.subD(1).BC0(obj.Fr2,size(Age,2)-1);
%             % BC (free slip)
%             [Abc2, Bbc2] = obj.subD(end).BChf(size(Age,2)-1);
%             A = [Age; Abc1; Abc2];
%             B = [Bge; Bbc1; Bbc2];
            % BC (truncated, exponential decay)
            [Abc2, Bbc2] = obj.subD(end).BChe(size(Age,2)-1);
            A = [Age; Abc1; Abc2];
            B = [Bge; Bbc1; Bbc2];    
        end
        function [A, B] = makeAB2(obj)
            % Governing equation
            [Age, Bge] = obj.subD(1).match(obj.subD);
            % BC (truncated, free surface)
            [Abc1, Bbc1] = obj.subD(1).BC0(obj.Fr2,size(Age,2)-1);
            % BC (free slip)
            [Abc2, Bbc2] = obj.subD(end).BChf(size(Age,2)-1);
            A = [Age; Abc1; Abc2];
            B = [Bge; Bbc1; Bbc2];
        end
    end
end