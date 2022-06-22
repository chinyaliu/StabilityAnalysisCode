classdef wMorland < handle
    properties
        ud = 2; delta = 1;        % Non-dimensional flow parameters
        k = 1;                    % Wavelength
        h = 6;                    % Truncation height
        N = 400;                  % Number of collocation points
    end
    properties (SetAccess = private)
        Re = 0;                 % Reynolds number
        zc;                       % Critical height
        subD; subDclass;          % Subdomains and its corresponding class function 
        method;                   % Numerical methods, set by @numMeth
        ord;                      % Order of GE, 2 for Rayleigh, 4 for Orr-Sommerfeld
        dm;                       % Differential matrix constructing method (function handle)
        baseflow;                 % Velocity profile (function handle)
        invbf;                    % Inverse function of base flow (function handle)
    end
    methods
        % Class constructor, set flow properties
        function obj = wMorland(N, h, ud, delta, lambda, meth, bf, Re)
            if (nargin >= 8)
                obj.N = N; obj.k = 2*pi/lambda; obj.h = h;
                obj.ud = ud; obj.delta = delta;
                obj.method = strings(1,2);
                obj.setprop('Re',Re);
                obj.numMeth(meth);
                obj.setbaseflow(bf);
            end
        end
        % Solve the stability equations of the flow
        [c, an, cA, errGEP, dob] = solver(obj, alg, des, bal, eigspec, funcN, addvar);
        % Solve the stability equations of the flow
        [c, an, cA, errGEP, dob] = solver2(obj, alg, des, bal, eigspec, funcN, addvar);
        % Select numerical methods and governing equations
        numMeth(obj,meth);
        % Return chosen properties of the flow
        function [out,varargout] = getprop(obj,prop,varargin)
            switch lower(prop)
            case 'u' % out:[z,U] 
                % Return the baseflow collocation points
                out = []; U = [];
                for i = 1:length(obj.subD)
                    out = [out;obj.subD(i).z];
                    U = [U;obj.subD(i).U];
                end
                varargout{1} = U;
            case 'cut' % out:cut positions
                % Return the domain decomposed locations
                out = nan(1,length(obj.subD)+1);
                for i = 1:length(obj.subD)
                    out(i) = obj.subD(i).ztop;
                end
                out(i+1) = obj.subD(end).zbot;
            case 'modeshape' % in:eigenvector, out:[z,phi] 
                an = varargin{1};
                num = 0; out = []; phi = [];
                for i = 1:length(obj.subD)
                    out = [out;obj.subD(i).z];
                    phi = [phi;obj.subD(i).modeshape(an(num+1:num+obj.subD(i).N+1))];
                    num = num + obj.subD(i).N + 1;
                end
                varargout{1} = phi;
            case 'lambda'
                out = 2*pi/obj.k;
            otherwise
                error('Invalid input for function getprop()');
            end
        end
        % Set properties of the flow
        function setprop(obj,varargin)
            if (mod(length(varargin),2) == 1 || isempty(varargin))
                error('Invalid number of inputs.');
            end
            nam = varargin(1:2:end);
            val = varargin(2:2:end);
            for i = 1:length(nam)
                switch nam{i}
                case 'Re'
                    if val{i}~=obj.Re
                        obj.Re = val{i};
                        if isinf(val{i})
                            obj.method(1) = "Ray";
                            obj.subDclass = @subRay;
                            obj.ord = 2;
                        else
                            obj.method(1) = "d4";
                            obj.subDclass = @subOrr;
                            obj.ord = 4;
                        end
                        obj.subD = obj.subDclass(); % initialize
                    end
                case 'method'
                    obj.numMeth(val{i});
                case 'baseflow'
                    obj.setbaseflow(val{i});
                case 'lambda'
                    obj.k = 2*pi/val{i};
                otherwise
                    if isprop(obj,nam{i})
                        obj.(nam{i}) = val{i};
                    else
                        error('Invalid property name.');
                    end
                end
            end
        end
        % Return the condition number of inverse for A, B
        [cA,cB,cAb,cBb] = calcond(obj, des, bal, funcN, addvar);
    end
    methods (Static)
        % Select the domain decompose strategy (meth, ord, dm)
        out = ddmtype(name);
    end
    methods (Access = private)
        % Set velocity profile used in this case (subDclass, invbf)
        setbaseflow(obj,bf);
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
        % Construct matrix A, B for the GEP
        function [A, B] = makeAB(obj)
            % Governing equation
            [Age, Bge] = obj.subD(1).match(obj.subD);
            % BC (truncated, free surface)
            [Abc1, Bbc1] = obj.subD(1).BC0(size(Age,2)-1);
%             % BC (free slip)
%             [Abc2, Bbc2] = obj.subD(end).BCh(size(Age,2)-1);
            % BC (truncated, exponential decay)
            [Abc2, Bbc2] = obj.subD(end).BCh3(size(Age,2)-1);
            A = [Age; Abc1; Abc2];
            B = [Bge; Bbc1; Bbc2];
        end
        % Construct matrix A, B for the GEP
        function [A, B] = makeAB2(obj)
            % Governing equation
            [Age, Bge] = obj.subD(1).match(obj.subD);
            % BC (truncated, free surface)
            [Abc1, Bbc1] = obj.subD(1).BC0(size(Age,2)-1);
             % BC (free slip)
            [Abc2, Bbc2] = obj.subD(end).BCh(size(Age,2)-1);
            A = [Age; Abc1; Abc2];
            B = [Bge; Bbc1; Bbc2];
        end
    end
end