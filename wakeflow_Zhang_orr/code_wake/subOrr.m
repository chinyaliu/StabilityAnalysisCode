classdef subOrr < subdomain
    properties (SetAccess = protected)
        Re = 1000;
    end
    methods
        function obj = subOrr(N,ztop,zbot,supObj,varargin)
            if (nargin >= 4)
                subargs = {N,ztop,zbot,supObj};
            else
                subargs = {};
            end
            obj@subdomain(subargs{:});
            if nargin >= 4
                obj.Re = supObj.Re;
            end
        end
        % Boundary conditions at free surface
        function [A, B] = BC0(obj,Fr2,N)
            Abc = [[obj.k*obj.D(1,:,1) obj.k*obj.U(1,1)];...
                [obj.D(1,:,3)+obj.k^2*obj.D(1,:,1) obj.U(1,3)];...
                [1i*obj.D(1,:,4)/obj.Re+(obj.k*obj.U(1,1)-3*1i*obj.k^2/obj.Re)*obj.D(1,:,2)-obj.k*obj.U(1,2)*obj.D(1,:,1) obj.k/Fr2]];
            A = [Abc(:,1:end-1) zeros(3,N-obj.N-1) Abc(:,end)];
            B = [zeros(1,N),1; zeros(1,N+1); obj.D(1,:,2),zeros(1,N-obj.N)];
        end
        % Boundary condition at bottom (truncated, free-slip)
        function [A, B] = BChf(obj,N)
            A = [zeros(2,N-obj.N-1), [obj.D(end,:,1); obj.D(end,:,3)], zeros(2,1)];
            B = zeros(2,N+1);
        end
        % Boundary condition at bottom (truncated, exponential decay)
        function [A,B] = BChe(obj,N)
            A = [zeros(2,N-obj.N-1), [obj.D(end,:,2)-obj.k*obj.D(end,:,1); obj.D(end,:,3)-obj.k*obj.D(end,:,2);], zeros(2,1)];
            B = zeros(2,N+1);
        end
        function makeAB(obj)
            % Matrix A, B (Orr-Sommerfeld equation)
            A_ge = 1i*obj.D(:,:,5)/obj.Re + (obj.k*obj.U(:,1) - 2*1i*obj.k^2/obj.Re).*obj.D(:,:,3) ...
                 + (1i*obj.k^4/obj.Re - obj.U(:,1)*obj.k^3 - obj.U(:,3)*obj.k).*obj.D(:,:,1);
            B_ge = obj.D(:,:,3) - obj.k^2*obj.D(:,:,1);
            obj.A = A_ge(2:end-1,:); obj.B = B_ge(2:end-1,:);
        end
    end
end