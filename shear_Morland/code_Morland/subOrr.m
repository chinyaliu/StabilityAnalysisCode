classdef subOrr < subdomain
    methods
        function obj = subOrr(N,ztop,zbot,supObj,varargin)
            if (nargin >= 4)
                subargs = {N,ztop,zbot,supObj};
            else
                subargs = {};
            end
            obj@subdomain(subargs{:});
        end
        % Boundary conditions at free surface
        function [A, B] = BC0(obj,N)
            k = obj.pFlow.k;
            Re = obj.pFlow.Re;
            c0 = obj.pFlow.k^2/4/pi+pi; % Non-dimensional constant in DFSBC
            Abc = [[k*obj.D(1,:,1) k*obj.U(1,1)];...
                [obj.D(1,:,3)+k^2*obj.D(1,:,1) obj.U(1,3)];...
                [1i*obj.D(1,:,4)/Re+(k*obj.U(1,1)-3*1i*k^2/Re)*obj.D(1,:,2)-k*obj.U(1,2)*obj.D(1,:,1) k*c0]];
            A = [Abc(:,1:end-1) zeros(3,N-obj.N-1) Abc(:,end)];
            B = [zeros(1,N),1; zeros(1,N+1); obj.D(1,:,2),zeros(1,N-obj.N)];
        end
        % Boundary condition at bottom (truncated, free-slip)
        function [A, B] = BCh(obj,N)
            A = [zeros(2,N-obj.N-1), [obj.D(end,:,1); obj.D(end,:,3)], zeros(2,1)];
            B = zeros(2,N+1);
        end
        % Boundary condition at bottom (truncated, exponential decay)
        function [A,B] = BCh3(obj,N)
            A = [zeros(2,N-obj.N-1), [obj.D(end,:,2)-obj.pFlow.k*obj.D(end,:,1); obj.D(end,:,3)-obj.pFlow.k*obj.D(end,:,2);], zeros(2,1)];
            B = zeros(2,N+1);
        end
        function makeAB(obj)
            % Matrix A, B (Orr-Sommerfeld equation)
            k = obj.pFlow.k;
            Re = obj.pFlow.Re;
            A_ge = 1i*obj.D(:,:,5)/Re + (k*obj.U(:,1) - 2*1i*k^2/Re).*obj.D(:,:,3) ...
                 + (1i*k^4/Re - obj.U(:,1)*k^3 - obj.U(:,3)*k).*obj.D(:,:,1);
            B_ge = obj.D(:,:,3) - k^2*obj.D(:,:,1);
            obj.A = A_ge(2:end-1,:); obj.B = B_ge(2:end-1,:);
        end
    end
end