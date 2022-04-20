classdef subRay < subdomain 
    methods
        % Call class constructor of @subdomain
        function obj = subRay(N,ztop,zbot,supObj,varargin)
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
            A = [[k*obj.D(1,:,1); (k*obj.U(1,1).*obj.D(1,:,2)-k*obj.U(1,2).*obj.D(1,:,1))],zeros(2,N-obj.N-1),[k*obj.U(1,1); k/obj.pFlow.Fr2]];
            B = [zeros(1,N),1; obj.D(1,:,2),zeros(1,N-obj.N-1),0];
       end
        % Boundary condition at bottom (truncated, free-slip)
        function [A,B] = BChf(obj,N)
            A = [zeros(1,N-obj.N-1), obj.D(end,:,1), 0];
            B = zeros(1,N+1);
        end
        % Boundary condition at bottom (truncated, exponential decay)
        function [A,B] = BChe(obj,N)
            A = [zeros(1,N-obj.N-1), obj.D(end,:,2)-obj.pFlow.k*obj.D(end,:,1), 0];
            B = zeros(1,N+1);
        end
       function makeAB(obj)
           % Matrix A, B (Rayleigh equation)
            k = obj.pFlow.k;
            A_ge = k*obj.U(:,1).*obj.D(:,:,3) - (obj.U(:,1)*k^3 + obj.U(:,3)*k).*obj.D(:,:,1);
            B_ge = obj.D(:,:,3) - k^2*obj.D(:,:,1);
            obj.A = A_ge(2:end-1,:); obj.B = B_ge(2:end-1,:);
       end
    end
end