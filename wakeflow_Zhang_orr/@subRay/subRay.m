classdef subRay < subdomain 
    methods
        function obj = subRay(N,ztop,zbot,dm,k,varargin)
            if (nargin >= 5)
                subargs = {N,ztop,zbot,dm,k};
            else
                subargs = {};
            end
            obj@subdomain(subargs{:});
        end
       function [A, B] = BC0(obj,Fr2,N)
           A = [[obj.k*obj.D(1,:,1); (obj.k*obj.U(1,1).*obj.D(1,:,2)-obj.k*obj.U(1,2).*obj.D(1,:,1))],zeros(2,N-obj.N-1),[obj.k*obj.U(1,1); obj.k/Fr2]];
           B = [zeros(1,N),1; obj.D(1,:,2),zeros(1,N-obj.N-1),0];
       end
       function [A,B] = BCh(obj,N)
%            % Exponential decay
%            A = [zeros(1,N-obj.N-1), obj.D(end,:,2)-obj.k*obj.D(end,:,1), 0];
           % Free slip
           A = [zeros(1,N-obj.N-1), obj.D(end,:,1), 0];
           B = zeros(1,N+1);
       end
       function makeAB(obj)
           % Matrix A, B (Rayleigh equation)
           A_ge = obj.k*obj.U(:,1).*obj.D(:,:,3) - (obj.U(:,1)*obj.k^3 + obj.U(:,3)*obj.k).*obj.D(:,:,1);
           B_ge = obj.D(:,:,3) - obj.k^2*obj.D(:,:,1);
           obj.A = A_ge(2:end-1,:); obj.B = B_ge(2:end-1,:);
       end
   end
end