classdef subRay < subdomain 
    methods
        function obj = subRay(N,ztop,zbot,dm,k,ud,delta,varargin)
            if (nargin >= 7)
                subargs = {N,ztop,zbot,dm,k,ud,delta};
            else
                subargs = {};
            end
            obj@subdomain(subargs{:});
        end
       function [A, B] = BC0(obj,N)
           c0 = obj.k^2/4/pi+pi;
           A = [[obj.D(1,:,1); (obj.U(1,1).*obj.D(1,:,2)-obj.U(1,2).*obj.D(1,:,1))],zeros(2,N-obj.N-1),[obj.U(1,1); c0]];
           B = [zeros(1,N),1; obj.D(1,:,2),zeros(1,N-obj.N-1),0];
       end
       function [A,B] = BCh(obj,N)
           A = [zeros(1,N-obj.N-1), obj.D(end,:,1), 0];
           B = zeros(1,N+1);
       end
       function [A,B] = BCh2(obj,N)
           A = [[zeros(1,N-obj.N-1), obj.D(end,:,1), 0, -exp(-obj.k*obj.zbot)];...
                [zeros(1,N-obj.N-1), obj.D(end,:,2)-obj.k*obj.D(end,:,1), 0, 0]];
           B = zeros(2,N+2);
       end
       function [A,B] = BCh3(obj,N)
           A = [zeros(1,N-obj.N-1), obj.D(end,:,2)-obj.k*obj.D(end,:,1), 0];
           B = zeros(1,N+1);
       end
       function makeAB(obj)
           % Matrix A, B (Rayleigh equation)
           A_ge = obj.U(:,1).*obj.D(:,:,3) - (obj.U(:,1)*obj.k^2 + obj.U(:,3)).*obj.D(:,:,1);
           B_ge = obj.D(:,:,3) - obj.k^2*obj.D(:,:,1);
           obj.A = A_ge(2:end-1,:); obj.B = B_ge(2:end-1,:);
       end
   end
    methods(Static)
        function Ur = baseflow(z,ud,delta)      
            % exponential profile
            Ur(:,1) = ud*exp(2*z/delta);
            Ur(:,2) = (2/delta)*Ur(:,1);
            Ur(:,3) = (2/delta)*Ur(:,2);
        end
    end
end