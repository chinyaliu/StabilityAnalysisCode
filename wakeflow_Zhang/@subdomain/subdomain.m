classdef subdomain < handle
    properties
        zeta; Din;
    end
    properties (SetAccess = private)
        z; ord; A; B; U; D;
        ztop; zbot; k; N; diffm;
    end
    methods
        function obj = subdomain(N,ztop,zbot,dm,k)
            if (nargin == 5)
                obj.N = N; obj.ztop = ztop; obj.k = k;
                obj.zbot = zbot;obj.diffm = dm;
                obj.diffm(obj);
                obj.ord = size(obj.Din, 3);
                obj.makeAB();
            end
        end
        function chgL(obj,N,ztop,zbot)
            if obj.N ~= N
                obj.N = N;
                obj.diffm(obj);
            end
            obj.ztop = ztop;
            obj.zbot = zbot;
            obj.makeAB();
        end
        function phi = modeshape(obj,an)
            phi = reshape(reshape(permute(obj.D(:,:,1:3),[1 3 2]),[],size(obj.D,2))*an,[],3);
        end
        function [A, Aq, B, Bq] = BC0(obj,Fr2)
            A = [obj.k*obj.D(1,:,1); obj.k*obj.U(1,1).*obj.D(1,:,2) - obj.k*obj.U(1,2).*obj.D(1,:,1)];
            B = [zeros(1,obj.N+1); obj.D(1,:,2)];
            Aq = [obj.k*obj.U(1,1); obj.k/Fr2];
            Bq = [1;0];
        end
        function A = BCh(obj)
            A = obj.D(end,:,1);
        end
    end
    methods (Access = private)
        function makeAB(obj)
            % Baseflow
            obj.z = 0.5*(obj.zbot-obj.ztop)*(obj.zeta-1)-obj.ztop;
            c1 = 0.9988; c2 = 0.8814;
            Ur(:,1) = (1-c1*cosh(c2*obj.z).^(-2));
            Ur(:,2) = 2*c1*c2*tanh(c2*obj.z).*(sech(c2*obj.z)).^2;
            Ur(:,3) = 2*c1*c2^2*(sech(c2*obj.z)).^2.*((sech(c2*obj.z)).^2-2*(tanh(c2*obj.z)).^2);
            obj.U = Ur;
            % Differential matrix
            w = (2/(obj.zbot-obj.ztop)).^(0:1:obj.ord-1);
            obj.D = reshape((reshape(obj.Din,[],obj.ord).*w),size(obj.Din,1),[],obj.ord);
            % Matrix A, B (Rayleigh)
            A_ge = obj.k*obj.U(:,1).*obj.D(:,:,3) - (obj.U(:,1)*obj.k^3 + obj.U(:,3)*obj.k).*obj.D(:,:,1);
            B_ge = obj.D(:,:,3) - obj.k^2*obj.D(:,:,1);
            obj.A = [permute(-obj.D(1,:,1:2),[3 2 1]); A_ge(2:end-1,:); permute(obj.D(end,:,1:2),[3 2 1])];
            obj.B = B_ge(2:end-1,:);
        end
    end
end