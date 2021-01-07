classdef subdomain < handle
    properties
        zeta; Din;
    end
    properties (SetAccess = private)
        z; ord; A; B; U; D;
        ztop; zbot; k; N; 
    end
    properties(SetObservable)
        diffm;
    end
    methods
        function obj = subdomain(N,ztop,zbot,dm,k)
            if (nargin == 5)
                obj.N = N; obj.ztop = ztop; obj.k = k;
                obj.zbot = zbot;obj.diffm = dm;
                dm(obj);
                obj.ord = size(obj.Din, 3);
                obj.makeAB();
                addlistener(obj,'diffm','PostSet',@obj.chgdm);
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
        function chgdm(obj, varargin)
            obj.diffm(obj);
        end
        function phi = modeshape(obj,an)
            phi = reshape(reshape(permute(obj.D(:,:,1:3),[1 3 2]),[],size(obj.D,2))*an,[],3);
        end
        function [A, B] = BC0(obj,Fr2,N)
            A = [[obj.k*obj.D(1,:,1); obj.k*obj.U(1,1).*obj.D(1,:,2) - obj.k*obj.U(1,2).*obj.D(1,:,1)],zeros(2,N-obj.N-1),[obj.k*obj.U(1,1); obj.k/Fr2]];
            B = [zeros(1,N),1; obj.D(1,:,2),zeros(1,N-obj.N-1),0];
        end
        function A = BCh(obj)
            A = obj.D(end,:,1);
        end
    end
    methods (Access = private)
        function makeAB(obj)
            % Baseflow
            zr = 0.5*(obj.zbot-obj.ztop)*(obj.zeta-1)-obj.ztop;
            c1 = 0.9988; c2 = 0.8814;
            Ur(:,1) = (1-c1*cosh(c2*zr).^(-2));
            Ur(:,2) = 2*c1*c2*tanh(c2*zr).*(sech(c2*zr)).^2;
            Ur(:,3) = 2*c1*c2^2*(sech(c2*zr)).^2.*((sech(c2*zr)).^2-2*(tanh(c2*zr)).^2);
            % Differential matrix
            w = (2/(obj.zbot-obj.ztop)).^(0:1:obj.ord-1);
            obj.D = reshape((reshape(obj.Din,[],obj.ord).*w),size(obj.Din,1),[],obj.ord);
            % Matrix A, B (Rayleigh)
            A_ge = obj.k*Ur(:,1).*obj.D(:,:,3) - (Ur(:,1)*obj.k^3 + Ur(:,3)*obj.k).*obj.D(:,:,1);
            B_ge = obj.D(:,:,3) - obj.k^2*obj.D(:,:,1);
            obj.A = A_ge(2:end-1,:); obj.B = B_ge(2:end-1,:);
            obj.U = Ur; obj.z = zr;
        end
    end
    methods (Static)
        function [A, B] = match(sub)
            ln = 2*length(sub);
            A = blkdiag(sub.A,zeros(ln+1,1));
            B = blkdiag(sub.B,zeros(ln+1,1));
            amat = cell(ln/2,1);
            for i = 1:length(sub)
                amat{i} = [zeros(2*(i-2)*(i>2),sub(i).N+1);permute(-sub(i).D(1,:,1:2),[3 2 1]);
                    permute(+sub(i).D(end,:,1:2),[3 2 1]);zeros((ln-2*i-2)*(2*i<ln-2),sub(i).N+1)];
            end
            amat{1} = amat{1}(3:end,:); amat{end} = amat{end}(1:end-2,:);
            A(end-ln:end-3,1:end-1) = horzcat(amat{:});
        end
    end
end