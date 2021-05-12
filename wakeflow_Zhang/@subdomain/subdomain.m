classdef subdomain < handle
    properties
        zeta; Din;
    end
    properties (SetAccess = protected)
        z; sz; A; B; U; D;
        ztop; zbot; k; N; 
    end
    properties(SetObservable)
        diffm;
    end
    methods
        function obj = subdomain(N,ztop,zbot,dm,k,varargin)
            if (nargin >= 5)
                obj.N = N; obj.ztop = ztop; obj.k = k;
                obj.zbot = zbot;obj.diffm = dm;
                dm(obj);
                obj.sz = size(obj.Din, 3);
                addlistener(obj,'diffm','PostSet',@obj.chgdm);
            end
        end
        function chgL(obj,N,ztop,zbot)
            if obj.N ~= N
                obj.N = N;
                obj.diffm(obj);
            end
            if ztop < zbot
                obj.ztop = ztop;
                obj.zbot = zbot;
                obj.makeAB();
            else
                error('Upper limit lower than lower limit.');
            end
        end
        function chgdm(obj, varargin)
            obj.diffm(obj);
        end
        function phi = modeshape(obj,an)
            phi = reshape(reshape(permute(obj.D(:,:,1:3),[1 3 2]),[],size(obj.D,2))*an,[],3);
        end
        function [A, B] = BC0(obj,Fr2,N)
            A = [[obj.k*obj.D(1,:,1); (obj.k*obj.U(1,1).*obj.D(1,:,2)-obj.k*obj.U(1,2).*obj.D(1,:,1))],zeros(2,N-obj.N-1),[obj.k*obj.U(1,1); obj.k/Fr2]];
            B = [zeros(1,N),1; obj.D(1,:,2),zeros(1,N-obj.N-1),0];
        end
        function [A,B] = BCh(obj,N)
            A = [zeros(1,N-obj.N-1), obj.D(end,:,2)-obj.k*obj.D(end,:,1), 0];
            B = zeros(1,N+1);
        end
        function makeAB(obj)
            % Baseflow
            zr = 0.5*(obj.zbot-obj.ztop)*(obj.zeta-1)-obj.ztop;
            c1 = 0.9988; c2 = 0.8814;
            Ur(:,1) = (1-c1*cosh(c2*zr).^(-2));
            Ur(:,2) = 2*c1*c2*tanh(c2*zr).*(sech(c2*zr)).^2;
            Ur(:,3) = 2*c1*c2^2*(sech(c2*zr)).^2.*((sech(c2*zr)).^2-2*(tanh(c2*zr)).^2);
            % Differential matrix
            w = (2/(obj.zbot-obj.ztop)).^(0:1:obj.sz-1);
            obj.D = reshape((reshape(obj.Din,[],obj.sz).*w),size(obj.Din,1),[],obj.sz);
            % Matrix A, B (Rayleigh)
            A_ge = obj.k*Ur(:,1).*obj.D(:,:,3) - (Ur(:,1)*obj.k^3 + Ur(:,3)*obj.k).*obj.D(:,:,1);
            B_ge = obj.D(:,:,3) - obj.k^2*obj.D(:,:,1);
            obj.A = A_ge(2:end-1,:); obj.B = B_ge(2:end-1,:);
            obj.U = Ur; obj.z = zr;
        end
    end
    methods (Static)
        function [A, B] = match(subd)
            ord = subd(1).sz - 1;
            ln = length(subd);
            A = blkdiag(subd.A);
            B = blkdiag(subd.B,zeros(ord*(ln-1),1));
            amat = cell(ln,1);
            for i = 1:length(subd)
                amat{i} = [zeros(ord*(i-2)*(i>1),subd(i).N+1);permute(-subd(i).D(1,:,1:ord),[3 2 1]);
                    permute(+subd(i).D(end,:,1:ord),[3 2 1]);zeros(ord*(ln-i-1)*(i<ln),subd(i).N+1)];
            end
            amat{1} = amat{1}(ord+1:end,:); amat{end} = amat{end}(1:end-ord,:);
            A = [[A; horzcat(amat{:})], zeros(size(B,1),1)];
        end
    end
end