classdef (Abstract) subdomain < handle
    properties (SetAccess = protected)
        z; A; B; U; D;
        ztop; zbot; k; N; 
    end
    methods
        function obj = subdomain(N,ztop,zbot,dm,k,ud,delta,varargin)
            if (nargin == 7)
                obj.N = N; obj.ztop = ztop; obj.k = k;
                obj.zbot = zbot;
                obj.setD(ztop,zbot,dm);
                obj.baseflow(ud,delta);
            end
        end
        function setD(obj,ztop,zbot,dm)
            [zeta,Din] = dm(obj);
            sz = size(Din,3);
            w = (2/(zbot-ztop)).^(0:1:sz-1);
            obj.D = reshape((reshape(Din,[],sz).*w),size(Din,1),[],sz);
            obj.z = 0.5*(zbot-ztop)*(zeta-1)-ztop;
        end
        function chgL(obj,N,ztop,zbot,dm,ud,delta)
            if ztop > zbot
                error('Upper limit lower than lower limit.');
            end
            if obj.N~=N
                obj.N = N;
                obj.setD(ztop,zbot,dm);
            else
                sz = size(obj.D,3);
                w = (((obj.zbot-obj.ztop)/(zbot-ztop)).^(0:1:sz-1));
                obj.D = reshape((reshape(obj.D,[],sz).*w),size(obj.D,1),[],sz);
                obj.z = (obj.z+obj.ztop)*(zbot-ztop)/(obj.zbot-obj.ztop)-ztop;
            end
            obj.ztop = ztop; obj.zbot = zbot;
            obj.baseflow(ud,delta);
            obj.makeAB();
        end
        function phi = modeshape(obj,an)
            phi = reshape(reshape(permute(obj.D(:,:,1:3),[1 3 2]),[],size(obj.D,2))*an,[],3);
        end
    end
    methods(Abstract)
        baseflow(obj,ud,delta)      
        [A, B] = BC0(obj,N)
        [A, B] = BCh(obj,N)
        makeAB(obj)
    end
    methods (Static)
        function [A, B] = match(subd)
            ord = size(subd(1).D,3) - 1;
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