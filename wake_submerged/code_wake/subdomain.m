classdef (Abstract) subdomain < handle
    properties (SetAccess = protected)
        ztop; zbot; % Upper and lower bounds of subdomain (both > 0)
        N = 0;      % Number of collocation points in this subdomain
        z;          % Collocation points defined using Gauss-Lobatto grids
        U;          % Base flow velocity up to nth derivative at z
        D;          % Chebyshev differential matrix up to nth deriavtive at z
        A; B;       % RHS/LHS of the GEP
        pFlow;      % Properties of the flow (class handle)
    end
    methods
        % Class constructor, set flow properties
        function obj = subdomain(N,ztop,zbot,supObj,varargin)
            if (nargin == 4)
                obj.pFlow = supObj;
                obj.chgL(N,ztop,zbot);
            end
        end
        % Set Chebyshev derivative matrix (z, D)
        function setD(obj,ztop,zbot)
            [zeta,Din] = obj.pFlow.dm(obj);
            sz = size(Din,3);
            w = (2/(zbot-ztop)).^(0:1:sz-1);
            obj.D = reshape((reshape(Din,[],sz).*w),size(Din,1),[],sz);
            obj.z = 0.5*(zbot-ztop)*(zeta-1)-ztop;
        end
        % Change subdomain upper/lower limit and collocation point 
        function chgL(obj,N,ztop,zbot)
            if ztop > zbot
                error('Upper limit lower than lower limit.');
            end
            if obj.N~=N % Call @setD only when N is modified
                obj.N = N;
                obj.setD(ztop,zbot);
            else
                sz = size(obj.D,3);
                w = (((obj.zbot-obj.ztop)/(zbot-ztop)).^(0:1:sz-1));
                obj.D = reshape((reshape(obj.D,[],sz).*w),size(obj.D,1),[],sz);
                obj.z = (obj.z+obj.ztop)*(zbot-ztop)/(obj.zbot-obj.ztop)-ztop;
            end
            obj.ztop = ztop; obj.zbot = zbot;
            obj.U = obj.pFlow.baseflow(obj.z);
            obj.makeAB(); % Update matrix A, B with new flow properties
        end
        % Return 0-2nd order derivative of input eigenvector
        function phi = modeshape(obj,an)
            phi = reshape(reshape(permute(obj.D(:,:,1:3),[1 3 2]),[],size(obj.D,2))*an,[],3);
        end
    end
    methods(Abstract)
        % Boundary condition at surface
        [A, B] = BC0(obj,N)
        % Boundary condition at bottom
        [A, B] = BChf(obj,N)
        % Boundary condition at bottom
        [A, B] = BChe(obj,N)
        % Construct matrix A, B for the GEP
        makeAB(obj)
    end
    methods (Static)
        % Combine matrix A, B of the input subdomains and set the matching conditions
        function [A, B] = match(subd)
            ord = size(subd(1).D,3) - 1; % match the eigenvector up to "ord"-order
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