function A = matAB(obj, D, U)
% D0 = D(:,:,1); D1 = D(:,:,2); D2 = D(:,:,3); D3 = D(:,:,4); D4 = D(:,:,5);
% U = U(:,1); Uz = U(:,2); Uzz = U(:,3); 
N = obj.N; o = obj.o;
switch obj.method(1)
    case 'ray'
    %% Rayleigh
        % Governing equation
        A0 = o*D(:,:,3);
        A1 = U(:,3).*D(:,:,1) - U(:,1).*D(:,:,3);
        A2 = -o.*D(:,:,1);
        A3 = U(:,1).*D(:,:,1);
        % FSBC
        A0(1,:) = D(1,:,2).*o^2;
        A1(1,:) = -2*D(1,:,2).*U(1,1)*o;
        A2(1,:) = D(1,:,2).*U(1,1).^2 - D(1,:,1)./obj.Fr2;
        A3(1,:) = 0;
        % z = -h
        A0(end,:) = D(end,:,1);
        A1(end,:) = 0;
        A2(end,:) = 0;
        A3(end,:) = 0;
        A = {A0,A1,A2,A3};
        % z = -h exponential decay
%         A0(end,:) = D(end,:,2);
%         A1(end,:) = -D(end,:,1);
%         A2(end,:) = 0;
%         A3(end,:) = 0;
    case 'ray2'
    %% Rayleigh
        A = cell(1,4);
        B = cell(1,4);
        BC1 = zeros(4,size(U,1)+1);
        BC2 = BC1; BC3 = BC1;
        % Governing 
        B{1} = o*D(:,:,3);
        B{2} = U(:,3).*D(:,:,1) - U(:,1).*D(:,:,3);
        B{3} = -o.*D(:,:,1);
        B{4} = U(:,1).*D(:,:,1);
        % FSBC
        BC1(1,end) = o;
        BC1(2,:) = [-D(1,:,1) -U(1,1)];
        BC2(1,:) = [o*D(1,:,1) 0];
        BC2(2,:) = [-U(1,1).*D(1,:,2) -1/obj.Fr2];
        % z = -h
        BC3(1,:) = [D(end,:,1) 0];
        for i = 1:4
            A{i} = [[B{i}(2:end-1,:),zeros(size(U,1)-2,1)];BC1(i,:);BC2(i,:);BC3(i,:)];
        end
    otherwise
        error('Invalid create AB matrix method name');
end
end