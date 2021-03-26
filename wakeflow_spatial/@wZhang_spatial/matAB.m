function [A0, A1, A2, A3] = matAB(obj, D, U)
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
        % z = -h exponential decay
%         A0(end,:) = D(end,:,2);
%         A1(end,:) = -D(end,:,1);
%         A2(end,:) = 0;
%         A3(end,:) = 0;
    otherwise
        error('Invalid create AB matrix method name');
end
end