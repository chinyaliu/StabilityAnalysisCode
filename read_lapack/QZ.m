function [w,xv] = QZ(A,B,algorithm,modenum,includenegative)
% do QZ alogarithm and filter the value we need

% use eig function
switch algorithm
    case 'qr'
        [V,D] = eig(A\B);
        ev = 1./diag(D);
    case 'qz'
        % same as LAPACK ZGGEV
        [AA,BB,~,~,V,~] = qz(A,B);
        ev = diag(AA)./diag(BB);
    case 'eig'
        [V,D] = eig(A,B,'qz');
        ev = diag(D);
    case 'eigs'
        % Arnoldi method
        [V,D] = eigs(A\B,length(A)-2);
        ev = diag(D);
        % eigenvector uncheck
    case 'polyeig'
        [V,e] = polyeig(A,-B);
        ev = e;
        % eigenvector uncheck
    case 'singgep'
        % MultiParEig package
        [ev,~,~,V,~,~,~,~,~] = singgep(A,B,0);
        % eigenvector uncheck
    case 'jdqz'
        options = [];
        options.disp = 1;
        ev = jdqz(A,B,10,'SM',options);
end
% figure(1), plot(ev,'+'), hold on

% filter the infinity value
e1 = ev(abs(ev)<1e+3 & abs(ev)>1e-5 & imag(ev)~=0);
e2 = e1(imag(e1)>0);
% if isempty(e2)~=1, Z(Z(:,2)==e2(1),3), else NaN, end
%
w = NaN(modenum,1)*(1+1i);
xv = NaN(length(A),1,modenum)*(1+1i);
ind = zeros(length(modenum),1);

if isempty(e1) == 0
    switch includenegative
        case 'y'
            wisort = sort(imag(e1),'descend');
            for i = 1:length(wisort)
                ind(i) = find(imag(ev) == wisort(i));
            end
            if length(ind) > modenum
                ii = modenum;
            else
                ii = length(ind);
            end
            w(1:ii) = ev(ind(1:ii));
            xv(:,:,1:ii) = V(:,ind(1:ii),:);
        case 'n'
            if isempty(e2) == 0
                wisort = sort(imag(e2),'descend');
                for i = 1:length(e2)
                    ind(i) = find(imag(ev) == wisort(i));
                end
                if length(ind) > modenum
                    ii = modenum;
                else
                    ii = length(ind);
                end
                w(1:ii) = ev(ind(1:ii));
                xv(:,:,1:ii) = V(:,ind(1:ii),:);
            end
    end
end
end