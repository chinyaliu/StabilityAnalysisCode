function [w,xv,err,cA] = solvePoly(As,modenum,algorithm)
switch algorithm
    case 'qz'
        [V,ev,cA] = polyeig(As{:});
end
% find unstable modes
a = 1:length(ev);
if (strcmpi(modenum,'all'))
    ev_ind = a(abs(ev)<1e+3);
    [~, ind] = sort(imag(ev(ev_ind)),'descend');
    ev_ind = ev_ind(ind);
    w = ev(ev_ind);
    xv = V(:,ev_ind);
else
    ev_ind = a(abs(ev)<1e+3 & abs(ev)>1e-5 & abs(imag(ev)) > 1e-10);
    if isempty(ev_ind)
        w = NaN*(1+1i);
        xv = NaN(length(A),1)*(1+1i);
    else
        [~, ind] = sort(imag(ev(ev_ind)),'descend');
        ev_ind = ev_ind(ind);
        w = ev(ev_ind);
        xv = V(:,ev_ind);
        if (strcmpi(modenum,'max')) || (w(1) < 0)
            w = w(1);
            xv = xv(:,1);
        else
            w = w(w>0);
            xv = xv(:,w>0);
        end
    end
end
if nargout > 2
    err = max(abs(A*xv-B*xv*diag(w)));
end
end