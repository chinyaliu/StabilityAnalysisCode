clear;
load('saveh5k004.mat');

% fr2_list = linspace(0.5,10,100).^2;
k_list = linspace(0.04,2.5,100);
cmodes = nan(4,length(k_list),length(fr2_list));
% h = 3*2*pi./real(k_list)+H;
zc = -0.74708299-H;

tic;
parfor ii = 1:length(fr2_list)
    p1 = wSubmerged(N,H,1,1,Re,fr2_list(ii));
    p1.numMeth(method);
    cs = css(:,ii);
    csall = findmode(cs,zc,H,p1,alg, de_singularize, do_balancing, eig_spectrum, f, Re);
    cmodes(:,:,ii) = csall;
end
toc;

%% plot
fr = fr2_list.^(0.5);
[X,Y] = meshgrid(k_list,fr);
oi4 = squeeze(imag(cmodes(4,:,:)));
contourf(X,Y,oi4);

function csall = findmode(cs,zc,H,p1,alg, de_singularize, do_balancing, eig_spectrum, f, Re)
    addvar = struct('zL1',zc);
    k_list = linspace(0.04,2.5,100);
    h = 3*2*pi./real(k_list)+H;
    csall = nan(length(cs),length(k_list));
    for i = 1:length(k_list)
        p1.k = k_list(i); p1.h = h(i);
        oall = p1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
        oall = oall(real(oall)>-50); % Remove the eigenvalues assigned by de-singularizing
        o = maxeig(oall);
        if real(o) > 0
            addvar.zL1 = p1.criticalH(real(o)/k_list(i));
        end

        call = oall/k_list(i);

        if isinf(Re)
            a = 1:length(call); 
            crange = ((real(call)-0.0012>-1e-5) & (real(call)-1<=1e-5));
            dis = ((real(call)-1).^2 +imag(call).^2)>1e-5;
            aa = a(crange&dis);
            abch = isoutlier(imag(call(aa)),'movmedian',20);
            aa = [a(dis&~crange) aa(abch)];
            cnind = call(aa);
        else
            if (i~=1) && (abs(k-0.8)>1e-5)
                cmatn = repmat(calln,1,length(call));
                cmat = repmat(call.',length(calln),1);
                [~,ind] = min(abs(cmatn-cmat),[],2);
                cnind = call(setdiff(1:end,ind));
                calln = call(ind);
            else
                cnind = call(imag(call)>-2);
                c1 = repmat(cs.',1,length(cnind));
                c2 = repmat(cnind.',length(cs),1);
                [~,ind] = min(abs(c1-c2),[],2);
                calln = cnind(setdiff(1:end,ind));
            end
        end

        % Plot selected eigenvalues
        for j = 1:length(cs)
            [~,ind] = min(abs(cnind-cs(j)));
            c = cnind(ind);
            [~,ind2] = min(abs(c-cs));
            if ind2 == j
                cs(j) = c;
                csall(j,i) = c;
            end
        end
    end
end

function o = maxeig(ev)
    a = 1:length(ev);
    ev_ind = a(abs(ev)<1e+3 & abs(ev)>1e-5 & abs(imag(ev)) > 1e-10);
    [~, ind] = sort(imag(ev(ev_ind)),'descend');
    ev_ind = ev_ind(ind);
    w = ev(ev_ind);
    if isempty(w)
        o = NaN;
    else
        o = w(1);
    end
end