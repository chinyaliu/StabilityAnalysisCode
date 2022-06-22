function f = pltneutral(ud,dltlist,varargin)
if ~isempty(varargin)
    f = varargin{1};
    ax = get(f,'CurrentAxes');
else
    f = figure;
    ax = gca;
end
hold(ax,'on');
for i = 1:length(dltlist)
    [lam1,lam2] = findneutral(ud,dltlist(i));
    plot(ax, dltlist(i),lam1,'k.','Markersize',4);
    plot(ax, dltlist(i),lam2,'k.','Markersize',4);
end
hold(ax,'off');
end