function c_chosen = filt(c)
    a = 1:length(c); 
    crange = ((real(c)-0.0012>-1e-5) & (real(c)-1<=1e-5));
    dis = ((real(c)-1).^2 +imag(c).^2)>1e-4;
    aa = a(crange&dis);
    abch = isoutlier(imag(c(aa)),'movmedian',50);
    cgood = abs(c)<10;
    aa = [a(dis&~crange&cgood) aa(abch)];
    c_chosen = c(aa);
    c_chosen = c_chosen(~isinf(c_chosen));
end