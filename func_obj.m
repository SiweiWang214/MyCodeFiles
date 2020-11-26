function obj = func_obj(y)

% The quadrupole strengths
x(1) = 16*y(1)-8;
x(2) = 16*y(2)-8;
x(3) = 16*y(3)-8;
x(4) = 16*y(4)-8;
x(5) = 16*y(5)-8;
x(6) = 16*y(6)-8;
x(7) = 16*y(7)-8;
% Drift and dipole length
x(8) = 2*y(8)+0.5; 
x(9) = 3*y(9)+1;
x(10) = 2*y(10)+0.5;
x(11) = 1*y(11)+0.3;
x(12) = 3*y(12)+1;
x(13) = 1*y(13)+0.3;
x(14) = 2*y(14)+1;
x(15) = 1*y(15)+0.5;
x(16) = 1*y(16)+0.5;

x(17) = 2*y(17)+1;
x(18) = 2*y(18)+1;
x(19) = y(19)-0.5;
x(20) = -0.5*y(20);

half_booster_v2;
global THERING
ring0 = THERING';

indqh01 = findcells(ring0,'FamName','qh01');
indqh02 = findcells(ring0,'FamName','qh02');
indqv01 = findcells(ring0,'FamName','qv01');
indqv02 = findcells(ring0,'FamName','qv02');
indqv03 = findcells(ring0,'FamName','qv03');
indqd01 = findcells(ring0,'FamName','qd01');
indqd02 = findcells(ring0,'FamName','qd02');

inddf = findcells(ring0,'FamName','df');
inddf0 = findcells(ring0,'FamName','df0');
inddf1 = findcells(ring0,'FamName','df1');
inddf2 = findcells(ring0,'FamName','df2');
inddf3 = findcells(ring0,'FamName','df3');
inddf4 = findcells(ring0,'FamName','df4');
inddf5 = findcells(ring0,'FamName','df5');
inddf6 = findcells(ring0,'FamName','df6');
inddf7 = findcells(ring0,'FamName','df7');

indbm05 = findcells(ring0,'FamName','bm05');
indbm10 = findcells(ring0,'FamName','bm10');

ring0 = atsetfieldvalues(ring0,indqh01,'PolynomB',[0 x(1)]);
ring0 = atsetfieldvalues(ring0,indqh02,'PolynomB',[0 x(2)]);
ring0 = atsetfieldvalues(ring0,indqv01,'PolynomB',[0 x(3)]);
ring0 = atsetfieldvalues(ring0,indqv02,'PolynomB',[0 x(4)]);
ring0 = atsetfieldvalues(ring0,indqv03,'PolynomB',[0 x(5)]);
ring0 = atsetfieldvalues(ring0,indqd01,'PolynomB',[0 x(6)]);
ring0 = atsetfieldvalues(ring0,indqd02,'PolynomB',[0 x(7)]);

ring0 = atsetfieldvalues(ring0,inddf,'Length',x(8));
ring0 = atsetfieldvalues(ring0,inddf0,'Length',x(9));
ring0 = atsetfieldvalues(ring0,inddf1,'Length',x(10));
ring0 = atsetfieldvalues(ring0,inddf2,'Length',x(11));
ring0 = atsetfieldvalues(ring0,inddf3,'Length',x(12));
ring0 = atsetfieldvalues(ring0,inddf4,'Length',x(13));
ring0 = atsetfieldvalues(ring0,inddf5,'Length',x(14));
ring0 = atsetfieldvalues(ring0,inddf6,'Length',x(15));
ring0 = atsetfieldvalues(ring0,inddf7,'Length',x(16));

ring0 = atsetfieldvalues(ring0,indbm10,'Length',x(17));
ring0 = atsetfieldvalues(ring0,indbm10,'PolynomB',[0 x(20)]);
ring0 = atsetfieldvalues(ring0,indbm05,'Length',x(18));
ring0 = atsetfieldvalues(ring0,indbm05,'PolynomB',[0 x(19)]);

% circum = findspos(ring0, length(ring0)+1);

[td,tune,chrom] = twissring(ring0,0,1:length(ring0)+1,'chrom');
betap = cat(1,td.beta);
displ = cat(2,td.Dispersion);

if isnan(tune(1)) || isnan(tune(2)) || isnan(betap(1,1)) || isnan(betap(1,2)) ||...
        ~isreal(betap(1,1)) || ~isreal(betap(1,2)) || ~isreal(tune(1)) || ~isreal(tune(2)) || any(displ(1,:)<0)
    obj = [1e8 1e8 1e8];
else
    [emitt,dampingtime] = calcdampingtime(ring0);
%     damp_max = 1e3*max(abs(dampingtime));
%     chrom_max = max(abs(chrom));
    d0 = displ(1,1);
%     dmax = max(displ(1,:));
    beta_m = max(max(betap));
%     circumf = findspos(ring0, length(ring0)+1);
    obj = [10*d0,1e9*abs(emitt),beta_m];
end