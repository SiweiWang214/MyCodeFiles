function [xx,zz,area]=atdynap_par(ring,nt,dpp,np,rfrac)
%ATDYNAP		Compute the dynamic aperture
%
%
%[XX,ZZ]=ATDYNAP(RING,NTURNS,DPP,RFRAC)
%
%XX,ZZ :	limit of the dynamic aperture (betatron amplitudes in m)
%RING :		Structure for tracking
%NTURNS:	Number of turns
%DPP :		Off-momentum value (default: 0)
%RFRAC :	Resolution of the grid for checking the stability
%			as a fraction of the maximum stable amplitude
%			(default: 0.02)
% Modified by S.W.Wang in Nov.14,2019, supporting parallel calculation.

rlist=0:0.001:0.1;
if nargin < 5, rfrac=0.02; end
if nargin < 4, np=5; end
if nargin < 3, dpp=0.0; end

if isnumeric(dpp)
clorb=[findorbit4(ring,dpp);dpp;0];
else
   clorb=findorbit6(ring);
end

t1=linspace(0,pi,2*np+3);
xpmax=ascan(ring,nt,clorb,0,rlist);
zmax=ascan(ring,nt,clorb,0.5*pi,rlist);
xmmax=ascan(ring,nt,clorb,pi,rlist);

% slist=[0:0.05:0.4,0.42:0.02:0.7,0.71:0.01:1];
slist = rfrac : rfrac : 1.5;
xx=NaN(2*np+3,1);
x0=zeros(2*np+3,1);
z0 =x0;
zz=xx;
for i=1:np+3
    x0(i) = xpmax*cos(t1(i));
    z0(i) = zmax*sin(t1(i));
   [xx(i),zz(i)]=bscan(ring,nt,clorb,xpmax*cos(t1(i))*slist,zmax*sin(t1(i))*slist);
end
for i=np+4:2*np+3
    x0(i) = xmmax*cos(t1(i));
    z0(i) = zmax*sin(t1(i));
   [xx(i),zz(i)]=bscan(ring,nt,clorb,xmmax*cos(t1(i))*slist,zmax*sin(t1(i))*slist);
end
area=0;
for i=1:length(xx)-1
    area=area+1/2*abs(xx(i)*zz(i+1)-zz(i)*xx(i+1));
end
disp(['DA area is: ',num2str(area)]);

% [xx1,indx1] = sort(xx);
% zz1 = zz(indx1);
% ind_0 = find(abs(xx1)<1e-5);
% for i = ind_0+1 : length(zz1)
%     if zz1(i) > zz1(i-1)
%         zz1(i) = zz1(i-1);
%     end
% end
% for i = ind_0-1 : -1 : 1
%     if zz1(i) > zz1(i+1)
%         zz1(i) = zz1(i+1);
%     end
% end
% 
% a_L = 0; a_R = 0;
% for i = 1:ind_0-1
%     a_L = a_L + 1/2*abs(xx1(i)*zz1(i+1)-zz1(i)*xx1(i+1));
% end
% for i = ind_0:length(xx1)-1
%     a_R = a_R + 1/2*abs(xx1(i)*zz1(i+1)-zz1(i)*xx1(i+1));
% end
% area_m = a_L+a_R-abs(a_L-a_R);

function rmax=ascan(ring,nt,clorb,theta,rlist)
% rmax = 0.0;
for i = 1:length(rlist)
   rr=rlist(i);
   rin=clorb+[rr*cos(theta);0;rr*sin(theta);0;0;0];
   [dummy,lost]=ringpass(ring,rin,nt); %#ok<ASGLU>
   ind_part(i) = lost;
end
x = find(ind_part==1);
% save('see')
if isempty(x)
    rmax = rlist(end);
elseif x(1)==1
    rmax = 0;
else
    rmax = rlist(x(1)-1);
end
fprintf('theta: %g, r: %g\n',theta,rmax);

function [xmax,zmax]=bscan(ring,nt,clorb,xlist,zlist)
% xmax = 0.0;
% zmax = 0.0;
for i=1:length(xlist)
   rin=clorb+[xlist(i);0;zlist(i);0;0;0];
   [dummy,lost]=ringpass(ring,rin,nt);
   ind_part(i) = lost;
end
x = find(ind_part==1);
if isempty(x) 
    xmax = xlist(end);
    zmax = zlist(end);
elseif x(1) == 1
    xmax = 0; zmax = 0;
else
    xmax=xlist(x(1)-1);
    zmax=zlist(x(1)-1);
end
fprintf('xm: %g, zm: %g\n',xmax,zmax);