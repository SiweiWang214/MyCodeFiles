function half_booster_v2
global THERING GLOBVAL
GLOBVAL.E0 = 2.2e9;
GLOBVAL.LatticeFile = 'half_booster_v2';

df = atdrift('df', 1.2983);
% dfs1 = atdrift('dfs1',0.17);
% dfs2 = atdrift('dfs2',0.6583);
% dfs3 = atdrift('dfs3',0.17);
df0 = atdrift('df0', 1.230000);
df1 = atdrift('df1', 0.520000);
% df1s1 = atdrift('df1s1',0.18);
% df1s2 = atdrift('df1s2',0.19);
df2 = atdrift('df2', 0.400000);
% df2s1 = atdrift('df2s1',0.12);
% df2s2 = atdrift('df2s2',0.13);

df3 = atdrift('df3', 1.3650);
% df3s1 = atdrift('df3s1',1.0550);
% df3s2 = atdrift('df3s2',0.1600);

df4 = atdrift('df4', 2.520000);
df5 = atdrift('df5', 2.520000);
df6 = atdrift('df6', 2.520000);
df7 = atdrift('df7', 2.520000);

qh01 = atquadrupole('qh01',0.30000, 1.1115,'StrMPoleSymplectic4Pass');
qh01.('NumIntSteps')=10; 
qh02 = atquadrupole('qh02',0.30000, 2.1193,'StrMPoleSymplectic4Pass');
qh02.('NumIntSteps')=10;

qv01 = atquadrupole('qv01',0.30000, -1.030000,'StrMPoleSymplectic4Pass');
qv01.('NumIntSteps')=10; 
qv02 = atquadrupole('qv02',0.300000, -0.0636,'StrMPoleSymplectic4Pass');
qv02.('NumIntSteps')=10; 
qv03 = atquadrupole('qv03',0.300000, -0.0636,'StrMPoleSymplectic4Pass');
qv03.('NumIntSteps')=10; 

qd01 = atquadrupole('qd01',0.30000, -1.030000,'StrMPoleSymplectic4Pass');
qd01.('NumIntSteps')=10; 
qd02 = atquadrupole('qd02',0.300000, -0.0636,'StrMPoleSymplectic4Pass');
qd02.('NumIntSteps')=10; 

bm05 = atsbend('bm05',0,0,0,'BndMPoleSymplectic4E2Pass');
bm05.('Length') = 2.5979;
bm05.('NumIntSteps')=10; 
bm05.('BendingAngle') = pi/32;
bm05.('K')=+0; 
bm05.('PolynomB')(2)=0; 
bm05.('EntranceAngle')=pi/64; 
bm05.('ExitAngle')=pi/64;
bm05.('MaxOrder') = 2;

bm10 = atsbend('bm10',0,0,0,'BndMPoleSymplectic4E2Pass');
bm10.('Length') = 2.5979;
bm10.('NumIntSteps')=10; 
bm10.('BendingAngle') = pi/16;
bm10.('K')=-0.229000; 
bm10.('PolynomB')(2)=-0.229000; 
bm10.('EntranceAngle')=pi/32; 
bm10.('ExitAngle')=pi/32;
bm10.('MaxOrder') = 2;

sb05 = atsextupole('sb05',0,0,'StrMPoleSymplectic4Pass');
sb05.('NumIntSteps')=10;
sb05.('Length')=0.15; 
sb05.('PolynomB')(3)=-0; 
sb05.('PolynomB')=sb05.('PolynomB')*1/2;
sb05.('MaxOrder') = 2;

sb10 = atsextupole('sb10',0,0,'StrMPoleSymplectic4Pass');
sb10.('NumIntSteps')=10;
sb10.('Length')=0.15;
sb10.('PolynomB')(3)=-0; 
sb10.('PolynomB')=sb10.('PolynomB')*1/2;
sb10.('MaxOrder') = 2;

sb03 = atsextupole('sb03',0,0,'StrMPoleSymplectic4Pass');
sb03.('NumIntSteps')=10;
sb03.('Length')=0.15;
sb03.('PolynomB')(3)=-0; 
sb03.('PolynomB')=sb03.('PolynomB')*1/2;
sb03.('MaxOrder') = 2;

sq01 = atsextupole('sq01',0,0,'StrMPoleSymplectic4Pass');
sq01.('NumIntSteps')=10;
sq01.('Length')=0.15;
sq01.('PolynomB')(3)=0; 
sq01.('PolynomB')=sq01.('PolynomB')*1/2;
sq01.('MaxOrder') = 2;

sq02 = atsextupole('sq02',0,0,'StrMPoleSymplectic4Pass');
sq02.('NumIntSteps')=10;
sq02.('Length')=0.15;
sq02.('PolynomB')(3)=0;
sq02.('PolynomB')=sq02.('PolynomB')*1/2;
sq02.('MaxOrder') = 2;

sq03 = atsextupole('sq03',0,0,'StrMPoleSymplectic4Pass');
sq03.('NumIntSteps')=10;
sq03.('Length')=0.15;
sq03.('PolynomB')(3)=0;
sq03.('PolynomB')=sq03.('PolynomB')*1/2;
sq03.('MaxOrder') = 2;

in = [{df0}, {qh01}, {df1}, {qv02}, {df2}, {bm05}, {df4},{qv01},{df3},{qd01},{df5},{qd02},{df3},{qv03},{df7},{bm05},{df6}, {qh02}];
out = [{df6},{bm05},{df7},{qv03}, {df3}, {qd02}, {df5},{qd01}, {df3},{qv01},{df4}, {bm05}, {df2}, {qv02}, {df1} {qh01}, {df0}];
fodo = [{df}, {bm10}, {df}, {qh02}];
% in = [{df0}, {qh01}, {df1s1}, {sq01}, {df1s2}, {qv02}, {df2s1}, {sb05},{df2s2}, {bm05}, {df3s1},{sb03},{df3s2}, {qh02}];
% out = [{df3s2},{sb03},{df3s1}, {bm05}, {df2s2}, {sb05},{df2s1}, {qv02}, {df1s2}, {sq01}, {df1s1}, {qh01}, {df0}];
% fodo = [{dfs1}, {sq02}, {dfs2},{sq03},{dfs3}, {bm10},{dfs3},{sq03}, {dfs2},{sb10},{dfs1}, {qh02}];

lat = [in, repmat(fodo,1,6), out];
THERING = repmat(lat,1,4);
THERING = setcellstruct(THERING, 'Energy', 1:length(THERING), GLOBVAL.E0);
evalin('caller','global THERING GLOBVAL')
disp('Finished loading booster')