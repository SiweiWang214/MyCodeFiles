function [emit,radiationDamping] = calcdampingtime(varargin)

% if isempty(varargin)
%     global THERING GLOBVAL;
% else
    THERING = varargin{1};
    GLOBVAL.E0 = THERING{1}.Energy;
    GLOBVAL.LatticeFile = '';  
% end

sum.e0            = GLOBVAL.E0*1e-9;
sum.circumference = findspos(THERING, length(THERING)+1);
sum.revTime       = sum.circumference / PhysConstant.speed_of_light_in_vacuum.value;
sum.revFreq       = PhysConstant.speed_of_light_in_vacuum.value / sum.circumference;
sum.gamma         = sum.e0 / PhysConstant.electron_mass_energy_equivalent_in_MeV.value*1e3; %0.51099906e-3; 
sum.beta          = sqrt(1 - 1/sum.gamma);

[TD, sum.tunes, sum.chromaticity] = twissring(THERING, 0, 1:length(THERING)+1, 'chrom', 1e-8);
sum.compactionFactor = mcf(THERING);

% For calculating the synchrotron integrals
temp  = cat(2,TD.Dispersion);
D_x   = temp(1,:)';
D_x_  = temp(2,:)';
beta  = cat(1, TD.beta);
alpha = cat(1, TD.alpha);
gamma = (1 + alpha.^2) ./ beta;
circ  = TD(length(THERING)+1).SPos;

% Synchrotron integral calculation
sum.integrals = zeros(1,6);

ii =0;
for i = 1:length(THERING)
    if isfield(THERING{i}, 'BendingAngle') && isfield(THERING{i}, 'EntranceAngle')
        ii = ii +1;
        rho = THERING{i}.Length/THERING{i}.BendingAngle;
        [dI1,dI2,dI3,dI4,dI5,curHavg1(ii), Dxavg(ii)] = ...
            calcRadInt(rho,THERING{i}.BendingAngle, ...
            alpha(i,1),beta(i,1),D_x(i),D_x_(i),...
            THERING{i}.PolynomB(2),THERING{i}.EntranceAngle,THERING{i}.ExitAngle);
        
        sum.integrals(1) = sum.integrals(1) + dI1;
        sum.integrals(2) = sum.integrals(2) + dI2;
        sum.integrals(3) = sum.integrals(3) + dI3;
        % For general wedge magnets
        sum.integrals(4) = sum.integrals(4) + dI4;
        sum.integrals(5) = sum.integrals(5) + dI5;
        %         sum.integrals(4) = sum.integrals(4) + 2*0.5*(D_x(i)+D_x(i+1))*THERING{i}.Length/rho^3;
        H1 = beta(i,1)*D_x_(i)*D_x_(i)+2*alpha(i)*D_x(i)*D_x_(i)+gamma(i)*D_x(i)*D_x(i);
        H0 = beta(i+1,1)*D_x_(i+1)*D_x_(i+1)+2*alpha(i+1)*D_x(i+1)*D_x_(i+1)+gamma(i+1)*D_x(i+1)*D_x(i+1);
        sum.integrals(6) = sum.integrals(6) + THERING{i}.PolynomB(2)^2*Dxavg(ii)^2*THERING{i}.Length;
    end
end

% Damping numbers
% Use Robinson's Theorem
sum.damping(1) = 1 - sum.integrals(4)/sum.integrals(2);
sum.damping(2) = 1;
sum.damping(3) = 2 + sum.integrals(4)/sum.integrals(2);

sum.radiation           = 8.846e-5*sum.e0.^4*sum.integrals(2)/(2*pi);
sum.naturalEnergySpread = sqrt(3.8319e-13*sum.gamma.^2*sum.integrals(3)/(2*sum.integrals(2) + sum.integrals(4)));
emit    = 3.8319e-13*(sum.e0*1e3/PhysConstant.electron_mass_energy_equivalent_in_MeV.value).^2*sum.integrals(5)/(sum.damping(1)*sum.integrals(2)); %% need to be replaced by constant ?

% Damping times
radiationDamping(1) = 1/(2113.1*sum.e0.^3*sum.integrals(2)*sum.damping(1)/circ);
radiationDamping(2) = 1/(2113.1*sum.e0.^3*sum.integrals(2)*sum.damping(2)/circ);
radiationDamping(3) = 1/(2113.1*sum.e0.^3*sum.integrals(2)*sum.damping(3)/circ);

% Slip factor
sum.etac = sum.gamma^(-2) - sum.compactionFactor;

function [dI1,dI2,dI3,dI4,dI5,curHavg, Dxavg] = calcRadInt(rho,theta, a0,b0,D0,D0p,K1,th1,th2)
%[dI1,dI2,dI3,dI4,dI5,curHavg] = calcRadInt(rho,theta, a0,b0,D0,D0p,K1)
%calculate the contribution to the radiation integrals of a dipole.
%  INPUTS
%  rho, theta, radius and angle of the dipole
%  a0, b0, are horizontal alpha and beta at the entrance of the dipole,
%  D0, D0p are dispersion at the entrance of the dipole
%  K1, the gradient parameter in AT's convention, i.e., positive for
%  horizontal focusing, K1=0 by default
%  th1, th2, the entrance and exit angle, respectively, th1=th2= 0 [theta/2] by
%  default.
%

% If not combine dipole
if nargin == 6
    K1=0;
end

%
if nargin<9
    th1 = 0; %theta/2.0;
    th2 = 0; %theta/2.0;
end

% Edge focusing
M21 = tan(th1)/rho;
D0p = M21*D0+D0p;
a0 = -M21*b0+a0;

% split the dipole in N pieces
N = 100;
th = (0:N)/N*theta;

% Compute Twiss parameters inside dipole
for ii=1:length(th)
    [Dx(ii), Dxp(ii)] = calcdisp(rho, th(ii), D0, D0p, K1);
    [ax, bx] = calctwiss(rho, th(ii), a0, b0, K1);
    curHavg1(ii) = (Dx(ii)^2+(ax*Dx(ii)+bx*Dxp(ii))^2)/bx;
end

% Edge focusing
M21 = tan(th2)/rho;
Dxp(end) =  M21*Dx(end)+Dxp(end);
ax  = -M21*bx+ax;
curHavg1(end) = (Dx(end)^2+(ax*Dx(end)+bx*Dxp(end))^2)/bx;

% Average data
curHavg = ((curHavg1(1)+curHavg1(end))/2.0 + sum(curHavg1(2:end-1)))/(length(th)-1);
Dxavg   = ((Dx(1)+Dx(end))/2.0 + sum(Dx(2:end-1)))/(length(th)-1);

dI1 = ((Dx(1) + Dx(end))/2.0 + sum(Dx(2:end-1)))*theta/N;
dI2 = abs(theta/rho);
dI3 = abs(theta/rho^2);
dI4 = (1/rho^2 + 2*K1)*dI1  - (Dx(1)/rho^2*tan(th1) + Dx(end)/rho^2*tan(th2));
dI5 = curHavg*abs(theta/rho^2);

function [Dx, Dxp] = calcdisp(rho, theta, D0, D0p, K1)
%calcdisp - calculate dispersion function inside a combined-function dipole
%  INPUTS
%  1. rho - curvature radius
%  2. theta - angle
%  3. D0 - Horizontal dispersion function at the entrance
%  4. D0p - DErivative of  Horizontal dispersion function at the entrance
%  5. K1 - Focusing
%
% Transfert matrix of A wedge dipole p58 Handbook af Accelerator Physics
s = rho*theta;
if K1>-1/rho^2 %horizontal focusing
    sqK = sqrt(1/rho^2+K1);
    Dx  =  D0*cos(sqK*s) + D0p/sqK*sin(sqK*s)+(1-cos(sqK*s))/rho/sqK^2;
    Dxp = -D0*sqK*sin(sqK*s)+D0p*cos(sqK*s)+sin(sqK*s)/rho/sqK;
else %horizontal defocusing
    sqK = sqrt(-(1/rho^2+K1));
    Dx =  D0*cosh(sqK*s) + D0p/sqK*sinh(sqK*s)+(-1+cosh(sqK*s))/rho/sqK^2;
    Dxp = D0*sqK*sinh(sqK*s)+D0p*cosh(sqK*s)+sinh(sqK*s)/rho/sqK;
end

function [ax, bx] = calctwiss(rho, theta, a0, b0, K1)
% calctwiss calculate twiss function inside a combined-function dipole manget
%  INPUTS
%  1. rho - curvature radius
%  2. theta - angle
%  3. a0 - Horizontal alpha function at the entrance
%  4. b0 - Horizontal beta function at the entrance
%  5. K1 - Focusing
%
%  [beta ] = [  Mx11^2        -2*MX11*Mx12         Mx12^2   ] [beta0 ]
%  [alpha] = [ -Mx11*Mx21 Mx11*Mx22 + Mx11*Mx21   -Mx12*Mx22] [alpha0]
%  [gamma] = [  Mx21^2        -2*MX21*Mx22         Mx22^2   ] [gamma0]

Mx = calcMx(rho, K1,theta);
g0 = (1+a0^2)/b0;
twx2 = [Mx(1,1)^2, -2*Mx(1,1)*Mx(1,2), Mx(1,2)^2;
    -Mx(1,1)*Mx(2,1), Mx(1,1)*Mx(2,2)+Mx(1,2)*Mx(2,1),-Mx(1,2)*Mx(2,2);
    Mx(2,1)^2, -2*Mx(2,1)*Mx(2,2),Mx(2,2)^2] * [b0, a0, g0]';
ax = twx2(2);
bx = twx2(1);

function Mx = calcMx(rho,K1,theta)
% calcMx calculate transfer matrice of a combined-function dipole manget

s = rho*theta;
if K1>-1/rho^2 %horizontal focusing
    sqK = sqrt(1/rho^2+K1);
    Mx = [cos(sqK*s), sin(sqK*s)/sqK; -sqK*sin(sqK*s), cos(sqK*s)];
else %horizontal defocusing
    sqK = sqrt(-(1/rho^2+K1));
    Mx = [cosh(sqK*s), sinh(sqK*s)/sqK; sqK*sinh(sqK*s), cosh(sqK*s)];
end