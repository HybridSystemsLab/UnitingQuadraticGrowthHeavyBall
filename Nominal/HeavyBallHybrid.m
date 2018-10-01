%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: HeavyBall.m
%--------------------------------------------------------------------------
% Project: Unifying local and global control, using the Heavy-Ball Method
% for convergence to a global minimum. Different parameters for lambda and
% gamma are used globally and locally. 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00

clear all

set(0,'defaultTextInterpreter','latex'); %trying to set the default

% global variables
global gamma_0 lambda_0 gamma_1 lambda_1 c_0 c_10 c_1 delta cTilde_0 cTilde_10 d_0 d_10 lambda gamma r0 r0Prime z1Star
%%%%%%%%%%%%%%%%%%%%%%%
% setting the globals %
%%%%%%%%%%%%%%%%%%%%%%%

% Heavy ball constants: Uniting system
lambda_0 = 11.5;%11.5; % Gravity. 
            % For gamma fixed, "large values of  lambda are seen to give rise to slowly converging 
            % solutions resembling the steepest descent’s while smaller values give 
            % rise to fast solutions with oscillations getting wilder as lambda decreases."
gamma_0 = 2/3;% Viscous friction to mass ratio.

lambda_1 = 1/5;

gamma_1 = 1/2; 

% Constants for the individual heavy ball controllers
lambda = lambda_0;
gamma = gamma_0;

% Convexity constant: Try to determine what these would be.
alpha_0 = 1;
alpha_10 = 1;

% Varepsilons
epsilon0 = 4; 
epsilon10 = 0.7; %0.7;

% Level set values for Lyapunov function V:
c_0 = 12.5; % \U_0
c_10 = 6.39; % \T_{1,0}

% Values for the new switching rules: 
cTilde_0 = alpha_0*epsilon0/2; 
cTilde_10 = alpha_10*epsilon10/2; 

d_0 = c_0 - gamma_0 * ((2*cTilde_0^2)/alpha_0 + CalculateLStar())
d_10 = c_10 - gamma_1 * ((2*cTilde_10^2)/alpha_10 + CalculateLStar())

r0 = gamma_0 * alpha2((2*cTilde_10)/alpha_10) + d_10
r0Prime = r0 + 1; % Can adjust this, as needed.

delta = 0.01;

%%%%%%%%%%%%%%%%%%%%%%
% setting the locals %
%%%%%%%%%%%%%%%%%%%%%%

timeToDelta = 0;
timeToDeltaIdx = 1;
z1delta = 0;
z2delta = 0;

timeToDeltaSlow = 0;
timeToDeltaIdxSlow = 1;
z1deltaSlow = 0;
z2deltaSlow = 0;

timeToDeltaOscillating = 0;
timeToDeltaIdxOscillating = 1;
z1deltaOscillating = 0;
z2deltaOscillating = 0;

% initial conditions
z1_0 = -10;
z2_0 = 0;
q_0 = 1;
% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0];

% simulation horizon
TSPAN=[0 400];
JSPAN = [0 10];

% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;

options = odeset('RelTol',1e-6,'MaxStep',.1);

% simulate the hybrid closed-loop system
[t,j,x] = HyEQsolver(@f,@g,@C,@D,...
    x0,TSPAN,JSPAN,rule,options);

% simulate slow heavy ball
[tSlow,jSlow,xSlow] = HyEQsolver(@f_Ind,@g_Ind,@C_Ind,@D_Ind,...
    x0,TSPAN,JSPAN,rule,options);

% Update parameters to those for fast, oscillatory controller
lambda = lambda_1;
gamma = gamma_1;

% simulate fast, oscillating heavy ball
[tOscillating,jOscillating,xOscillating] = HyEQsolver(@f_Ind,@g_Ind,@C_Ind,@D_Ind,...
    x0,TSPAN,JSPAN,rule,options);


%Add code here to look through x1, and see the last time it dips below
%delta
for i=2:length(x(:,1))
    if (distance(x(i,1)) <= delta) && (distance(x(i-1,1)) > delta)
        timeToDeltaIdx = i;
        z1delta = x(i,1);
    end
end

z2delta = x(timeToDeltaIdx,2);
timeToDelta = t(timeToDeltaIdx,1);

%Add code here to look through x1, and see the last time it dips below
%delta
for i=2:length(xSlow(:,1))
    if (distance(xSlow(i,1)) <= delta) && (distance(xSlow(i-1,1)) > delta)
        timeToDeltaIdxSlow = i;
        z1deltaSlow = xSlow(i,1);
    end
end

z2deltaSlow = xSlow(timeToDeltaIdxSlow,2);
timeToDeltaSlow = tSlow(timeToDeltaIdxSlow,1);

%Add code here to look through x1, and see the last time it dips below
%delta
for i=2:length(xOscillating(:,1))
    if (distance(xOscillating(i,1)) <= delta) && (distance(xOscillating(i-1,1)) > delta)
        timeToDeltaIdxOscillating = i;
        z1deltaOscillating = xOscillating(i,1);
    end
end

z2deltaOscillating = xOscillating(timeToDeltaIdxOscillating,2);
timeToDeltaOscillating = tOscillating(timeToDeltaIdxOscillating,1);

% Prepare data to plot on same figure
minarc = min([length(x),length(xSlow),length(xOscillating)]);
ta = [t(1:minarc),tSlow(1:minarc),tOscillating(1:minarc)];
ja = [j(1:minarc),jSlow(1:minarc),jOscillating(1:minarc)];
xa = [x(1:minarc,1),xSlow(1:minarc,1),xOscillating(1:minarc,1)];
xb = [x(1:minarc,2),xSlow(1:minarc,2),xOscillating(1:minarc,2)];

% plot solution
figure(1) 
clf
modificatorF{1} = '';
modificatorF{2} = 'LineWidth';
modificatorF{3} = 1.5;
modificatorJ{1} = '*--';
modificatorJ{2} = 'LineWidth';
modificatorJ{3} = 1.5;
subplot(2,1,1), plotHarc(ta,ja,xa,[],modificatorF,modificatorJ);
hold on
plot(timeToDelta,z1delta,'k.','MarkerSize',14)
strDelta = [num2str(timeToDelta), 's'];
text(timeToDelta,z1delta,strDelta,'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',12);
plot(timeToDeltaSlow,z1deltaSlow,'k.','MarkerSize', 14)
strDeltaSlow = [num2str(timeToDeltaSlow), 's'];
text(timeToDeltaSlow,z1deltaSlow,strDeltaSlow,'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',12);

plot(timeToDeltaOscillating,z1deltaOscillating,'k.','MarkerSize', 14)
strDeltaOscillating = [num2str(timeToDeltaOscillating), 's'];
text(timeToDeltaOscillating,z1deltaOscillating,strDeltaOscillating,'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',12);
grid on
ylabel('$\mathrm{z_1}$','FontSize',16)
axis([0 250 -10 6]);
subplot(2,1,2), plotHarc(ta,ja,xb,[],modificatorF,modificatorJ);
hold on
plot(timeToDelta,z2delta,'k.','MarkerSize', 14)
plot(timeToDeltaSlow,z2deltaSlow,'k.','MarkerSize', 14)
plot(timeToDeltaOscillating,z2deltaOscillating,'k.','MarkerSize', 14)
grid on
ylabel('$\mathrm{z_2}$','FontSize',16)
axis([0 250 -3 4]);
saveas(gcf,'Plots\MultiSubPlots','epsc')


% Generating contours for c_0 and c_10
x1 = -11:0.005:10;
x2 = -4:0.005:6;
[Z1,Z2] = meshgrid(x1,x2);
V0_0 = gamma_0*(0.25*(Z1-z1Star).^2 - CalculateLStar()) + (1/2)*Z2.^2;
V1_1 = gamma_1*(0.25*(Z1-z1Star).^2 - CalculateLStar()) + (1/2)*Z2.^2;

% Generating contours for cTilde_0 and cTilde_10, and
% d_0 and d_10
for i=1:length(x1)
    cTilde0Vec(i) = cTilde_0;
    cTilde10Vec(i) = cTilde_10;
end

for i=1:length(x2)
    d0Vec(i) = d_0; 
    d10Vec(i) = d_10;
end

% plot phase plane
figure(2) % position
clf
pos1 = [0.1 0.6 0.35 0.35];
subplot('Position', pos1), plotHarcColor(xSlow(:,1),jSlow,xSlow(:,2),tSlow);
xlabel('$\mathrm{z_1}$','FontSize',16)
ylabel('$\mathrm{z_2}$','FontSize',16)
contour(Z1,Z2,V0_0,[c_0 c_0],'-r','ShowText','on')
contour(Z1,Z2,V1_1,[c_10 c_10],'-g','ShowText','on')
text(0,5,'$c_0$','FontSize',12,'VerticalAlignment','bottom')
text(0,3.6,'$c_{1,0}$','FontSize',12,'VerticalAlignment','bottom')
plot(x1,cTilde0Vec, '-r')
text(3,cTilde_0,'$\tilde{c}_0$','FontSize',12,'VerticalAlignment','bottom')
plot(x1,cTilde10Vec, '-g')
text(3,cTilde_10,'$\tilde{c}_{1,0}$','FontSize',12,'VerticalAlignment','bottom')
plot(d0Vec,x2,'-r')
text(d_0,5,'$d_0$','FontSize',12,'HorizontalAlignment','left')
plot(d10Vec,x2,'-g')
text(d_10,5,'$d_{1,0}$','FontSize',12,'HorizontalAlignment','right')
plot(z1deltaSlow,z2deltaSlow,'r.','MarkerSize', 14)
strDeltaSlow = [num2str(timeToDeltaSlow), 's'];
text(z1deltaSlow,z2deltaSlow,strDeltaSlow,'HorizontalAlignment','right','VerticalAlignment','bottom');
axis([-11 10 -0.5 6]);
grid on

pos2 = [0.55 0.6 0.35 0.35];
subplot('Position', pos2), plotHarcColor(xOscillating(:,1),jOscillating,xOscillating(:,2),tOscillating);
xlabel('$\mathrm{z_1}$','FontSize',16)
ylabel('$\mathrm{z_2}$','FontSize',16)
contour(Z1,Z2,V0_0,[c_0 c_0],'-r','ShowText','on')
contour(Z1,Z2,V1_1,[c_10 c_10],'-g','ShowText','on')
text(0,5,'$c_0$','FontSize',12,'VerticalAlignment','bottom')
text(0,3.6,'$c_{1,0}$','FontSize',12,'VerticalAlignment','bottom')
plot(x1,cTilde0Vec, '-r')
text(3,cTilde_0,'$\tilde{c}_0$','FontSize',12,'VerticalAlignment','top')
text(3,cTilde_10,'$\tilde{c}_{1,0}$','FontSize',12,'VerticalAlignment','top')
plot(x1,cTilde10Vec, '-g')
plot(d0Vec,x2,'-r')
plot(d10Vec,x2,'-g')
text(d_0,5,'$d_0$','FontSize',12,'HorizontalAlignment','left')
text(d_10,5,'$d_{1,0}$','FontSize',12,'HorizontalAlignment','right')
plot(z1deltaOscillating,z2deltaOscillating,'r.','MarkerSize', 14)
strDeltaOscillating = [num2str(timeToDeltaOscillating), 's'];
text(z1deltaOscillating,z2deltaOscillating,strDeltaOscillating,'HorizontalAlignment','right','VerticalAlignment','top');
axis([-11 10 -4 6]);
grid on

pos3 = [0.33 0.15 0.35 0.35];
subplot('Position', pos3), plotHarcColor(x(:,1),j,x(:,2),t);
xlabel('$\mathrm{z_1}$','FontSize',16)
ylabel('$\mathrm{z_2}$','FontSize',16)
contour(Z1,Z2,V0_0,[c_0 c_0],'-r','ShowText','on')
contour(Z1,Z2,V1_1,[c_10 c_10],'-g','ShowText','on')
text(0,5,'$c_0$','FontSize',12,'VerticalAlignment','bottom')
text(0,3.6,'$c_{1,0}$','FontSize',12,'VerticalAlignment','bottom')
plot(x1,cTilde0Vec, '-r')
text(3,cTilde_0,'$\tilde{c}_0$','FontSize',12,'VerticalAlignment','bottom')
text(3,cTilde_10,'$\tilde{c}_{1,0}$','FontSize',12,'VerticalAlignment','bottom')
plot(x1,cTilde10Vec, '-g')
plot(d0Vec,x2,'-r')
plot(d10Vec,x2,'-g')
text(d_0,5,'$d_0$','FontSize',12,'HorizontalAlignment','left')
text(d_10,5,'$d_{1,0}$','FontSize',12,'HorizontalAlignment','right')
plot(z1delta,z2delta,'r.','MarkerSize', 14)
strDelta = [num2str(timeToDelta), 's'];
text(z1delta,z2delta,strDelta,'HorizontalAlignment','right','VerticalAlignment','bottom');
axis([-11 10 -0.5 6]);
grid on
saveas(gcf,'Plots\MultiPlane','epsc')