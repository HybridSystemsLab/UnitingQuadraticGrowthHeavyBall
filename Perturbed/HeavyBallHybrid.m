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
global gamma_0 lambda_0 gamma_1 lambda_1 c_0 c_10 delta z1Star cTilde_0 cTilde_10 d_0 d_10 lambda gamma r0 r0Prime randomsInterp randomsIndex sigma
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
deltaNoise = 0.01;



%%%%%%%%%%%%%%%%%%%%%%
% setting the locals %
%%%%%%%%%%%%%%%%%%%%%%

timeToDelta = 0;
timeToDeltaIdx = 1;
z1delta = 0;
z2delta = 0;
jumpsx1 = [];
jumpsx2 = [];
jumpst = [];
jumpIndex = 1;

timeToDeltab = 0;
timeToDeltaIdxb = 1;
z1deltab = 0;
z2deltab = 0;
jumpsx1b = [];
jumpsx2b = [];
jumpstb = [];
jumpIndexb = 1;

timeToDeltac = 0;
timeToDeltaIdxc = 1;
z1deltac = 0;
z2deltac = 0;
jumpsx1c = [];
jumpsx2c = [];
jumpstc = [];
jumpIndexc = 1;

timeToDeltaNom = 0;
timeToDeltaIdxNom = 1;
z1deltaNom = 0;
z2deltaNom = 0;

% initial conditions
z1_0 = -10;
z2_0 = 0;
q_0 = 1;
zHat1_0 = -10; 
zHat2_0 = 0;
% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0];

% simulation horizon
TSPAN=[0 200];
JSPAN = [0 10];

% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;

options = odeset('RelTol',1e-6,'MaxStep',.1);

% simulate the nominal system
[tNom,jNom,xNom] = HyEQsolver(@fNom,@g,@C,@D,...
    x0,TSPAN,JSPAN,rule,options);

nomLen = 1000*length(xNom);
sample = 0.1;

% simulate the perturbed system: make this bigger if increase sim time!
randoms = randn(1,nomLen);
randInd = 1:nomLen;
randomsInterpIndex = 1:sample:nomLen;
randomsInterp = interp1(randInd,randoms,randomsInterpIndex);
randomsIndex = 1;
sigma = 0.1;
[ta,ja,xa] = HyEQsolver(@f,@g,@C,@D,...
    x0,TSPAN,JSPAN,rule,options);

randomsIndex = 1;
sigma = 0.5;
% simulate the perturbed system
[tb,jb,xb] = HyEQsolver(@f,@g,@C,@D,...
    x0,TSPAN,JSPAN,rule,options);

randomsIndex = 1;
sigma = 1;
% simulate the perturbed system
[tc,jc,xc] = HyEQsolver(@f,@g,@C,@D,...
    x0,TSPAN,JSPAN,rule,options);

% Finding time of convergence for the perturbed system with sigma = 0.01
for i=2:length(xa(:,1))
    if (distance(xa(i,1)) <= deltaNoise) && (distance(xa(i-1,1)) > deltaNoise)
        timeToDeltaIdx = i;
        z1delta = xa(i,1);
    end
end

z2delta = xa(timeToDeltaIdx,2);
timeToDelta = ta(timeToDeltaIdx,1);

% Find the jumps for sigma = 0.1:
for i=2:length(ja)
    if(ja(i,1) ~= ja(i-1,1))
        jumpsx1(jumpIndex) = xa(i,1);
        jumpsx2(jumpIndex) = xa(i,2);
        jumpst(jumpIndex) = ta(i,1);
        jumpIndex = jumpIndex + 1;
    end
end

% Finding time of convergence for the perturbed system with sigma = 0.05
for i=2:length(xb(:,1))
    if (distance(xb(i,1)) <= deltaNoise) && (distance(xb(i-1,1)) > deltaNoise)
        timeToDeltaIdxb = i;
        z1deltab = xb(i,1);
    end
end

z2deltab = xb(timeToDeltaIdxb,2);
timeToDeltab = tb(timeToDeltaIdxb,1);

% Find the jumps for sigma = 0.5:
for i=2:length(jb)
    if(jb(i,1) ~= jb(i-1,1))
        jumpsx1b(jumpIndexb) = xb(i,1);
        jumpsx2b(jumpIndexb) = xb(i,2);
        jumpstb(jumpIndexb) = tb(i,1);
        jumpIndexb = jumpIndexb + 1;
    end
end

% Finding time of convergence for the perturbed system with sigma = 0.1
for i=2:length(xc(:,1))
    if (distance(xc(i,1)) <= deltaNoise) && (distance(xc(i-1,1)) > deltaNoise)
        timeToDeltaIdxc = i;
        z1deltac = xc(i,1);
    end
end

z2deltac = xc(timeToDeltaIdxc,2);
timeToDeltac = tc(timeToDeltaIdxc,1);

% Find the jumps for sigma = 0.5:
for i=2:length(jc)
    if(jc(i,1) ~= jc(i-1,1))
        jumpsx1c(jumpIndexc) = xc(i,1);
        jumpsx2c(jumpIndexc) = xc(i,2);
        jumpstc(jumpIndexc) = tc(i,1);
        jumpIndexc = jumpIndexc + 1;
    end
end


% Finding time of convergence for the nominal system
for i=2:length(xNom(:,1))
    if (distance(xNom(i,1)) <= delta) && (distance(xNom(i-1,1)) > delta)
        timeToDeltaIdxNom= i;
        z1deltaNom = xNom(i,1);
    end
end

z2deltaNom = xNom(timeToDeltaIdxNom,2);
timeToDeltaNom = tNom(timeToDeltaIdxNom,1);

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
figure(1) % position
clf
pos1 = [0.1 0.6 0.35 0.35];
subplot('Position', pos1), plotHarc(xa(:,1),ja,xa(:,2));
hold on
xlabel('$\mathrm{z_1}$','FontSize',16)
ylabel('$\mathrm{z_2}$','FontSize',16)
contour(Z1,Z2,V0_0,[c_0 c_0],'-r','ShowText','on')
contour(Z1,Z2,V1_1,[c_10 c_10],'-g','ShowText','on')
text(0,5,'$c_0$','FontSize',14,'VerticalAlignment','bottom','HorizontalAlignment','right')
text(0,3.6,'$c_{1,0}$','FontSize',14,'VerticalAlignment','bottom','HorizontalAlignment','right')
plot(x1,cTilde0Vec, '-r')
plot(x1,cTilde10Vec, '-g')
text(-7,cTilde_0,'$\tilde{c}_0$','FontSize',14,'VerticalAlignment','bottom')
text(-7,cTilde_10,'$\tilde{c}_{1,0}$','FontSize',14,'VerticalAlignment','bottom')
plot(d0Vec,x2,'-r')
plot(d10Vec,x2,'-g')
text(d_0,5,'$d_0$','FontSize',14,'HorizontalAlignment','left')
text(d_10,5,'$d_{1,0}$','FontSize',14,'HorizontalAlignment','right')

plot(z1delta,z2delta,'r.','MarkerSize', 14)
strDelta = [num2str(timeToDelta), 's'];
text(z1delta,z2delta,strDelta,'HorizontalAlignment','left','VerticalAlignment','bottom');
axis([-11 10 -0.5 6]);
grid on

pos2 = [0.55 0.6 0.35 0.35];
subplot('Position', pos2), plotHarc(xb(:,1),jb,xb(:,2));
hold on
xlabel('$\mathrm{z_1}$','FontSize',16)
ylabel('$\mathrm{z_2}$','FontSize',16)
contour(Z1,Z2,V0_0,[c_0 c_0],'-r','ShowText','on')
contour(Z1,Z2,V1_1,[c_10 c_10],'-g','ShowText','on')
text(0,5,'$c_0$','FontSize',14,'VerticalAlignment','bottom','HorizontalAlignment','right')
text(0,3.6,'$c_{1,0}$','FontSize',14,'VerticalAlignment','bottom','HorizontalAlignment','right')
plot(x1,cTilde0Vec, '-r')
plot(x1,cTilde10Vec, '-g')
text(-7,cTilde_0,'$\tilde{c}_0$','FontSize',14,'VerticalAlignment','bottom')
text(-7,cTilde_10,'$\tilde{c}_{1,0}$','FontSize',14,'VerticalAlignment','bottom')
plot(d0Vec,x2,'-r')
plot(d10Vec,x2,'-g')
text(d_0,5,'$d_0$','FontSize',14,'HorizontalAlignment','left')
text(d_10,5,'$d_{1,0}$','FontSize',14,'HorizontalAlignment','right')

plot(z1deltab,z2deltab,'r.','MarkerSize', 14)
strDeltab = [num2str(timeToDeltab), 's'];
text(z1deltab,z2deltab,strDeltab,'HorizontalAlignment','left','VerticalAlignment','bottom');
axis([-11 10 -0.5 6]);
grid on

pos3 = [0.1 0.15 0.35 0.35];
subplot('Position', pos3), plotHarc(xc(:,1),jc,xc(:,2));
hold on
xlabel('$\mathrm{z_1}$','FontSize',16)
ylabel('$\mathrm{z_2}$','FontSize',16)
contour(Z1,Z2,V0_0,[c_0 c_0],'-r','ShowText','on')
contour(Z1,Z2,V1_1,[c_10 c_10],'-g','ShowText','on')
text(0,5,'$c_0$','FontSize',14,'VerticalAlignment','bottom','HorizontalAlignment','right')
text(0,3.6,'$c_{1,0}$','FontSize',14,'VerticalAlignment','bottom','HorizontalAlignment','right')
plot(x1,cTilde0Vec, '-r')
plot(x1,cTilde10Vec, '-g')
text(-7,cTilde_0,'$\tilde{c}_0$','FontSize',14,'VerticalAlignment','bottom')
text(-7,cTilde_10,'$\tilde{c}_{1,0}$','FontSize',14,'VerticalAlignment','bottom')
plot(d0Vec,x2,'-r')
plot(d10Vec,x2,'-g')
text(d_0,5,'$d_0$','FontSize',14,'HorizontalAlignment','left')
text(d_10,5,'$d_{1,0}$','FontSize',14,'HorizontalAlignment','right')

plot(z1deltac,z2deltac,'r.','MarkerSize', 14)
strDeltac = [num2str(timeToDeltac), 's'];
text(z1deltac,z2deltac,strDeltac,'HorizontalAlignment','left','VerticalAlignment','bottom');
axis([-11 10 -0.5 6]);
grid on

pos4 = [0.55 0.15 0.35 0.35];
subplot('Position', pos4), plotHarc(xNom(:,1),jNom,xNom(:,2));
hold on
xlabel('$\mathrm{z_1}$','FontSize',16)
ylabel('$\mathrm{z_2}$','FontSize',16)
contour(Z1,Z2,V0_0,[c_0 c_0],'-r','ShowText','on')
contour(Z1,Z2,V1_1,[c_10 c_10],'-g','ShowText','on')
text(0,5,'$c_0$','FontSize',14,'VerticalAlignment','bottom','HorizontalAlignment','right')
text(0,3.6,'$c_{1,0}$','FontSize',14,'VerticalAlignment','bottom','HorizontalAlignment','right')
plot(x1,cTilde0Vec, '-r')
plot(x1,cTilde10Vec, '-g')
text(-7,cTilde_0,'$\tilde{c}_0$','FontSize',14,'VerticalAlignment','bottom')
text(-7,cTilde_10,'$\tilde{c}_{1,0}$','FontSize',14,'VerticalAlignment','bottom')
plot(d0Vec,x2,'-r')
plot(d10Vec,x2,'-g')
text(d_0,5,'$d_0$','FontSize',14,'HorizontalAlignment','left')
text(d_10,5,'$d_{1,0}$','FontSize',14,'HorizontalAlignment','right')

plot(z1deltaNom,z2deltaNom,'r.','MarkerSize', 14)
strDeltaNom = [num2str(timeToDeltaNom),'s'];
text(z1deltaNom,z2deltaNom,strDeltaNom,'HorizontalAlignment','left','VerticalAlignment','bottom');
axis([-11 10 -0.5 6]);
grid on
saveas(gcf,'Plots\PhasePlaneMulti','epsc')

% Plot all the examples on the same plot
figure(2)
clf
pos1 = [0.1 0.15 0.45 0.75];
subplot('Position', pos1),plotHarc(xNom(:,1),jNom,xNom(:,2));
hold on
plot(xa(:,1),xa(:,2),'Color',[0 .5 .0]);
plot(jumpsx1(1),jumpsx2(1),'*','Color',[0 .5 .0]);

plot(xb(:,1),xb(:,2),'-r');
plot(jumpsx1b(1),jumpsx2b(1),'*r');

plot(xc(:,1),xc(:,2),'Color',[.44 0 .44]);
plot(jumpsx1c(1),jumpsx2c(1),'*','Color',[.44 0 .44]);
xlabel('$\mathrm{z_1}$','FontSize',16)
ylabel('$\mathrm{z_2}$','FontSize',16)
contour(Z1,Z2,V0_0,[c_0 c_0],'-r','ShowText','on')
contour(Z1,Z2,V1_1,[c_10 c_10],'-g','ShowText','on')
plot(z1deltaNom,z2deltaNom,'r.','MarkerSize', 14)
strDeltaNom = [num2str(timeToDeltaNom),'s'];
text(z1deltaNom,z2deltaNom,strDeltaNom,'HorizontalAlignment','left','VerticalAlignment','bottom');
plot(z1delta,z2delta,'r.','MarkerSize', 14)
strDelta = [num2str(timeToDelta), 's'];
text(z1delta,z2delta,strDelta,'HorizontalAlignment','left','VerticalAlignment','top');
plot(z1deltab,z2deltab,'r.','MarkerSize', 14)
strDeltab = [num2str(timeToDeltab), 's'];
text(z1deltab,z2deltab,strDeltab,'HorizontalAlignment','right','VerticalAlignment','bottom');
plot(z1deltac,z2deltac,'r.','MarkerSize', 14)
strDeltac = [num2str(timeToDeltac), 's'];
text(z1deltac,z2deltac,strDeltac,'HorizontalAlignment','right','VerticalAlignment','top');
axis([-11 2 -1 5.5]);
grid on
pos2 = [0.65 0.15 0.25 0.75];
subplot('Position', pos2),plotHarc(xNom(:,1),jNom,xNom(:,2));
hold on
plot(xa(:,1),xa(:,2),'Color',[0 .5 .0]);
plot(jumpsx1(1),jumpsx2(1),'*','Color',[0 .5 .0]);

plot(xb(:,1),xb(:,2),'-r');
plot(jumpsx1b(1),jumpsx2b(1),'*r');

plot(xc(:,1),xc(:,2),'Color',[.44 0 .44]);
plot(jumpsx1c(1),jumpsx2c(1),'*','Color',[.44 0 .44]);
xlabel('$\mathrm{z_1}$','FontSize',16)
ylabel('$\mathrm{z_2}$','FontSize',16)
contour(Z1,Z2,V0_0,[c_0 c_0],'-r','ShowText','on')
contour(Z1,Z2,V1_1,[c_10 c_10],'-g','ShowText','on')
plot(z1deltaNom,z2deltaNom,'r.','MarkerSize', 14)
strDeltaNom = [num2str(timeToDeltaNom),'s'];
text(z1deltaNom,z2deltaNom,strDeltaNom,'HorizontalAlignment','left','VerticalAlignment','bottom');
plot(z1delta,z2delta,'r.','MarkerSize', 14)
strDelta = [num2str(timeToDelta), 's'];
text(z1delta,z2delta,strDelta,'HorizontalAlignment','left','VerticalAlignment','top');
plot(z1deltab,z2deltab,'r.','MarkerSize', 14)
strDeltab = [num2str(timeToDeltab), 's'];
text(z1deltab,z2deltab,strDeltab,'HorizontalAlignment','right','VerticalAlignment','bottom');
plot(z1deltac,z2deltac,'r.','MarkerSize', 14)
strDeltac = [num2str(timeToDeltac), 's'];
text(z1deltac,z2deltac,strDeltac,'HorizontalAlignment','right','VerticalAlignment','top');
axis([-0.5 0.25 -0.5 4]);
grid on
saveas(gcf,'Plots\PhasePlaneNomAllCloseup','epsc')