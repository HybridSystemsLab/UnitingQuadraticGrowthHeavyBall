% figure(1) % position
% clf
% pos1 = [0.1 0.6 0.35 0.35];
% subplot('Position', pos1), plotHarc(xa(:,1),ja,xa(:,2));
% hold on
% xlabel('$\mathrm{z_1}$','FontSize',16)
% ylabel('$\mathrm{z_2}$','FontSize',16)
% contour(Z1,Z2,V0_0,[c_0 c_0],'-r','ShowText','on')
% contour(Z1,Z2,V1_1,[c_10 c_10],'-g','ShowText','on')
% text(0,5,'$c_0$','FontSize',14,'VerticalAlignment','bottom','HorizontalAlignment','right')
% text(0,3.6,'$c_{1,0}$','FontSize',14,'VerticalAlignment','bottom','HorizontalAlignment','right')
% plot(x1,cTilde0Vec, '-r')
% plot(x1,cTilde10Vec, '-g')
% text(-7,cTilde_0,'$\tilde{c}_0$','FontSize',14,'VerticalAlignment','bottom')
% text(-7,cTilde_10,'$\tilde{c}_{1,0}$','FontSize',14,'VerticalAlignment','bottom')
% plot(d0Vec,x2,'-r')
% plot(d10Vec,x2,'-g')
% text(d_0,5,'$d_0$','FontSize',14,'HorizontalAlignment','left')
% text(d_10,5,'$d_{1,0}$','FontSize',14,'HorizontalAlignment','right')
% 
% plot(z1delta,z2delta,'r.','MarkerSize', 14)
% strDelta = [num2str(timeToDelta), 's'];
% text(z1delta,z2delta,strDelta,'HorizontalAlignment','left','VerticalAlignment','bottom');
% axis([-11 10 -0.5 6]);
% grid on
% 
% pos2 = [0.55 0.6 0.35 0.35];
% subplot('Position', pos2), plotHarc(xb(:,1),jb,xb(:,2));
% hold on
% xlabel('$\mathrm{z_1}$','FontSize',16)
% ylabel('$\mathrm{z_2}$','FontSize',16)
% contour(Z1,Z2,V0_0,[c_0 c_0],'-r','ShowText','on')
% contour(Z1,Z2,V1_1,[c_10 c_10],'-g','ShowText','on')
% text(0,5,'$c_0$','FontSize',14,'VerticalAlignment','bottom','HorizontalAlignment','right')
% text(0,3.6,'$c_{1,0}$','FontSize',14,'VerticalAlignment','bottom','HorizontalAlignment','right')
% plot(x1,cTilde0Vec, '-r')
% plot(x1,cTilde10Vec, '-g')
% text(-7,cTilde_0,'$\tilde{c}_0$','FontSize',14,'VerticalAlignment','bottom')
% text(-7,cTilde_10,'$\tilde{c}_{1,0}$','FontSize',14,'VerticalAlignment','bottom')
% plot(d0Vec,x2,'-r')
% plot(d10Vec,x2,'-g')
% text(d_0,5,'$d_0$','FontSize',14,'HorizontalAlignment','left')
% text(d_10,5,'$d_{1,0}$','FontSize',14,'HorizontalAlignment','right')
% 
% plot(z1deltab,z2deltab,'r.','MarkerSize', 14)
% strDeltab = [num2str(timeToDeltab), 's'];
% text(z1deltab,z2deltab,strDeltab,'HorizontalAlignment','left','VerticalAlignment','bottom');
% axis([-11 10 -0.5 6]);
% grid on
% 
% pos3 = [0.1 0.15 0.35 0.35];
% subplot('Position', pos3), plotHarc(xc(:,1),jc,xc(:,2));
% hold on
% xlabel('$\mathrm{z_1}$','FontSize',16)
% ylabel('$\mathrm{z_2}$','FontSize',16)
% contour(Z1,Z2,V0_0,[c_0 c_0],'-r','ShowText','on')
% contour(Z1,Z2,V1_1,[c_10 c_10],'-g','ShowText','on')
% text(0,5,'$c_0$','FontSize',14,'VerticalAlignment','bottom','HorizontalAlignment','right')
% text(0,3.6,'$c_{1,0}$','FontSize',14,'VerticalAlignment','bottom','HorizontalAlignment','right')
% plot(x1,cTilde0Vec, '-r')
% plot(x1,cTilde10Vec, '-g')
% text(-7,cTilde_0,'$\tilde{c}_0$','FontSize',14,'VerticalAlignment','bottom')
% text(-7,cTilde_10,'$\tilde{c}_{1,0}$','FontSize',14,'VerticalAlignment','bottom')
% plot(d0Vec,x2,'-r')
% plot(d10Vec,x2,'-g')
% text(d_0,5,'$d_0$','FontSize',14,'HorizontalAlignment','left')
% text(d_10,5,'$d_{1,0}$','FontSize',14,'HorizontalAlignment','right')
% 
% plot(z1deltac,z2deltac,'r.','MarkerSize', 14)
% strDeltac = [num2str(timeToDeltac), 's'];
% text(z1deltac,z2deltac,strDeltac,'HorizontalAlignment','left','VerticalAlignment','bottom');
% axis([-11 10 -0.5 6]);
% grid on
% 
% pos4 = [0.55 0.15 0.35 0.35];
% subplot('Position', pos4), plotHarc(xNom(:,1),jNom,xNom(:,2));
% hold on
% xlabel('$\mathrm{z_1}$','FontSize',16)
% ylabel('$\mathrm{z_2}$','FontSize',16)
% contour(Z1,Z2,V0_0,[c_0 c_0],'-r','ShowText','on')
% contour(Z1,Z2,V1_1,[c_10 c_10],'-g','ShowText','on')
% text(0,5,'$c_0$','FontSize',14,'VerticalAlignment','bottom','HorizontalAlignment','right')
% text(0,3.6,'$c_{1,0}$','FontSize',14,'VerticalAlignment','bottom','HorizontalAlignment','right')
% plot(x1,cTilde0Vec, '-r')
% plot(x1,cTilde10Vec, '-g')
% text(-7,cTilde_0,'$\tilde{c}_0$','FontSize',14,'VerticalAlignment','bottom')
% text(-7,cTilde_10,'$\tilde{c}_{1,0}$','FontSize',14,'VerticalAlignment','bottom')
% plot(d0Vec,x2,'-r')
% plot(d10Vec,x2,'-g')
% text(d_0,5,'$d_0$','FontSize',14,'HorizontalAlignment','left')
% text(d_10,5,'$d_{1,0}$','FontSize',14,'HorizontalAlignment','right')
% 
% plot(z1deltaNom,z2deltaNom,'r.','MarkerSize', 14)
% strDeltaNom = [num2str(timeToDeltaNom),'s'];
% text(z1deltaNom,z2deltaNom,strDeltaNom,'HorizontalAlignment','left','VerticalAlignment','bottom');
% axis([-11 10 -0.5 6]);
% grid on
% saveas(gcf,'Plots\PhasePlaneMulti','epsc')

% Plot all the examples on the same plot
figure(2)
clf
modificatorF{1} = '';
modificatorF{2} = 'LineWidth';
modificatorF{3} = 1.5;
modificatorJ{1} = '*--';
modificatorJ{2} = 'LineWidth';
modificatorJ{3} = 1.5;

pos1 = [0.1 0.15 0.45 0.75];
subplot('Position', pos1),plotHarc(xNom(:,1),jNom,xNom(:,2),[],modificatorF,modificatorJ);
hold on
plot(xa(:,1),xa(:,2),'Color',[0 .5 .0],'Linewidth',1.5);
plot(jumpsx1(1),jumpsx2(1),'*','Color',[0 .5 .0]);

plot(xb(:,1),xb(:,2),'-r','Linewidth',1.5);
plot(jumpsx1b(1),jumpsx2b(1),'*r');

plot(xc(:,1),xc(:,2),'Color',[.44 0 .44],'Linewidth',1.5);
plot(jumpsx1c(1),jumpsx2c(1),'*','Color',[.44 0 .44]);
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
subplot('Position', pos2),plotHarc(xNom(:,1),jNom,xNom(:,2),[],modificatorF,modificatorJ);
hold on
plot(xa(:,1),xa(:,2),'Color',[0 .5 .0],'Linewidth',1.5);
plot(jumpsx1(1),jumpsx2(1),'*','Color',[0 .5 .0]);

plot(xb(:,1),xb(:,2),'-r','Linewidth',1.5);
plot(jumpsx1b(1),jumpsx2b(1),'*r');

plot(xc(:,1),xc(:,2),'Color',[.44 0 .44],'Linewidth',1.5);
plot(jumpsx1c(1),jumpsx2c(1),'*','Color',[.44 0 .44]);
xlabel('$\mathrm{z_1}$','FontSize',16)
ylabel('$\mathrm{z_2}$','FontSize',16)
contour(Z1,Z2,V0_0,[c_0 c_0],'-r','ShowText','on')
contour(Z1,Z2,V1_1,[c_10 c_10],'-g','ShowText','on')
text(0,5,'$c_0$','FontSize',14,'VerticalAlignment','bottom','HorizontalAlignment','right')
text(0,3.6,'$c_{1,0}$','FontSize',14,'VerticalAlignment','bottom','HorizontalAlignment','right')
plot(x1,cTilde0Vec, '-r')
plot(x1,cTilde10Vec, '-g')
text(-0.6,cTilde_0,'$\tilde{c}_0$','FontSize',14,'VerticalAlignment','bottom')
text(-0.6,cTilde_10,'$\tilde{c}_{1,0}$','FontSize',14,'VerticalAlignment','bottom')
plot(d0Vec,x2,'-r')
plot(d10Vec,x2,'-g')
text(d_0,5,'$d_0$','FontSize',14,'HorizontalAlignment','left')
text(d_10,5,'$d_{1,0}$','FontSize',14,'HorizontalAlignment','right')

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
axis([-0.75 0.25 -0.5 4]);
grid on
saveas(gcf,'Plots\PhasePlaneNomAllCloseup','epsc')