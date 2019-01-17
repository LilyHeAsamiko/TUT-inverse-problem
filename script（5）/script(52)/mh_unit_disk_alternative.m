n = 3; 
sigma = 0.2;
iteraatiot = 10000;
burn_in = 500;
gamma = 0.1;

particle = [0.5 ; 0.5];
particle = particle(1)*[cos(particle(2)) ; sin(particle(2))];

sensors_vali = 2*pi/n;
sensors_phi = [0 : sensors_vali : 2*pi-sensors_vali];
sensors = [ ones(1, n) ; sensors_phi];
sensors = [cos(sensors(2,:)) ; sin(sensors(2,:))];

etaisyydet = sqrt(sum((sensors - particle(:,ones(1,n))).^2));
v = 1./etaisyydet;
%v = v + sigma*randn(size(v));

alkupiste = [0 ; 0];

x = alkupiste;
Pisteet = zeros(2, iteraatiot);
Pisteet(:,1) = [x];
acceptance_ind = 0;

 etaisyydet = sqrt(sum((sensors - x(:, ones(1,n))).^2));
 pi_x = exp( -sum((v - 1./etaisyydet).^2)/(2*sigma^2));

for k = 1 : iteraatiot

  Fy = x/(1 - x'*x) + gamma*randn(size(x));
 y = (-1 + sqrt(1 + 4*Fy'*Fy))*Fy/(2*Fy'*Fy);

 if norm(y) <= 1
  etaisyydet  = sqrt(sum((sensors - y(:,ones(1, n))).^2)); 
  pi_y = exp( -sum((v - 1./etaisyydet).^2)/(2*sigma^2));
 else
  pi_y = 0;
 end
 
 acceptance_ratio = min([ 1 ; (pi_y*(1 + x'*x)*(1-y'*y)^3)/...
 (pi_x*(1+y'*y)*(1-x'*x)^3)]);
 if (acceptance_ratio >= rand(1))
  x = y;
  pi_x = pi_y;
  acceptance_ind = acceptance_ind + 1;
 end
 
 Pisteet(:,k) = x;

 end
 
 Pisteet = Pisteet(:, burn_in + 1 : end);
 acceptance_pros = acceptance_ind/iteraatiot;

res = 200;
r = [0 : 1/(res-1) : 1];
phi = [0 : 2*pi/(res-1) : 2*pi];
[R, Phi] = meshgrid(r, phi);

posterior_values = zeros(size(r,2),size(phi,2));
for i = 1 : size(r, 2)
for j = 1 : size(phi, 2)
etaisyydet = sqrt(sum((sensors - r(i)*repmat([ cos(phi(j)) ;...
  sin(phi(j))], 1, n)).^2));
posterior_values(i, j) = exp(-sum((v - 1./etaisyydet).^2)...
/(2*sigma^2) );
end 
end

posterior_values = 100*posterior_values/max(posterior_values(:));

cm_point = mean(Pisteet');

figure(1); clf;
set(gcf, 'Renderer', 'painters'); 
pcolor(R.*cos(Phi), R.*sin(Phi), posterior_values');
shading flat;
c_map = colormap('gray');
colormap(flipud(c_map));
hold on;
set(gca,'ylim',[-1.1 1.1]);
set(gca,'xlim',[-1.1 1.1]);
plot(cos([0:0.002:2*pi]), sin([0:0.002:2*pi]),'linewidth',10);
plot( sensors(1, :), sensors(2, :), 'ro','markersize',24, ...
    'linewidth',10)
plot( particle(1), particle(2),'mo', ...
    'markersize',24,'linewidth',5);
set(gca,'visible','off');
h_colorbar = colorbar('vertical','fontsize',35);
set(h_colorbar,'ylim',[0 100]);
set(h_colorbar,'ytick',[0  50 100]);
set(h_colorbar,'yticklabel',{'0','50','100'});
axis equal;
hold off;

figure(2); clf;
set(gcf, 'Renderer', 'painters'); 
pcolor(R.*cos(Phi), R.*sin(Phi), posterior_values');
shading flat;
c_map = colormap('gray');
colormap(flipud(c_map));
hold on;
set(gca,'ylim',[-1.1 1.1]);
set(gca,'xlim',[-1.1 1.1]);
plot(cos([0:0.002:2*pi]), sin([0:0.002:2*pi]),'linewidth',10);

plot( sensors(1, :), sensors(2, :), 'ro','markersize',24, ...
    'linewidth',10)
plot( particle(1), particle(2),'mo', ...
    'markersize',24,'linewidth',10);
plot(Pisteet(1, :), Pisteet(2, :), 'b.'); 
plot( particle(1), particle(2),'mo', ...
    'markersize',24,'linewidth',5);
plot( cm_point(1), cm_point(2),'go', ...
    'markersize',24,'linewidth',5);
set(gca,'visible','off');
h_colorbar = colorbar('vertical','fontsize',35);
set(h_colorbar,'ylim',[0 100]);
set(h_colorbar,'ytick',[0  50 100]);
set(h_colorbar,'yticklabel',{'0','50','100'});
axis equal;
text(-0.58,-1.2,['AR = ' num2str(acceptance_pros) ' %'],'fontsize',35);
hold off;

figure(3); clf;
set(gcf, 'Renderer', 'painters'); 
plot(Pisteet(1,:)','linewidth',2);
hold on;
std_val = std(Pisteet(1,:));
h_line = line([1 size(Pisteet,2)], [ particle(1) particle(1)]);
set(h_line,'linewidth',5,'color','g','linestyle','--');
h_line = line([1 size(Pisteet,2)], [ cm_point(1) cm_point(1)]);
set(h_line,'linewidth',5,'color','m','linestyle','--');
h_line = line([1 size(Pisteet,2)], [ cm_point(1) cm_point(1)] + 2*std_val);
set(h_line,'linewidth',5,'color','k','linestyle','--');
h_line = line([1 size(Pisteet,2)], [ cm_point(1) cm_point(1)] - 2*std_val);
set(h_line,'linewidth',5,'color','k','linestyle','--');
set(gca,'xlim',[1 size(Pisteet,2)]);
set(gca,'xtick',[1 round(size(Pisteet,2)/2) size(Pisteet,2)]);
set(gca,'xticklabel',{'1',int2str(round(size(Pisteet,2)/2)),int2str(size(Pisteet,2))});
xlabel(['STD = ' num2str(std_val)],'fontsize',35);
set(gca,'ylim',[-1 1]);
set(gca,'linewidth',5);
set(gca,'fontsize',35);
pbaspect([2 1 1]); 
hold off;

figure(4); clf;
plot(Pisteet(2,:)','linewidth',2);
hold on;
std_val = std(Pisteet(2,:));
h_line = line([1 size(Pisteet,2)], [ particle(2) particle(2)]);
set(h_line,'linewidth',5,'color','g','linestyle','--');
h_line = line([1 size(Pisteet,2)], [ cm_point(2) cm_point(2)]);
set(h_line,'linewidth',5,'color','m','linestyle','--');
h_line = line([1 size(Pisteet,2)], [ cm_point(2) cm_point(2)] + 2*std_val);
set(h_line,'linewidth',5,'color','k','linestyle','--');
h_line = line([1 size(Pisteet,2)], [ cm_point(2) cm_point(2)] - 2*std_val);
set(h_line,'linewidth',5,'color','k','linestyle','--');
set(gca,'xlim',[1 size(Pisteet,2)]);
set(gca,'xtick',[1 round(size(Pisteet,2)/2) size(Pisteet,2)]);
set(gca,'xticklabel',{'1',int2str(round(size(Pisteet,2)/2)),int2str(size(Pisteet,2))});
xlabel(['STD = ' num2str(std_val)],'fontsize',35);
set(gca,'ylim',[-1 1]);
set(gca,'linewidth',5);
set(gca,'fontsize',36);
pbaspect([2 1 1]); 
hold off;

figure(5); clf;
set(gcf, 'Renderer', 'painters'); 
pbaspect([2 1 1]); 
hold on;
std_val = std(Pisteet(2,:));
res_aux = 30;
t_vec = [-1 : 2/(res_aux-1) : 1];
[hist_vec_1 hist_vec_2] = hist(Pisteet(2,:),t_vec);
hist_vec_1 = hist_vec_1/sum(hist_vec_1*2/(res_aux-1));
h_bar = bar(hist_vec_2,hist_vec_1,1);
set(h_bar,'facecolor','b','linewidth',5);
res_aux = 200;
t_vec = [-1 : 2/(res_aux-1) : 1];
marginal_vec = zeros(res,1);
for i = 1 : res_aux
for j = 1 : res_aux
etaisyydet = sqrt(sum((sensors - repmat([ t_vec(j) ;...
  t_vec(i) ], 1, n)).^2));
marginal_vec(i) = marginal_vec(i) + exp(-sum((v - 1./etaisyydet).^2)...
/(2*sigma^2) );
end 
end
marginal_vec = marginal_vec/sum(marginal_vec*2/(res_aux-1));
plot(t_vec,marginal_vec,'linewidth',5,'color','r','linestyle','--');
set(gca,'fontsize',36);
set(gca,'xtick',[-1 0 1]);
set(gca,'xlim',[-1 1]);
y_lim = get(gca,'ylim');
set(gca,'ytick',[0 round(5*y_lim(2))/10 y_lim(2)]);
h_line = line([particle(2) particle(2)],[0 y_lim(2)]);
set(h_line,'linewidth',5,'color','g','linestyle','--');
h_line = line([ cm_point(2) cm_point(2)],[0 y_lim(2)]);
set(h_line,'linewidth',5,'color','m','linestyle','--');
h_line = line([ cm_point(2) cm_point(2)] + 2*std_val,[0 y_lim(2)]);
set(h_line,'linewidth',5,'color','k','linestyle','--');
h_line = line([ cm_point(2) cm_point(2)] - 2*std_val,[0 y_lim(2)]);
set(h_line,'linewidth',5,'color','k','linestyle','--');
set(gca,'linewidth',5);
set(gca,'box','on');
hold off;

figure(6); clf;
set(gcf, 'Renderer', 'painters'); 
pbaspect([2 1 1]); 
hold on;
std_val = std(Pisteet(1,:));
res_aux = 30;
t_vec = [-1 : 2/(res_aux-1) : 1];
[hist_vec_1 hist_vec_2] = hist(Pisteet(1,:),t_vec);
hist_vec_1 = hist_vec_1/sum(hist_vec_1*2/(res_aux-1));
h_bar = bar(hist_vec_2,hist_vec_1,1);
set(h_bar,'facecolor','b','linewidth',5);
res_aux = 200;
t_vec = [-1 : 2/(res_aux-1) : 1];
marginal_vec = zeros(res,1);
for i = 1 : res_aux
for j = 1 : res_aux
etaisyydet = sqrt(sum((sensors - repmat([ t_vec(i) ;...
  t_vec(j) ], 1, n)).^2));
marginal_vec(i) = marginal_vec(i) + exp(-sum((v - 1./etaisyydet).^2)...
/(2*sigma^2) );
end 
end
marginal_vec = marginal_vec/sum(marginal_vec*2/(res_aux-1));
plot(t_vec,marginal_vec,'linewidth',5,'color','r','linestyle','--');
set(gca,'fontsize',36);
set(gca,'xtick',[-1 0 1]);
set(gca,'xlim',[-1 1]);
y_lim = get(gca,'ylim');
set(gca,'ytick',[0 round(5*y_lim(2))/10 y_lim(2)]);
h_line = line([particle(1) particle(1)],[0 y_lim(2)]);
set(h_line,'linewidth',5,'color','g','linestyle','--');
h_line = line([ cm_point(1) cm_point(1)],[0 y_lim(2)]);
set(h_line,'linewidth',5,'color','m','linestyle','--');
h_line = line([ cm_point(1) cm_point(1)] + 2*std_val,[0 y_lim(2)]);
set(h_line,'linewidth',5,'color','k','linestyle','--');
h_line = line([ cm_point(1) cm_point(1)] - 2*std_val,[0 y_lim(2)]);
set(h_line,'linewidth',5,'color','k','linestyle','--');
set(gca,'linewidth',5);
set(gca,'box','on');
hold off;

aux_length = 250;

figure(7); clf;
set(gcf, 'Renderer', 'painters'); 
plot(Pisteet(1,1:aux_length)','linewidth',5);
hold on;
std_val = std(Pisteet(1,:));
h_line = line([1 size(Pisteet,2)], [ particle(1) particle(1)]);
set(h_line,'linewidth',5,'color','g','linestyle','--');
h_line = line([1 size(Pisteet,2)], [ cm_point(1) cm_point(1)]);
set(h_line,'linewidth',5,'color','m','linestyle','--');
h_line = line([1 size(Pisteet,2)], [ cm_point(1) cm_point(1)] + 2*std_val);
set(h_line,'linewidth',5,'color','k','linestyle','--');
h_line = line([1 size(Pisteet,2)], [ cm_point(1) cm_point(1)] - 2*std_val);
set(h_line,'linewidth',5,'color','k','linestyle','--');
set(gca,'xlim',[1 aux_length]);
set(gca,'xtick',[1 round(aux_length/2) aux_length]);
set(gca,'xticklabel',{'1',int2str(round(aux_length/2)),int2str(aux_length)});
set(gca,'ylim',[-1 1]);
set(gca,'linewidth',5);
set(gca,'fontsize',35);
pbaspect([2 1 1]); 
hold off;


figure(8); clf;
plot(Pisteet(2,:)','linewidth',5);
hold on;
std_val = std(Pisteet(2,1:aux_length));
h_line = line([1 size(Pisteet,2)], [ particle(2) particle(2)]);
set(h_line,'linewidth',5,'color','g','linestyle','--');
h_line = line([1 size(Pisteet,2)], [ cm_point(2) cm_point(2)]);
set(h_line,'linewidth',5,'color','m','linestyle','--');
h_line = line([1 size(Pisteet,2)], [ cm_point(2) cm_point(2)] + 2*std_val);
set(h_line,'linewidth',5,'color','k','linestyle','--');
h_line = line([1 size(Pisteet,2)], [ cm_point(2) cm_point(2)] - 2*std_val);
set(h_line,'linewidth',5,'color','k','linestyle','--');
set(gca,'xlim',[1 aux_length]);
set(gca,'xtick',[1 round(aux_length/2) aux_length]);
set(gca,'xticklabel',{'1',int2str(round(aux_length/2)),int2str(aux_length)});
set(gca,'ylim',[-1 1]);
set(gca,'linewidth',5);
set(gca,'fontsize',36);
pbaspect([2 1 1]); 
hold off;

print(1,'-r300','-dpng','pic_1.png')
print(2,'-r300','-dpng','pic_2.png')
print(3,'-r300','-dpng','pic_3.png')
print(4,'-r300','-dpng','pic_4.png')
print(5,'-r300','-dpng','pic_5.png')
print(6,'-r300','-dpng','pic_6.png')
print(7,'-r300','-dpng','pic_7.png')
print(8,'-r300','-dpng','pic_8.png')
