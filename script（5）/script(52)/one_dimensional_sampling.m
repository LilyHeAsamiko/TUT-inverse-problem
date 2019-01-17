
sample_size = 100000;
hist_res = 20;
t = [-1 : 0.01 : 1];
t_hist = [-1:2/(hist_res - 1):1];


for i = 1 : 4

if i == 1    
y = cos((pi/2)*t);
elseif i == 2  
y = (t + 1).^2/4;
elseif i == 3
y = (t + 1)/2;
elseif i == 4
y = log(t + 2);
end
    

z = cumsum(y);
z = z/z(end);
    

x = zeros(1,sample_size);

for k  = 1 : sample_size
    
aux_ind = find(z > rand(1),1);
x(k) = t(aux_ind);

end

figure(1); clf; cla;
set(gcf, 'Renderer', 'painters'); 
hold on;
[a,b] = hist(x,t_hist);
a = a*max(y)/max(a);
h_bar = bar(b,a,1);
set(h_bar,'facecolor','b','linewidth',7);
plot(t,y,'r','linewidth',7);
axis equal;
set(gca,'xlim',[-1 1]);
set(gca,'xtick',[-1 0 1]);
set(gca,'ylim',[0 1.2]);
set(gca,'ytick',[0 1]);
set(gca,'fontsize',35);
set(gca,'linewidth',7);
set(gca,'box','on');
hold off;
print(1,'-r300','-dpng',['pic_' int2str(i) '.png']);
end