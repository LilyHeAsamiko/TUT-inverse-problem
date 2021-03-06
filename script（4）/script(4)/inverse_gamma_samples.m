%Example 7. Sampsa Pursiainen 2012.
beta = 1.5;
theta0 = 0.00125;

for i = 1 : 4

figure(1); clf;    
imagesc(1./gamrnd(beta*ones(64),1./(theta0*ones(64))));
colormap('gray'); colorbar('vertical','fontsize',20); axis equal;
set(gca,'visible','off');
print(1,'-r300','-dpng',['sample_' int2str(i) '.png']);

end