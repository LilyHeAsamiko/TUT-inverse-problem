%This is a script for Example 5. @ Sampsa Pursiainen, 2012.
N = 64;
sigma = 30; 
alpha = 100;

load ground_layers_matrix.mat A;
A = full(A);

for k_1 = 1 : 2

    if k_1 == 1 
     fn = 'layers_a';
    elseif k_1 == 2
     fn = 'layers_b';
    end
    
    for k_2 = 1 : 2
       
    W = eye(N^2,N^2);    
    if k_2 == 2 
     Ind_mat = reshape([1 : N^2]',N,N); 
     if k_1 == 1
      phi = pi/9;
     elseif k_2 == 2
      phi = pi/6;
     end
     W = cos(phi)*W + sin(phi)*W;
     for i = 2 : N
      for j = 1 : N-1 
       W(Ind_mat(i,j),Ind_mat(i,j+1)) = W(Ind_mat(i,j),Ind_mat(i,j+1)) - cos(phi);
       W(Ind_mat(i,j),Ind_mat(i-1,j)) = W(Ind_mat(i,j),Ind_mat(i-1,j)) - sin(phi);
      end
     end
    end
        
image_exact = double(imread([fn '.bmp'],'bmp')); 
figure(1); clf; imagesc(image_exact); colormap('gray');
axis equal; set(gca,'visible','off');
drawnow; 

y = A*image_exact(:);
y = y + sigma*randn(size(y));


A = A*inv(W); 
x = W\(A'*(( A * A' + (sigma^2/alpha^2)*eye(size(A,1),size(A,1)))\y));

image_x = reshape(x,N,N);
figure(2); clf; imagesc(image_x); colormap('gray');
axis equal; set(gca,'visible','off');
drawnow;
print(2,'-r300','-dpng',[fn '_' int2str(k_2) '.png'])

    end
end
