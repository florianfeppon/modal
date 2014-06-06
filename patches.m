clf;
getd = @(p)path(p,path);

getd('toolbox_signal/');
getd('toolbox_general/');

% Size of the image : n*n
n = 108;
c = [100 200];

% Loading the image


%f0 = load_image('input_small.png');
f0=imread('input_small.png');
f0=double(f0(1:108,1:108));

%f0 = rescale( crop(f0,n, c) );


% Adding noise to the image
sigma = .04;
f = f0 + randn(n,n)*sigma;

%Size of the patches : 2*w+1
w = 2;
w1 = 2*w+1;

% location of pixels
[Y,X] = meshgrid(1:n,1:n);
% offsets
[dY,dX] = meshgrid(-w:w,-w:w);
% location of pixels to extract
dX = reshape(dX, [1 1 w1 w1]);
dY = reshape(dY, [1 1 w1 w1]);
X = repmat(X, [1 1 w1 w1]) + repmat(dX, [n n 1 1]);
Y = repmat(Y, [1 1 w1 w1]) + repmat(dY, [n n 1 1]);

% Boundary condition handled by reflection
X(X<1) = 2-X(X<1); Y(Y<1) = 2-Y(Y<1);
X(X>n) = 2*n-X(X>n); Y(Y>n) = 2*n-Y(Y>n);

% Patch extractor operator
patch = @(f0)f0(X + (Y-1)*n);

P = patch(f0);

 clf;
 %imageplot( squeeze(P(1,1,:,:)), 'P0', 1,2,1 );
% 
% % Setting parameters 
 s=3;
 neighbors=9;
% 
% % Looking for the nearest neighbors
 Pprim = reshape(P,[n*n w1*w1]);
 disp('Calcul des voisins');
 [IDX,D]=knnsearch(Pprim,Pprim,'K',neighbors);
 disp('Fin du calcul');
 % Size of the patches : 2*w+1
 
   fHR=imresize(f0,s,'nearest');
% % % 
% % % 
%   w1HR = w1*s;
%   wHR=(w1HR-1)/2;
% % % A
% % % % location of pixels
%   [YHR,XHR] = meshgrid(1:n*s,1:n*s);
% % 
% % % % offsets
%   [dYHR,dXHR] = meshgrid(-wHR:wHR,-wHR:wHR);
% % % % location of pixels to extract
%   dXHR = reshape(dXHR, [1 1 w1HR w1HR]);
%   dYHR = reshape(dYHR, [1 1 w1HR w1HR]);
%   XHR = repmat(XHR, [1 1 w1HR w1HR]) + repmat(dXHR, [n*s n*s 1 1]);
%  YHR = repmat(YHR, [1 1 w1HR w1HR]) + repmat(dYHR, [n*s n*s 1 1]);
%  
% % % Boundary condition handled by reflection
%   XHR(XHR<1) = 2-XHR(XHR<1); YHR(YHR<1) = 2-YHR(YHR<1);
%   XHR(XHR>n*s) = 2*n*s-XHR(XHR>n*s); YHR(YHR>n*s) = 2*n*s-YHR(YHR>n*s);
% % % 
% % % % Patch extractor operator
%   patchHR = @(fHR)fHR(XHR + (YHR-1)*n*s);
%  % 
%   PHR = patchHR(fHR);

%PHR(s*x,s*y)<->P(x,y)
%imageplot(squeeze(P(65,65,:,:)),'',1,2,1);
%imageplot(squeeze(PHR(65*3-1,65*3-1,:,:)),'',1,2,2);

Gy=[1,2,1;0,0,0;-1,-2,-1];
GradientY = conv2(f0, Gy, 'same');

Gx=[1,0,-1;2,0,-2;1,0,-1];
GradientX=conv2(f0,Gx,'same');
%  A=zeros(n,n,neighbors);
%  B=zeros(n,n,neighbors);
gx=patch(GradientX);
gy=patch(GradientY);
disp('Début des misalignments');

% for i=1:n
%     for j=1:n
%         for k=2:neighbors
%             index=(j-1)*n+i;
%             gt=-Pprim(index,:)+Pprim(IDX(index,k),:);
%             gxx=reshape(squeeze(gx(i,j,:,:)),[1,w1*w1]);
%             gyy=reshape(squeeze(gy(i,j,:,:)),[1,w1*w1]);
%             det=sum(sum(gxx.*gxx))*sum(sum(gyy.*gyy))-sum(sum(gxx.*gyy))^2;
%             A(i,j,k)=(sum(sum(gxx.*gt))*sum(sum(gyy.*gyy))-sum(sum(gxx.*gyy))*sum(sum(gt.*gyy)))/det;
%             B(i,j,k)=(sum(sum(gxx.*gxx))*sum(sum(gyy.*gt))-sum(sum(gxx.*gt))*sum(sum(gxx.*gyy)))/det;
%         end
%     end
% end
% disp('Fin des misalignments');
fHR2=255*ones(s*n,n*s);
fHR2=fHR;
disp('Début du processing');
for nb=1:1
    %fHR = fHR2;
    for i=1:n
        for j=1:n
            for k=2:neighbors
                for l=-w:w
                    for m=-w:w
                if abs(A(i,j,k))<1 && abs(B(i,j,k))<1
                    index1=s*(i-1)+1+round(A(i,j,k)*s);
                    index2=s*(j-1)+1+round(B(i,j,k)*s);
                    PprimPrim=reshape(Pprim(IDX((j-1)*n+i,k),:),[w1,w1]);
                    if index1+s*(l+w)<=n*s && index2+s*(m+w)<=n*s
                    fHR2(index1+s*(l+w),index2+s*(m+w))= fHR(index1+s*(l+w),index2+s*(m+w))+1/4*(PprimPrim(w+l+1,w+m+1)-fHR(index1+s*(l+w),index2+s*(m+w)));
                    end
                end
                    end
                end
            end
        end
    end
end
disp('Fin du processing');
imageplot(fHR2,'SR',1,2,1);
imageplot(fHR,'LR',1,2,2);
%Fprim=image_offset(f0,s,3,2);

%disp('Calcul des offsets')
%   A=zeros(n,n,s,s);
%   for l1=1:s
%       for l2=1:s
%           A(:,:,l1,l2)=image_offset(f0,s,l1-1,l2-1);
%       end
%   end
% disp('Fin du calcul')
% 
% %Compute offset
% disp('debut des misalignments')
%  
%   B = zeros(s, s, n, n, w1, w1);
%   for l1=1:s
%      for l2=1:s
%           B(l1, l2, :, :, :, :) = patch(A(:, :, l1, l2));
%      end
%   end
%   
%   subpixelalignment_l1 = zeros(n,n,neighbors);
%   subpixelalignment_l2 = zeros(n,n,neighbors);
%   
%   C=zeros(n,n,neighbors,2);
%   for i=1:n
%       for j=1:n
%           for k=1:neighbors
%               C(i,j,k,1)=max([n mod(IDX(n*(j-1)+i, k)-1, n)+1]);
%               C(i,j,k,2)=max([floor((IDX(n*(j-1)+i, k)-1) / n)+1 n]);
%           end
%       end
%   end
%   
%   disp('début knn');
%   voisins=reshape(B,[n*n*s*s, w1*w1]);
%    (l1-1)*s*n*n+(l2-1)*n*n+(i-1)*n+(j-1)
%   [meilleurs_voisins,D]=knnsearch(voisins,voisins,'K',neighbors);
%   disp('fin');
%   
%    for i=1:n
%        for j=1:n
%            for k=1:neighbors
%                min = -1;
%                for l1=1:s
%                    for l2=1:s
%                       %disp(n*(i-1)+j)
%                       temp=distanceBetweenPatches(squeeze(P(i,j,:,:)),squeeze(B(l1, l2, C(i,j,k,1), C(i,j,k,2), :, :)));
%                       
%                       if min==-1 || temp < min
%                           min = temp;
%                           subpixelalignment_l1(i,j,k) = l1-1;
%                           subpixelalignment_l2(i,j,k) = l2-1;
%                       end
%                   end
%               end
%    %        supbixelalignment_l1(i,j,k)=
%          end
%       end
%   end
% disp('Fin des misalignments')
% 
% % 
%  sommevaleurspossibles = zeros(n*s, n*s);
%  nombrevaleurspossibles = zeros(n*s, n*s);
% moyennevaleurspossibles = -ones(n*s, n*s);
% % 
%  for i=1:n
%      for j=1:n
%          for k=1:neighbors
%              for i_p=-w:w
%                  for j_p=-w:w
%                      i_k = mod(IDX(n*(j-1)+i, k)-1, n)+1;
%                      j_k = floor((IDX(n*(j-1)+i, k)-1) / n)+1;
%                      if i+i_p > 0 && j+j_p > 0 && i+i_p < n && j+j_p < n
%                          sommevaleurspossibles((i+i_p-1)*s+subpixelalignment_l1(i,j,k)+1, (j+j_p-1)*s+subpixelalignment_l2(i,j,k)+1) = sommevaleurspossibles((i+i_p-1)*s+subpixelalignment_l1(i,j,k)+1, (j+j_p-1)*s+subpixelalignment_l2(i,j,k)+1)+P(i_k,j_k,i_p+w+1,j_p+w+1);
%                          nombrevaleurspossibles((i+i_p-1)*s+subpixelalignment_l1(i,j,k)+1, (j+j_p-1)*s+subpixelalignment_l2(i,j,k)+1) = nombrevaleurspossibles((i+i_p-1)*s+subpixelalignment_l1(i,j,k)+1, (j+j_p-1)*s+subpixelalignment_l2(i,j,k)+1)+1;
%                      end
%                  end
%              end
%          end
%      end
%  end
% disp('fin du calcul des valeurs possibles')
% 
% 
% for i=1:n*s
%     for j=1:n*s
%         if nombrevaleurspossibles(i,j) ~= 0
%             moyennevaleurspossibles(i,j) = sommevaleurspossibles(i,j)/nombrevaleurspossibles(i,j);
%         end
%     end
% end
% 
% fenetre=2;
% for i=1:n*s
%     for j=1:n*s
%         somme=0;nombre=0;
%         if moyennevaleurspossibles(i,j) == -1
%             for l=-fenetre:fenetre
%                 for m=-fenetre:fenetre
%                     if 0<i+l && i+l<=n*s && 0<j+m && j+m<=n*s
%                         if moyennevaleurspossibles(i+l,j+m) ~= -1
%                %   for i=1:n
%        for j=1:n
%            for k=1:neighbors
%                min = -1;
%                for l1=1:s
%                    for l2=1:s
%                       %disp(n*(i-1)+j)
%                       temp=distanceBetweenPatches(squeeze(P(i,j,:,:)),squeeze(B(l1, l2, C(i,j,k,1), C(i,j,k,2), :, :)));
%                       
%                       if min==-1 || temp < min
%                           min = temp;
%                           subpixelalignment_l1(i,j,k) = l1-1;
%                           subpixelalignment_l2(i,j,k) = l2-1;
%                       end
%                   end
%              end
%           end
%        end
%    end
%          somme= somme+moyennevaleurspossibles(i+l,j+m);
%                         nombre=nombre+1;
%                         end
%                     end
%                 end
%             end
%             moyennevaleurspossibles(i,j)=double(somme)/double(nombre);
%         end
%     end
% end
% 
% 
% 
% imageplot(f0, '', 1, 2, 1);
% imageplot(moyennevaleurspossibles, '', 1, 2, 2);
% disp('fin du programme')
% 
% imageplot(squeeze(P(4, 4, :, :)), '', 1, 2, 1);
% imageplot(f0(2:6, 2:6), '', 1, 2, 2);
    
% 
% imageplot(Fprim,'fprim',1,2,2);
% imageplot(f,'f',1,2,1);
