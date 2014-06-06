%Read image
f0=double(imread('input_small.png'));
[height,width]=size(f0);

imageplot(f0,'',1,1,1);

%Get Patches
w=2;
Patches=transpose(im2col(f0,[5,5]));

disp('Calcul des plus proches voisins');
[IDX,D]=knnsearch(Patches,Patches,'K',9);
disp('Fin du calcul');
 

i=42;j=18;
imageplot(f0(i-w:i+w,j-w:j+w),'',1,3,1);
index=getIndex(i,j,w,height);
[a,b]=getIJ(index,w,height);
P=reshape(Patches(index,:),[5 5]);
imageplot(P,'',1,3,2);
imageplot(f0(a-w:a+w,b-w:b+w),'',1,3,3);

