clear all ;

%   Lecture de l'image et importation des données relevées
   
im = imread('20170201_155908.jpg');
im_625 = imread('20170201_160625.jpg');
xyz = importdata('xyz.txt',',');
uv = importdata('uv_908.txt',',');
uv_625 = importdata('uv_625.txt',',');


%   Création de la matrice C 64*12 avec toutes les valeurs à 0.

C = zeros(64,12);
C2 = zeros(64,12);

%   On calcul la matrice 64*12

for i =1:32
    C(2*i-1,:) = [xyz(i,:) 1 0 0 0 0 -uv(i,1)*xyz(i,1) -uv(i,1)*xyz(i,2) -uv(i,1)*xyz(i,3) -uv(i,1)];
    C(2*i,:) = [0 0 0 0 xyz(i,:) 1 -uv(i,2)*xyz(i,1) -uv(1,2)*xyz(i,2) -uv(i,2)*xyz(i,3) -uv(i,2)];
    C2(2*i-1,:) = [xyz(i,:) 1 0 0 0 0 -uv_625(i,1)*xyz(i,1) -uv_625(i,1)*xyz(i,2) -uv_625(i,1)*xyz(i,3) -uv_625(i,1)];
    C2(2*i,:) = [0 0 0 0 xyz(i,:) 1 -uv_625(i,2)*xyz(i,1) -uv_625(1,2)*xyz(i,2) -uv_625(i,2)*xyz(i,3) -uv_625(i,2)];
end;

% 

B=C'*C ;
B2=C2'*C2 ;

%

[eigenVectors, eigenValues] = eig(B);
MinIndice = 1 ;

[eigenVectors2, eigenValues2] = eig(B2);
MinIndice2 = 1 ;

%

for i = 2:12
    if((eigenValues(i,i) < eigenValues(MinIndice, MinIndice)) && not(eigenValues(i,i)==0))
        MinIndice = i ;
    end;
    if((eigenValues2(i,i) < eigenValues2(MinIndice2, MinIndice2)) && not(eigenValues2(i,i)==0))
        MinIndice2 = i ;
    end;
end;

M = eigenVectors(:,MinIndice)

M2 = eigenVectors2(:,MinIndice2)

xt = 0;
yt = 59;
zt = 151;

xt2 = 0;
yt2 = 59;
zt2 = 151;

% Calibrations
ut = (M(1)*xt + M(2)*yt + M(3)*zt + M(4))/(M(9)*xt+M(10)*yt+M(11)*zt+M(12));
vt = (M(5)*xt + M(6)*yt + M(7)*zt + M(8))/(M(9)*xt+M(10)*yt+M(11)*zt+M(12));

ut2 = (M2(1)*xt2 + M2(2)*yt2 + M2(3)*zt2 + M2(4))/(M2(9)*xt2+M2(10)*yt2+M2(11)*zt2+M2(12));
vt2 = (M2(5)*xt2 + M2(6)*yt2 + M2(7)*zt2 + M2(8))/(M2(9)*xt2+M2(10)*yt2+M2(11)*zt2+M2(12));

% Stéréovision
ud = 2083 ;
vd = 424 ;
ug = 1850 ;
vg = 501 ;

% AX=E
A=[M2(1)-ug*M2(9) M2(2)-ug*M2(10) M2(3)-ug*M2(11);M2(5)-vg*M2(9) M2(6)-vg*M2(10) M2(7)-vg*M2(11);M(1)-ud*M(9) M(2)-ud*M(10) M(3)-ud*M(11);M(5)-vd*M(9) M(6)-vd*M(10) M(7)-vd*M(11)];
%X=[x;y;z];
E=[-(M2(4)-ug*M2(12));-(M2(8)-vg*M2(12));-(M(4)-ud*M(12));-(M(8)-vd*M(12))];

X=pinv(A)*E

image(im);
%image(im_625);

hold on ;

plot(uv(:,1), uv(:,2),'s');
plot(ut,vt,'md');
%plot(uv_625(:,1), uv_625(:,2),'s');
%plot(ut2,vt2,'md');
