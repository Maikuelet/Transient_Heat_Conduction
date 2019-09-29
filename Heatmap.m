%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% PROCESSING 2D HT  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
close all;
clc

%Read .txt data
A=dlmread('nodeMATRIX.txt');
sizeX = size(A,1);
sizeY = size(A,2);

%Generate plot vectors
x=linspace(0,1,sizeX);
y=linspace(0,1,sizeY);


%PLOTTING
%-------------------------------------------------
figure(1);
h = heatmap(A);  %heatmap
colormap(jet);
colorbar;
title('Temperature (Steady State)');
xlabel('x'), ylabel('y');

figure(2);
contour(A);     %mostra les isotermes
colormap(jet);
colorbar;
title('Temperature isoterms (Steady State)');
xlabel('x'), ylabel('y');

figure(3);
contourf(x,y,A);     %mostra les isotermes
colormap(jet);
colorbar;
title('Temperature isoterms (Steady State)');
xlabel('x'), ylabel('y');

figure(4);
pcolor(x,y,A);
shading interp; 
colormap(jet);
colorbar;
title('Temperature (Steady State)');
xlabel('x'), ylabel('y');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% POSTPROCESSING  %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T=zeros(sizeX,1);

for j=1:sizeY
    i = round(sizeX/2);
    T(j) = A(i,j);
end

figure(5)
yy= zeros(sizeY)+  max(T);
p=plot(x,T,x,yy,'--k','LineWidth',1);
title('Temperature over middle nodes in Y direction');
xlabel('Domain size (x)'), ylabel('Temperature (ºC)');

ylim([min(T)-10 max(T)+10]);
p(1).LineWidth = 2;


T=zeros(sizeY,1);

for i=1:sizeX
    j = round(sizeY/2);
    T(i) = A(i,j);
end

figure(6)
yy= zeros(sizeY)+  max(T);
p=plot(x,T,x,yy,'--k','LineWidth',1);
title('Temperature over middle nodes in X direction');
xlabel('Domain size (x)'), ylabel('Temperature (ºC)');

ylim([min(T)-10 max(T)+10]);
p(1).LineWidth = 2;

