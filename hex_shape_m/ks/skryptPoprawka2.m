function b = skryptPoprawka2(nhex,nkap,r1,r2,chaos,gauss)
tic;
%close all;
value = 1; % zmiana: parametr opisujacy fluktuacje intensywnosci
sigmaintens=0.2;
sigmareject=0.002;

if (r1-r2)<=3
    r2 = r1-3;
end
wekt = Kolo3poprawka(nkap,r1,r2,value,chaos,gauss,sigmaintens,sigmareject);

[wektsizey, wektsizex] = size(wekt);
dsizex = 2*wektsizex;
rozmiar = nhex*dsizex+1;
xhcpsr = rozmiar;
yhcpsr = rozmiar;
dhx = wektsizex;
dhy = wektsizey;
points = zeros(2*rozmiar);
%%%%%
points(yhcpsr+dhy,xhcpsr) = -1;
for i=1:(nhex-1)
    
    yy = yhcpsr + ceil((i+1)*dhy);
    yyy = yhcpsr - ceil((i-1)*dhy);
    
    xx = xhcpsr + round(3*dhx*i/4)+2*i*(r1-r2);
    xxx = xhcpsr - round(3*dhx*i/4)-2*i*(r1-r2);
    
    points(yy,xhcpsr) = -1;
    points(yyy,xhcpsr) = -1;
    
    for j=1:(2*nhex-1-i)
        dy = (-ceil((2*nhex-i-2)*dhy)/2) + j*dhy;
        dyy = (ceil((2*nhex-i-2)*dhy)/2) - j*dhy;
        yy = ceil(yhcpsr+dy);
        yyy = ceil(yhcpsr - dyy);
        points(yy,xx) = -1;
        points(yyy,xxx) = -1;
    end
end


for i=1:size(points,1)
    i
    for j=1:size(points,2)
        if points(i,j) == -1
           points(i,j) = 0;
            wekt = Kolo3poprawka(nkap,r1,r2,value,chaos,gauss,sigmaintens,sigmareject);
            [wektsizey, wektsizex] = size(wekt);
            for k=1:wektsizey
                for m=1:wektsizex
                     if points(i-k,j-m) == 0
                    points(i-k,j-m) =points(i-k,j-m) + wekt(k,m);
                     end
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%
disp('Here!')
 %figure;imagesc(points),colormap('gray'),axis('square');
[x,y] = find(points);
colymin = y(1);
colymax = y(end);
rowxmin = min(x);
rowxmax = max(x);
%% wiersze,kolumny
angleData = points(rowxmin:rowxmax,colymin:colymax);
clear points;
%figure;imagesc(angleData),colormap('gray'),axis('equal');
[fy fx] = size(angleData);
if fy>=fx
    fs = fy;
    fd = fx;
else
    fs = fx;
    fd = fy;
end
liczbaC = floor((fs-fd)/2);
addCol = zeros(fs,1);
% disp(size(angleData));
for i=0:(liczbaC-1)
angleData = cat(2,angleData,addCol);  % Insert as column    
angleData = cat(2,addCol,angleData);
end
% disp(size(angleData));
disp('Here2!')

imagesc(angleData),colormap('gray'),axis('equal');
fsx = nearestpow2(fs);
fsx = 2^fsx;
% disp([fs,fsx]);
 b = matrixResize(angleData,fsx);
 clear angleData;
 imagesc(b),colormap('gray'),axis('equal');
%  disp(size(b));
 toc;
% 
