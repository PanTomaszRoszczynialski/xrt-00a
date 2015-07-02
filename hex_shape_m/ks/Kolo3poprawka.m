function angleData = Kolo3poprawka(n,r1,r2,value,chaos,gauss,sigmaintens,sigmareject)

if (r1-r2)<=3
    r2 = r1-3;
end
%% wybieramy liczbe hexagonalna
hn = 2*n^2-n; %% hn - nta liczba hexagonalna
% disp(hn);
%% srodek siatki HCP
d = 2*r1;
rozmiar = (n-1)*d+1;
xhcpsr = rozmiar+r1;
yhcpsr = rozmiar+r1;

points=zeros(3*rozmiar);
clear rozmiar;
% disp(rozmiar);
rozn = r1-r2-2;
points(yhcpsr,xhcpsr+d) =(1-sigmaintens*abs(randn(1))).*(rand(1)>sigmareject);;
for i=1:(n-1)
    if chaos==0
        chx = 0; chy = 0;
        chxx = 0; chyy = 0;
    else
        chx = randi([-rozn,rozn],1);
        chy = randi([-rozn,rozn],1);
        if ((abs(chx) + r2) >= r1) || ((abs(chy) + r2) >= r1)
            chx = 0; chy = 0;
        end
        chxx = randi([-rozn,rozn],1);
        chyy = randi([-rozn,rozn],1);
        if ((abs(chxx) + r2) >= r1) || ((abs(chyy) + r2) >= r1)
            chxx = 0; chyy = 0;
        end
        
        
    end
    
    yy = yhcpsr + ceil(sqrt(3)*d*i/2);
    yyy = yhcpsr - ceil(sqrt(3)*d*i/2);
    xx = xhcpsr + (i+1)*d;
    xxx = xhcpsr - (i-1)*d;
    
    
    
    points(yhcpsr+chy,xx+chx) = (1-sigmaintens*abs(randn(1))).*(rand(1)>sigmareject);;
    points(yhcpsr+chyy,xxx+chxx) = (1-sigmaintens*abs(randn(1))).*(rand(1)>sigmareject);;
    
    for j=1:(2*n-1-i)
        if chaos==0
            chx = 0; chy = 0;
            chxx = 0; chyy = 0;
        else
            chx = randi([-rozn,rozn],1);
            chy = randi([-rozn,rozn],1);
            if ((abs(chx) + r2) >= r1) || ((abs(chy) + r2) >= r1)
                chx = 0; chy = 0;
            end
            chxx = randi([-rozn,rozn],1);
            chyy = randi([-rozn,rozn],1);
            if ((abs(chxx) + r2) >= r1) || ((abs(chyy) + r2) >= r1)
                chxx = 0; chyy = 0;
            end
            
        end
        dx = (-(2*n-i-2)*d)/2 + j*d;
        xx = xhcpsr+dx;
        if(i == (n-1)) && (j == 1)
            chx = 0; chy = 0; chxx = 0; chyy = 0;
            points(yy+chy,xx+chx) = (1-sigmaintens*abs(randn(1))).*(rand(1)>sigmareject);;
            points(yyy+chyy,xx+chxx) = (1-sigmaintens*abs(randn(1))).*(rand(1)>sigmareject);;
        elseif (i == (n-1)) && (abs(chy) >= rozn/2)
            points(yy+chy,xx+chx) = 0;
            points(yyy+chyy,xx+chxx) = (1-sigmaintens*abs(randn(1))).*(rand(1)>sigmareject);;
        elseif (i == (n-1)) && (abs(chyy) >= rozn/2)
            points(yy+chy,xx+chx) = (1-sigmaintens*abs(randn(1))).*(rand(1)>sigmareject);;
            points(yyy+chyy,xx+chxx) = 0;
        elseif (j == 1) && (chx < 0)
            chy = abs(chx);
            points(yy+chy,xx+chx) = (1-sigmaintens*abs(randn(1))).*(rand(1)>sigmareject);;
        elseif (j == 1) && (chxx < 0)
            chyy = abs(chxx);
            points(yyy+chyy,xx+chxx) = (1-sigmaintens*abs(randn(1))).*(rand(1)>sigmareject);;
        elseif (j == (2*n-1-i)) && (chx > 0)
            chy = abs(chx);
            points(yy+chy,xx+chx) = (1-sigmaintens*abs(randn(1))).*(rand(1)>sigmareject);;
        elseif (j == (2*n-1-i)) && (chxx > 0)
            chyy = abs(chxx);
            points(yyy+chyy,xx+chxx) =(1-sigmaintens*abs(randn(1))).*(rand(1)>sigmareject);;
        else
            points(yy+chy,xx+chx) = (1-sigmaintens*abs(randn(1))).*(rand(1)>sigmareject);;
            points(yyy+chyy,xx+chxx) = (1-sigmaintens*abs(randn(1))).*(rand(1)>sigmareject);;
        end
    end
end
% imagesc(points),colormap('gray'),axis('square');
if gauss == 0
    wekt = Kolo1(r1,r2,value);
elseif gauss == 1
    wekt = Kolo1Gaus(r1,r2,value);
end
pointsconv2 = conv2(points,wekt,'same');
clear points;
% figure;imagesc(pointsconv2),colormap('gray'),axis('equal');
% imagesc(pointsconv2),colormap('gray'),axis('square');
%% usuwamy puste kolumny i wiersze
[x,y] = find(pointsconv2);
colymin = y(1);
colymax = y(end);
rowxmin = min(x);
rowxmax = max(x)+2*(r1-r2);
%% wiersze,kolumny
angleData = pointsconv2(rowxmin:rowxmax,colymin:colymax);
clear pointsconv2;
[dhy, dhx] = size(angleData);
str=sprintf('Szesciokat zlozony z kapilar, Chaos: %d w x i %d w y',chx,chy);
% figure;imagesc(angleData),colormap('gray'),title(str),axis('equal');

