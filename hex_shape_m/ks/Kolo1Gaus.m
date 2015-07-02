function i = Kolo1Gaus(r1,r2,value)
%r1 - promien ochronki
%r2 - promien srodka kapilary
%dx,dy - przesuniecie, chaos

m = 2*r1+1;
n = 2*r1+1;
xc = m/2;
yc = n/2;
i = zeros(r1);

%chaos


%srodek kapilary
xc = int16(xc);
yc = int16(yc);
x = int16(0);
y = int16(r2);
d = int16(1-r2);
i(xc, yc+y) = value;
i(xc, yc-y) = value;
i(xc+y, yc) = value;
i(xc-y, yc) = value;

while ( x < y - 1 )
    x = x + 1;
    if ( d < 0 )
        d = d + x + x + 10;
    else
        y = y - 1;
        a = x - y + 1;
        d = d + a + a;
    end
    
    i( x+xc,  y+yc) = value;
    i( y+xc,  x+yc) = value;
    i( y+xc, -x+yc) = value;
    i( x+xc, -y+yc) = value;
    i(-x+xc, -y+yc) = value;
    i(-y+xc, -x+yc) = value;
    i(-y+xc,  x+yc) = value;
    i(-x+xc,  y+yc) = value;
    
end
%wypelnienie pixelami
radius2 = r2;
for s =1:(xc+radius2)
    for w = 1:(yc+radius2)
        dx = s - xc;
        dy = w - yc;
        if( (dx*dx + dy*dy) <= (radius2*radius2))
            i(s,w) = value;
        end
    end
end



% ochronka
xc = int16(m/2);
yc = int16(n/2);
x = int16(0);
y = int16(r1);
value = 0;
i(xc, yc+y) = value;
i(xc, yc-y) = value;
i(xc+y, yc) = value;
i(xc-y, yc) = value;
while ( x < y - 1 )
    x = x + 1;
    if ( d < 0 )
        d = d + x + x + 10;
    else
        y = y - 1;
        a = x - y + 1;
        d = d + a + a;
    end
    i( x+xc,  y+yc) = value;
    i( y+xc,  x+yc) = value;
    i( y+xc, -x+yc) = value;
    i( x+xc, -y+yc) = value;
    i(-x+xc, -y+yc) = value;
    i(-y+xc, -x+yc) = value;
    i(-y+xc,  x+yc) = value;
    i(-x+xc,  y+yc) = value;
end

%imagesc(i), axis('square'),title('Pojedyncza kapilara - funkcja Kolo1'), colormap('gray');

%Gauss
%imagesc(i),colormap('gray');
pause(1);
mu = [0 0];
Sigma = [1 .10; .10 1];
n = 0.1*2*r1;
x1 = -n:.2:n; x2 = -n:.2:n;

[X1,X2] = meshgrid(x1,x2);
F = floor(100*mvnpdf([X1(:) X2(:)],mu,Sigma));
F = reshape(F,length(x2),length(x1));
%figure;imagesc(F),colormap('gray');
i = i + F;
%figure;
%str = sprintf('Pojedyncza kapilara  - Gauss');
%imagesc(i),colormap('gray'),title(str),axis('equal');




