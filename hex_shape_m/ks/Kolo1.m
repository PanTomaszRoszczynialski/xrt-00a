function i = Kolo1(r1,r2,value)
%r1 - promien ochronki
%r2 - promien srodka kapilary
% value - wartosc piksela
% 
% if (r1-r2)<2
%    r2 = r1-2;
% end

%value=-abs(randn(1));
%value=1+value.*(value>1); % fluktuacje intensywnosci
%rr=rand(1);
%value=value+(rr>=0.002)+eps.*(rr<0.002); % kompletne usuniecie

m = 2*r1+1;
n = 2*r1+1;
xc = m/2;
yc = n/2;
i = zeros(r2);

%% srodek kapilary
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
%% wypelnienie pixelami
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

%% ochronka - wyznaczane sa 4 punkty, to jest wystarczajace o ile pixele ochronki nie musza miec konkretnej wartosci (bo domyslnie maja 0)
%% ale przypisujemy wartosc = 0, zeby zachowac odpowiednie odleglosci miedzy kapilarami
xc = int16(m/2);
yc = int16(n/2);
x = int16(0);
y = int16(r1);
value = 0;
i(xc, yc+y) = value;
i(xc, yc-y) = value;
i(xc+y, yc) = value;
i(xc-y, yc) = value;
%% gdyby pixele ochronki musialyby miec jakas konkretna wartosc rozna od 0:
% while ( x < y - 1 )
%     x = x + 1;
%     if ( d < 0 )
%         d = d + x + x + 10;
%     else
%         y = y - 1;
%         a = x - y + 1;
%         d = d + a + a;
%     end
%     i( x+xc,  y+yc) = value;
%     i( y+xc,  x+yc) = value;
%     i( y+xc, -x+yc) = value;
%     i( x+xc, -y+yc) = value;
%     i(-x+xc, -y+yc) = value;
%     i(-y+xc, -x+yc) = value;
%     i(-y+xc,  x+yc) = value;
%     i(-x+xc,  y+yc) = value;
% end
%  disp(size(i));
%imagesc(i), axis('square'),title('Pojedyncza kapilara - funkcja Kolo1'), colormap('gray');




