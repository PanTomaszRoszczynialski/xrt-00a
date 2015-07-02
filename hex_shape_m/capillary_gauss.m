function [z]=capillary_circle(X,Y,x0,y0,capillary_diameter,channel_diameter)

R2=(X-x0).^2+(Y-y0).^2;
%z=(R2<=(capillary_diameter/2).^2)/2+(R2<=(channel_diameter/2).^2)/2;
z=exp(-0.5*R2/(channel_diameter/4)^2);


