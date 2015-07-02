function [xci,yci]=capillary_bundle(xbundle,ybundle,capillary_diameter,channel_diameter,nx_capillary,sigma)

indeks=0;
nxpol_capillary=(nx_capillary-1)/2;

atemp=capillary_diameter*nxpol_capillary;


for ix=-2*nxpol_capillary:2*nxpol_capillary
     
    for iy=-2*nxpol_capillary:2*nxpol_capillary       
 
        
        x0=capillary_diameter*ix + capillary_diameter/2*iy+sigma*(rand(1,1)-0.5)*(capillary_diameter-channel_diameter);
            y0=sqrt(3)/2*capillary_diameter*iy+sigma*(rand(1,1)-0.5)*(capillary_diameter-channel_diameter);
             
             
          
            
            [in_bundle]=isinhexagon(x0,y0,capillary_diameter*nxpol_capillary);
             x0=capillary_diameter*ix + capillary_diameter/2*iy+sigma*(rand(1,1)-0.5)*(capillary_diameter-channel_diameter);
            y0=sqrt(3)/2*capillary_diameter*iy+sigma*(rand(1,1)-0.5)*(capillary_diameter-channel_diameter);
             
            
            if in_bundle
            
            x=xbundle+x0;
            y=ybundle+y0;
            indeks=indeks+1;
            xci(indeks)=x;
            yci(indeks)=y;
            
                   
            end
        
    end
end
 
% binning by Ken Smith @ MATLAB Central
%{
c=reshape(sum(sum(reshape(c,[m n m n]),1),3),[n n]);

or more 

c=reshape(c,[m n m n]);
c=sum(c,1);
c=sum(c,3);
c=reshape(c,[n n]); 
%}

