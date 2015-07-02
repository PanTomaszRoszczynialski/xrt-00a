clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Size of S in pixels - exit surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nS=8192;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Capillary parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Outer diameter (touching)
capillary_diameter=0.0025;
% Inner diameter (x-ray microsource)
channel_diameter=capillary_diameter*0.5;

% number of capillaries in bundle along horizontal direction (must be odd)
nx_capillary=21 ;
if ~odd(nx_capillary)
    error('Number of capillaries must be odd!');
end


% randomness
sigma_position=0.1;
sigma_intensity=0.2;
reject_ratio=0.002;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bundle parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spacing between capilaries (touching)
minbundlespacing=(nx_capillary-1)*capillary_diameter*sqrt(3)/2+capillary_diameter;

% or pherhaps nontouching
bundlespacing=minbundlespacing;


% number of bundles along vertical direction (must be odd)
ny_bundle=21;
if ~odd(ny_bundle)
    error('Number of bundles must be odd!');
end



[Sxbundle,Sybundle,Sx,Sy]=capillary_lens_xy(ny_bundle,bundlespacing,capillary_diameter,channel_diameter,nx_capillary,sigma_position);



% Reject random capiilaries

r=rand(size(Sx));
gdzie=r>reject_ratio;
Sx=Sx(gdzie);
Sy=Sy(gdzie);




alpha=0*pi/180;
Sxrot=Sx*cos(alpha)-Sy*sin(alpha);
Syrot=Sx*sin(alpha)+Sy*cos(alpha);
Sx=Sxrot;
Sy=Syrot;

maxxy=max([max2(Sx) max2(Sy)]);
minxy=min([min2(Sx) min2(Sy)]);



Syscaled=(Sy-minxy)*(nS-1)/(maxxy-minxy)+1;
Sxscaled=(Sx-minxy)*(nS-1)/(maxxy-minxy)+1;

Sxbundlescaled=(Sxbundle-minxy)*(nS-1)/(maxxy-minxy)+1;
Sybundlescaled=(Sybundle-minxy)*(nS-1)/(maxxy-minxy)+1;

Ssparse=sparse(round(Syscaled),round(Sxscaled),ones(size(Sxscaled)),nS,nS);
Sfull=full(Ssparse);


% Add intensity fluctuations

Sfull=full(Ssparse);
maksimum=max2(Sfull);
szum=-maksimum*abs(randn(size(Sfull)))*sigma_intensity;
szum=szum;
Sfull=Sfull+szum.*(Sfull>0);



%Ssparsebundle=sparse(round(Sybundlescaled),round(Sxbundlescaled),ones(size(Sybundlescaled)),4096*2,4096*2);
%Sfullbundle=full(Ssparsebundle);

%f=fft2(Sfull);
%fbundle=fft2(Sfullbundle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convolve with gaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rcap=0.6;
nR=20;


x=linspace(-1,1,nR);
[X,Y]=meshgrid(x,x);
g=(X.^2+Y.^2)<Rcap;
g=exp(-(X.^2+Y.^2)/0.5^2);
Sfullwithcapillaryshape=conv_fft2(Sfull,g,'same');
S=bin2d(Sfullwithcapillaryshape,2,4096);

%x=linspace(-1,1,4096);
%[X,Y]=meshgrid(x,x);
%G=gauss2d(X,Y,1/0.4);
%S=S.*G;

save S_discrete_after_convolution_gaussian_d15 S