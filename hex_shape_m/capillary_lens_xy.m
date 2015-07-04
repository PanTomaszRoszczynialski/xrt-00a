function [xi,yi,xci,yci] = capillary_lens_xy(nY_bundle,bundlespacing,capillary_diameter,channel_diameter,nx_capillary,sigma)

    %
    xci=[];
    yci=[];
    nypol_bundle=(nY_bundle-1)/2;


    index=0;
    for ix=-2*nypol_bundle:2*nypol_bundle
        ix
        for iy=-2*nypol_bundle:2*nypol_bundle       
     
            x0 = sqrt(3)/2*bundlespacing*iy+0*sigma* ... 
                (rand(1,1)-0.5)* ... 
                (capillary_diameter-channel_diameter);

            y0 = bundlespacing*ix + bundlespacing/2*iy+ ... 
                0*sigma*(rand(1,1)-0.5)* ... 
                (capillary_diameter-channel_diameter);
            
            
            [in_lens] = isinhexagon(y0,x0,bundlespacing*nypol_bundle);
           
            
            if in_lens
                [xci0,yci0]= capillary_bundle_xy(x0,y0,... 
                    capillary_diameter,channel_diameter,... 
                    nx_capillary,sigma);
                index=index+1;
                xi(index)=x0;
                yi(index)=y0;
                xci=[xci xci0];
                yci=[yci yci0];
            end
        end
    end
end
