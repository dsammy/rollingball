%This script is intended for background correction and phase fraction measurements in greyscale BSE images. 
%It is currently written for analysis of two-phase materials whose phase of interest is the darker one.
%The first step for a successful segmentation is a good sample preparation.
%In the input parameters section it is possible to define the file names. 
%For the sake of simplicity it is recommended to use the file name system as follows image1.tif, image2.tif etc.
%------------------------------------------------------------------
%Written by Mariana Rodrigues (mariana.rodrigues@alumni.ubc.ca)
%V1.0
%2019.10.11.

clc;
close all;

%----------------------------------------------------------------------
%INPUT
%----------------------------------------------------------------------

% Here define the path where the images are located
myFolder = 'C:\Users\Mariana\Matlab'; 

% Read images1.tif through m.tif from "myFolder" directory
for k = 2:2
    tifFilename = sprintf('Ti5553_650C_14min%d.tif', k);
    fullFileName = fullfile(myFolder, tifFilename);
    if exist(fullFileName, 'file')
        im_original = imread(fullFileName);
        
        % If image is stored in RGB, convert it to 8-bit
        if size(im_original,3) == 4
            im_original(:,:,4) = [];
            im_original = rgb2gray(im_original);
        end
        
%----------------------------------------------------------------------
%PRE-PROCESSING AND BACKGROUND CORRECTION
%----------------------------------------------------------------------
        
        % Reading the size of the image in pixels
        [a,b] = size(im_original);
        
        % Image normalization using linear stretching
        im_normalz = imadjust(im_original);
        
        % Apply median filter using 3x3 pixels to remove the noise before applying rolling ball
        im_filtered = medfilt2(im_normalz,[3 3]);
        
        % Checking different ball radii
        for radius = 1:20
            SE = strel('sphere',radius);
            
            % Calculating the image background
            im_bckg1 = imclose(im_filtered,SE);
            x = 255-im_bckg1;
            
            % Correcting image background
            im_closing = imadd(im_filtered,x);
            im_sharpen = imsharpen(im_closing);
            
            % Calculating the sum and area fraction of pixels whose intensities are equal to 255
            Sum255 = sum(im_sharpen(:) == 255);
            Areafrac255(radius) = Sum255/(a*b);
            
            % Calculating the first derivative with respect to ball radius
            Derivative = diff(Areafrac255);
            
            % Selecting the best ball radius
            Flat = abs(Derivative) > 0.01; %results in a logical vector where 1 is assigned to all the derivative values above 0.01
            best_radius = find(Flat,1,'last') + 1;% find the corresponding ball radius
            sprintf('Ball radius = %d, Area fraction of pixels with 255 intensity = %d', radius, Areafrac255(radius))
        end
        
        % Returning the best ball radius and applying the background correction using this radius
        sprintf('The best rolling ball radius is = %d pixels', best_radius)
        SE_best = strel('sphere',best_radius);
        im_bckgbest = imclose(im_filtered,SE_best);
        x_best = 255-im_bckgbest;
        im_closingbest = imadd(im_filtered,x_best);
        im_sharpenbest = imsharpen(im_closingbest);
        
        % Plotting the 3-D structural element
        figure;
        isosurface(SE_best.Neighborhood);
        ax = gca;
        axis([0 25 0 25 0 25])
        ax.FontSize = 16;
        ax.XTick = 0:5:25;
        ax.YTick = 0:5:25;
        ax.ZTick = 0:5:25;
        xlabel('X (pixels)','FontSize',18)
        ylabel('Y (pixels)','FontSize',18)
        zlh = zlabel('Z (pixels)','FontSize',18);
        set(zlh, 'position', get(zlh,'position')-[0.7,0,0]);
        xh = get(gca,'xlabel'); 
        set(xh, 'Units', 'Normalized')
        pos = get(xh, 'Position');
        set(xh, 'Position',pos.*[1,0.6,1],'Rotation',15)
        yh = get(gca,'ylabel'); 
        set(yh, 'Units', 'Normalized')
        pos = get(yh, 'Position');
        set(yh, 'Position',pos.*[1,0.3,1],'Rotation',-25)
        colormap(gray)
        set(gcf,'color','w');
        grid on
        
        % Plotting the 3-D grey profile
        x = 0:size(im_bckgbest,2)-1;
        y = 0:size(im_bckgbest,1)-1;
        [X,Y] = meshgrid(x,y);
        figure;
        mesh(X, Y, im_bckgbest);
        ax = gca;
        axis([0 1024 0 768 0 255])
        ax.FontSize = 17;
        ax.XTick = 0:200:1024;
        ax.YTick = 0:200:768;
        ax.ZTick = 0:50:255;
        xlabel('X (pixels)','FontSize',19)
        ylabel('Y (pixels)','FontSize',19)
        zlabel('Grey intensity','FontSize',19)
        xh = get(gca,'xlabel');
        set(xh, 'Units', 'Normalized')
        pos = get(xh, 'Position');
        set(xh, 'Position',pos.*[1,0.6,1],'Rotation',15)
        yh = get(gca,'ylabel');
        set(yh, 'Units', 'Normalized')
        pos = get(yh, 'Position');
        set(yh, 'Position',pos.*[1,0.3,1],'Rotation',-25)
        set(gca,'YDir','reverse')
        set(gcf,'color',[0.95 0.95 0.95]);
        colormap(gray)
        
%----------------------------------------------------------------------
%IMAGE THRESHOLDING
%----------------------------------------------------------------------
        
        % Getting the threshold grey value
        prompt = 'Is the threshold from EBSD map known? [Y/N]: \n';
        str = input(prompt,'s');
        if str == 'Y'
            prompt2 = 'What is the threshold value? [integer from 1 to 255]\n';
            threshold = input(prompt2);
            
% ------- Analysing BSE image which is correspondent to the EBSD image
        elseif str == 'N'
            disp('Calculating the threshold...')
            
            % Here define the path where the image is located
            myFolder = 'C:\Users\Mariana\Matlab';
            for w = 1:1 %Here the last number must be chanegd according to the total number of images
                tifFilename2 = sprintf('Ti5553_650C_3h_ebsd_map0%d.tif', w);
                fullFileName2 = fullfile(myFolder, tifFilename2);
                if exist(fullFileName2, 'file')
                    im_EBSD1 = imread(fullFileName2);
                    if size(im_EBSD1,3)==4 
                        im_EBSD1(:,:,4) = []; 
                        im_EBSD1 = rgb2gray(im_EBSD1); 
                    end
                    % Image information
                    minGL = min(im_EBSD1(:));
                    maxGL = max(im_EBSD1(:));
                    meanGL = mean(im_EBSD1(:));
                    imEBSD_normalz = imadjust(im_EBSD1); 
                    [o,p] = size(imEBSD_normalz); 
                    message = sprintf('The min grey level = %d\nThe max grey level = %d\nThe mean grey level = %d\nThe image size is %dx%d', ...
                              minGL, maxGL, meanGL, p, o);
                    uiwait(msgbox(message, 'Info'));
                    
                    % Checking best rolling ball radius for the BSE image corresponding to the EBSD map
                    for radius_EBSD = 1:20
                        SE_EBSD = strel('sphere',radius_EBSD);
                        imEBSD_filtered = medfilt2(imEBSD_normalz,[3 3]);
                        
                        imEBSD_bckg = imclose(imEBSD_filtered,SE_EBSD);
                        h = 255-imEBSD_bckg;
                        imEBSD_closed = imadd(imEBSD_filtered,h);
                        imEBSD_sharpen = imsharpen(imEBSD_closed);
                        
                        imEBSD_Sum255 = sum(imEBSD_sharpen(:) == 255);
                        imEBSD_Areafrac255(radius_EBSD) = imEBSD_Sum255/(o*p);
                        Derivative_EBSD = diff(imEBSD_Areafrac255);
                        sprintf('Ball radius = %d, Area frac. of 255 intensity pixels. = %d', radius_EBSD, imEBSD_Areafrac255(radius_EBSD));
                        Flat_EBSD = abs(Derivative_EBSD) > 0.01;
                        bestradius_EBSD =  find(Flat_EBSD,1,'last') + 1;
                    end
                    
                    sprintf('The best rolling ball radius is = %d pixels', bestradius_EBSD)
                    SE_EBSD2 = strel('sphere',bestradius_EBSD);
                    
                    imEBSD_bckg2 = imclose(imEBSD_filtered,SE_EBSD2);
                    h2 = 255-imEBSD_bckg2;
                    imEBSD_closed2 = imadd(imEBSD_filtered,h2);
                    imEBSD_sharpen2 = imsharpen(imEBSD_closed2);
                    
                    % Enter fraction obtained from EBSD analysis
                    prompt = {'Enter the EBSD fraction in %:'};
                    title = 'EBSD fraction';
                    dims = [1 40];
                    definput = {'22.5'};
                    answer = inputdlg(prompt,title,dims,definput);
                    user_fracEBSD = str2num(answer{:});
                    Sumpixels = (user_fracEBSD/100)*(o*p);
                    
                    % Calculating the cumulative distribution function and the threshold grey value
                    hist_EBSD2 = imhist(imEBSD_sharpen2);
                    cdf_EBSD2 = cumsum(hist_EBSD2); %this array goes from 1 to 256
                    xIndex = find(cdf_EBSD2 > Sumpixels,1,'first'); %finds the first 1 value above Sumpixels
                    thresh_index = xIndex-1; %minus 1 to change to 0-255 grey level
                    sprintf('File Ti5553_650C_3h_ebsd_map: --> The threshold grey value is = %d', w, thresh_index)
                    
                else
                    warningMessage = sprintf('Image file does not exist!\n%s', fullFileName2);
                    uiwait(warndlg(warningMessage));
                end
            end
            message = sprintf('The threshold value is %d\n', thresh_index);
            uiwait(msgbox(message));
            threshold = thresh_index;
            
% ------- Finish analysing BSE image which is correspondent to the EBSD image and back to the main cycle         
        else
            f = msgbox('Invalid answer!','Error','error');
            break
        end
        
        % Checking each pixel in the image, finding the pixels below and above the threshold and changing the image to binary
        im_binary = zeros(a,b);
        for i = 1:a
            for j = 1:b
                if(im_sharpenbest(i,j) < threshold)
                    im_binary(i,j) = 0;
                else
                    im_binary(i,j) = 1;
                end
            end
        end
        
%----------------------------------------------------------------------
%PHASE FRACTION CALCULATION
%----------------------------------------------------------------------
        
        % Calculating the fraction of black pixels
        sum_black = sum(im_binary(:) == 0);
        alphafrac = sum_black/(a*b)*100;
        sprintf('File Ti5553_650C_14min%d: --> The best rolling ball radius is = %d pixels, The alpha fraction is = %d', k, best_radius, alphafrac)
        
        % Calculating the total phase boundary length (perimeter)
        invertthesh = imcomplement(im_binary);
        notinvertthesh = true(size(invertthesh)+2);
        notinvertthesh(2:end-1, 2:end-1) = ~invertthesh;
        topedges = invertthesh & notinvertthesh(1:end-2, 2:end-1);
        leftedges = invertthesh & notinvertthesh(2:end-1, 1:end-2);
        bottomedges = invertthesh & notinvertthesh(3:end, 2:end-1);
        rightedges = invertthesh & notinvertthesh(2:end-1, 3:end);
        perim = sum(topedges(:)) + sum(leftedges(:)) + ...
            + sum(bottomedges(:)) + sum(rightedges(:));
        
        % Printing the total boundary length
        sprintf('The total boundary length is %d\n', perim)
        
        % Saving the binary image
        imwrite(im_binary,sprintf('Binary_image_Ti5553_650C_14min%d.tif', k));
        
    else
        warningMessage = sprintf('Image file does not exist!\n%s', fullFileName);
        uiwait(warndlg(warningMessage));
    end
end
