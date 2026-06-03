% Performs Weisskoff analysis on EPI images.
% [N, data_all] = Weisskoff_analysis(EPIdata, slice, R, detrend_order, ...
%   discard, center)
%
% It take a EPIdata, a matrix of size n by n by slice (top to bottom) by
% time (x,y,z,t) and performs Weisskoff analysis with the specified
% parameters. Pass empty var, i.e. [], if you would like to use the default
% value.
%
% slice: the slice to analize. Default is number of slice / 2.
%
% R: size of the square to use, as defined by the number of voxels per
% side. Default is 80% of the largest square that will fit the phantom.
%
% detrend_order: the order of polynomial detrending to use. Default is
% none.
%
% discard: the number of volumes at the beginning to discard. Default value
% is 0.
%
% centers: manually defines the center of the image. Default is the center
% of the phantom.
%
% Output is the Radius of Decorrelation (RDC), and the CV vs N curve

function [RDC, N, CV] = Weisskoff_simple(EPIdata, slice, R,...
    detrend_order, discard, centers)

figure;

temp = EPIdata(:,:,slice,1);  % take the first volume
% fit a circle to the signal region with min radius of 6 up to the entire
% slices
[centers_default, radii] = imfindcircles(temp,[6 min(size(temp))],...
    'ObjectPolarity','bright');

% overwrite slice if specified
if isempty(slice)
    slice = size(EPIdata, 3)/2;
end

% find out what size square to go up to or use the specified size unless
% specified
if isempty(R)
    % find the largest square that will fit, and do 80% of that
    R = round((radii*0.8)/sqrt(2))*2;
end

% overwrite center of ROI if specified
if isempty(centers)
    squared_centers = round(centers_default);
else
    squared_centers = centers;
end

% overwrite volume discard if specified
if ~isempty(discard)
    EPIdata = EPIdata(:,:,:,discard+1:end);
end

% generate best case curve
Nideal  = (1:R);
squared_mask = false(size(temp));   % start a blank mask

% it is not the most efficient computing growing the square in a for
% loop... but it does show the square ROI growing and show the progress
for kk = 1:R
    ro2 = fix(kk/2);
    x1 = squared_centers(2) - ro2;
    x2 = x1 + kk - 1;
    y1 = squared_centers(1) - ro2;
    y2 = y1 + kk - 1;
    squared_mask(x1:x2,y1:y2)=true;
    for ll = 1:size(EPIdata,4)
        temp = EPIdata(:,:,slice,ll);
        data(ll) = mean(temp(squared_mask));
    end

    if isempty(detrend_order)
        data2 = data;
    else
        data2  = detrend(data,detrend_order);
    end

    CV(kk) = std(data2)/mean(data)*100;  % in percent
    N(kk)  = kk;

    subplot(2,1,1);
    imshow(mat2gray(temp))
    hold on
    viscircles(squared_centers,radii);
    visboundaries(squared_mask, 'color', 'b');
    hold off

    subplot(2,1,2);
    loglog(Nideal,CV(1)./Nideal,'g-');
    hold on
    loglog(N,CV, 'bo');
    hold off
    xlabel('ROI width (voxels)')
    ylabel('Coefficient of Variation (%)');

    pause(1/12); % pause of update the graph
end

RDC = CV(1)/CV(kk);
