function [] = generate_ishihara(image_name, varargin)
%% This program creates an Ishihara test for Color Blindness
% This generates a drawing of randomly placed colored dots with
% distribution of sizes and colors on a mesh grid. The color of these circles
% is dictated by the strucutre of an image that is loaded in the function
% call by the input 'image_name' This image is binarized 
%
% The image in binaraized by a simple threshold, and the dots that fall
% inside the areas of the binarized image above the threshold are given a
% color specified by the input from varargin
%
%
% *** Example:
% generate_ishihara('TemplateImage.png',[0.8 0.2 0.7], [0 0.4 0.8], [0.05 0.05 0.05])
% will generate an image determined by the file 'TemplateImage, with areas
% of high intensity having colors around the RGB Triplet [0.8 0.2 0.7], 
% areas of low intensity being colored the RGB triple [0 0.4 0.8]
% The three channels of these colors are given a random gausian noise whose 
% variance is determined by the final imput [0.05 0.05 0.05]
%
% Because the radii of the dots is random, it is computationally expensive
% to determin where sphere may be placed.
% This program has sacrificed speed for the aesthetic of a classic
% Ishihara-test
%
% Ben Lucas - 2/29/16
tic
%% Define Color Settings from varargin
if isempty(varargin)
    on_color = [0.8 0.2 0]; % mean value of each color channel
    off_color = [0.2 0.8 0];
    color_variations = [0.01 0.01 0.01]; % standard deviation of each color channel
elseif length(varargin) == 2
    on_color = varargin{1};
    off_color = varargin{2};
    color_variations = [0.01 0.01 0.01];
else
    on_color = varargin{1};
    off_color = varargin{2};
    color_variations = varargin{3};
end
dot_alpha = 1; % transparency of dots

    
%% Define a meshgrid of points
xmin = 0; xmax = 1; ymin = 0; ymax = 1;
mesh_spacing = 0.001;
[x,y] = meshgrid(xmin:mesh_spacing:xmax, ymin:mesh_spacing:ymax);
t = 0:0.01:(2*pi);  % grid spacing for defining parametric plots of graphics


%% Load Image for Mask
char_image = imread(image_name);
if size(char_image)>2
    char_image = sum(char_image,3); %Sum intensity for Binarization
end
mask_thresh = 1; % for binarization of mask image
char_image = fliplr(char_image'); % Dont really know why but images needed rotating...

[x_mask, y_mask] = make_mesh_mask(char_image, x, y, mask_thresh);


%% flatten the grid for application of boolean arrays
x_flat = x(:);
x_mask = x_mask(:);

y_flat = y(:);
y_mask = y_mask(:);





%% Define distribution of dots parameters of dots: (these settings are arbitrary)
num_dots = 8000;
min_separation = 0.01;  % the smallest distance between two centers;
dot_spacing = 0.45;     % under half distance
centers= rand(num_dots,2);
dists = squareform(pdist(centers));

good_pts = true(1,num_dots);
for i = 1:num_dots;
    if good_pts(i)
        too_close = dists(i,:) < min_separation;
        too_close(i) = false;
        good_pts = good_pts & ~too_close;
    end
end
centers = centers(good_pts,:);

% take the radius of each point to be half the distance to the nearest pt
dists = squareform(pdist(centers));
dists(dists==0) = 10*(xmax-xmin);
radius = dot_spacing*min(dists,[],1);
 

%% Plot Dots
figure; hold on;
for i = 1:length(radius)
    circ_x = radius(i)*cos(t)+centers(i,1);
    circ_y = radius(i)*sin(t)+centers(i,2);

    center_x_inds = find(abs(x_flat-centers(i,1)) == min(abs(x_flat-centers(i,1))));
    center_y_inds = find(abs(y_flat-centers(i,2)) == min(abs(y_flat-centers(i,2))));
    center_inds = intersect(center_x_inds, center_y_inds);
    center_index = center_inds(1);

    if x_mask(center_index) & y_mask(center_index)
        color = on_color + color_variations.*randn(1,3);
        color = min([1 1 1],max([0 0 0], color)); % truncate triplets beyond [0 1].
        patch(circ_x, circ_y, color,'EdgeColor', 'none', 'FaceAlpha', dot_alpha)
    else
        color = off_color + color_variations.*randn(1,3);
        color = min([1 1 1],max([0 0 0], color));
        patch(circ_x, circ_y, color,'EdgeColor', 'none', 'FaceAlpha', dot_alpha)
    end

end

toc
end






%% Helper Functions


function [x_mask, y_mask] = make_mesh_mask(char_image, x, y, mask_thresh);
%% This function takes the x and y matricies from a grid and creates boolean matricies
% that are true when isometric points on char_image are above mask_thresh
x_mask = false(size(x));
y_mask = false(size(y));
% for each point in the mesh find the value of the closest index in charimage
for i = 1:size(x)
    for j = 1:size(y)
        row_ind = round(x(i,j)*(size(char_image,1)-1))+1;
        col_ind = round(y(i,j)*(size(char_image,2)-1))+1;
        x_mask(i,j) = (char_image(row_ind, col_ind) > mask_thresh);
        y_mask(i,j) = (char_image(row_ind, col_ind) > mask_thresh);
    end
end


end


function [in_circle] = get_interior(x,y, center_x, center_y, radius)
% function takes a mesh grid (x,y) as well as the center and radius of a
% circle and returns a boolean vector in_circle that is true is the
% corresponding to the points of x and y are on the interior of the
% circle prescirbed by center and radius
in_circle = false(size(x));
center_dists = sqrt((x-repmat(center_x,size(x))).^2 + (y-repmat(center_y,size(y))).^2);
in_circle = center_dists < radius;
end
