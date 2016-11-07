function [] = generate_ishihara(image_name, varargin)
%% This program creates an Ishihara test by generating a drawing a dots with a random
% distribution of sizes on a mesh grid. The color of these circles
% is dictated by the strucutre of an image that is loaded in the function
% call by the input 'image_name' This image is binarized 
%
% The image in binaraized by a simple threshold, and the dots that fall
% inside the areas of the binarized image above the threshold are given a
% color specified by the input from varargin
%
%
% ***

%% Define Color Settings from varargin
if isempty(varargin)
    on_color = [0.8 0.2 0]; % mean value of each color channel
    off_color = [0.2 0.8 0];
    color_variations = [0.1 0.1 0.1]; % standard deviation of each color channel
elseif length(varargin) == 2
    on_color = varargin{1};
    off_color = varargin{2};
    color_variations = [0.1 0.1 0.1];
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
x_mask = x_mask(:);


%% flatten the grid for application of boolean arrays
x_flat = x(:);
y_flat = y(:);
y_mask = y_mask(:);





%% Define size and frequency parameters of dots: (these settings are arbitrary)
n = 5;                          % number of dot sizes
dot_sizes = 0.015*1./sqrt(1+1:n); % sizes of the various random dots
dot_probs = cumsum(1./dot_sizes); % probability that each of the above sizes occur (these must be same dimension)
dot_probs = dot_probs/max(dot_probs);

%% Define a vector of allowed x,y coordinates
filled_space = ~get_interior(x,y, (xmax-xmin)/2, (ymax-ymin)/2, (xmax-xmin)/2);
is_overlap = false; % parameter used to determine if the previous dot was overlapping

%% Loop through and add patches
fill_fraction = sum(filled_space(:))/length(filled_space(:));
fill_thresh = 0.55; % Sets the limit on how much space can be filled in
i = 0; max_iterations = 5000;
%figure; hold on;
while (fill_fraction < fill_thresh) & (i < max_iterations)
    i = i+1;
    %% Get Radius
    event_radius = rand; % rand from 0 to 1
    rad_prob = min(dot_probs(dot_probs>event_radius));
    radius = dot_sizes(find(dot_probs == rad_prob,1)); % Randomly select dot radius
    spacing = 0.95; % buffer to prevent dots from touching
    
    % Get X and Y Coordinates of the center
    if ~is_overlap % eliminates steps if the previous dot was blocked for overlap (is_clear must be initialized FALSE)
        is_allowed = ~filled_space(:);
    
        allowed_x = x_flat(is_allowed);
        allowed_y = y_flat(is_allowed);
    
        allowed_mask_y = y_mask(is_allowed);
        allowed_mask_x = x_mask(is_allowed);
    end
    
    center_index = randperm(length(allowed_x),1); % randomly select index for center coords
    center_x = allowed_x(center_index);
    center_y = allowed_y(center_index);
    
    %% Determine if dot has overlap with existing dots
    is_overlap = sum(get_interior(x,y, center_x, center_y, radius) & filled_space) == 0;

    if is_overlap
        circ_x = spacing*radius*cos(t)+center_x;
        circ_y = spacing*radius*sin(t)+center_y;
        
        
        if allowed_mask_x(center_index) & allowed_mask_y(center_index)
            color = on_color + color_variations.*randn(1,3);
            color = min([1 1 1],max([0 0 0], color)); % truncate triplets beyond [0 1].
            patch(circ_x, circ_y, color,'EdgeColor', 'none', 'FaceAlpha', dot_alpha)
        else
            color = off_color + color_variations.*randn(1,3);
            color = min([1 1 1],max([0 0 0], color));
            patch(circ_x, circ_y, color,'EdgeColor', 'none', 'FaceAlpha', dot_alpha)
        end

        filled_space = filled_space | get_interior(x,y, center_x, center_y, radius);
        fill_fraction = sum(filled_space(:))/length(filled_space(:));
    end


end


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
