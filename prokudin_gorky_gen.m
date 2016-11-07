function [] = prokudin_gorky_gen();
%% This function takes three images (files must be loaked


%% Load image (needs customization)
ims = 120:122;

data = [];
for i = 1:3
    im_name = ['IMG_0' num2str(ims(i)) '.jpg'];
    img = imread(im_name);
    data = cat(3, data, squeeze(img(:,:,i)));
end

figure; imshow(data) 
axis image off

%% Realign image (assumes x-y shift, with no rotation)

buffer = 50; % amount of size
base = squeeze(data(buffer:(end-buffer),buffer:(end-buffer),1));

%% define a shift function
shift_xy = @ (x,x_shift, y_shift) x((buffer-x_shift):(end-buffer-x_shift), (buffer-y_shift):(end-buffer-y_shift));
sum_squares = @(x,y) sum(sum((x-y).^2));
% flatten = @(x) x(:);
data_aligned = base;
for j = 2:3
    test = squeeze(data(:,:,j));
    %% GRID SEARCH!!!
    [x_shifts, y_shifts] = meshgrid(-(buffer-1):(buffer-1), -(buffer-1):(buffer-1));
    x_shifts = x_shifts(:);
    y_shifts = y_shifts(:);
    ssq_x = zeros(1,length(x_shifts));
    for i = 1:length(x_shifts)
        test_shift = shift_xy(test, x_shifts(i), y_shifts(i));
        ssq_x(i) = sum_squares(test_shift, base);
    end
    [~,opt_shift] = min(ssq_x);
    
    data_aligned = cat(3, data_aligned, shift_xy(test, x_shifts(opt_shift), y_shifts(opt_shift)));
end

figure; imshow(data_aligned) 
axis image off


end