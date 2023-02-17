function plot_3d_data_by_2d_surface(gridpos, value, slice_direction, slice_index)
% PLOT_3D_DATA_BY_2D_SURFACE plots the 2 dimensional slice (surface) of the 
%    3 dimensional data.
%
%    plot_3d_data_by_2d_surface(gridpos, value, slice_direction, slice_index)
%    plots the 2d slice of 3d data pair (gridpos, value) according to the
%    direction of the slice slice_direction (1, 2, 3) and index of the 
%    slice slice_index. slice_index is user-provided with default value the 
%    middle position of that directions.
%
%    NOTE: gridpos should be a cell with 3 elements containing position of 
%    grid in each directions. value can be a three dimensional tensor or 
%    an one dimension vector which was reshaped by the three dimensional 
%    tensor.  
%
%    See also plot_3d_data_by_1d_curve.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


numGrid = zeros(1, 3);
for d = 1 : 3
    numGrid(d) = length(gridpos{d});
end

if nargin == 3
    slice_index = floor(numGrid(slice_direction) / 2);
end 

if size(value, 2) == 1
    val3D = reshape(value, numGrid);
else
    val3D = value;
end

gridx = gridpos{1};
gridy = gridpos{2};
gridz = gridpos{3};

if slice_direction == 3
    [X, Y] = ndgrid(gridx, gridy);
    Z = squeeze(val3D(:, :, slice_index));
    surf(X, Y, Z, 'EdgeColor','none', 'FaceColor', 'interp');
    xlabel('x');
    ylabel('y');
    title("z index: " + num2str(slice_index) + " / " + num2str(numGrid(3)));
elseif slice_direction == 2
    [Z, X] = ndgrid(gridz, gridx);
    Y = squeeze(val3D(:, slice_index, :));
    surf(Z, X, Y, 'EdgeColor','none', 'FaceColor', 'interp');
    xlabel('z');
    ylabel('x');
    title("y index: " + num2str(slice_index) + " / " + num2str(numGrid(2)));
elseif slice_direction == 1
    [Y, Z] = ndgrid(gridy, gridz);
    X = squeeze(val3D(slice_index, :, :));
    surf(Y, Z, X, 'EdgeColor', 'none', 'FaceColor', 'interp');
    xlabel('y');
    ylabel('z');
    title("x index: " + num2str(slice_index) + " / " + num2str(numGrid(1)));    
else
    error('wrong input of slice direction');
end

colormap hsv

end