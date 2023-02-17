function plot_3d_data_by_1d_curve(gridpos, value, curve_dim, idxX, idxY, idxZ)
% PLOT_3D_DATA_BY_1D_CURVE plots the 1 dimensional slice (curve) of the 3
%    dimensional data.
%
%    plot_3d_data_by_1d_curve(gridpos, value, curve_dim, idxX, idxY, idxZ)
%    plots the 1d slice of 3d data pair (gridpos, value) according to the
%    dimension of the slice curve_dim and indices of other two dimensions.
%    For example, when curve_dim=1, idxX is [], idxY and idxZ are user
%    provided with default value the middle position of that directions.
%
%    NOTE: gridpos should be a cell with 3 elements containing position of 
%    grid in each directions. value can be a three dimensional tensor or 
%    an one dimension vector which was reshaped by the three dimensional 
%    tensor.  
%
%    See also plot_3d_data_by_2d_surface.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


numGrid = zeros(1, 3);
for d = 1 : 3
    numGrid(d) = length(gridpos{d});
end

if nargin == 3
    if curve_dim == 1
        idxY = floor(numGrid(2) / 2);
        idxZ = floor(numGrid(3) / 2);
    elseif curve_dim == 2
        idxX = floor(numGrid(1) / 2);
        idxZ = floor(numGrid(3) / 2);
    elseif curve_dim == 3
        idxX = floor(numGrid(1) / 2);
        idxY = floor(numGrid(2) / 2);
    end
end 

if size(value, 2) == 1
    val3D = reshape(value, numGrid);
else
    val3D = value;
end

if curve_dim == 1
    grid = gridpos{1};
    val = squeeze(val3D(:, idxY, idxZ));
    plot(grid, val, 'b-o');
    title("x direction with y index: " + num2str(idxY) + "/" + num2str(numGrid(2)) ...
        + ", z index: " + num2str(idxZ) + "/" + num2str(numGrid(3)));
elseif curve_dim == 2
    grid = gridpos{2};
    val = squeeze(val3D(idxX, :, idxZ));
    plot(grid, val, 'b-o');
    title("y direction with x index: " + num2str(idxX) + "/" + num2str(numGrid(1)) ...
        + ", z index: " + num2str(idxZ) + "/" + num2str(numGrid(3)));
elseif curve_dim == 3
    grid = gridpos{3};
    val = squeeze(val3D(idxX, idxY, :));
    plot(grid, val, 'b-o');
    title("z direction with x index: " + num2str(idxX) + "/" + num2str(numGrid(1)) ...
        + ", y index: " + num2str(idxY) + "/" + num2str(numGrid(2)));
else
    error('wrong input of dimension of curve');
end
box on

end