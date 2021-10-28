function h = plotCellNumbers3D(G, varargin)
%
% SYNOPSIS:
%   h = plotNodeNumbers3D(G)
%
% DESCRIPTION: 
%   Plots cells' numbers of the 3D grid
%    
% PARAMETERS:
%   G       - Grid structure as described by grid_structure.
%   cells   - cells, whos numbers to be displayed. 
%   
%   
% RETURNS:
%   h - handle to plot
%   
%
% EXAMPLE:
%
% SEE ALSO:
% 
% AUTHOR:
%   Aidan Kubeyev  

%{

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}


if mod(nargin, 2)
    c = 1:g.cells.num;
else
    c = varargin{1};
    varargin = varargin(2:end);
end

if G.griddim == 3 

    markertext = G.cells.indexMap; 

    h=text(G.cells.centroids(c, 1), G.cells.centroids(c, 2), G.cells.centroids(c, 3), ...
        num2str(markertext(c)), 'FontSize', 9, 'Color', 'blue');

else
    
    error('not suitable for 2D, use plotCellNumbers function'); 
    
end

end

