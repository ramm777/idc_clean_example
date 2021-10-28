function coords =  cellCoords(G, cell)
% 
% SYNOPSIS:
%   coords =  cellCoords(G, cell)
%
% DESCRIPTION: 
%   - Extract coordinate of a specific cell from grid G. Based on extractSubgrid.m using only part of it to save comput time/memory 
%   - Doesn't work in mirrowed intracells due to penetration

% PARAMETERS:
%   G       - Grid structure as described by grid_structure.
%   cells   - Cells to investigate
% 
% RETURNS:
%   coords  - Cell coords
% 
% EXAMPLE:
%
% SEE ALSO:
% 
% AUTHOR:
%  Aidan Kubeyev

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


    c = cell;
    ix = @(p,i) mcolon(double(p(i)), double(p(i+1)) - 1);

    
    cells      = false([G.cells.num + 1, 1]); % Sort the subset cells, c, in ascending order.  Add one for outside.
    cells(c+1) = true;
    c = find(cells)-1;                        % replace c by sorted version of c  

        
    faces = false([G.faces.num, 1]);          % Extract faces connected to 'c'. 
    faces(G.cells.faces(ix(G.cells.facePos, c))) = true;

    
    nodes = false([G.nodes.num, 1]);          % Extract nodes connected to 'faces' (in cell subset). 
    ff    = find(faces);
    nodes(G.faces.nodes(ix(G.faces.nodePos, ff))) = true;


    coords    = G.nodes.coords(nodes, :);

    
end

