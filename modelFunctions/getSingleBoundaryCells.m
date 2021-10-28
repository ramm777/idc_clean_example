function bc = getSingleBoundaryCells(G, bc, globalBoundary)


% SYNOPSIS:
%   bc = getSingleBoundaryCells(G, bc, 'Left'); 
%
% DESCRIPTION: Find and insert to bc structure, a global boundary cells of the selected
%              boundary from 1) 'Left' 2) 'Right', 3) 'Lower' 4) 'Upper' 
%                                   
% 
% PARAMETERS:
%   G               - Grid structure
%   bc              - boundary condition structure usually created by findBoundaryFacesAndCreateBC()
%   globalBoundary  - 1) 'Left' 2) 'Right', 3) 'Lower' 4) 'Upper' 
% 
% RETURNS:
%   bc       - updated bc structure, where the required boundary cells are inserted
%   
%
% EXAMPLE:
%
% SEE ALSO:
%
% AUTHOR:
%   Aidan Kubeyev 
%
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

    globalBoundaries = {'Left', 'Right', 'Lower', 'Upper'}; 
    boundary_id = find(contains(globalBoundaries, globalBoundary)); 
        
    
    G.cells.faces = [rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .', G.cells.faces]; % Restore G.cells.faces mapping
    [~, ind] = ismember(bc{boundary_id}.face, G.cells.faces(:, 2));                                  % Find which cells have faces that belongs to outer bc{1}
    bc{boundary_id}.cells = [G.cells.faces(ind, 1)];                                                 % Cells that is on outer bc 
    %G.cells.faces(:, 1) = []; % Delete 1st raw or this may cause error later


    
end

