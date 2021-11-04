function faces = findCellFaces(G, cell)
%
% SYNOPSIS:
%   faces = findCellFaces(G, cell)
%
% DESCRIPTION:
%   Function finds all faces of a specified cell(s) 
%
% PARAMETERS:
%   G          - Grid structure as described by grid_structure. 
%   cells      - cells of whose faces will be accumulated.
%   
% RETURNS:
%   faces      - unique faces of a specified cells.
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


% Reconstruct the cells columns in G.cells.faces
c = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';

% Indices of cell 'i' in G.cells.faces
[LiA, ~] = ismember(c, cell);
ind = find(LiA);

% Collect faces of cell 'i'
faces = G.cells.faces(ind, 1);
faces = unique(faces); 


end