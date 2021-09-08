function prop = intracellProp(G, intracells, varargin) 
%
% SYNOPSIS:
%   prop = intracellProp(G, intracells, varargin) 
%
% DESCRIPTION:
%   Find properties of the intracells that represent fractures, such as cells that are neighbours to intracells, faces of fracture (discontinuity) etc.
%
% PARAMETERS:
%   intracells      - cells whose props to find 
% 
% RETURNS:
%   neiintracells   - neighbours to intracells (for the use in ConDetect2(...) function)
%   fneiintracells  - faces of neighbouring intracells
%   sfaces          - surface faces of the fracture
% 
% EXAMPLE:
%
% SEE ALSO:
% 
% AUTHOR:
%   Aman Kubeyev  

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


    % -------------------------------------------------------------------------
    % Find neighbouring cells to the intracells
    
    prop = struct(); 
    cells = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .'; % Reconstruct the cells columns in G.cells.faces

    
    [LiA, ~] = ismember(cells, intracells); % Indices of intracells in G.cells.faces
    ind = find(LiA);

    
    intfaces = G.cells.faces(ind, 1);
    intfaces = unique(intfaces); 
    clear LiA ind

    
    [LiA, ~] = ismember(G.cells.faces(:, 1), intfaces);
    ind = find(LiA);
    cls = cells(ind);
    cls = unique(cls); 
    neiintracells = setdiff(cls, intracells);

    
    % -------------------------------------------------------------------------
    % Find surface faces of the fracture and faces of neighbouring intracells
    
    [LiA, ~] = ismember(cells, neiintracells);
    ind = find(LiA);

    % Collect faces of neiintracells
    fneiintracells = G.cells.faces(ind, 1);
    fneiintracells = unique(fneiintracells); 
    clear LiA ind

    % Find fracture (discontinuity) surface faces
    [~, ind] = ismember(intfaces, fneiintracells);
    ind(ind == 0) = []; 
    sfaces = fneiintracells(ind); 

    
    prop.neiintracells  = neiintracells; 
    prop.fneiintracells = fneiintracells; 
    prop.sfaces         = sfaces; 

    
end