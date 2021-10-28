function bc = findBoundaryNodesAndCreateBC(G, bc)


% SYNOPSIS:
%   bc = findBoundaryNodesAndCreateBC(G, bc)
%
% DESCRIPTION: Find the node of the different outer boundary sides and 
%              set up the elastisity boundary conditions
%                             
% 
% PARAMETERS:
%   G         - Grid structure
%   bc        - boundary condition structure usually created by findBoundaryFacesAndCreateBC()
% 
% 
% RETURNS:
%   bc       - updated bc structure, where in each of 4 bcs, disp_bc structure is added and force_bc empty added
%   
%
% EXAMPLE:
%
% SEE ALSO:
%
% AUTHOR:
%   Aidan Kubeyev, based on MRST existing routines

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


    for i = 1 : 2*G.griddim
        inodes = mcolon(G.faces.nodePos(bc{i}.face), G.faces.nodePos(bc{i}.face+1)-1);
        nodes = unique(G.faces.nodes(inodes));
        disp_bc = struct('nodes'  , nodes, ...
                         'uu'     , [], ...
                         'faces'  , bc{i}.face, ...
                         'uu_face', [], ...
                         'mask'   , true(numel(nodes), G.griddim));
        bc{i}.el_bc = struct('disp_bc', disp_bc, 'force_bc', []);
        
    end


end