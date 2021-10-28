function disp_bc = collectUBCs(bc_el_sides)

% SYNOPSIS:
%   bc = findBoundaryNodesAndCreateBC(G, bc)
%
% DESCRIPTION: Collect the pre-defined displacement boundary conditions into a structure
%              
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
    
    nodes = [];
    faces = [];
    mask  = [];
    disp_node  = [];
    disp_faces = [];

    for i = 1 : numel(bc_el_sides)
        if(~isempty(bc_el_sides{i}))
            nodes = [nodes; bc_el_sides{i}.el_bc.disp_bc.nodes]; %#ok
            faces = [faces; bc_el_sides{i}.el_bc.disp_bc.faces]; %#ok
            mask  = [mask; bc_el_sides{i}.el_bc.disp_bc.mask];   %#ok
            disp_node = [disp_node; bc_el_sides{i}.el_bc.disp_bc.uu];
            disp_faces = [disp_faces; bc_el_sides{i}.el_bc.disp_bc.uu_face];
        end
    end
    
    
    disp_bc = struct('nodes', nodes, 'uu', disp_node, 'faces', faces, 'uu_face', disp_faces, 'mask', mask); 


end