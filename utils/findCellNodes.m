function nodes = findCellNodes(G, cells) 
%
% SYNOPSIS:
%   nodes = findCellNodes(G, cells) 
%
% DESCRIPTION:
%   Returns a unique sets of nodes of the specified cells. Grouped to LHS(-) and RHS(+) of the cell. Grouping is valid only if if cells are MRST convention nodes
%   placement => upper_left, lower_left, lower_right_upper_right
%    
% PARAMETERS:
%   G       - Grid structure as described by grid_structure.
%   cells   - cells, whos nodes to be output
%   
%   
% RETURNS:
%   nodes   - unique nodes of the cell
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



    
    cn = cellNodes(G);                     % Mapping from cells to nodes (vertices) matrix, in correct format
    rows = find(ismember(cn(:,1), cells)); % find rows of all cells in fracture
    nodes = cn(rows, 3);                   % find nodes of all these cells
    nodes = reshape(nodes, 2, [])';        % Group LHS(-) and RHS(+) of the cell
    
    
    nodes = unique(nodes, 'rows');         % Remove duplicates, sorted on rows 

    
    % Sort fnodes from the bottom, then ascending
    d = sqrt(sum(bsxfun(@minus, repmat([0, 0], numel(nodes(:, 1)), 1), ...
                                            G.nodes.coords(nodes(:, 1), :)).^2, 2)); 
    [~, i] = sort(d);
    nodes = nodes(i, :); 

    
    
    
end