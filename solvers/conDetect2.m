function [closestpts, incontact] = conDetect2(G, demcells)
%
% SYNOPSIS:
%   incontact = conDetect2(G, demcells)
%
% DESCRIPTION:
%   - Find demcells that are in contact with nodes of other cells based on inpolygon(). Also finds closest points to the
%     each contacting point
%   - Node-to-node Contact Detection Algorithms here.
%    
% PARAMETERS:
%   G          - Grid structure as described by grid_structure.
%   demcells   - Cells (or Demcells) to investigate for contact
%   
% RETURNS:
%   incontact  - contacting pair of demcell vs node of global names
%   closestpts - closest points to each of the contacting point 
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


%% Data preparation


% -------------------------------------------------------------------------
% Matrix of demcells coords for same format as in inpolygon()

x = zeros(numel(demcells), 4); 
y = zeros(numel(demcells), 4);

for z = 1:numel(demcells)
    
    coord = cellCoords(G, demcells(z)); % For z demcell extract coords
    x(z, :) = coord(:, 1)';             % Place x and y coords as 2 separate matrices 
    y(z, :) = coord(:, 2)'; 
    
end


% Sort for function inpolygon() to work, otherwise it is not a polygon
sx = [x(:, 1), x(:, 2), x(:, 4), x(:, 3)]; % Sorted x
sy = [y(:, 1), y(:, 2), y(:, 4), y(:, 3)]; % Sorted y 
cn = cellNodes(G);                         % Mapping from cells to nodes (vertices) matrix, in correct format


clear coord

% -------------------------------------------------------------------------
% Collect demcells' nodes only

dcrows = [];
for j = 1 : numel(demcells) 
    dcrows  =  [dcrows; find(cn(:, 1) == demcells(j))];
end
dcnodes = cn(dcrows, 3);
dcnodes = unique(dcnodes, 'stable');

clear j dcrows


% -------------------------------------------------------------------------
% Get coords from the Grid for demcells nodes only, rest is NaN

ncoords = nan(size(G.nodes.coords));
ncoords(dcnodes, :) = G.nodes.coords(dcnodes, :);


%% Loop over all demcells to check if incontact with nodes (vertices)


incontact  = [];
closestpts = []; 
k = 0;


for i = demcells' 
   k = k + 1;
   
   % Find nodes of the current demcell
   rows = find(cn(:,1) == i); % find rows of the currect demcell in cn mapping
   nodes = cn(rows, 3);       % find nodes of the current demcell in cn mapping
   
   
   in = inpolygon(ncoords(:, 1), ncoords(:, 2), sx(k, :), sy(k, :)); % Find demcells' nodes incontact including current demcell (self-contact)
    
   
   in(nodes, :) = 0; % Exclude current demcell-i nodes from contact 
   
   
   % Helper figure to visualise
%    figure(1);
%    ex_ncoords = ncoords; 
%    ex_ncoords(nodes, :) = []; % where excluded self-contacting nodes
%    '*' - Nodes excluding polygon nodes
%    plot(sx(k, :), sy(k, :), ex_ncoords(:, 1), ex_ncoords(:, 2), '*')
%    Add node numbers as the markers 
%    markertext = (1 : numel(ncoords(:, 1)))';
%    text(ncoords(:, 1), ncoords(:, 2), num2str(markertext), 'FontSize', 10, 'Color','red')

%    % Plot cell numbers and Grid
%    figure(4); 
%    plotGrid(G, 'FaceColor', 'none');
%    plotCellNumbers(G, demcells);
%    markertext = (1 : numel(ncoords(:, 1)))';
%    text(ncoords(:, 1), ncoords(:, 2), num2str(markertext), 'FontSize', 10, 'Color','red')
   

   %% In case of a contact, extract information
   
   if any(in) 

       
       pt = find(in); % Find which one point is in contact

       % Loop over points that are incontact (for several points incontact with one cell)
       for m = pt'
             
             incontact = [incontact; [i, m]]; % Collect demcells vs contacting points matrix

             % Connect closest nodes if not previously connected in closestpts
             if ~ismember(m, closestpts) 

                 % Create nan matrix similar to ncoords
                 ex_ncoords2 = nan(size(ncoords));

                 % Find i-s demcell coords
                 ex_ncoords2(nodes, :) = ncoords(nodes, :);

                 % Remove coords that's in closestpts or incontact
                 ex_ncoords2(closestpts, :) = nan;

                 % Find distances btwn contacting point and all other points excluding former
                 distances = sqrt(sum(bsxfun(@minus, ncoords(m, :), ex_ncoords2).^2, 2));

                 % Find coords of the closest point
                 %closest_coords = ncoords(distances==min(distances), :);

                 % Find closest point, 1st point that identified if 2 
                 closest_pt = find(distances==min(distances), 1);

                 % Collect contacting points vs closest point nearby matrix
                 closestpts = [closestpts; [m, closest_pt]];
             end

       end 

   end
   
end




% APPENDIX
% Sort on a row direction to remove duplicates
% closestpts = sort(closestpts, 2); 

% Remove duplicates of nodes connecttions
% closestpts = unique(closestpts, 'rows'); 