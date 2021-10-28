function h = plotNodeNumbers(G, varargin)
%
% SYNOPSIS:
%   h = plotNodeNumbers(G)
%
% DESCRIPTION: 
%   Plots nodes' numbers of the grid, can plot selected nodes only
%
% PARAMETERS:
%   G       - Grid structure as described by grid_structure.
%   nodes   - (optional) nodes which to be plotted only    
%   
% RETURNS:
%   h - handle to plot
%
% EXAMPLE:
%
% SEE ALSO: gridCellNo()
%
% AUTHOR:
%   Amanzhol Kubeyev  

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
if G.griddim == 2
        
    % If no specific nodes to be plotted are specified, plot all nodes 
    if  mod(nargin, 2)
        ns = [1 : G.nodes.num]'; 
    else
       ns = varargin{1}; % nodes to be plotted
       varargin = varargin(2:end);
        
       % Make sure ns is a vector, if not make it.  
       ns = ns(:); 
    
    end

    
    % Names of nodes
    markertext = [repmat('  ', length(ns), 1), num2str(ns)]; 
    
    % Plot names
    h = text(G.nodes.coords(ns, 1), G.nodes.coords(ns, 2), markertext,...
        'FontSize', 8, 'Color', 'r', 'Rotation', -45);
    hold on

    
    % Plot ticks
    h = plot(G.nodes.coords(ns, 1), G.nodes.coords(ns, 2), '.', 'Color','r');
    hold off
    
        
% -------------------------------------------------------------------------
% Plot specific cells

%         % Below doesn't work for the re-meshed grid with Hanging Nodes 
%        
%         % If no specific cells to be plotted are specified, plot all 
%         if mod(nargin, 2)
%              c = 1:G.cells.num;
%         % If specific cells to be plotted are specified 
%         else
%             c = varargin{1}; % cells to be plotted
%             varargin = varargin(2:end);
%             
%             % Select nodes of specified cells to plot
%             warning('if only selected cells to be plot, this works only for quads');
%             cellnodes = reshape(G.cells.nodes, 4, [])';
%             ns = cellnodes(c, :);
%             
%             % Reshape back to the original format
%             ns = reshape(ns', [], 1);
%             ns = unique(ns); 
%             
%             % Find names of nodes
%             markertext = [repmat(' ', numel(G.nodes.coords(:, 1)), 1),...
%                 num2str((1 : numel(G.nodes.coords(:, 1)))')]; 
%             
%             
%             % Plot
%             h = text(G.nodes.coords(ns, 1), G.nodes.coords(ns, 2), markertext(ns, :),...
%                 'FontSize', 10, 'Color','red');  % Names
%             hold on
%             h = plot(G.nodes.coords(ns, 1), G.nodes.coords(ns, 2), 'x'); % Ticks
%             hold off
%                        
%         end
        
end
    

% ---------------------------------------------------------------------
% WARNING: must be updated in same way as for 2D above

% if G.griddim == 3 
% 
%     warning('If only selected cells are to be plotted, 3D doesnt work yet, need update');
% 
%     markertext = [repmat('  ', numel(G.nodes.coords(:, 1)), 1), num2str((1 : numel(G.nodes.coords(:, 1)))')]; 
% 
%     h = text(G.nodes.coords(:, 1), G.nodes.coords(:, 2), G.nodes.coords(:, 3), num2str(markertext),...
%         'FontSize', 10, 'Color','red');
%     hold on
% 
%     h = plot(G.nodes.coords(:, 1), G.nodes.coords(:, 2), G.nodes.coords(:, 3), 'x');
% 
%     hold off
% end
    

end