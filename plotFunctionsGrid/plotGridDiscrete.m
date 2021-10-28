function h = plotGridDiscrete(G, rhsCells, lhsCells, varargin)


% SYNOPSIS:
%   h = plotGridDiscrete(G, u, rhsCells, lhsCells, varargin)
%
% DESCRIPTION: Plot discrete grid model that consists of two bodies lhs body and rhs body described by the cell indices
%
% PARAMETERS:
%   G        - Grid structure
%   rhsCells - cells that defines the rhs body
%   lhsCells - cells that defines the rhs body 
%   'plotCellNumbers' - (optional) if true, then plot cell numbers (indices) superimposed on the grid
%   'plotFaceNumbers' - (optional) if true, then plot face indices superimposed on the grid
%   'plotNodeNumbers' - (optional) if true, then plot node indices superimposed on the grid
% 
% RETURNS:
%   h - handle to plot
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
    
    % -------------------------------------------------------------------------
    % Defaults
    
    opt     = struct('plotCellNumbers',      false, ...
                     'plotFaceNumbers',      false, ...
                     'plotNodeNumbers',      false ); 
                 
    % Overwrite defaults if specified
    opt     = merge_options(opt, varargin{:});
    
    
    % -------------------------------------------------------------------------
    % Plot
    
    h = plotGrid(G, rhsCells, 'FaceColor', 'yellow', 'EdgeAlpha', 0.1); hold on
    h = plotGrid(G, lhsCells, 'FaceColor', 'red',    'EdgeAlpha', 0.1);   
    
    
    intracells = setdiff(G.cells.indexMap, [lhsCells; rhsCells]); 
    
    h = plotGrid(G, intracells, 'EdgeAlpha', 0.01, 'FaceColor', 'none');
    
    title('Discrete Initial Grid');
    axis equal;
    
    % -------------------------------------------------------------------------
    % Plot additional grid info if requested
    
    if nargin > 3
        if opt.plotCellNumbers == true
            h = plotCellNumbers(G,'FontSize', 7);
        end
        if opt.plotFaceNumbers == true
            h = plotFaceNumbers(G,'FontSize', 6, 'Color','red');
        end
        if opt.plotNodeNumbers == true
            plotNodeNumbers(G);
        end
    end
    
    % -------------------------------------------------------------------------
    
    hold off
    legend('RHS body', 'LHS body');  
    xlabel('x-axis, m');
    ylabel('y-axis, m');

end