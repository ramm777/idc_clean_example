function h = plotGcd(Gcd, lhsbody, rhsbody, varargin) 


% SYNOPSIS:
%   h = plotGcd(Gcd, lhsbody, rhsbody) 
%
% DESCRIPTION: Plot contact detection grid Gcd with the lhsbody and rhsbody. 
%  
% PARAMETERS:
%   Gcd      - Contact detection grid
%   lhsbody  - left hand side body
%   rhsbody  - right hand side body
%   closestpts - nodes pairs that detected in contact detection algorithm
%
% RETURNS:
%   h  - handle to plot
% 
% EXAMPLE:
%
% SEE ALSO: plotGrid()
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

% Screensize for future figures 
ss = get(groot, 'Screensize'); 

if nargin > 3
    closestpts = varargin{1}; 
else 
    closestpts = []; 
end

% -------------------------------------------------------------------------
% Plot 


%plotGrid(G2, 'FaceColor', 'none');
h = plotGrid(Gcd, rhsbody, 'FaceColor', 'none', 'EdgeAlpha', 0.05);
h = plotGrid(Gcd, lhsbody, 'FaceColor', 'cyan', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.05);
%plotCellNumbers(G2, demcells);
%plotNodeNumbers(G2);
axis tight; % zoom(10);
hold on

if ~isempty(closestpts)
    % Plot connodes 1st column
    p3 = plot(Gcd.nodes.coords(closestpts(:, 1), 1), Gcd.nodes.coords(closestpts(:, 1), 2), 'square');
    hold on
    % Plot connodes 2nd column
    p4 = plot(Gcd.nodes.coords(closestpts(:, 2), 1), Gcd.nodes.coords(closestpts(:, 2), 2), 'o');
    hold off
    
    legend([p3, p4], '=> Must be on LHSbody', '=> Must be on RHSbody'); 
    
end

title('Contact Detection Grid: Gcd'); 

end