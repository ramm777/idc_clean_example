function [closestpts, lagnum, g] = dataForContact(G, ginfo, closestpts, gapr)
%
% SYNOPSIS:
%   [closestpts, lagnum, g] = dataForContact(G, ginfo, closestpts, gapr)
%
% DESCRIPTION:
%   1) Sorts connodes (closestpts) correclty for Contact Problem 
%   2) Calculates gap
%   3) Calculates lagrange multipliers number lagnum     
%   Note: Includes residual gap inside gap
%
% PARAMETERS:
%   G           - grid
%   ginfo       - grid information (additional)
%   closestpts  - connecting nodes for the Contact Problem
%   gapr        - 'single' calculate a single kf based on one mean Vm
%       
%   
%  
% RETURNS:
%   closestpts  - sorted correctly connecting nodes 
%   lagnum      - number of lagrange multipliers
%   g           - gaps for each connodes 
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




    rhsbody    = ginfo.rhsbody;
    lhsbody    = ginfo.lhsbody;

    %--------------------------------------------------------------------------
    % Prepare data to Solve the Contact. This is important for accurate solution. 

    
    closestpts = sort(closestpts, 2, 'ascend');  % Connecting nodes must always be sorted: [LHS body,  RHS body] or [LOWER, UPPER]
    closestpts = unique(closestpts, 'rows');     % Remove duplicates of nodes connecttions, and sort
    lagnum = 2*numel(closestpts(:, 1));          % Number of LMs = number of unique points incontact*2 (as 2D) 
    
    
    rhsnodes = findCellNodes(G, rhsbody);  
    lhsnodes = findCellNodes(G, lhsbody);
    
    
    [Lia, ~] = ismember(closestpts(:, 1), rhsnodes(:)); 
    if any(Lia)
        closestpts(Lia, :) = flip(closestpts(Lia, :), 2); 
    end
    
    
    [Lia, ~] = ismember(closestpts(:, 2), lhsnodes(:));
    assert(~any(Lia), 'Connecting nodes must be sorted: [LHS body,  RHS body] and gap positive for equations to be correct'); 
    
         
    %----------------------------------------------------------------------
    % Calculate gap (distances) from G

    
    A = [G.nodes.coords(closestpts(:, 1), 1), G.nodes.coords(closestpts(:, 1), 2)]; % xy coords of the first contacting (closest) points columns
    B = [G.nodes.coords(closestpts(:, 2), 1), G.nodes.coords(closestpts(:, 2), 2)]; % xy coords of the second contacting (closest) points columns 
    horzdist = B(:, 1) - A(:, 1);
   
    
    g = horzdist; 
    if any(g < gapr)
        g(g < gapr) = 0; 
    end
    g(g >= gapr) = g(g >= gapr) - gapr;
       
    
    assert(~any(g < 0), 'Error: gaps must be positive, fix');
    
    
end