function G = moveBodiesCloser(G, gapr, rhsCells, lhsCells, fnodes)


% SYNOPSIS:
%   G = moveBodiesCloser(G, gapr, rhsCells, lhsCells)
%
% DESCRIPTION: To save simulation time, this function moves bodies close to each other. The 
%              The distance between them is small and equal to gapr (residual gap)              
% 
% PARAMETERS:
%   G        - Grid structure
%   gapr     - residual gap, the distance that will be present between 2 bodies  
%   rhsCells - cells that defines the rhs body
%   lhsCells - cells that defines the rhs body 
%   fnodes   - fracture nodes, must be sorted correctly [LHS body,  RHS body]
% 
% RETURNS:
%   G - updated grid where rhsbody and lhsbody are closer to each other. 
%
% EXAMPLE:
%
% SEE ALSO:
%
% AUTHOR:
%   Amanzhol Kubeyev  

%{
Copyright 2009-2021 SINTEF ICT, Applied Mathematics.

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


    rhsnodes = findCellNodes(G, rhsCells);  
    lhsnodes = findCellNodes(G, lhsCells);
    
    % -------------------------------------------------------------------------
    % Assertions
    
    [Lia, ~] = ismember(fnodes(:, 1), rhsnodes(:)); 
    if any(Lia)
        fnodes(Lia, :) = flip(fnodes(Lia, :), 2); 
    end

    [Lia, ~] = ismember(fnodes(:, 2), lhsnodes(:));
    assert(~any(Lia), 'fnodes must be sorted: [LHS body,  RHS body]'); 
    
    % -------------------------------------------------------------------------
    % Find the smallest distance and move bodies 
    
    A = [G.nodes.coords(fnodes(:, 1), 1), G.nodes.coords(fnodes(:, 1), 2)]; % xy coords of the first points columns
    B = [G.nodes.coords(fnodes(:, 2), 1), G.nodes.coords(fnodes(:, 2), 2)]; % xy coords of the second points columns
    horzdist = B(:, 1) - A(:, 1);    
    xmov =  min(horzdist) - 10*gapr; % Find smallest distance btwn LHSbody and RHSbody, then move LHSbody for [min - 10*gapr]
        
    
    % -------------------------------------------------------------------------
    % Move bodies closer 
    
    G.nodes.coords(lhsnodes(:), 1) = G.nodes.coords(lhsnodes(:), 1) + xmov; 
    G = createAugmentedGrid(G);
    G = computeGeometry(G);


end