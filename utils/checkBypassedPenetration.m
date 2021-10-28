function y = checkBypassedPenetration(G, closestpts, fraccells)  
%
% SYNOPSIS:
%   y = checkBypassedPenetration(G, closestpts, fraccells)  
%
% DESCRIPTION: 
%   Check for the 'Bypassed Penetration'. If the any node in closestpts does not belong to fraccells then 
%   => it's a 'Bypassed Penetration', error.
%    
% PARAMETERS:
%   G          - Grid structure as described by grid_structure.
%   closestpts - Connecting nodes (in node-to-node it's closestpts)
%   fraccells  - Cells that represent fracture (sometimes intracells) 
%   
% RETURNS:
%   error if there is a bypassed penetration
%   
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



    
    fracnodes = findCellNodes(G, fraccells);  
    
    [lia] = ismember(closestpts(:), fracnodes); 
    assert(all(lia), 'IDC error: Bypassed Penetration, fix by decreasing contact detection step or in other way');        
    
    
    
end