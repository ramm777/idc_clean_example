function C = Enu2C_TransIso(E, nus, G, G2, varargin)
% 
% SYNOPSIS:
%   function C = Enu2C_TransIso(E, nus, G, G2, varargin)
%
% DESCRIPTION: 
%   For each cell, construct the 3x3 (in 2D) matrix representing the elasticity tensor of the Transverse Isotropy (Anisotropy). Default is mode is: 'mode', 'plane_strain',
%   can be changed to: 'mode', 'plane_stress'. Equations:
%
%   n = E1 / E2
%   m = G2 / E2
%
%   In 2D, the matrix of the elasticity tensor for a given cell is written as 
%   Plane Strain
%                E2                | n*(1 - n*nu2^2), n*nu2*(1 + nu1),              0                   |            
%  ------------------------------x | n*nu2*(1 + nu1), (1 - nu1^2)    ,              0                   |     
%  (1 + nu1)(1 - nu1 - 2*n*nu2^2)  |        0       ,      0         ,   m(1 + nu1)(1 - nu1 - 2n*nu2^2) |            
%
%
%   Plane Stress
%
%         E2           | n     ,      n*nu2  ,      0           |            
%  ------------------x | n*nu2 ,      1      ,      0           |     
%     (1 - n*nu2^2)    | 0     ,      0      ,   m*(1 - n*nu2^2)|            
%
%   Transverse Isotropic material can be described by 5 independant constants: E1,E2,nu1,nu2,G2 
%   For more infomration, see Zienkievicz, chapter 4
%
% PARAMETERS:
%   E   - Young's modulus on 2 directions (two entry per cell)
%   nus - Poisson ratio on 2 directions (two entry per cell)
%   G   - Grid
%   G2  - Shear Modulus (one entry per cell)
% 
% RETURNS:
%   C - (k,n) matrix, where k=3^2 (2D), and n is the number of cells. 
%   Each row thus represents the entries of the elasticity tensor for a specific cell.
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



% Assert that it's a 2D
assert(G.griddim == 2, 'Error: this function is for 2D only');

% Default mode: plane_strain
opt = struct('mode', 'plane_strain');

% Override default control options if 'mode' is 'plane_stress' in varargin
opt = merge_options(opt, varargin{:});

% Extract axial elasticity parameters
E1 = E(:, 1); % Ex
E2 = E(:, 2); % Ey

nu1 = nus(:, 1); % nu1
nu2 = nus(:, 2); % nu2

n = E1./E2;
m = G2./E2; 

if opt.mode == 'plane_strain'
% Construct matrix 
    z = zeros(numel(nu1), 1);
    C = [reshape([n.*(1 - n.*(nu2.^2)), n.*(nu2.*(1 + nu1)), z                                    ]', [], 1), ...
       reshape([n.*nu2.*(1 + nu1)   ,(1 - nu1.^2)          , z                                    ]', [], 1), ...
       reshape([z                   ,z                     , m.*(1 + nu1).*(1 - nu1 - 2*n.*(nu2.^2)) ]', [], 1) ];

    nlin = 3;

    fac = (E2 ./ ((1 + nu1) .* (1 - nu1 - 2 *n.*(nu2.^2))));
    % Row is each cell's matrix
    C   = reshape(C', nlin * nlin, [])';
    C   = bsxfun(@times, C, fac);

elseif opt.mode == 'plane_stress'
    
    z = zeros(numel(nu1), 1);
    x = ones(numel(nu1), 1);
    C = [reshape([n      ,  n.*nu2,         z              ]', [], 1), ...
       reshape([n.*nu2 ,  x     ,         z              ]', [], 1), ...
       reshape([z      ,  z     ,  m.*(1 - n.*(nu2.^2))  ]', [], 1) ];

    nlin = 3;
    
    fac = E2./(1 - n.*(nu2.^2));
    % Row is each cell's matrix
    C   = reshape(C', nlin * nlin, [])';
    C   = bsxfun(@times, C, fac);

else
    warning('No such mode. Select mode: plane_strain (default if empty) OR plane_stress')

end
