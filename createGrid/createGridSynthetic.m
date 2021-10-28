function [G, ginfo] = createGridSynthetic(grid_type, varargin)
%
% SYNOPSIS:
%   [G, ginfo] = createGridSynthetic(grid_type, varargin)
%
% DESCRIPTION:
%   - Create a valid numerical grid for the IDC simulations. These grids are always hard-coded in science.  
%   - Here, only one out of ~100 created grids is shown.  
% 
% PARAMETERS:
%   grid_type  - Name of the pre-constructed grid, e.g. 'GaussGrid2'
%
%
% RETURNS:
%   G          - grid with all the structure
%   ginfo      - grid information such as intracells, lhsbody etc.. 
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
% 


%--------------------------------------------------------------------------------------
% CONTROL

mrstModule add vemmech mrst-gui coarsegrid % incomp coarsegrid - for plotCellNumbers(G)
opt = struct('getGridInfo', true);
opt = merge_options(opt, varargin{:});


% Pre-allocate grid information
lhsbody    = []; 
rhsbody    = [];
intracells = [];
demcells   = [];
ellipcells = [];
fraccells  = [];
ellength   = []; 


% -------------------------------------------------------------------------

switch grid_type
          
    
    case('Synthetic1')      
        error('IDC: This grid is temporary abscent'); 

    case('GaussGrid2')
    % Gauss distribution grid 2, it's fixed and doesn't change every time you run it. 
    % If you want to change the roughness distribution change cl, cr in the code below.      
                
        
        nx = 20;  ny = 20;          % Grid dimensions input
        Lx = 2;  Ly = 2;
        x = linspace(0, Lx, nx+1);  % Coords of boundaries
        y = linspace(0, Ly, ny+1);

        
        block_size_av = (((Lx/nx)+(Ly/ny))*0.5);     % Thickness of intra cells
        delta = (0.5)*block_size_av;                 % 1/million of average block face
        x = [x(1:11), x(11)+delta, x(12:end)+delta]; % Add intracell to boundary coords

       
        G = tensorGrid(x, y);
        G.type    = [G.type, { mfilename }];
        G = createAugmentedGrid(G);
        G = computeGeometry(G);
        isintracell = zeros(G.cells.num, 1);
        

        % Collect isintracell (logical). Specify cell number based on whose volume the intracells will be collected 
        benchmarkcell = 1;
        tol = 0.7; 
        isintracell(G.cells.volumes < tol*G.cells.volumes(benchmarkcell)) = 1;
        clear benchmarkcell tol ind

        
        demcells = G.cells.indexMap(~isintracell);  
        intracells = G.cells.indexMap(~~isintracell);
        fraccells = intracells; % Here, intracells are fraccells 

        
        % Collect LHS cells
        clear lhsbody
        low  = intracells - 10;
        high = intracells - 1 ;
        for i = 1:numel(intracells)
            lhsbodyo(i, :) = low(i):1:high(i);
        end
        lhsbody = reshape(lhsbodyo', [], 1);
        rhsbody = setdiff(G.cells.indexMap, [intracells; lhsbody]);

       
    
    %% Collect data of intracells for twisting
        
        %--------------------------------------------------------------------------
        % Collect intracells' nodes only
        
        cn = cellNodes(G);
        introws = [];
        for j = 1 : numel(intracells) 
            introws  =  [introws; find(cn(:, 1) == intracells(j))];
        end

      
        intnodes = cn(introws, 3);
        intnodes = unique(intnodes, 'stable');
        IntCoords = G.nodes.coords(intnodes, :); 

        
        %--------------------------------------------------------------------------
        % Twist coords of intracells

        
        TwIntCoords = IntCoords; 

        % Pre-defined distribution
        cl = [1.00001743125225;0.998937441802593;0.962705741189451;1.00871758484014;0.967113476200512;0.965210795400351;1.00119900565448;0.985772785246442;1.00617235932135;1.01015466708526;1.01286598817808;0.989632613119255;1.00674066434750;1.00150950025473;1.01239104997705;1.00804235619889;1.01346832638978;0.998020931981131;0.997791978157731;1.01511660107958;0.968145168063764]; 
        cr = [0.992431203917290;0.980941083252870;0.994261227959385;1.00973018893073;1.01238590723863;0.984775845359798;0.992933951309753;1.00205537311195;0.995622049363696;1.00452727832892;1.00599896414434;0.986050576615898;0.997347546011062;0.968018581012577;1.01718042565778;0.990563638585088;0.981942250389967;0.996190829748658;0.978570297030519;0.999687135734476;0.991590025052077];
        
        
        % Twist coords of the LHS-x and RHS-x and modify
        TwIntCoords(find(IntCoords(:, 1) == 1), 1) = cl.*IntCoords(find(IntCoords(:, 1) == 1), 1);
        TwIntCoords(find(IntCoords(:, 1) == 1 + delta), 1) = cr.*IntCoords(find(IntCoords(:, 1) == 1 + delta), 1);
        G.nodes.coords(intnodes, :) = TwIntCoords;

        
        G = createAugmentedGrid(G);
        G = computeGeometry(G);

                        
     otherwise
        error('Select grid');
    
end

    if  opt.getGridInfo  % Collect grid information   
        
        ginfo.lhsbody    = lhsbody; 
        ginfo.rhsbody    = rhsbody;
        ginfo.intracells = intracells;
        ginfo.demcells   = demcells;
        ginfo.ellipcells = ellipcells;
        ginfo.fraccells  = fraccells;  
        ginfo.ellength   = ellength ; 
        
    end

end
        
        
        
        
        
        
        
        
        
        
        
        