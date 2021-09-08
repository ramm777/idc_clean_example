function [uu, extra] = VEM_linElast_LMd(G, C, el_bc, load, varargin)
% 
%
% SYNOPSIS:
%   function [uu, extra] = VEM_linElast(G, C, el_bc, load, varargin)
%
% DESCRIPTION: Directional Lagrange Multipliers connecting DOFs on either 'x' or 'y' in 2D. Contacting nodes pairs removed if any of the node is Dirichlet.
%              Important: connodes must follow [LHSbody, RHSbody] or [Lower, Upper]. 
%              
%
%
% PARAMETERS:
%   G        - Grid structure as described by full_grid_structure
%   C        - Elasticity tensor
%   el_bc    - Elastic boundary condition structure. It contains the fields
%             'disp_bc'  : displacement boundary condition. It contains the
%                          fields
%                  'nodes'    : nodes where the displacement condition is applied
%                  'uu'       : value for the displacement
%
%                  The two following fields are not used in the VEM implementation
%                  but becomes relevant for other methods such as MPSA, see paper
%
%                  'faces'    : faces displacement
%                  'uu_faces' : value for the displacement
%
%                  'mask'     : if false then displacement values that are
%                               imposed in given Cartesian directions are in
%                               fact ignored.
%             'force_bc'  : force boundary condition applied on faces. It contains the
%                           fields
%                  'faces' : faces where the force is applied
%                  'force' : value of the force that is applied
%
%   load        - loading term
%   lagnum      - number of lagrange multipliers to be placed (depends on the conDetect2())
%   connodes    - connecting (contacting) nodes matrix, where each row is one contacting node pair.
%                 must always be sorted: [LHS body, RHS body].   
%   g           - gap (taken before the penetration and all equations are in consistency)
%   direction   - direction in which the Contact Problem will be applied: 'x' or 'y'
%   ginfo       - grid information structure 
%
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%  'linsolve'             - Linear solver
%  'blocksize'            - block size used in the assembly
%  'add_operators'        - Add operator in output
%  'force_method'         - Method for computing the loading term, see function calculateVolumeTerm below for
%                           a list of the possible alternatives.
%  'alpha_scaling'        - Coefficient of the stabilisation term (default 1)
%  'S'                    - Stabilization matrix to use (used only in very special cases, experimental feature)
%  'experimental_scaling' - Scaling proposed in [Andersen et al: http://arxiv.org/abs/1606.09508v1]
%  'pressure'             - Pressure field, used at loading term
%
% RETURNS:
%   uu    - Displacement field
%   extra - Extra outputs
%
% EXAMPLE:
%
% SEE ALSO:
%
% AUTHOR:
%   Aman Kubeyev (based on SINTEF's VEM_linElast() original code)

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

    opt = struct('linsolve'     , @mldivide                       , ...
                 'blocksize'    , 30000 - (G.griddim - 2) * 20000 , ...
                 'add_operators', true                            , ...
                 'force_method' , 'dual_grad_type'                , ...
                 'alpha_scaling', 1                               , ...
                 'S'            , []                              , ...
                 'experimental_scaling', false                    , ...
                 'pressure'     , []                              , ...
                 'no_solve'     , false                           , ...
                 'lagnum'       , 0                               , ...
                 'connodes'   , []                              , ...
                 'gap'          , []                              , ...
                 'direction'    , []                              , ...
                 'ginfo'        , []    );
    opt = merge_options(opt, varargin{:}); 
    opt.add_operators = opt.add_operators && nargout>1;
    if(opt.add_operators)
        [S, extra] = VEM_assemble(G, C, ...
                                  'blocksize'            , opt.blocksize, ...
                                  'alpha_scaling'        , opt.alpha_scaling, ...
                                  'S'                    , opt.S, ...
                                  'experimental_scaling' , opt.experimental_scaling);
    else
        S = VEM_assemble(G, C, ...
                         'blocksize'           , opt.blocksize, ...
                         'alpha_scaling'       , opt.alpha_scaling, ...
                         'S'                   , opt.S, ...
                         'experimental_scaling', opt.experimental_scaling);
    end

    % Recalculate "weights" (all are calculated in the assembly, they could in fact only be calculated once).
    if (G.griddim == 2)
        qf_all = G.faces.areas / 2;
        [qc_all, qcvol] = calculateQF(G);
    else
        [qc_all, qf_all, qcvol] = calculateQC(G);
    end

    % Apply Diriclet boundary conditions
    % NB: Only rolling conditions in the Cartesian directions is allowed for now.
    % If duplicated  values exist, remove them.
    bc = el_bc.disp_bc;
    
    
    %----------------------------------------------------------------------
    % MODIFIED (fixed a bug) by Amanzhol Kubeyev, how u_bc is collected
    
    [bcnodes, ~, i] = unique(bc.nodes); % i will be used in the original mask assignment
      
    
    % Sum all the non-unique bc to the unique for each column [x y] 
    for z = 1 : G.griddim
        u_bc(:, z) = accumarray(i, bc.uu(:, z));  % make size of unique number of nodes?
    end
    clear z
    
      
    if numel(find(bc.uu(:, 1))) ~= numel(find(u_bc(:, 1))) || numel(find(bc.uu(:, 2))) ~= numel(find(u_bc(:, 2)))
        error('disp_bc has more that 1 definiction per DOF')
    end
    
        
    mask = zeros(numel(bcnodes), G.griddim);
    for k = 1:G.griddim
        mask(:, k) = accumarray(i, bc.mask(:, k), [numel(bcnodes), 1]);
    end
    mask = mask > 0;

    
    %----------------------------------------------------------------------     
     
    % Find the logical vector to remove Diriclet boundary conditions
    ndof = G.griddim * G.nodes.num;
    if(all(mask(:) == true))
        dirdofs = mcolon(G.griddim * (bcnodes - 1) + 1, G.griddim * (bcnodes - 1) + G.griddim)';
        u_bc = reshape(u_bc', [], 1);
    else
        dirdofs = mcolon(G.griddim * (bcnodes - 1) + 1, G.griddim * (bcnodes - 1) + G.griddim)';
        dirdofs = reshape(dirdofs, G.griddim, []);
        mm = mask';
        ind = find(mm);
        u_bc = reshape(u_bc', [], 1);
        dirdofs = reshape(dirdofs, [], 1);
        dirdofs = dirdofs(ind);
        u_bc = u_bc(ind);
    end
    isdirdofs = false(ndof, 1);
    isdirdofs(dirdofs) = true;


    
    V_dir = nan(ndof, 1);      % Calculate the boundary conditions
    V_dir(isdirdofs) = u_bc(:);
    V_dir(~isdirdofs) = 0;
    rhso = -S * V_dir;

    % Calculate load terms
    f = calculateVolumeTerm(G, load, qc_all, qcvol, opt);
    if(~isempty(opt.pressure))
        div = VEM_div(G);
        f_pressure = div'*opt.pressure;
    else
        f_pressure = zeros(size(f));
    end

    
    rhs = rhso + f + f_pressure;
    forces = f + f_pressure;

        
    % Add boundary forces
    if(~isempty(el_bc.force_bc))
        if(numel(el_bc.force_bc.faces) > 0)
            bc       = el_bc;
            bc_force = bc.force_bc;
            bcf      = bc_force;

            % Find weights
            faces    = bcf.faces;
            assert(numel(faces) == numel(unique(faces))); % (A) returns the same data as in A, but with no repetitions
            lnn      = G.faces.nodePos(faces + 1) - G.faces.nodePos(faces);
            inodes   = mcolon(G.faces.nodePos(bcf.faces), G.faces.nodePos(bcf.faces + 1) - 1);
            nodes    = G.faces.nodes(inodes);
            lfacenum = rldecode([1:numel(faces)]', lnn);
            assert(all(sum(G.faces.neighbors(faces, :) == 0, 2) == 1));

            if(G.griddim == 2)
                qf = qf_all(rldecode(faces(:), lnn), :);
            else
                qf = qf_all(inodes);
            end

            % Find nodes corresponding to forces
            dofs  = mcolon(G.griddim * (nodes - 1) + 1, G.griddim * (nodes - 1) + G.griddim);
            force = bsxfun(@times, bc_force.force(lfacenum, :), qf);
            % Transform to degree of freedom indexing
            force = force';
            fbc   = accumarray(dofs', force(:)', [ndof, 1]);
            rhs   = rhs + fbc;
            
            % Collect forces
            forces = forces + fbc;
        end
    end
    
    %% CONTACT INTERACTION (modified by Amanzhol Kubeyev)
    
    
    connodes = opt.connodes;
    lagnum = opt.lagnum / 2; % Number of lambdas-x only
    indof = ndof;            % Save size(x) of initial S for later = ndof
    
    
    % Enlarge matrices depending on lagnum
    S = [S, zeros(size(S, 2), lagnum)];
    S = [S; zeros(lagnum, size(S, 2))];
    rhs = [rhs; zeros(lagnum,1)];
    isdirdofs = [isdirdofs; false(lagnum, 1)];
    V_dir = [V_dir; zeros(lagnum, 1)];
    ndof = ndof + lagnum; 
   
    % Assert connodes sortion, based on lhsnodes/rhsnodes and connodes(:, 2) 
    rhsnodes = findCellNodes(G, opt.ginfo.rhsbody);  % rhsbody's nodes
    lhsnodes = findCellNodes(G, opt.ginfo.lhsbody);  % lhsbody's nodes
    [Lia, ~] = ismember(connodes(:, 2), lhsnodes(:));
    assert(~any(Lia), 'Connecting nodes must be sorted: [LHS body,  RHS body] and gap positive for equations to be correct'); 
    clear Lia
    [Lia, ~] = ismember(connodes(:, 1), rhsnodes(:));
    assert(~any(Lia), 'Connecting nodes must be sorted: [LHS body,  RHS body] and gap positive for equations to be correct'); 
    connodes(any(connodes == 0, 2), :) = []; % If any row contains 0, then make this row empty
    
    
    % Construct LM matrix
    xy = sort(reshape(connodes, [], 1)); 
    DOFx = xy*2 - 1;
    DOFy = xy*2;  
    DOF = [DOFx, DOFy]; 
    
    
    % Remove contacting pair where any node's DOF == Dirichlet
    for x = 1 : numel(connodes)
         nd  = connodes(x);      % Select one node 
         ds = DOF(xy == nd, :);  % Get [DOFx, DOFy] of that node 
         if any(ds == dirdofs)   % Check if any of [DOFx, DOFy] are same as dirdofs  
              connodes(x) = 0;   % If yes set this node to 0 
         end

    end
    
    
    %---------------------------------------------------------------------- 
    % Insert 'contact problem' equations by modifying the matrices
    
    switch opt.direction
        
        case('x')
           
            for z = 1 : numel(connodes(:, 1))
                
               % Take each node of contacting nodes pairs to connect them
               n1 = connodes(z, 1); 
               n2 = connodes(z, 2);

               % Find DOFs (x,y) of each node
               dx1 = DOF(xy == n1, 1); % node1 DOF-x
               dx2 = DOF(xy == n2, 1); % node2 DOF-x

               % Insert lambdas to S for contacting nodes 1 and 2 
               S(dx1, indof + z) =  1; % lambda1 horiz
               S(dx2, indof + z) = -1; % lambda1 horiz

               % Insert u = u + g to the equation to S 
               S(indof + z, dx1) =  1; % u horiz
               S(indof + z, dx2) = -1; % u horiz

               % Insert gaps to the rhs of the equation
               rhs(indof + z) = opt.gap(z); % u horiz
       
            end
            
        case('y')
                      
            for z = 1 : numel(connodes(:, 1))
                
               % Take each node of contacting nodes pairs to connect them
               n1 = connodes(z, 1); 
               n2 = connodes(z, 2);

               % Find DOFs (x,y) of each node
               dy1 = DOF(xy == n1, 2); % node1 DOF-y
               dy2 = DOF(xy == n2, 2); % node2 DOF-y

               % Insert lambdas to S for contacting nodes 1 and 2 
               S(dy1, indof + z)   =  1; % lambda2 vert
               S(dy2, indof + z)   = -1; % lambda2 vert

               % Insert u = u + g to the equation to S 
               S(indof + z, dy1)   =  1; % u vert
               S(indof + z, dy2)   = -1; % u vert

               % Insert gaps to the rhs of the equation
               rhs(indof + z) = opt.gap(z); % u vert

            end   
        
        otherwise
            disp('Specify either x or y DOFs to connect in LM')
    end
       
    
    B = S(:, indof + 1 : end); % Output B matrix    
    
    
    %% SOLVE
    
    
    rhs = rhs(~isdirdofs);          % Reduce the degrees of freedom
    A   = S(~isdirdofs, ~isdirdofs);  
    if(opt.no_solve)
        x = nan(sum(~isdirdofs),1);
    else
        x = opt.linsolve(A, rhs);
    end
    
    
    %% OUTPUT
    
    
    lagmult = x(end - (lagnum-1): end);      % Extract Lagrange Multipliers
    x = x(1 : end-lagnum);                   % Cut 'x' from bottom to lagnum are lagmults
    isdirdofs = isdirdofs(1 : end - lagnum); % Cut isdirdofs and V_dir, note: logical
    V_dir = V_dir(1 : end - lagnum);
    ndof = ndof - lagnum;                    % Decrease number of DOFs
    
    
    u   = nan(ndof, 1); 
    u(isdirdofs)  = V_dir(isdirdofs);        % u is 0 where Dirichlet
    u(~isdirdofs) = x;                       % u is solution where no Dirichlet
    uu = reshape(u, G.griddim, [])';

    
    if(nargout == 2) 
        extra.A    = A;
        extra.S    = S;
        extra.rhs  = rhs;
        vdiv       = VEM_div(G);
        extra.disc = struct('A'         , A                        , ...
                            'S'         , S                        , ...
                            'isdirdofs' , isdirdofs                , ...
                            'rhs'       , rhs                      , ...
                            'V_dir'     , V_dir                    , ...
                            'ovol_div'  , vdiv                     , ...
                            'gradP'     , vdiv(:    , ~isdirdofs)' , ...
                            'div'       , vdiv(:    , ~isdirdofs)  , ...
                            'divrhs'    , vdiv * V_dir             , ... 
                            'forces'    , forces                   , ...
                            'lagmult'   , lagmult                  , ...
                            'B'         , B);

    disp('SOLVED: Directional Lagrange Multipliers (LMd)'); 
    
    end

end


function f = calculateVolumeTerm(G, load, qc_all, qcvol, opt)

    cells  = [1:G.cells.num]';
    inodes = mcolon(G.cells.nodePos(cells), G.cells.nodePos(cells + 1) - 1');
    nodes  = G.cells.nodes(inodes);

    switch opt.force_method

      case 'node_force'
        % Evaluate forces at nodes. The result is weighted, using adjacent
        % cell volume contributions, see paper [Gain et al: doi:10.1016/j.cma.2014.05.005]

        X = G.nodes.coords(nodes, :);
        w = qcvol;
        ll = bsxfun(@times, load(X), w)';

      case  'cell_force_baric'
        %
        % Evaluate the forces at nodes. Then, from the nodal values, compute the (exact)
        % L^2 projection on each cell and, using this cell values, assemble
        % the load term in term of the degrees of freedom.
        %
        % The VEM theory tells us that there exist a virtual basis such that the
        % two steps above can be done exactly. See Ahmad et al (doi:10.1016/j.camwa.2013.05.015)
        %

        nlc = G.cells.nodePos(cells + 1) - G.cells.nodePos(cells);
        X = rldecode(G.cells.centroids(cells, :), nlc);
        lcellnum = rldecode(cells, nlc);
        BB = nan(numel(cells), 1);
        for i = 1:G.griddim
            BB(:, i) = accumarray(lcellnum, G.nodes.coords(nodes, i), [numel(cells), 1]);
        end
        fac  = accumarray(lcellnum, 1, [numel(cells), 1]);
        BB   = bsxfun(@rdivide, BB, fac);
        XB   = X - rldecode(BB, nlc);
        vols = rldecode(G.cells.volumes(cells, :), nlc);

        % Weights to calculate the volume force term from the nodal values.
        if(G.griddim == 3)
            w = (vols ./ rldecode(nlc, nlc) + sum(qc_all(inodes, :) .* (XB), 2));
        else
            assert(G.griddim == 2)
            w = (vols ./ rldecode(nlc, nlc) + sum(qc_all(inodes, :) .* (XB), 2));
        end
        ll   = bsxfun(@times, load(X), w)';

      case 'cell_force'
        % Evaluate the force at the cell centroids. Then, for each node, sum up each
        % adjacent cell contributions after weighting them with volume
        % contribution.

        nlc = G.cells.nodePos(cells + 1) - G.cells.nodePos(cells);
        X   = rldecode(G.cells.centroids(cells, :), nlc);
        w   = qcvol;
        ll  = bsxfun(@times, load(X), w)';

      case 'dual_grad_type'
        % For the virtual basis, there exists a natural divergence operator (from node
        % values to cell values) that can be computed exactly. This operator
        % does not depend on the stabilisation term which is chosen for the
        % stiffness matrix
        %
        % By duality we can define a gradient operator (from cell-values to
        % node values). When the force can be expressed as a gradient, this
        % gives us a way to compute the load term which is implemented here.
        %
        % Such computation has appeared to be more stable, see [Andersen et al: http://arxiv.org/abs/1606.09508v1].
        %
        %
        nlc     = G.cells.nodePos(cells + 1) - G.cells.nodePos(cells);
        X       = rldecode(G.cells.centroids(cells, :), nlc);
        rel_vec = -(X-G.nodes.coords(nodes, :));
        ll      = bsxfun(@times, load(X), qc_all.*rel_vec)';

      otherwise
        error('No such force  calculation')
    end

    ndof = G.griddim * G.nodes.num;
    dofs = mcolon(G.griddim * (nodes - 1) + 1, G.griddim * (nodes - 1) + G.griddim)';
    f    = accumarray(dofs(:), ll(:), [ndof, 1]);

end