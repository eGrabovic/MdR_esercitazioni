classdef kinematicTree < handle
    %
    % obj = kinematicTree(p, s, k, mu, ni, gT, Xu, h)
    %
    properties
        
        Nb         % number of bodies
        p          % predecessor array (parent nodes of i-th joint with polarity)
        s          % successors array (children nodes of i-th joint with polarity)
        lambda     % parent array
        k          % cell array of joints from base to node ; i -> element of the array
        mu         % cell array of direct childrens of node i
        ni         % cell array of all the nodes cointained in the sub-tree with node i as root
        leaves     % array of all the nodes without childrens
        Goff       % array of offset homogenous matrices from joint i-1 to joint i (both with respect to body lambda(i) )
        h          % array of helical lead of each joint (0 -> revolute; Inf-> prismatic; otherwise helical joint
        I          % inertias
        G_local    % local transform matrices for each body
        G_localJ   % local transform matrices for each joint
        G_global   % global transform matrices for each body
        G_globalJ  % global transform matrices for each joint
        V          % twists
        Vdot       % twist derivative
        xdot_sym   % state space dynamics in casadi symbolic form
        Cbar_fun
        pinvC
        lam
        CCRdq
        Xu % unitary twists (body) of the joints
        
        % ***!!!  mu{i} and ni{i} are associated to the (i-1)-th node !!!***
    end
    
    properties (Hidden = true)
        
        graphics
        adjoints
        
    end
    
    methods (Static, Hidden = true)
        
        function [x, y, z] = createCylinder(r, len)
            
            theta = linspace(0, 2*pi, 15);
            x = r.*cos(theta);
            y = r.*sin(theta);
            z = len/2.*ones(1,15);
            
        end
        
        function [points, faces] = createParallelepiped(edge, len)
            
            a = edge;
            if length(len) == 2
                h1 = len(1);
                h2 = len(2);
            else
                h1 = 0;
                h2 = len;
            end
            P1 = [a/2;-a/2;-h1]; P2 = [a/2;a/2;-h1]; P3 = [a/2;a/2;h2]; P4 = [a/2;-a/2;h2];
            P5 = [-a/2;a/2;-h1]; P6 = [-a/2;a/2;h2]; P7 = [-a/2;-a/2;-h1]; P8 = [-a/2;-a/2;h2];
            points = [P1,P2,P3,P4,P5,P6,P7,P8];
            faces = [1,2,3,4;2,5,6,3; 3,6,8,4; 8,6,5,7; 1,4,8,7; 1,7,5,2];
            
        end
        
    end
    
    methods
        
        function obj = kinematicTree(p, s, gT, X, h)
            
            obj.Xu = X;
            obj.p = p;
            obj.s = s;
            obj.lambda = min(p, s);
            obj.Goff = gT;
            obj.h = h;
            obj.Nb = length(s);
            %             obj.Nj = length(h);
            
            % compute the i-th body sub tree from node to body
            for i = 1:obj.Nb
                c = p(i);
                obj.k{i} = i;
                while c ~= 0 % until we arrive to root
                    obj.k{i} = [c, obj.k{i}];
                    if c == p(c) % or until the predecessor is the node itself
                        break
                    end
                    c = p(c);
                end
            end
            
            % compute the i-th body direct childrens
            for i = 1:obj.Nb+1
                obj.mu{i} = find(obj.lambda == i-1);
            end
            
            for i = 1:obj.Nb+1 % compute i-th body sub tree with all the childerns
                index = i;
                obj.ni{i} = i-1;
                while ~isempty(index)
                    obj.ni{i} = [obj.ni{i} obj.mu{index}];
                    index = [obj.mu{index}] + 1;
                end
            end
            
            obj.leaves = nan(1, obj.Nb+1);
            cont = 0;
            for i = 1:length(obj.mu)
                if isempty(obj.mu{i})
                    cont = cont+1;
                    obj.leaves(cont) = i-1;
                end
            end
            
        end
        
        function FWKin(this, q)
            %
            % FWKin(this, q)
            % forward kinematics computed with the local POE formulation
            % all the intermediate transform matrices are stored
            for i = 1:this.Nb
                this.G_localJ{i} = expTw_helix(this.Xu(:, i), q(i), this.h(i)); % twist exponential joint formulation
                this.G_local{i} = this.Goff(:,:,i)*this.G_localJ{i};
                if this.lambda(i) == 0
                    this.G_global{i} = this.G_local{i};
                    this.G_globalJ{i} = this.Goff(:,:,i);
                else
                    this.G_global{i} = this.G_global{this.lambda(i)}*this.G_local{i}; % property that lambda(i) < i
                    this.G_globalJ{i} = this.G_global{this.lambda(i)}*this.Goff(:,:,i);
                end
            end
        end
        
        function [Jb, Jtot] = bodyJacobian(this, numNode)
            % body jacobian of numNode w.r.t. the root in frame numNode(and thus his pole)
            % for how we computed the FWKin the B_0,j matrix is just
            % G_globalJ{i}
            % we do not actually compute the expTw(twist, q) here but
            % borrow from the FWKin computations.
            % the jacobian uses the FWKin joints values
            
            subTree = this.k{numNode};
            n = length(subTree);
            g = eye(4);
            Jb = zeros(6, n, class(this.G_localJ{1}));
            
            for j = n : -1 : 1
                %hl = this.h(subTree(j)); % helical lead to identify joint type
                X = this.Xu(:, j);
                g = this.G_localJ{subTree(j)}*g;
                Jb(:, j) = adjoint(rigidInverse(g))*X; % adjoint(g_j,n^-1)*Y_i
                g = this.Goff(:,:,subTree(j))*g;
                
            end
            
            Jtot = zeros(6, this.Nb, class(this.G_localJ{1}));
            nn = 1;
            for kk = 1:1:this.Nb
                if ismember(kk, subTree)
                    Jtot(:, kk) = Jb(:, nn);
                    nn = nn+1;
                end
            end
        end
        
        function plotInit(this, options)
            %
            %
            %
            arguments
                this
                options.parent = [];
                options.EEoffset = [];
            end
            if isempty(options.parent)
                figure('color', 'w'); hold on; axis equal
                this.graphics.axes = gca;
            else
                this.graphics.axes = options.parent;
            end
            
            if isempty(this.G_local{end})
                this.FWKin(zeros(this.Nb, 1));
            end
            % precomputing end effector graphics points
            S(:,1) = [0;0;0];
            S(:,2) = [-1;0;0];
            S(:,3) = [0;0;0];
            S(:,4) = [0;0;0];
            
            % precomputing revolute joint gaphics
            [x, y, z] = kinematicTree.createCylinder(0.1, 0.5);
            ptsR = [x;y;z];
            % precomputing prismatic joint gaphics
            [ptsP, faces] = kinematicTree.createParallelepiped(0.1, [0.25 0.25]);
            ptsPslide = kinematicTree.createParallelepiped(0.05, [0.4 0.4]);
            nn = 1;
            
            for i = 1:this.Nb % joint plots
                % initi the hgtr handles
                this.graphics.homogTransfJoint{i} = hgtransform(Parent = this.graphics.axes);
                this.graphics.homogTransfBody{i} = hgtransform(Parent = this.graphics.axes);
                this.graphics.homogTransfJointNext{i} = hgtransform(Parent = this.graphics.axes);
                
                % update the hgtr handles
                this.graphics.homogTransfJoint{i}.Matrix = this.G_globalJ{i};
                this.graphics.homogTransfJointNext{i}.Matrix = this.G_global{i};
                if this.lambda(i) == 0
                    this.graphics.homogTransfBody{i}.Matrix = eye(4);
                else
                    this.graphics.homogTransfBody{i}.Matrix = this.G_global{this.lambda(i)};
                end
                
                if this.h(i) ~= Inf % revolute joint or helical
                    
                    surf([ptsR(1,:); ptsR(1,:)], [ptsR(2,:); ptsR(2,:)], [ptsR(3,:); -ptsR(3,:)], 'facecolor', 'r', 'edgecolor', 'none', 'Parent', this.graphics.homogTransfJoint{i});
                    fill3(ptsR(1,:), ptsR(2,:), ptsR(3,:), 'r', 'Parent', this.graphics.homogTransfJoint{i});
                    fill3(ptsR(1,:), ptsR(2,:), -ptsR(3,:), 'r', 'Parent', this.graphics.homogTransfJoint{i});
                else   % prismatic joint
                    
                    patch('faces', faces, 'vertices', ptsP.',  'facecolor', 'green', 'Parent', this.graphics.homogTransfJoint{i});
                    patch('faces', faces, 'vertices', ptsPslide.',  'facecolor', 'green', 'Parent', this.graphics.homogTransfJointNext{i});
                end
                
                P = this.Goff(:,:,i);
                P = P(1:3, 4);
                L(:, 1) = [0;0;0];
                %                 L(:, 2) = [P(1);0;0];
                %                 L(:, 3) = [P(1);P(2);0];
                L(:, 4) = [P(1);P(2);P(3)];
                line(L(1,:), L(2,:), L(3,:), 'color', 'k', 'linewidth', 1.4, 'Parent', this.graphics.homogTransfBody{i});
                
                if isempty(this.mu{i+1}) && ~isempty(options.EEoffset)% plot also a placeholder end effector if the i-th body has no children
                    S = options.EEoffset(:,:,nn);
                    S = S(1:3, 4);
                    S = [[0;0;0], S];
                    line(S(1,:), S(2,:), S(3,:), 'color', 'b', 'linewidth', 1.4, 'Parent', this.graphics.homogTransfJointNext{i});
                    nn = nn+1;
                elseif isempty(this.mu{i+1}) && isempty(options.EEoffset)
                    line(S(1,:), S(2,:), S(3,:), 'color', 'b', 'linewidth', 1.4, 'Parent', this.graphics.homogTransfJointNext{i});
                end
            end
            
        end
        
        function updatePlot(this)
            %
            %
            %
            
            % check if plot initialization still exists
            if isempty(this.graphics) || ~isvalid(this.graphics.axes)
                error('no current plot to update exists, initialize the plot first with "treeObject.plotInit" ');
            end
            
            % updating the kinematics
            %             FWKin(this, q);
            
            for i = 1:this.Nb % joint plots
                % just update the hgtransform handles and we're done
                this.graphics.homogTransfJoint{i}.Matrix = full(this.G_globalJ{i});
                this.graphics.homogTransfJointNext{i}.Matrix = full(this.G_global{i});
                if this.lambda(i) == 0
                    this.graphics.homogTransfBody{i}.Matrix = eye(4);
                else
                    this.graphics.homogTransfBody{i}.Matrix = full(this.G_global{this.lambda(i)});
                end
                
            end
        end
        
        function [tau, F] = RNEAdyn(this, q, qd, qdd, V0, V0d, Fext)
            %
            %
            %
            
            %% forward recursive computation
            % posture computation
            this.FWKin(q);
            
            % twist computation
            twists = nan(6, this.Nb+1);
            twists(:, 1) = V0;
            twistsd(:, 1) = V0d;
            X = this.Xu;
            
            for i = 1:this.Nb
                %                 hl = this.h(i); % helical lead to identify joint type
                ADjG = adjointInv(this.G_local{i});
                ADjGV = ADjG*twists(:,this.lambda(i) + 1);
                Xd = X(:,i).*qd(i);
                twists(:, i+1) = ADjGV + Xd;
                twistsd(:, i+1) = X(:,i).*qdd(i) + ADjG*twistsd(:,this.lambda(i) + 1) + ad(ADjGV)*Xd;
            end
            %% backward recursive computation
            
            % generalized forces initialization
            F = nan(6, this.Nb);
            % actuator component
            tau = nan(1, this.Nb);
            
            for i = this.Nb:-1:1
                % computing the forces on all the direct childrens
                Fch = zeros(6,1);
                m = this.mu{i+1};
                if ~isempty(m)
                    for j = 1:length(m)
                        Fch = Fch + adjointStar(this.G_local{m(j)})*F(:, m(j));
                    end
                end
                II = this.I{i}*twistsd(:, i+1) + adStar(twists(:,i+1))*this.I{i}*twists(:, i+1); % inertial force
                F(:, i) = Fch - Fext(:, i) + II;  % equilibrium
                tau(i) = X(:, i).'*F(:, i);       % extracting the joint component
            end
        end
        
        function qdd = ABAdyn(this, q, qd, tau, V0, V0d, Fext)
            %
            % q = istantaneous joint configuration
            % qd = istantaneous joint velocity
            % V0 = istantaneous base velocity
            % V0d = istantaneous base acceleration
            % tau = istantaneous joint action (forces/moments)
            % Fext = external wrenches applied to each body (6xN matrix)
            
            cl = class(q);
            % compute posture
            this.FWKin(q);
            
            % twist computation
            twists = nan(6, this.Nb+1, cl);
            twists(:, 1) = V0;
            X = this.Xu;
            AdGinv = nan(6, 6, this.Nb, cl);
            a      = nan(6, 6, this.Nb, cl);
            % forward propagation of twists
            for i = 1:this.Nb
                
                AdGinv(:,:,i) = adjointInv(this.G_local{i});
                twists(:, i+1) = AdGinv(:,:,i)*twists(:,this.lambda(i) + 1) + X(:,i).*qd(i);
                a(:,:,i) = ad(X(:,i))*qd(i);
            end
            
            % backward propagation of projected inertia and biases
            
            % intialization of inertias
            Mtilde = nan(6, 6, this.Nb, cl);
            Mbar   = nan(6, 6, this.Nb, cl);
            btilde = nan(6, this.Nb, cl);
            b_bar  = nan(6, this.Nb, cl);
            bi     = nan(6, this.Nb, cl);
            mscalar= nan(this.Nb, cl);
            
            for i = this.Nb:-1:1
                % see literature for documentation
                % (meccanica dei robot Gabiccini ABA)
                Mtilde(:,:,i) = this.I{i}; % Mtilde = M_i as initialization
                bi(:,i) = adStar(twists(:,i+1))*this.I{i}*twists(:,i+1);
                btilde(:,i) = bi(:,i) - Fext(:, i);
                m = this.mu{i+1};
                for j = 1:length(m)
                    children = m(j);
                    ADstar = AdGinv(:,:,children).';
                    Mtilde(:,:,i) =  Mtilde(:,:,i) + ADstar*Mbar(:,:,children)*AdGinv(:,:,children);
                    btilde(:,i) =  btilde(:,i) + ADstar*( b_bar(:,children) - Mbar(:,:,children)*a(:,:,children)*AdGinv(:,:,children)*twists(:, i+1)+...
                        (1./(mscalar(children)))*Mtilde(:,:,children)*X(:,children)*tau(children) );
                end
                
                mscalar(i) = (X(:,i)).'*Mtilde(:,:,i)*X(:,i);
                P = ( eye(6) - Mtilde(:,:,i)*X(:,i)*X(:,i).'./mscalar(i) );     % projector
                Mbar(:,:,i) = P*Mtilde(:,:,i);                                      % projected inertia
                b_bar(:,i) = P*btilde(:,i);                                 % projected bias
            end
            
            % forward propagation of accelerations
            twistsD = nan(6, this.Nb+1, cl);
            twistsD(:,1) = V0d;
            qdd = nan(1,this.Nb, cl);
            for i = 1:this.Nb
                Aji = AdGinv(:,:,i);
                Vlambda = Aji*twists(:,this.lambda(i)+1); % i+1-th twist actually refers to the i-th body
                VlambdaD = Aji*twistsD(:,this.lambda(i)+1);
                qdd(i) = (tau(i) - X(:,i).'*(Mtilde(:,:,i)*(VlambdaD - a(:,:,i)*Vlambda) + btilde(:,i)))./(X(:,i).'*Mtilde(:,:,i)*X(:,i));
                twistsD(:, i+1) = VlambdaD + X(:,i).*qdd(i) - a(:,:,i)*Vlambda;
            end
        end
        
        function qdd = ABAdynCasadi(this, q, qd, tau, V0, V0d, Fext, options)
            %
            % q = istantaneous joint configuration
            % qd = istantaneous joint velocity
            % V0 = istantaneous base velocity
            % V0d = istantaneous base acceleration
            % tau = istantaneous joint action (forces/moments)
            % Fext = external wrenches applied to each body (6xN matrix)
            % C = constraint matrix
            arguments
                this
                q
                qd
                tau
                V0
                V0d
                Fext
                options.constraint = [];
            end
            
            cl = class(q);
            % compute posture
            this.FWKin(q);
            C = options.constraint.'; % lets get directly A^T
            
            % twist computation
            twists = nan(6, this.Nb+1, cl);
            twists(:, 1) = V0;
            X = this.Xu;
            AdGinv = cell(1, this.Nb);
            a      = cell(1, this.Nb);
            
            % forward propagation of twists
            for i = 1:this.Nb
                AdGinv{i} = adjointInv(this.G_local{i});
                twists(:, i+1) = AdGinv{i}*twists(:,this.lambda(i) + 1) + X(:,i).*qd(i);
                a{i} = ad(X(:,i))*qd(i);
            end
            
            % backward propagation of projected inertia and biases
            
            % intialization of inertias
            Mtilde = cell(1, this.Nb);
            mscalar = nan(1, this.Nb, cl);
            Mbar = cell(1, this.Nb);
            btilde = nan(6, this.Nb, cl);
            b_bar = nan(6, this.Nb, cl);
            bi     = nan(6, this.Nb, cl);
            
            for i = this.Nb:-1:1
                % heavy math computations, see literature for documentation
                % (meccanica dei robot Gabiccini ABA)
                Mtilde{i} = this.I{i}; % Mtilde = M_i as initialization
                bi(:,i) = adStar(twists(:,i+1))*this.I{i}*twists(:,i+1);
                btilde(:,i) = bi(:,i) - Fext(:, i);
                m = this.mu{i+1};
                
                for j = 1:length(m)
                    children = m(j);
                    ADstar = AdGinv{children}.';
                    
                    Mtilde{i} =  Mtilde{i} + ADstar*Mbar{children}*AdGinv{children};
                    
                    btilde(:,i) =  btilde(:,i) + ADstar*( b_bar(:,children) - Mbar{children}*a{children}*AdGinv{children}*twists(:, i+1)+...
                        (1./(mscalar(children)))*Mtilde{children}*X(:,children)*tau(children) );
                end
                
                mscalar(i) = (X(:,i)).'*Mtilde{i}*X(:,i);
                
                P = ( eye(6) - Mtilde{i}*X(:,i)*X(:,i).'./mscalar(i) ); % projector
                Mbar{i} = P*Mtilde{i};
                b_bar(:,i) = P*btilde(:,i);
            end
            
            % forward propagation of accelerations
            twistsD = nan(6, this.Nb+1, cl);
            twistsD(:,1) = V0d;
            qdd = nan(1,this.Nb, cl);
            for i = 1:this.Nb
                Aji = AdGinv{i};
                Vlambda = Aji*twists(:,this.lambda(i)+1); % i+1-th twist actually refers to the i-th body
                VlambdaD = Aji*twistsD(:,this.lambda(i)+1);
                qdd(i) = (tau(i) - X(:,i).'*(Mtilde{i}*(VlambdaD - a{i}*Vlambda) + btilde(:,i)))./(X(:,i).'*Mtilde{i}*X(:,i));
                twistsD(:, i+1) = VlambdaD + X(:,i).*qdd(i) - a{i}*Vlambda;
            end
        end
        
        function computeCasadi_state_space_dyn(this, options)
            %
            % returns dyn eq. in state space form with ABA
            %
            arguments
                this
                options.constraint = []
                options.error = [];
            end
            
            import casadi.*
            
            q = SX.sym('q', this.Nb, 1);
            qd = SX.sym('qd', this.Nb, 1);
            tau = SX.sym('tau', this.Nb, 1);
            V0 = SX.sym('V0', 6, 1);
            V0d = SX.sym('V0d', 6, 1);
            Fext = SX.sym('Fe', 6, this.Nb);
            if ~isempty(options.constraint)
                C = options.constraint(q);
            else
                C = [];
            end
            [ABA, Cbar] = this.ABAdynCasadi(q, qd, tau, V0, V0d, Fext, 'constraint', C);
            ABA = ABA.';
            qdd = ABA;
            delta = options.error(q);
            if ~isempty(options.constraint)
                CA = Cbar*pinv(C*Cbar);
                Cdq = jacobian(C,q);
                Cdq = Cdq(:,1).*qd(1) + Cdq(:,2).*qd(2) + Cdq(:,3).*qd(3) + Cdq(:,4).*qd(4);
                Cdq = reshape(Cdq, size(C));
                qdd = (eye(this.Nb) - CA*C)*ABA - CA*(Cdq*qd);
                %                 lamb = -pinv(Cbar)*(pinv(C)*Cdq*qd + ABA);
            end
            xdot = [qd; qdd];
            this.xdot_sym = Function('xdot', {[q; qd], [tau; V0; V0d; Fext(:)]}, {xdot});
            %             this.Cbar_fun = Function('C', {q}, {Cbar});
            %             this.pinvC= Function('C', {q}, {pinv(C)});
            %             this.lam = Function('lam', {[q;qd], [tau; V0; V0d; Fext(:)]}, {lamb});
            %             this.CCRdq = Function('CCR', {[q;qd], [tau; V0; V0d]}, {CCR*pinv(C)*(Cdq*qd)});
        end
        
        function compileDyn(this)
            
            opts = struct('main', false,...
                'mex', true);
            this.xdot_sym.generate('treeDyn_mex.c',opts);
            mex treeDyn_mex.c -largeArrayDims
            
        end
        
    end
    
    
end