%% closed-loop ABA

%% GRAPH BUILD
p = [0,0,1,2]; % predecessor array
s = [1,2,3,4]; % successor array

% create the graph (Remember matlab indexing origin starts from 1)
G = digraph(p+1, s+1);

jointLabels = cell(length(s), 1);
nodeLabels = cell(length(s), 1);
nodeLabels{1} = '0';
for i = 1:length(s)
    jointLabels{i} = ['arc' num2str(i)];
    nodeLabels{i+1} = num2str(i);
end

% plot the kinematic graph (useful to remember the joint and link
% enumeration)
figure;
kinematicTreeGraph = plot(G, 'NodeLabel', nodeLabels, 'EdgeLabel', jointLabels);
set(gca, 'YDir','reverse');

%% LINK/FRAMES OFFSET TRANSFORMATIONS

% start defining 2D offset displacement (x and y coords) between each
% reference frame (LOCAL POE param.)
dx = [0, 0, -1,  1, 2, 0];
dy = [0, 0,  1,  1, 1, 1];

% build the transformation matrices from the displacement values.
% NOTE! there are 6 offset transformations (there is an additional end-eff. transformation
% for each leaf)

gT = nan(4,4,length(dx));
gT(:,:,1) = eye(4);
for i = 2:length(dx)
    
    gT(:,:,i) = Ttx(dx(i))*Tty(dy(i));
    
end

% joint twists
X = [[0;0;0;0;0;1],[0;0;0;0;0;1],[0;0;0;0;0;1],[0;0;0;0;0;1]]; 

% helical lead
h(1:4) = 0; % all revolute

% kinematic tree object
tree = kinematicTree(p, s, gT, X, h);

% assign some inertia to links
m = 5;
for i = 1:tree.Nb
    tree.I{i} = eye(6).*m;
    tree.I{i}(4:5) = 1/2.*tree.I{i}(4:5).*0.01;
    tree.I{i}(6) = 1./12*tree.I{i}(6).*4;
    tree.I{i} = adjointStar(gT(:,:,i+2))*tree.I{i}*adjointInv(gT(:,:,i+2));
end

GEE3 = gT(:,:,end-1);
GEE4 = gT(:,:,end);

%% CASADI VARIABLES
import casadi.*
q_sym = SX.sym('q', tree.Nb, 1);
qd_sym = SX.sym('qd', tree.Nb, 1);
V0_sym = SX.sym('V0', 6, 1);
V0d_sym = SX.sym('V0d', 6, 1);
Fext_sym = SX.sym('F', 6, tree.Nb);
tau_sym = SX.sym('tau', tree.Nb, 1);

%% forward kinematic
tree.FWKin(q_sym);

%% KINEAMTIC CONSTRAINTS
% body jacobians:
% root to third link
[J3, J3tot] = tree.bodyJacobian(3);
% root to fourth link
[~, J4tot] = tree.bodyJacobian(4);

% change pole to extremity of e-e links
J3totP = twistPole(-[dx(5);dy(5);0])*J3tot; %
J4totP = twistPole(-[dx(6);dy(6);0])*J4tot; %

% change to spatial coordinates
G3 = tree.G_global{3};
G4 = tree.G_global{4};
R03 = tree.G_global{3}(1:3, 1:3);
R04 = tree.G_global{4}(1:3, 1:3);

L03 = twistTransf(R03);
L04 = twistTransf(R04);

% pfaffian matrix
A = L03*J3totP - L04*J4totP;

% extract only the first two rows (2D linear vel. of the e-e's)
A = A([1,2],:);
Afun = Function('A', {q_sym}, {A});

% build the state-space dynamic system with the augmented Lagrangian method
% and with the kinematic constraint
tree.computeCasadi_state_space_dyn('constraint', Afun);

% compile the dynamic eq. function
tree.compileDyn();

%% SIMULATION

% intial conditions
q0 = zeros(tree.Nb, 1);
q0d = zeros(tree.Nb, 1);
Fext = zeros(6,tree.Nb);
V0 = zeros(6,1);
V0d = [0;9.81;0;0;0;0];
tau = zeros(tree.Nb, 1);
tree.FWKin(q0)
tree.plotInit('EEoffset', gT(:,:,end-1:end));
axis off
xlim([-4 4])
ylim([-4 4])
x = [q0;q0d];

nFrames = 10000;
% vidObj = VideoWriter('quadrilatero.avi');
% vidObj.Quality = 100;
% vidObj.FrameRate = 60;
% open(vidObj);
for i = 1:nFrames
    % uodate the posture
    tree.FWKin(x(1:tree.Nb));
    % update graphics
    tree.updatePlot();
    % update plot with either pause(seconds) or drawnow()
    pause(0.01)
    % add some damping
    tau = -5.*x(tree.Nb+1:end);
    % compute an integration step
    x = full(RK4_step(x, [tau;V0;V0d;Fext(:)], @treeDyn_mex, 0.005));
%    writeVideo(vidObj, getframe(gca));
    

end

% close(vidObj);