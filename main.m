addpath('C:\Program Files\Mosek\9.3\toolbox\R2015a')
addpath(genpath('YALMIP'))

%%
clc

g = 9.81;
Mc = 1;
l = 1;
dMp = ureal('dMp',0,'Range',[-0.9,0.9]);

% Bring the plant in LFR form
A = [0, 1, 0, 0;
     0, 0, -1.1*g/Mc, 0;
     0, 0, 0, 1;
     0, 0, g*(Mc+1.1)/(Mc*l), 0];
B = [0,0;
     -g/Mc,1/Mc;
     0,0;
     g/(Mc*l),-1/Mc];
C = [0, 0, 1, 0;
     eye(4)];
D = zeros(5,2);

G = ss(A,B,C,D);
G.OutputName = {'z','x1','x2','x3','x4'};
G.InputName = {'v','u'};

smb = sumblk('u = uK + p');

K = ss([7.40 14.96 125.82 27.73]);
K.InputName = {'x1','x2','x3','x4'};
K.OutputName = {'uK'};

Gcl = connect(G,K,smb,{'v','p'},{'z','x1','x2','x3','x4'})

%%
yalmip('clear')

iqc1 = IQC_RRP(0.9,1,9,-10);
delta = Delta(1,1,1);
delta = delta.attachIQC(iqc1,1);
delta = delta.compile();

iqcp = IQC_L2G(0.1416,4,1);
deltap = Delta(4,1,1);
deltap = deltap.attachIQC(iqcp,1);
deltap = deltap.compile();

problem = IQCProblem(Gcl(1,1),delta,Gcl(2:5,1),Gcl(1,2),Gcl(2:5,2),deltap);
problem = problem.compile();
problem.solve()

%% Experimental valid.

Gcl_ = lft(dMp, Gcl);
N = 1e3;
L2norms = hinfnorm(usample(Gcl_,N));
mean(L2norms)
max(L2norms)
