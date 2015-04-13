% This script runs all 3 demos which can each be run
% individually with different inputs after editing.

clear; close all;

RandStream.setDefaultStream(RandStream('mt19937ar','seed',31415));

disp('Solver: 1 = LMaFit, 2 = IALM, 1:2 both')
%Solver = input(' Solvers = ');
% Solver = 1:2;
Solver = 1;

t0 = cputime;

disp('+++ Results from demo_chkb.m +++')
m = 256; impulse = .75;
demo_chkb(Solver,m,impulse); 
set(gcf,'Position',[0 400 1200 400]); drawnow

disp('+++ Results from demo_spim.m +++')
rankD = 8; impulse = .08;
demo_spim(Solver,rankD,impulse);

disp('+++ Results from demo_rand.m +++')
m = 400; n = m; rankD = 60; impulse = 0.15; tol = 1e-6;
demo_rand(Solver,m,n,rankD,impulse,tol);

fprintf('Total CPU: %g seconds\n\n\n',cputime-t0)