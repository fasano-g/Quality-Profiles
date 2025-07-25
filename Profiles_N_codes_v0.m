%=========================================================================
% This Matlab code can be adopted to benchmark solvers for optimization
% problems. In particular, this code can be used to accomplish the 
% following three specific purposes:
%
%   1) Computing PERFORMANCE PROFILES (E.D.Dolan, J.More',
%      "Benchmarking optimization software with performance profiles",
%      Mathematical Programming 91, pp. 201-213, 2002)
%
%   2) Computing single/multiple DATA PROFILES (J.More', S.Wild,
%      Benchmarking Derivative-Free Optimization Algorithms, SIAM
%      J. Optimization 20, pp. 172-191, 2009)
%
%   3) Computing QUALITY PROFILES (G.Fasano, C.Piermarini, M.Roma,
%      Benchmarking Optimization Algorithms with Quality Profiles,
%      submitted to Mathemtical Programming Computations, 2025)
%
% All the above mentioned profiles can be plotted using either linear or
% logarithmic scale on the abscissa axis. Hereafter we describe the steps 
% for each of the above three tasks, including settings of the required
% parameters in the code.
%
%------------------------------------------------------------------
% TO RUN PERFORMANCE PROFILES the next quantities must be set, all
% the remaining quantities in the code can be ignored. The code will
% automatically detect both the number of test problems and solvers
% from the provided data:
%------------------------------------------------------------------
%
%   profiles = must be set equal to the value 1
%
%   limit = represents the largest value on the abscissa axis
%           used to plot the performance profiles
%   
%   logscaling = (true) for a logarithmic scale on the abscissa axis
%                (false) for a linear scale on the abscissa axis
%
%   MM = represents a PxS real matrix, being P the number of test
%        problems and S the number of solvers to be compared. Each entry 
%        of MM represents the performance of the solver 's' on the test 
%        problem 'p'. In case the solver 's' has a failure on the test 
%        problem 'p', then provide a very large value in the 
%        corresponding position of the matrix
%
% The code provides warnings and error messages in case data is
% inconsistent (e.g. some entries of MM are missing or too small, so that 
% performance profiles cannot be computed)
%
%
%------------------------------------------------------------------
% TO RUN DATA PROFILES the next quantities must be set, all the 
% remaining quantities in the code can be ignored. The code will
% automatically detect both the number of test problems and solvers
% from the provided data:
%------------------------------------------------------------------
%
%   profiles = must be set equal to the value 2
%
%   logscaling = (true) for a logarithmic scale on the abscissa axis
%                (false) for a linear scale on the abscissa axis
%
%   FF = represents a PxSxT real matrix, being P the number of test
%        problems, S the number of solvers to be compared and T the
%        levels of precision (namely the entries of the vector 'tau'. 
%        Each entry of FF represents the number of function evaluations 
%        necessary to the solver 's' to find a solution of the test 
%        problem 'p', within the precision level 't'. In case the solver 
%        's' has a failure on the test problem 'p' and for the precision 
%        level 't', then provide a very large value in the corresponding 
%        position of the matrix
%
%  tau = represents the precision associated to data profile(s). It can 
%        be a value in (0,1], or a vector (even number of entries, being 
%        each entry in (0,1])
%
%   TT = represents a Px(1+S) matrix. The first column (real number >= 1) 
%        gives the dimension (number of unknowns) of each of the P test 
%        problems, the second column (any real number) is the value 
%        of the objective function at the initial iterate x_0, while each
%        of the remaining S columns reports the optimal value (any real 
%        number) of the objective function found by the solver 's' on the 
%        test problem 'p'
%  
% The code provides warnings and error messages in case data is
% inconsistent (e.g. the solution found by a solver is larger than the
% initial objective function value, so that data profiles cannot be 
% computed)
%
%
%------------------------------------------------------------------
% TO RUN QUALITY PROFILES the next quantities must be set, all the 
% remaining quantities in the code can be ignored. The code will
% automatically detect both the number of test problems and solvers
% from the provided data:
%------------------------------------------------------------------
%
%   profiles = must be set equal to the value 3
%
%   logscaling = (true) for a logarithmic scale on the abscissa axis
%                (false) for a linear scale on the abscissa axis
%
%   TT = represents a Px(2+S) matrix. The first column gives the 
%        dimension (number of unknowns) of each of the P test problems.
%        It is ignored when building Quality Profiles, but is used to 
%        build Data Profiles. The second column (any real number) is the 
%        value of the objective function at the initial iterate x_0, 
%        while each of the remaining S columns reports the optimal value 
%        (any real number) of the objective function found by the solver 
%        's' on the test problem 'p'
% 
%   r1 = represents a zooming parameter for the 'abscisa' axis in a 
%        Quality Profile.
%        r1 = 1 implies no zoomin
%        r1 < 1 implies zooming nearby the abscissa value 1
%        r1 > 1 implies zooming nearby the abscissa value 0 
%
%   r2 = represents a zooming parameter for the 'ordinate' axis in a 
%        Quality Profile.
%        r2 = 1 implies no zoomin
%        r2 < 1 implies zooming nearby the ordinate value 1
%        r2 > 1 implies zooming nearby the ordinate value 0 
%
%   strong = represents the choice between weak and strong quality 
%            profiles, i.e.
%        strong = true implies that strong quality profiles are generated
%        strong = false implies that weak quality profiles are generated
%  
% The code provides warnings and error messages in case data is
% inconsistent (e.g. the solution found by a solver is larger than the
% initial objective function value, so that quality profiles cannot be 
% computed)


%==================================================================
%--------------------------------
% Clear screen and all variables
%--------------------------------
clc; clf;
close all;
clear all;

%-----------------------
% Start to compute time
%-----------------------
tic;

%-----------------------------
% Select the type of profiles
%-----------------------------
%profiles = 1;   % Performance profiles
%profiles = 2;   % Data Profiles
profiles = 3;   % Quality Profiles

% ------------------------------------------------------------------------
% Set the 'limit' for the abscissa parameter of Performance Profiles
% Set the required precision values for Data Profiles
% Set the required scaling parameters 'r1' and 'r2' for Quality Profiles
% ------------------------------------------------------------------------
limit = 15;         % Upper bound for the abscissa value of the 
                    % Performance Profiles (user's choice)
                   
%tau = [0.01];      % Precision level in case only one data profile 
                    % is sought: the user can select any number
                    % in (0,1]
tau = [0.1 
       0.01 
       0.001 
       0.0001 
       0.00001 
       0.000001];   % Precision levels in case multiple (even number) 
                    % of data profiles are sought: any entry must be
                    % selected by the user in (0,1], and will be
                    % associated to a different data profile

r1 = 1; r2=1;       % Scaling parameters for Quality Profiles. Each of
                    % these parameters must be a real positive number       

%strong = true;     % The user creates 'Strong' Quality Profiles
strong = false;     % The user creates 'Weak' Quality Profiles

%-----------------------------------------
% Decide the scaling of the abscissa axis
% for all the types of profiles
%-----------------------------------------
logscaling = true;      % log10 scaling for the abscissa axis
%logscaling = false;      % no log scaling for the abscissa axis


% -----------------------------------------------------------------------
% Matrix used to plot Performance Profiles. Its dimensions are PxS, being 
%       P = # test Problems 
%       S = # of Solvers
% Each entry of MM represents the 'performance' of the solve 's' on the
% problem 'p'. Hence, in case of the Performance Profiles it represents 
% (e.g.) the number of iterations, the number of function evaluations, 
% etc. In case the solver 's' has a failure on the test  problem 'p', 
% then provide a very large value in the corresponding position of the 
% matrix
% -----------------------------------------------------------------------
% MM = [];
 MM = 10 + (1000-10)*rand(100,3);     % Just an example for matrix MM
                                      % including 100 test problems
                                      % and 3 solvers

% -----------------------------------------------------------------------
% Matrix used to plot Data Profiles. Its dimensions are PxSxT, being 
%       P = # test Problems 
%       S = # of Solvers
%       T = # of precision levels
% Each entry of FF represents the 'performance' of the solve 's' on the
% problem 'p'. Hence, in case of the Data Profiles it represents 
% the number of function evaluations necessary to solve the problem 'p'
% by the solver 's' at the given level of precision 't'. In case the 
% solver 's' has a failure on the test problem 'p' when the precision 
% level 't' is sought, then provide a very large value in the 
% corresponding position of the matrix
% -----------------------------------------------------------------------
% FF = [];
 FF = zeros(100,3,6);
 FF(:,:,1) = 1 + (10-1)*rand(100,3);   % Just an example for matrix FF
 FF(:,:,2) = 1 + (100-1)*rand(100,3);   % including 100 test problems, 
 FF(:,:,3) = 10 + (1000-10)*rand(100,3);   % 3 solvers and 6 levels of
 FF(:,:,4) = 10 + (5000-10)*rand(100,3);   % precision
 FF(:,:,5) = 100 + (10000-100)*rand(100,3);
 FF(:,:,6) = 100 + (20000-100)*rand(100,3);

% -----------------------------------------------------------------
% Provide the next matrix TT whose columns are (respectively): 
%   1st column  = the dimensions of the P test problems 
%   2nd column  = the objective function value at the initial iterate
%   s-th column = (s > 2) the objective function value at the optimal 
%                 point detected by the solver 's'
% This data will be uniquely used to plot either Data Profiles (just the
% 1st column is used) or Quality Profiles (all the columns apart from 
% the 1st one are used), and will be ignored if Performance Profiles 
% are sought
% -----------------------------------------------------------------
% TT = [];
 TT = [randi([1,50],100,1)  (-1 + 2*rand(100,1)) ...
       (-2 + 2*rand(100,1)) (-2.01 + 2*rand(100,1)) (-3.0 + 2*rand(100,1))];
       % Just an example for matrix TT including at most 50 unknowns, 
       % 100 test problems and 3 solvers

%--------------------------------------------
% Assign the dimensions of the test problems 
%--------------------------------------------
n_p = TT(:,1);

%------------------------------------------------------------------
% Assign the value of the objective function at the first iterate, 
% for each test problem 
%------------------------------------------------------------------
f_0 = TT(:,2);

%------------------------------------
% Print a summary of the user's data
%------------------------------------
if (profiles ==1)
    [P,S] = size(MM);
else
    [P,S] = size(TT);
    S = S - 2;
end
fprintf('\n====================================================== \n');
fprintf('  Number of TEST PROBLEMS ------>  %i  \n',P);
fprintf('  Number of SOLVERS       ------>  %i  \n',S);
if (profiles ==1)
    fprintf('  Performance Profiles: user abscissa limit = %8.3f  \n',limit);
elseif (profiles ==2)
    fprintf('  Data Profiles: precision value tau = %d  \n',tau);
elseif(profiles ==3)
    fprintf('  Quality Profiles: scaling parameter (x) r1 =%6.2f  \n',r1);
    fprintf('  Quality Profiles: scaling parameter (y) r2 =%6.2f  \n',r2);
end
fprintf('====================================================== \n \n');

% ---------------------------
% Check for data consistency
% ---------------------------
[dim_f_0,dummy1] = size(f_0);
[dim_n_p,dummy2] = size(n_p);
if ((profiles == 2 || profiles == 3) && (min(n_p) < 1))
    fprintf('  ERROR: Wrong Data (problems dimension) !!!');
    return;
end
if (profiles == 2)
    [a,b] = size(tau);
    if (b ~= 1) && (mod(b,2) ~= 0) 
        fprintf('  ERROR: odd number of required Data Profiles !!!');
        return;
    end
end

% --------------------------------------------------------
% In case of "Quality Profiles" then
%       - select the stepsize for the parameter 
%         representing the precision 'tau'
%       - update the vector f_0
% --------------------------------------------------------
if (profiles ==3)
    intervals = floor(max([1000 1000*r1 1000/r1]));
    stepsize = 1/intervals;
    tau = zeros(intervals+1,1);
    if (strong == true)
        for i=1:P
            f_0(i) = max(TT(i,3:(3+S-1)));
        end
    end
end

% -----------------------------------------------------------------
% Compute the best value among the solvers, for each test problem. 
% The best value is only used for Quality Profiles
% -----------------------------------------------------------------
if (profiles == 3)
    f_star = zeros(P,1);
    for p=1:P
        f_star(p) = min(TT(p,3:(3+S-1)));
    end
end

% ------------------------------
% Compute the required profiles
% ------------------------------
% ----------------------
%  PERFORMANCE PROFILES
% ----------------------
if (profiles == 1)  
    % Here we create the Performance Profiles "RHO_s(tau)"
    % for any solver 's'. The parameter 'tau' has no dimension and
    % represents the abscissa parameter in Performance Profiles
    RHO = zeros(S,1000);
    tps = MM;
    rps = zeros(P,S);
    r_M = 0;        % initialize the maximum performance ratio
    min_tps = zeros(P,1);
    for p=1:P
        min_tps(p) =  min(MM(p,:));
        for s=1:S
           if (abs(tps(p,s)) <= 1.0e-10)
                fprintf('WARNING: very small data in Performance Profiles !!!');
                pause;
           end
           rps(p,s) =  tps(p,s) / min_tps(p);
           r_M = max(r_M , rps(p,s));
        end    
    end
    r_M = 1.1*min(r_M,limit);   % compute an upper bound to 'r_M' so that 
                                % the plot of the Performance Profiles
                                % will correctly fit the user's interval 
    tau = linspace(1, r_M, 1000);
    for t=1:1000
        for s=1:S
            counter = 0;
            for p=1:P
                if (rps(p,s) <= tau(t)) 
                    counter = counter + 1;
                end 
            end
            RHO(s,t) = counter/P;
        end
    end

% ---------------
%  DATA PROFILES
% ---------------
elseif (profiles == 2)  
    % Here we create the Data Profiles "d_s(alpha)" for any 
    % solver 's'. The parameter 'alpha' has no dimension and
    % represents a multiple of the cost of one simplex gradient 
    % evaluation, i.e. n_p + 1, being 'n_p' the number of 
    % variables of the test problem 'p'
    [entries,dummy] = size(tau);
    d = zeros(S,5000,entries);
    tps = FF;
    i_M = ones(entries,1);
    
    % Compute the maximum performance ratio
    delta = zeros(entries);
    r_max = zeros(entries);
    for t=1:entries
        for p=1:P
             for s=1:S
                 r_max(t) = max(r_max(t) , tps(p,s,t)/(n_p(p)+1));
             end
        end 
    end
    for t=1:entries
        r_max(t) = 1.05*r_max(t);     % Compute an upper bound on 'r_M'
        delta(t) = r_max(t)/5000;     % Compute an elementary interval
    end
    for t=1:entries         % Compute Data Profiles for each precision value
        for s=1:S
            for i=1:5000    % Increase the abscissa parameter 'alpha'
                counter = 0;
                for p=1:P
                    if (tps(p,s,t)/(n_p(p)+1) <= i*delta(t))
                         counter = counter + 1;
                    end
                end
                d(s,i,t) = counter/P;
            end    
        end
    end
    
    % Compute the largest abscissa such that, for a given "tau", the plot
    % of "data profiles" for the solvers does not change (this is useful in
    % order to scale the abscissa axis of the different data profiles, 
    % depending on "tau")
    for t=1:entries
        for i=1:4999
            i_found = false;
            for s=1:S
                if ( abs(d(s,5000-i+1,t) - d(s,5000-i,t)) > 1.0e-5 )
                    i_M(t) = floor( min((5000-i+1)*1.1,5000) );
                    i_found = true;
                    break
                end
            end
            if i_found
                break
            end
        end
    end
   
% ------------------
%  QUALITY PROFILES
% ------------------
else
    % Here we create the Quality Profiles "Qs(tau)". Each row of the 
    % matrix "Q" is associated with a different precision level.
    % Each column of the matrix "Q" is associated with a different solver
    Q = zeros(intervals+1,S);
    for s=1:S
        for t=1:(intervals+1)
            tau(t) = ((t-1)*stepsize)^(r1);
            counter = 0;
            for w=1:P
                if ((TT(w,2+s) - f_star(w)) <= tau(t)*(f_0(w) - f_star(w)))
                    counter = counter + 1;
                end
            end
            % Scaling, through the parameter "r2", of the Quality Profiles 
            Q(t,s) = (counter/P)^(r2);
        end
    end
end

%----------------------------%
% Plot the required profiles %
%----------------------------%
%---------------------------%
% Plot PERFORMANCE PROFILES %
%---------------------------%
if (profiles ==1)     
    figure1 = figure;
        % colors = [linspace(1, 0, S)', rand(S,1), linspace(0, 1, S)'];
        colors = [linspace(1, 0, S)', zeros(S,1), linspace(0, 1, S)'];
        if (logscaling)
            for s=1:S
                semilogx(tau,RHO(s,:),':','Color',colors(s,:),'LineWidth',1.5);
                hold on
            end
        else
            for s=1:S
                plot(tau,RHO(s,:),':','Color',colors(s,:),'LineWidth',1.5);
                hold on
            end
        end
        title(['\it {\color{red}Performance Profiles:} # solvers = ',num2str(S), ...
               '; # test problems = ',num2str(P), ...
               '; \tau--range = [1 , ',num2str(r_M),']'],'FontName','Baskerville Old Face', ...
               'Fontsize',35);
        grid on
        xlabel('\it \tau','FontWeight','bold','FontSize',35); 
        ylabel(' \it \rho_s(\tau)','FontWeight','bold','FontSize',30);
        set(gca, 'FontSize', 15);
        % Here you can specify the name of the solvers in the vector 
        % "SolversNames": the following example is provided 
        SolversNames = strseq(' Solver - ',1:S);
        legend1 = legend(SolversNames);
        set(legend1,'Location','southeast','FontSize',20);
        hold off
        
% --------------------
%  Plot DATA PROFILES
% --------------------
elseif (profiles ==2)  
    figure2 = figure;
        % colors = [linspace(1, 0, S)', rand(S,1), linspace(0, 1, S)'];
        colors = [linspace(1, 0, S)', zeros(S,1), linspace(0, 1, S)'];
        % ---------------------------
        %  create 'entries' subplots
        % ---------------------------
        if (entries ==1)
                if (logscaling)
                    for s=1:S
                        semilogx((1:i_M(1))*delta(1),d(s,1:i_M(1),1),':','Color',colors(s,:),'LineWidth',1.5);
                        hold on
                    end
                else
                    for s=1:S
                        plot((1:i_M(1))*delta(1),d(s,1:i_M(1),1),':','Color',colors(s,:),'LineWidth',1.5);
                        hold on
                    end
                end
                title(['{\color{red}Data Profile:} |S| = ', ...
                        num2str(S),';  |P| = ', ...
                        num2str(P),';  \tau = ',num2str(tau(1))],'FontName','Baskerville Old Face', ...
                        'Fontsize',35);
                grid on
                xlabel('\it \alpha','FontWeight','bold','FontSize',35); 
                ylabel(' \it d_s(\alpha)','FontWeight','bold','FontSize',30);
                set(gca, 'FontSize', 15);
                % Here you can specify the name of the solvers in the 
                % vector "SolversNames": the following example is provided
                SolversNames = strseq(' Solver - ',1:S);
                legend2 = legend(SolversNames);
                set(legend2,'Location','southeast','FontSize',10);
        else
            for j=1:entries
                subplot(entries/2,2,j)
                if (logscaling)
                    for s=1:S
                        semilogx((1:i_M(j))*delta(j),d(s,1:i_M(j),j),':','Color',colors(s,:),'LineWidth',1.5);
                        hold on
                    end
                else
                    for s=1:S
                        plot((1:i_M(j))*delta(j),d(s,1:i_M(j),j),':','Color',colors(s,:),'LineWidth',1.5);
                        hold on
                    end
                end
                title(['{\color{red}Data Profile:} |S| = ', ...
                        num2str(S),';  |P| = ', ...
                        num2str(P),';  \tau = ',num2str(tau(j))],'FontName','Baskerville Old Face', ...
                        'Fontsize',20);
                grid on
                xlabel('\it \alpha','FontWeight','bold','FontSize',20); 
                ylabel(' \it d_s(\alpha)','FontWeight','bold','FontSize',20);
                set(gca, 'FontSize', 13);
                % Here you can specify the name of the solvers in the 
                % vector "SolversNames": the following example is provided
                SolversNames = strseq(' Solver - ',1:S);
                legend2 = legend(SolversNames);
                set(legend2,'Location','southeast','FontSize',10);
            end
        end

% -----------------------
%  Plot QUALITY PROFILES
% -----------------------
else
    figure3 = figure;
        % colors = [linspace(1, 0, S)', rand(S,1), linspace(0, 1, S)'];
        colors = [linspace(1, 0, S)', zeros(S,1), linspace(0, 1, S)'];
        if (logscaling)
            for s=1:S
                semilogx(tau(1:(intervals+1)),Q(:,s),':','Color',colors(s,:),'LineWidth',2.5);
                %loglog(tau(1:(intervals+1)),Q(:,s),':','Color',colors(s,:),'LineWidth',2.5);
                hold on
            end
        else
            for s=1:S
                plot(tau(1:intervals+1),Q(:,s),':','Color',colors(s,:),'LineWidth',2.5);
                hold on
            end
            [rows,cols] = size((0:stepsize:1)');
            plot(0.001*ones(rows,1),(0:stepsize:1)','-k',...
                 0.01*ones(rows,1),(0:stepsize:1)','-k',...
                 0.1*ones(rows,1),(0:stepsize:1)','-k','LineWidth',0.5);
        end
        if (strong)
            title([' {\color{red}(Strong) Quality Profiles:} |S| = ', ...
                            num2str(S),';  |P| = ', ...
                            num2str(P),';  x-scaling (r_1 = ',num2str(r1), ...
                   ');  y-scaling (r_2 = ',num2str(r2),')'],'FontName','Baskerville Old Face', ...
                   'Fontsize',35);
        else
            title([' {\color{red}(Weak) Quality Profiles:} |S| = ', ...
                            num2str(S),';  |P| = ', ...
                            num2str(P),';  x-scaling (r_1 = ',num2str(r1), ...
                   ');  y-scaling (r_2 = ',num2str(r2),')'],'FontName','Baskerville Old Face', ...
                   'Fontsize',35);
        end
        grid on
        xlabel('0 \leq \tau^{r_1} \leq 1','FontWeight','bold','FontSize',35); 
        ylabel(' [ Q_s(\tau) ]^{r_2}','FontWeight','bold','FontSize',35);
        set(gca, 'FontSize', 15);
        % Here you can specify the name of the solvers in the 
        % vector "SolversNames": the following example is provided
        SolversNames = strseq(' Solver - ',1:S);
        legend3 = legend(SolversNames);
        set(legend3,'Location','southeast','FontSize',15);
        hold off
end

   
%----------------------
% Stop to compute time
%----------------------
toc;
