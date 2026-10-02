clear all
close all

%% Parameters
mu = 0.1;
b=0.14286;
a=0.28571;

% Vary lambda
lambda_values = [linspace(0,10,200) linspace(10,12,500) linspace(12,15,200)];

n_l = length(lambda_values);
%% Initial guesses
xGuess = linspace(0,60,10);
yGuess = linspace(0,2,10);

[X0,Y0] = meshgrid(xGuess,yGuess);
initialGuesses = [X0(:),Y0(:)];

%% Solver options
unique_solutions1 = [];
unique_solutions2 = [];
unique_solutions4 = [];
unique_solutions8 = [];

R0_fp_stable  = NaN(n_l,1); 
R0_2_stable  = NaN(n_l,1);
R0_4_stable  = NaN(n_l,1);
R0_8_stable = NaN(n_l,1);

R0_fp_unstable  = NaN(n_l,1); 
R0_2_unstable  = NaN(n_l,1);
R0_4_unstable  = NaN(n_l,1);
R0_8_unstable = NaN(n_l,1);

options = optimoptions('fsolve',...
    'Display','off',...
    'FunctionTolerance',1e-12,...
    'StepTolerance',1e-12);

for k = 1:length(lambda_values)

    lambda = lambda_values(k);
    k
    
    solutions1 = [];
    solutions2 = [];
    solutions4 = [];
    solutions8 = [];

    for j = 1:size(initialGuesses,1)
        guess = initialGuesses(j,:)';
        %% Find equilibria
    
        [x_sol1, fval1, exitflag1] = fsolve(@(X) mySystemParasitoid(X,b,a,lambda,mu),guess,options);
        [x_sol2, fval2, exitflag2] = fsolve(@(X) mySystemParasitoid2(X,b,a,lambda,mu),guess,options);
        [x_sol4, fval4, exitflag4] = fsolve(@(X) mySystemParasitoid4(X,b,a,lambda,mu),guess,options);
        [x_sol8, fval8, exitflag8] = fsolve(@(X) mySystemParasitoid8(X,b,a,lambda,mu),guess,options);
    
       % Only keep valid solutions
        if exitflag1 > 0  && all(x_sol1>=0) && norm(fval1) < 1e-4 % exitflag>0 means fsovle converged; all(x_sol>=0) keeps only non-negative solutions, norm(fval)<1e-6 ensures small residual
            solutions1 = [solutions1; x_sol1(:)']; % stores solultions  
        end
        if exitflag2 > 0  && all(x_sol2>=0) && norm(fval2) < 1e-8
            solutions2 = [solutions2; x_sol2(:)']; 
    
            x1 = x_sol2;
            x2 = parasitoidMap(x1,b,a,lambda,mu);
    
            solutions2 = [solutions2; x2(:)'];
        end
    
        if exitflag4 > 0  && all(x_sol4>=0) && norm(fval4) < 1e-8
            solutions4 = [solutions4; x_sol4(:)']; 
    
            x1 = x_sol4;
            x2 = parasitoidMap(x1,b,a,lambda,mu);
            x3 = parasitoidMap(x2,b,a,lambda,mu);
            x4 = parasitoidMap(x3,b,a,lambda,mu);
    
            solutions4 = [solutions4; x4(:)'];
        end
    
        if exitflag8 > 0  && all(x_sol8>=0) && norm(fval8) < 1e-8
            solutions8 = [solutions8; x_sol8(:)']; 
    
            x1 = x_sol8;
            x2 = parasitoidMap(x1,b,a,lambda,mu);
            x3 = parasitoidMap(x2,b,a,lambda,mu);
            x4 = parasitoidMap(x3,b,a,lambda,mu);
            x5 = parasitoidMap(x4,b,a,lambda,mu);
            x6 = parasitoidMap(x5,b,a,lambda,mu);
            x7 = parasitoidMap(x6,b,a,lambda,mu);
            x8 = parasitoidMap(x7,b,a,lambda,mu);
    
            solutions8 = [solutions8; x8(:)'];
        end
    
    end
    %% Unique solutions
    unique_solutions1 = unique(round(solutions1, 1), 'rows');      
    unique_solutions2 = unique(round(solutions2, 2), 'rows');
    unique_solutions4 = unique(round(solutions4, 2), 'rows');
    unique_solutions8 = unique(round(solutions8, 2), 'rows');

    % Remove duplicates

    if ~isempty(unique_solutions1)
        ref_fp = unique_solutions1;
    else
        ref_fp = [];                  % no reference solutions
    end

    % Now compare system 2 to system 1
    if ~isempty(unique_solutions2)

        if isempty(ref_fp)
            % Nothing to compare: keep all system 2 solutions
            nonmatching_fp2{k} = unique_solutions2;

        else
            % Compute differences between each system 2 solution and each system 1 solution
            % Using pdist2 so all pairwise distances are computed
            D = pdist2(unique_solutions2, ref_fp);

            % Keep system 2 points that do NOT match ANY system 1 point
            nonmatch_mask = min(D, [], 2) > 1e-2;

            nonmatching_fp2{k} = unique_solutions2(nonmatch_mask, :);

        end

    else
        nonmatching_fp2{k} = [];
    end

    % Build reference set from system 1 and system 2
    ref_all = [];
    
    if ~isempty(unique_solutions1)
        ref_all = [ref_all; unique_solutions1];
    end
    
    if ~isempty(unique_solutions2)
        ref_all = [ref_all; unique_solutions2];
    end
    
    % Remove duplicate reference points
    if ~isempty(ref_all)
        ref_all = uniquetol(ref_all,1e-10,'ByRows',true);
    end
    
    % Compare system 4 against both system 1 and system 2
    if ~isempty(unique_solutions4)
    
        if isempty(ref_all)
            nonmatching_fp4{k} = unique_solutions4;
        else
            D = pdist2(unique_solutions4, ref_all);
            nonmatch_mask = min(D,[],2) > 1e-2;
            nonmatching_fp4{k} = unique_solutions4(nonmatch_mask,:);
        end
    
    else
        nonmatching_fp4{k} = [];
    end

    % Build reference set from system 1, system 2, and system 4
    ref_all2 = [];
    
    if ~isempty(unique_solutions1)
        ref_all2 = [ref_all2; unique_solutions1];
    end
    
    if ~isempty(unique_solutions2)
        ref_all2 = [ref_all2; unique_solutions2];
    end

    if ~isempty(unique_solutions4)
        ref_all2 = [ref_all2; unique_solutions4];
    end
    
    % Remove duplicate reference points
    if ~isempty(ref_all2)
        ref_all2 = uniquetol(ref_all2,1e-10,'ByRows',true);
    end
    
    % Compare system 8 against system 1, system 2, and system 4
    if ~isempty(unique_solutions8)
    
        if isempty(ref_all2)
            nonmatching_fp8{k} = unique_solutions8;
        else
            D = pdist2(unique_solutions8, ref_all2);
            nonmatch_mask = min(D,[],2) > 1e-2;
            nonmatching_fp8{k} = unique_solutions8(nonmatch_mask,:);
        end
    
    else
        nonmatching_fp8{k} = [];
    end

    stable1 = []; unstable1 = [];
    stable2 = []; unstable2 = [];
    stable4 = []; unstable4 = [];
    stable8 = []; unstable8 = [];

    %% Stability

    % System 1 stability
    for s1 = 1:size(unique_solutions1,1)
        X = unique_solutions1(s1,:)';

        J_fun = makeJacobianParasitoid(b,a,lambda,mu);

        J = J_fun(X);
        eigvals = eig(J);
        if max(abs(eigvals))<1
            stable1 = [stable1; X'];
        else
            unstable1 = [unstable1; X'];
        end
    end

    % System 2 stability (ONLY nonmatching points)
    these_fp2 = nonmatching_fp2{k};       % <- filtered list

    J_fun = makeJacobianParasitoid(b,a,lambda,mu);
    for s = 1:size(these_fp2,1)
    
        X1 = these_fp2(s,:)';
        X2 = parasitoidMap(X1,b,a,lambda,mu);
    
        J1 = J_fun(X1);
        J2 = J_fun(X2);
        M = J2*J1;
        eigvals = eig(M);
    
        if max(abs(eigvals)) < 1
            stable2 = [stable2; X1'];
        else
            unstable2 = [unstable2; X1'];
        end
    
    end

    % System 4 stability (ONLY nonmatching points)
    these_fp4 = nonmatching_fp4{k};       % <- filtered list

    J_fun = makeJacobianParasitoid(b,a,lambda,mu);
    for s = 1:size(these_fp4,1)
    
        X1 = these_fp4(s,:)';
        X2 = parasitoidMap(X1,b,a,lambda,mu);
        X3 = parasitoidMap(X2,b,a,lambda,mu);
        X4 = parasitoidMap(X3,b,a,lambda,mu);
    
        J1 = J_fun(X1);
        J2 = J_fun(X2);
        J3 = J_fun(X3);
        J4 = J_fun(X4);

        M = J4*J3*J2*J1;
        eigvals = eig(M);
    
        if max(abs(eigvals)) < 1
            stable4 = [stable4; X1'];
        else
            unstable4 = [unstable4; X1'];
        end
    
    end

    % System 8 stability (ONLY nonmatching points)
    these_fp8 = nonmatching_fp8{k};       % <- filtered list

    J_fun = makeJacobianParasitoid(b,a,lambda,mu);
    for s = 1:size(these_fp8,1)
    
        X1 = these_fp8(s,:)';
        X2 = parasitoidMap(X1,b,a,lambda,mu);
        X3 = parasitoidMap(X2,b,a,lambda,mu);
        X4 = parasitoidMap(X3,b,a,lambda,mu);
        X5 = parasitoidMap(X4,b,a,lambda,mu);
        X6 = parasitoidMap(X5,b,a,lambda,mu);
        X7 = parasitoidMap(X6,b,a,lambda,mu);
        X8 = parasitoidMap(X7,b,a,lambda,mu);
    
        J1 = J_fun(X1);
        J2 = J_fun(X2);
        J3 = J_fun(X3);
        J4 = J_fun(X4);
        J5 = J_fun(X5);
        J6 = J_fun(X6);
        J7 = J_fun(X7);
        J8 = J_fun(X8);

        M = J8*J7*J6*J5*J4*J3*J2*J1;
        eigvals = eig(M);
    
        if max(abs(eigvals)) < 1
            stable8 = [stable8; X1'];
        else
            unstable8 = [unstable8; X1'];
        end
    
    end

    stable_fp1{k}   = stable1;
    unstable_fp1{k} = unstable1;
    stable_fp2{k}   = stable2;
    unstable_fp2{k} = unstable2;
    stable_fp4{k}   = stable4;
    unstable_fp4{k} = unstable4;
    stable_fp8{k}   = stable8;
    unstable_fp8{k} = unstable8;
end

%% --- Plot results (2D) ---
varNames = {'Host', 'Parasitoid'}; % variable labels
nVars = numel(varNames);

for var = 1:nVars
    figure; hold on;

    hStable = plot(nan, nan, 'b.', 'MarkerFaceColor','b','DisplayName', 'Stable');
    hUnstable = plot(nan, nan, 'r.', 'MarkerFaceColor','r','DisplayName', 'Unstable');

    xlabel('lambda', 'FontSize', 14);
    ylabel(varNames{var}, 'FontSize', 14);
    grid on;

    for k = 1:n_l
        lambda = lambda_values(k);

        % Plot unstable system 8 fps
            if ~isempty(unstable_fp8{k})
                plot(lambda*ones(size(unstable_fp8{k},1),1), unstable_fp8{k}(:,var), 'r.', 'MarkerFaceColor', 'r');
            end

        % Plot unstable system 4 fps
            if ~isempty(unstable_fp4{k})
                plot(lambda*ones(size(unstable_fp4{k},1),1), unstable_fp4{k}(:,var), 'r.', 'MarkerFaceColor', 'r');
            end
        
        % Plot unstable system 2 fps
            if ~isempty(unstable_fp2{k})
                plot(lambda*ones(size(unstable_fp2{k},1),1), unstable_fp2{k}(:,var), 'r.', 'MarkerFaceColor', 'r');
            end

            % Plot unstable system 1 fps
            if ~isempty(unstable_fp1{k})
                plot(lambda*ones(size(unstable_fp1{k},1),1), unstable_fp1{k}(:,var), 'r.', 'MarkerFaceColor', 'r');
            end
                       
            % Plot stable system 8 fps
            if ~isempty(stable_fp8{k})
                plot(lambda*ones(size(stable_fp8{k},1),1),stable_fp8{k}(:,var), 'b.', 'MarkerFaceColor', 'b');
            end

            % Plot stable system 4 fps
            if ~isempty(stable_fp4{k})
                plot(lambda*ones(size(stable_fp4{k},1),1),stable_fp4{k}(:,var), 'b.', 'MarkerFaceColor', 'b');
            end

            % Plot stable system 2 fps
            if ~isempty(stable_fp2{k})
                plot(lambda*ones(size(stable_fp2{k},1),1),stable_fp2{k}(:,var), 'b.', 'MarkerFaceColor', 'b');
            end

            % Plot stable system 1 fps
            if ~isempty(stable_fp1{k})
                plot(lambda*ones(size(stable_fp1{k},1),1),stable_fp1{k}(:,var), 'b.', 'MarkerFaceColor', 'b');
            end
    end
end

%% --- Plot results (3D) ---
figure; hold on; grid on;
xlabel('$\lambda$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Hosts','FontSize',14);
zlabel('Parasitoids','FontSize',14);
view(3);
hUnstable = plot3(nan, nan, nan,'r', 'LineWidth', 1.5);
hStable = plot3(nan, nan, nan,'b', 'LineWidth', 1.5);
hold on

    for k = 1:n_l
        lambda = lambda_values(k);

        % Plot unstable system 8 fps
            if ~isempty(unstable_fp8{k})
                plot3(lambda*ones(size(unstable_fp8{k},1),1), unstable_fp8{k}(:,1),unstable_fp8{k}(:,2), 'r.', 'MarkerFaceColor', 'r');
            end

        % Plot unstable system 4 fps
            if ~isempty(unstable_fp4{k})
                plot3(lambda*ones(size(unstable_fp4{k},1),1), unstable_fp4{k}(:,1), unstable_fp4{k}(:,2),'r.', 'MarkerFaceColor', 'r');
            end
        
        % Plot unstable system 2 fps
            if ~isempty(unstable_fp2{k})
                plot3(lambda*ones(size(unstable_fp2{k},1),1), unstable_fp2{k}(:,1), unstable_fp2{k}(:,2), 'r.', 'MarkerFaceColor', 'r');
            end

            % Plot unstable system 1 fps
            if ~isempty(unstable_fp1{k})
                plot3(lambda*ones(size(unstable_fp1{k},1),1), unstable_fp1{k}(:,1), unstable_fp1{k}(:,2), 'r.', 'MarkerFaceColor', 'r');
            end
                       
            % Plot stable system 8 fps
            if ~isempty(stable_fp8{k})
                plot3(lambda*ones(size(stable_fp8{k},1),1),stable_fp8{k}(:,1),stable_fp8{k}(:,2), 'b.', 'MarkerFaceColor', 'b');
            end

            % Plot stable system 4 fps
            if ~isempty(stable_fp4{k})
                plot3(lambda*ones(size(stable_fp4{k},1),1),stable_fp4{k}(:,1),stable_fp4{k}(:,2), 'b.', 'MarkerFaceColor', 'b');
            end

            % Plot stable system 2 fps
            if ~isempty(stable_fp2{k})
                plot3(lambda*ones(size(stable_fp2{k},1),1),stable_fp2{k}(:,1),stable_fp2{k}(:,2), 'b.', 'MarkerFaceColor', 'b');
            end

            % Plot stable system 1 fps
            if ~isempty(stable_fp1{k})
                plot3(lambda*ones(size(stable_fp1{k},1),1),stable_fp1{k}(:,1),stable_fp1{k}(:,2), 'b.', 'MarkerFaceColor', 'b');
            end
    end


%% Clean plotting

% close all;
% Plotter for cycles
n = length(lambda_values);

% Fixed point
fp_branch_stable = nan(n,2);
fp_branch_unstable_zeros = nan(n,2);
fp_branch_unstable_host = nan(n,2);
fp_branch_unstable = nan(n,2);

% 2-cycle
lower_branch_stable = nan(n,2);
upper_branch_stable = nan(n,2);

lower_branch_unstable_sync = nan(n,2);
upper_branch_unstable_sync = nan(n,2);
lower_branch_unstable_async = nan(n,2);
upper_branch_unstable_async = nan(n,2);

% 4-cycle
lower_branch_stable_4cyc_1 = nan(n,2);
lower_branch_stable_4cyc_2 = nan(n,2);
upper_branch_stable_4cyc_1 = nan(n,2);
upper_branch_stable_4cyc_2 = nan(n,2);

lower_branch_unstable_4cyc_1 = nan(n,2);
lower_branch_unstable_4cyc_2 = nan(n,2);
upper_branch_unstable_4cyc_1 = nan(n,2);
upper_branch_unstable_4cyc_2 = nan(n,2);

% 8-cycle
lower_branch_stable_8cyc_1 = nan(n,2);
lower_branch_stable_8cyc_2 = nan(n,2);
lower_branch_stable_8cyc_3 = nan(n,2);
lower_branch_stable_8cyc_4 = nan(n,2);
upper_branch_stable_8cyc_1 = nan(n,2);
upper_branch_stable_8cyc_2 = nan(n,2);
upper_branch_stable_8cyc_3 = nan(n,2);
upper_branch_stable_8cyc_4 = nan(n,2);

lower_branch_unstable_8cyc_1 = nan(n,2);
lower_branch_unstable_8cyc_2 = nan(n,2);
lower_branch_unstable_8cyc_3 = nan(n,2);
lower_branch_unstable_8cyc_4 = nan(n,2);
upper_branch_unstable_8cyc_1 = nan(n,2);
upper_branch_unstable_8cyc_2 = nan(n,2);
upper_branch_unstable_8cyc_3 = nan(n,2);
upper_branch_unstable_8cyc_4 = nan(n,2);

for i = 1:n
    if ~isempty(stable_fp1{i})
        fp_branch_stable(i,:) = stable_fp1{i}(1,:);
    end
    if ~isempty(unstable_fp1{i})
    if size(unstable_fp1{i},1)==1
        fp_branch_unstable_zeros(i,:) = unstable_fp1{i}(1,:);
    elseif size(unstable_fp1{i},1)==2
        if unstable_fp1{i}(2,2)>0
            fp_branch_unstable_zeros(i,:) = unstable_fp1{i}(1,:);
            fp_branch_unstable(i,:) = unstable_fp1{i}(2,:);
        else
            fp_branch_unstable_zeros(i,:) = unstable_fp1{i}(1,:);
            fp_branch_unstable_host(i,:) = unstable_fp1{i}(2,:);
        end
    elseif size(unstable_fp1{i},1)==3
        if unstable_fp1{i}(2,2) ==0
            fp_branch_unstable_zeros(i,:) = unstable_fp1{i}(1,:);
            fp_branch_unstable_host(i,:) = unstable_fp1{i}(2,:);
        elseif unstable_fp1{i}(3,2)>0
            fp_branch_unstable_zeros(i,:) = unstable_fp1{i}(1,:);
            fp_branch_unstable(i,:) = unstable_fp1{i}(2,:);
        else
            fp_branch_unstable_zeros(i,:) = unstable_fp1{i}(1,:);
            fp_branch_unstable(i,:) = unstable_fp1{i}(2,:);
            fp_branch_unstable_host(i,:) = unstable_fp1{i}(3,:);
        end
    else
        warning('Unstable_fp1: Unexpected number of rows in cell %d',i)
    end
    end

    if ~isempty(stable_fp2{i})
    if size(stable_fp2{i},1) == 2
        lower_branch_stable(i,:) = stable_fp2{i}(1,:);
        upper_branch_stable(i,:) = stable_fp2{i}(2,:);
    else
        warning('Stable_fp2: Unexpected number of rows in cell %d',i);
    end
    end
end

for i = 1:n
    if ~isempty(unstable_fp2{i})
    if size(unstable_fp2{i},1) == 2
        if unstable_fp2{i}(1,2)==0
            lower_branch_unstable_sync(i,:) = unstable_fp2{i}(1,:);
            upper_branch_unstable_sync(i,:) = unstable_fp2{i}(2,:);
        else
            lower_branch_unstable_async(i,:) = unstable_fp2{i}(1,:);
            upper_branch_unstable_async(i,:) = unstable_fp2{i}(2,:);
        end

    elseif size(unstable_fp2{i},1) == 4
        if all(unstable_fp2{i}(:,2)==0)
            lower_branch_unstable_sync(i,:) = unstable_fp2{i}(1,:);
            upper_branch_unstable_sync(i,:) = unstable_fp2{i}(4,:);
        else
            lower_branch_unstable_sync(i,:) = unstable_fp2{i}(1,:);
            lower_branch_unstable_async(i,:) = unstable_fp2{i}(2,:);
            upper_branch_unstable_async(i,:) = unstable_fp2{i}(3,:);
            upper_branch_unstable_sync(i,:) = unstable_fp2{i}(4,:);
        end
    elseif size(unstable_fp2{i},1) == 1
        x1 = unstable_fp2{i}(1,:);
        x2 = map(x1,b,a,lambda_values(i),mu);
        lower_branch_unstable_sync(i,:) = x2;
        upper_branch_unstable_sync(i,:) = x1;
    elseif size(unstable_fp2{i},1) == 3
    x1 = unstable_fp2{i}(1,:);
    x2 = map(x1,b,a,lambda_values(i),mu);
    lower_branch_unstable_async(i,:) = x1;
    upper_branch_unstable_async(i,:) = x2;
    else
        warning('Unstable_fp2: Unexpected number of rows in cell %d',i);
    end
    end
end

% 8-cycle
for i = 1:n
    if ~isempty(stable_fp8{i})
    cycle = zeros(8,size(stable_fp8{i},2));
    cycle(1,:) = stable_fp8{i}(1,:);
    for j = 2:8
        cycle(j,:) = map(cycle(j-1,:),b,a,lambda_values(i),mu);
    end
    [~,idx] = sort(cycle(:,1));
    cycle = cycle(idx,:);

    lower_branch_stable_8cyc_1(i,:) = cycle(1,:);
    lower_branch_stable_8cyc_2(i,:) = cycle(2,:);
    lower_branch_stable_8cyc_3(i,:) = cycle(3,:);
    lower_branch_stable_8cyc_4(i,:) = cycle(4,:);
    upper_branch_stable_8cyc_1(i,:) = cycle(5,:);
    upper_branch_stable_8cyc_2(i,:) = cycle(6,:);
    upper_branch_stable_8cyc_3(i,:) = cycle(7,:);
    upper_branch_stable_8cyc_4(i,:) = cycle(8,:);
    end
    if ~isempty(unstable_fp8{i})
    cycle = zeros(8,size(unstable_fp8{i},2));
    cycle(1,:) = unstable_fp8{i}(1,:);
    for j = 2:8
        cycle(j,:) = map(cycle(j-1,:),b,a,lambda_values(i),mu);
    end
    [~,idx] = sort(cycle(:,1));
    cycle = cycle(idx,:);

    lower_branch_unstable_8cyc_1(i,:) = cycle(1,:);
    lower_branch_unstable_8cyc_2(i,:) = cycle(2,:);
    lower_branch_unstable_8cyc_3(i,:) = cycle(3,:);
    lower_branch_unstable_8cyc_4(i,:) = cycle(4,:);
    upper_branch_unstable_8cyc_1(i,:) = cycle(5,:);
    upper_branch_unstable_8cyc_2(i,:) = cycle(6,:);
    upper_branch_unstable_8cyc_3(i,:) = cycle(7,:);
    upper_branch_unstable_8cyc_4(i,:) = cycle(8,:);
    end
end

% Correct as needed
% for i = 1:423
% lower_branch_unstable_sync(i,:) = nan(1,2);
% upper_branch_unstable_sync(i,:) = nan(1,2);
% end
% 
% for i = 735:n
% lower_branch_unstable_sync(i,:) = nan(1,2);
% upper_branch_unstable_sync(i,:) = nan(1,2);
% end
% 
% for i = 1:849
% lower_branch_unstable_4cyc_1(i,:) = nan(1,2);
% lower_branch_unstable_4cyc_2(i,:) = nan(1,2);
% upper_branch_unstable_4cyc_1(i,:) = nan(1,2);
% upper_branch_unstable_4cyc_2(i,:) = nan(1,2);
% end
% 
% for i = 1:734
% lower_branch_stable_4cyc_1(i,:) = nan(1,2);
% lower_branch_stable_4cyc_2(i,:) = nan(1,2);
% upper_branch_stable_4cyc_1(i,:) = nan(1,2);
% upper_branch_stable_4cyc_2(i,:) = nan(1,2);
% end
% 
% for i = 1:875
% lower_branch_unstable_8cyc_1(i,:) = nan(1,2);
% lower_branch_unstable_8cyc_2(i,:) = nan(1,2);
% lower_branch_unstable_8cyc_3(i,:) = nan(1,2);
% lower_branch_unstable_8cyc_4(i,:) = nan(1,2);
% upper_branch_unstable_8cyc_1(i,:) = nan(1,2);
% upper_branch_unstable_8cyc_2(i,:) = nan(1,2);
% upper_branch_unstable_8cyc_3(i,:) = nan(1,2);
% upper_branch_unstable_8cyc_4(i,:) = nan(1,2);
% end
% 
% for i = 1:849
% lower_branch_stable_8cyc_1(i,:) = nan(1,2);
% lower_branch_stable_8cyc_2(i,:) = nan(1,2);
% lower_branch_stable_8cyc_3(i,:) = nan(1,2);
% lower_branch_stable_8cyc_4(i,:) = nan(1,2);
% upper_branch_stable_8cyc_1(i,:) = nan(1,2);
% upper_branch_stable_8cyc_2(i,:) = nan(1,2);
% upper_branch_stable_8cyc_3(i,:) = nan(1,2);
% upper_branch_stable_8cyc_4(i,:) = nan(1,2);
% end
% 
% for i = 878:n
% lower_branch_stable_8cyc_1(i,:) = nan(1,2);
% lower_branch_stable_8cyc_2(i,:) = nan(1,2);
% lower_branch_stable_8cyc_3(i,:) = nan(1,2);
% lower_branch_stable_8cyc_4(i,:) = nan(1,2);
% upper_branch_stable_8cyc_1(i,:) = nan(1,2);
% upper_branch_stable_8cyc_2(i,:) = nan(1,2);
% upper_branch_stable_8cyc_3(i,:) = nan(1,2);
% upper_branch_stable_8cyc_4(i,:) = nan(1,2);
% end

%% Plotting
varNames = {'Hosts','Parasitoids'};
for col = 1:2

    figure; hold on;
    xlabel('$\lambda$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel(varNames{col},'FontSize',14);
    hUnstable = plot(nan, nan, 'r', 'LineWidth', 1.5);
    hStable = plot(nan, nan, 'b', 'LineWidth', 1.5);
    hold on
    % Zeros
    plot(lambda_values, zeros(size(lambda_values)),'r','LineWidth',1.5);

    % 8-cycle unstable
    y1 = lower_branch_unstable_8cyc_1(:,col);
    mask = ~isnan(y1);
    plot(lambda_values(mask), y1(mask),'r:','LineWidth',1.5)
    
    y1 = lower_branch_unstable_8cyc_2(:,col);
    mask = ~isnan(y1);
    plot(lambda_values(mask), y1(mask),'r:','LineWidth',1.5)

    y1 = lower_branch_unstable_8cyc_3(:,col);
    mask = ~isnan(y1);
    plot(lambda_values(mask), y1(mask),'r:','LineWidth',1.5)

    y1 = lower_branch_unstable_8cyc_4(:,col);
    mask = ~isnan(y1);
    plot(lambda_values(mask), y1(mask),'r:','LineWidth',1.5)
    
    y1 = upper_branch_unstable_8cyc_1(:,col);
    mask = ~isnan(y1);
    plot(lambda_values(mask), y1(mask),'r:','LineWidth',1.5)

    y1 = upper_branch_unstable_8cyc_2(:,col);
    mask = ~isnan(y1);
    plot(lambda_values(mask), y1(mask),'r:','LineWidth',1.5)

    y1 = upper_branch_unstable_8cyc_3(:,col);
    mask = ~isnan(y1);
    plot(lambda_values(mask), y1(mask),'r:','LineWidth',1.5)

    y1 = upper_branch_unstable_8cyc_4(:,col);
    mask = ~isnan(y1);
    plot(lambda_values(mask), y1(mask),'r:','LineWidth',1.5)

    % 4-cycle unstable
    y1 = lower_branch_unstable_4cyc_1(:,col);
    mask = ~isnan(y1);
    plot(lambda_values, y1,'r--','LineWidth',1.5)

    
    y1 = lower_branch_unstable_4cyc_2(:,col);
    mask = ~isnan(y1);
    plot(lambda_values, y1,'r--','LineWidth',1.5)
    
    y1 = upper_branch_unstable_4cyc_1(:,col);
    mask = ~isnan(y1);
    plot(lambda_values, y1,'r--','LineWidth',1.5)
    
    y1 = upper_branch_unstable_4cyc_2(:,col);
    mask = ~isnan(y1);
    plot(lambda_values, y1,'r--','LineWidth',1.5)

    % 2-cycle unstable
    y1 = lower_branch_unstable_sync(:,col);
    mask = ~isnan(y1);
    plot(lambda_values, y1,'r--','LineWidth',1.5)
    
    y1 = upper_branch_unstable_sync(:,col);
    mask = ~isnan(y1);
    plot(lambda_values, y1,'r--','LineWidth',1.5)
    
    y1 = lower_branch_unstable_async(:,col);
    mask = ~isnan(y1);
    plot(lambda_values(mask), y1(mask),'r-.','LineWidth',1.5)
    
    y1 = upper_branch_unstable_async(:,col);
    mask = ~isnan(y1);
    plot(lambda_values(mask), y1(mask),'r-.','LineWidth',1.5)

    % 8-cycle stable
    y1 = lower_branch_stable_8cyc_1(:,col);
    mask = ~isnan(y1);
    plot(lambda_values(mask), y1(mask),'b:','LineWidth',1.5)
    
    y1 = lower_branch_stable_8cyc_2(:,col);
    mask = ~isnan(y1);
    plot(lambda_values(mask), y1(mask),'b:','LineWidth',1.5)

    y1 = lower_branch_stable_8cyc_3(:,col);
    mask = ~isnan(y1);
    plot(lambda_values(mask), y1(mask),'b:','LineWidth',1.5)

    y1 = lower_branch_stable_8cyc_4(:,col);
    mask = ~isnan(y1);
    plot(lambda_values(mask), y1(mask),'b:','LineWidth',1.5)
    
    y1 = upper_branch_stable_8cyc_1(:,col);
    mask = ~isnan(y1);
    plot(lambda_values(mask), y1(mask),'b:','LineWidth',1.5)
    
    y1 = upper_branch_stable_8cyc_2(:,col);
    mask = ~isnan(y1);
    plot(lambda_values(mask), y1(mask),'b:','LineWidth',1.5)

    y1 = upper_branch_stable_8cyc_3(:,col);
    mask = ~isnan(y1);
    plot(lambda_values(mask), y1(mask),'b:','LineWidth',1.5)

    y1 = upper_branch_stable_8cyc_4(:,col);
    mask = ~isnan(y1);
    plot(lambda_values(mask), y1(mask),'b:','LineWidth',1.5)

    % 4-cycle stable
    y1 = lower_branch_stable_4cyc_1(:,col);
    mask = ~isnan(y1);
    plot(lambda_values(mask), y1(mask),'b--','LineWidth',1.5)
    
    y1 = lower_branch_stable_4cyc_2(:,col);
    mask = ~isnan(y1);
    plot(lambda_values(mask), y1(mask),'b--','LineWidth',1.5)
    
    y1 = upper_branch_stable_4cyc_1(:,col);
    mask = ~isnan(y1);
    plot(lambda_values(mask), y1(mask),'b--','LineWidth',1.5)
    
    y1 = upper_branch_stable_4cyc_2(:,col);
    mask = ~isnan(y1);
    plot(lambda_values(mask), y1(mask),'b--','LineWidth',1.5)

    % 2-cycle stable
    y1 = lower_branch_stable(:,col);
    mask = ~isnan(y1);
    plot(lambda_values(mask), y1(mask),'b-.','LineWidth',1.5)
    
    y1 = upper_branch_stable(:,col);
    mask = ~isnan(y1);
    plot(lambda_values(mask), y1(mask),'b-.','LineWidth',1.5)

    % Fixed points
    y1 = fp_branch_unstable_host(:,col);
    mask = ~isnan(y1);
    plot(lambda_values(mask), y1(mask),'r','LineWidth',1.5);
    y1 = fp_branch_unstable(:,col);
    mask = ~isnan(y1);
    plot(lambda_values(mask), y1(mask),'r','LineWidth',1.5);
    plot(lambda_values, fp_branch_stable(:,col),'b','LineWidth',1.5);
    
    legend([hStable, hUnstable], {'Stable', 'Unstable'}, 'FontSize',14,'Location','best');
    
xticks([      0,     1,     5,   7.3891,  10, 10.889, 10.9218, 11.5882, 12.5126, 14.2462, 14.6683])
xticklabels({'0', '$B_1$', '5', '$B_2$', '10','', '',   '$B_5$', '$B_6$',   '$B_7$',       ''})
    ax = gca; 
    ax.TickLabelInterpreter = 'latex';

xlim([0 15])

ax = gca;
ax.TickLabelInterpreter = 'latex';

h1 = annotation('textarrow', [0.665 0.69], [0.05 0.1], ...
     'String', '$B_3$', 'Interpreter','LaTex');
 
h2 = annotation('textarrow', [0.72 0.7], [0.05 0.1], ...
    'String', '$B_4$', 'Interpreter','LaTex');
 
h3 = annotation('textarrow', [0.888 0.888], [0.05 0.1], ...
    'String', '$B_8$', 'Interpreter','LaTex');


h1.LineWidth = 1;
h1.HeadWidth = 4;
h1.HeadLength = 4;

h2.LineWidth = 1;
h2.HeadWidth = 4;
h2.HeadLength = 4;

h3.LineWidth = 1;
h3.HeadWidth = 4;
h3.HeadLength = 4;

end

%% Plot 3D
figure; hold on; grid on;
xlabel('$\lambda$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Hosts','FontSize',14);
zlabel('Parasitoids','FontSize',14);
view(3);
hUnstable = plot3(nan, nan, nan,'r', 'LineWidth', 1.5);
hStable = plot3(nan, nan, nan,'b', 'LineWidth', 1.5);
hold on
% Zeros
plot3(lambda_values, zeros(size(lambda_values)),zeros(size(lambda_values)),'r--','LineWidth',1.5);
% 8-cycle unstable
y1 = lower_branch_unstable_8cyc_1(:,1);
y2 = lower_branch_unstable_8cyc_1(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'r:','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = lower_branch_unstable_8cyc_2(:,1);
y2 = lower_branch_unstable_8cyc_2(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'r:','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = lower_branch_unstable_8cyc_3(:,1);
y2 = lower_branch_unstable_8cyc_3(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'r:','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = lower_branch_unstable_8cyc_4(:,1);
y2 = lower_branch_unstable_8cyc_4(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'r:','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = upper_branch_unstable_8cyc_1(:,1);
y2 = upper_branch_unstable_8cyc_1(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'r:','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = upper_branch_unstable_8cyc_2(:,1);
y2 = upper_branch_unstable_8cyc_2(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'r:','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = upper_branch_unstable_8cyc_3(:,1);
y2 = upper_branch_unstable_8cyc_3(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'r:','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = upper_branch_unstable_8cyc_4(:,1);
y2 = upper_branch_unstable_8cyc_4(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'r:','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)
% 4-cycle unstable
y1 = lower_branch_unstable_4cyc_1(:,1);
y2 = lower_branch_unstable_4cyc_1(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'r:','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = lower_branch_unstable_4cyc_2(:,1);
y2 = lower_branch_unstable_4cyc_2(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'r:','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = upper_branch_unstable_4cyc_1(:,1);
y2 = upper_branch_unstable_4cyc_1(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'r:','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = upper_branch_unstable_4cyc_2(:,1);
y2 = upper_branch_unstable_4cyc_2(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'r:','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)
% 2-cycle unstable
y1 = lower_branch_unstable_sync(:,1);
y2 = lower_branch_unstable_sync(:,2);
mask = ~isnan(y1);
plot3(lambda_values, y1,y2,'r:','LineWidth',1.5)
plot3(lambda_values, y1,y2,'w','LineWidth',1)

y1 = upper_branch_unstable_sync(:,1);
y2 = upper_branch_unstable_sync(:,2);
mask = ~isnan(y1);
plot3(lambda_values, y1,y2,'r:','LineWidth',1.5)
plot3(lambda_values, y1,y2,'w','LineWidth',1)

y1 = lower_branch_unstable_async(:,1);
y2 = lower_branch_unstable_async(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'r--','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = upper_branch_unstable_async(:,1);
y2 = upper_branch_unstable_async(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'r--','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)
% 8-cycle stable
y1 = lower_branch_stable_8cyc_1(:,1);
y2 = lower_branch_stable_8cyc_1(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'b:','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = lower_branch_stable_8cyc_2(:,1);
y2 = lower_branch_stable_8cyc_2(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'b:','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = lower_branch_stable_8cyc_3(:,1);
y2 = lower_branch_stable_8cyc_3(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'b:','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = lower_branch_stable_8cyc_4(:,1);
y2 = lower_branch_stable_8cyc_4(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'b:','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = upper_branch_stable_8cyc_1(:,1);
y2 = upper_branch_stable_8cyc_1(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'b:','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = upper_branch_stable_8cyc_2(:,1);
y2 = upper_branch_stable_8cyc_2(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'b:','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = upper_branch_stable_8cyc_3(:,1);
y2 = upper_branch_stable_8cyc_3(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'b:','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = upper_branch_stable_8cyc_4(:,1);
y2 = upper_branch_stable_8cyc_4(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'b:','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)
% 4-cycle stable
y1 = lower_branch_stable_4cyc_1(:,1);
y2 = lower_branch_stable_4cyc_1(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'b:','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = lower_branch_stable_4cyc_2(:,1);
y2 = lower_branch_stable_4cyc_2(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'b:','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = upper_branch_stable_4cyc_1(:,1);
y2 = upper_branch_stable_4cyc_1(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'b:','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = upper_branch_stable_4cyc_2(:,1);
y2 = upper_branch_stable_4cyc_2(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'b:','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)
% 2-cycle stable
y1 = lower_branch_stable(:,1);
y2 = lower_branch_stable(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'b:','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

y1 = upper_branch_stable(:,1);
y2 = upper_branch_stable(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'b:','LineWidth',1.5)
plot3(lambda_values(mask), y1(mask),y2(mask),'w','LineWidth',1)

% Fixed points
y1 = fp_branch_unstable_host(:,1);
y2 = fp_branch_unstable_host(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'r','LineWidth',1.5);

y1 = fp_branch_unstable(:,1);
y2 = fp_branch_unstable(:,2);
mask = ~isnan(y1);
plot3(lambda_values(mask), y1(mask),y2(mask),'r','LineWidth',1.5);

plot3(lambda_values, fp_branch_stable(:,1),fp_branch_stable(:,2),'b','LineWidth',1.5);

legend([hStable, hUnstable], {'Stable', 'Unstable'}, 'FontSize',14,'Location','best');

xticks([      0,     1,     5,   7.3891,  10, 10.889, 10.9218, 11.5882, 12.5126, 14.2462, 14.6683])
xticklabels({'0', '$B_1$', '5', '$B_2$', '10','', '',   '$B_5$', '$B_6$',   '$B_7$',       ''})
ax = gca; 
ax.TickLabelInterpreter = 'latex';


xlim([0 15])

ax = gca;
ax.TickLabelInterpreter = 'latex';

h1 = annotation('textarrow', [0.81 0.782], [0.13 0.18], ...
    'String', '$B_4$', 'Interpreter','LaTex');

h2 = annotation('textarrow', [0.95 0.91], [0.19 0.2], ...
    'String', '$B_8$', 'Interpreter','LaTex');

h3 = annotation('textarrow', [0.79 0.775], [0.11 0.18], ...
    'String', '$B_3$', 'Interpreter','LaTex');


h1.LineWidth = 1;
h1.HeadWidth = 4;
h1.HeadLength = 4;

h2.LineWidth = 1;
h2.HeadWidth = 4;
h2.HeadLength = 4;

h3.LineWidth = 1;
h3.HeadWidth = 4;
h3.HeadLength = 4;