function print_coefs(opti_coefficients,new_time, policy_dates,SmoothVar,tests)

    Npop = 59240329; % Total Population of Italy
    customColors2 = {   [1, 0.6, 0.2],
                        [1, 0.3, 0.05],
                        [1, 0.2, 0.05],
                        [0.8, 0.2, 0.1],
                        [1, 0.2, 0.05],            
                        [0.8, 0.2, 0.1],
                        [1, 0.3, 0.05],
                        [1, 0.6, 0.2]
                    };
    
    for ii = 1:length(policy_dates)-1
    
        area.x(ii, :) = [policy_dates(ii) policy_dates(ii) policy_dates(ii+1) policy_dates(ii+1)];
        area.y_alpha(ii, :) = [0 max(opti_coefficients.alpha)*1.1 max(opti_coefficients.alpha)*1.1 0];
        area.y_beta(ii, :) = [0 max(opti_coefficients.beta)*1.1 max(opti_coefficients.beta)*1.1 0];
        area.y_gamma(ii, :) = [0 max(opti_coefficients.gamma)*1.1 max(opti_coefficients.gamma)*1.1 0];
        
    end

    % Alpha and Variants
    figure()
    subplot(2,1,1)
    yyaxis right;
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.SC2 zeros(1,57)], [0 0.4470 0.7410], 'LineStyle', '-','FaceAlpha',.1,'EdgeColor', [0 0.4470 0.7410])
    hold on
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.alpha zeros(1,57)], [0.8500 0.3250 0.0980], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.8500 0.3250 0.0980])
    hold on 
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.gamma zeros(1,57)], [0.9290 0.6940 0.1250], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.9290 0.6940 0.1250])
    hold on
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.delta zeros(1,57)], [0.4940 0.1840 0.5560], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.4940 0.1840 0.5560])
    ax2 = gca;
    ax2.YColor = 'k';
    yyaxis left;
    plot(new_time, opti_coefficients.alpha, 'k','LineWidth', 2);
    ylabel('$\beta$ Coefficient Values','Interpreter','latex')
    ax1 = gca;
    ax1.YColor = 'k';
    title('\textbf{}','Interpreter','latex')
    title('\textbf{SARS-CoV-2 Variants - Coefficient $\alpha$}','Interpreter','latex')
    legend('$\alpha$','SARS-CoV-2', 'B.1.1.7 - ALPHA VARIANT', 'P.1 - GAMMA VARIANT', ' B.1.617.2 - DELTA VARIANT', 'Interpreter','latex','Location','northwest')
    grid on
    xlim([SmoothVar.date(1), SmoothVar.date(end)])
    ylim([0, max(opti_coefficients.alpha)*1.1])
    subplot(2,1,2)
    for ii = 1:length(policy_dates)-1
        fill(area.x(ii, :) ,area.y_alpha(1, :), customColors2{ii,1} ,'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', 'off')
        hold on
        xline(policy_dates(ii),":",'HandleVisibility', 'off')
        hold on 
    end
    hold on
    plot(new_time, opti_coefficients.alpha,'k','LineWidth',1.5, 'DisplayName', '$\alpha$')
    ylabel('Coefficients Values','Interpreter','latex')
    title('\textbf{$\alpha$ coefficient}','Interpreter','latex')
    grid on
    legend('Interpreter','latex','location','southeast')
    ylim([0, max(opti_coefficients.alpha)*1.1])
    xlim([new_time(1), new_time(end)])
    set(gca, 'TickLabelInterpreter', 'Latex')
    
    
    % Beta and Variants
    figure()
    subplot(2,1,1)
    yyaxis right;
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.SC2 zeros(1,57)], [0 0.4470 0.7410], 'LineStyle', '-','FaceAlpha',.1,'EdgeColor', [0 0.4470 0.7410])
    hold on
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.alpha zeros(1,57)], [0.8500 0.3250 0.0980], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.8500 0.3250 0.0980])
    hold on 
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.gamma zeros(1,57)], [0.9290 0.6940 0.1250], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.9290 0.6940 0.1250])
    hold on
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.delta zeros(1,57)], [0.4940 0.1840 0.5560], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.4940 0.1840 0.5560])
    ax2 = gca;
    ax2.YColor = 'k';
    yyaxis left;
    plot(new_time, opti_coefficients.beta, 'k','LineWidth', 2);
    ylabel('$\beta$ Coefficient Values','Interpreter','latex')
    ax1 = gca;
    ax1.YColor = 'k';
    title('\textbf{}','Interpreter','latex')
    title('\textbf{SARS-CoV-2 Variants - Coefficient $\beta$}','Interpreter','latex')
    legend('$\beta$','SARS-CoV-2', 'B.1.1.7 - ALPHA VARIANT', 'P.1 - GAMMA VARIANT', ' B.1.617.2 - DELTA VARIANT', 'Interpreter','latex','Location','northwest')
    grid on
    xlim([SmoothVar.date(1), SmoothVar.date(end)])
    ylim([0, max(opti_coefficients.beta)*1.1])
    subplot(2,1,2)
    for ii = 1:length(policy_dates)-1
        fill(area.x(ii, :) ,area.y_beta(1, :), customColors2{ii,1} ,'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', 'off')
        hold on
        xline(policy_dates(ii),":",'HandleVisibility', 'off')
        hold on 
    end
    hold on
    plot(new_time, opti_coefficients.beta,'k','LineWidth',1.5, 'DisplayName', '$\beta$')
    ylabel('Coefficients Values','Interpreter','latex')
    title('\textbf{$\beta$ coefficient}','Interpreter','latex')
    grid on
    legend('Interpreter','latex')
    xlim([new_time(1), new_time(end)])
    ylim([0, max(opti_coefficients.beta)*1.1])
    set(gca, 'TickLabelInterpreter', 'Latex')


    % Gamma coefficient
    figure()
    yyaxis left;
    plot(new_time, opti_coefficients.gamma,'LineWidth', 1.5);
    ylabel('$\gamma$ Coefficients Values','Interpreter','latex')
    ax1 = gca;
    ax1.YColor = 'k';
    yyaxis right;
    plot(new_time, tests.data/Npop,'-.','LineWidth',1.5,'DisplayName', 'Testing Activity')
    ylabel('$\char"0023$ of Tests','Interpreter','latex')
    ax2 = gca;
    ax2.YColor = 'k';
    title('\textbf{Comparison Testing Activity and $\gamma$ Coefficient}','Interpreter','latex')
    legend(' $\gamma$', 'Tests','Interpreter','latex')
    grid on;
    xlim([new_time(1), new_time(end)])


    % Deltas and Variants
    figure()
    yyaxis right;
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.SC2 zeros(1,57)], [0 0.4470 0.7410], 'LineStyle', '-','FaceAlpha',.1,'EdgeColor', [0 0.4470 0.7410])
    hold on
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.alpha zeros(1,57)], [0.8500 0.3250 0.0980], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.8500 0.3250 0.0980])
    hold on 
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.gamma zeros(1,57)], [0.9290 0.6940 0.1250], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.9290 0.6940 0.1250])
    hold on
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.delta zeros(1,57)], [0.4940 0.1840 0.5560], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.4940 0.1840 0.5560])
    ax2 = gca;
    ax2.YColor = 'k';
    yyaxis left;
    plot(new_time, opti_coefficients.delta1, 'k-', LineWidth=2)
    hold on
    plot(new_time, opti_coefficients.delta2, 'r-', LineWidth=2)
    ylabel('$\delta$ Coefficient Values','Interpreter','latex')
    ax1 = gca;
    ax1.YColor = 'k';
    title('\textbf{}','Interpreter','latex')
    title('\textbf{SARS-CoV-2 Variants - Coefficient $\delta$}','Interpreter','latex')
    legend('$\delta_1$','$\delta_2$','SARS-CoV-2', 'B.1.1.7 - ALPHA VARIANT', 'P.1 - GAMMA VARIANT', ' B.1.617.2 - DELTA VARIANT', 'Interpreter','latex','Location','northwest')
    grid on
    xlim([SmoothVar.date(1), SmoothVar.date(end)])
    ylim([0, max(opti_coefficients.delta1)*1.1])

    % Sigmas and Variants
    figure()
    yyaxis right;
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.SC2 zeros(1,57)], [0 0.4470 0.7410], 'LineStyle', '-','FaceAlpha',.1,'EdgeColor', [0 0.4470 0.7410])
    hold on
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.alpha zeros(1,57)], [0.8500 0.3250 0.0980], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.8500 0.3250 0.0980])
    hold on 
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.gamma zeros(1,57)], [0.9290 0.6940 0.1250], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.9290 0.6940 0.1250])
    hold on
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.delta zeros(1,57)], [0.4940 0.1840 0.5560], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.4940 0.1840 0.5560])
    ax2 = gca;
    ax2.YColor = 'k';
    yyaxis left;
    plot(new_time, opti_coefficients.sigma1, 'k-', LineWidth=2)
    hold on
    plot(new_time, opti_coefficients.sigma2, 'r-', LineWidth=2)
    ylabel('$\sigma$ Coefficient Values','Interpreter','latex')
    ax1 = gca;
    ax1.YColor = 'k';
    title('\textbf{}','Interpreter','latex')
    title('\textbf{SARS-CoV-2 Variants - Coefficient $\sigma$}','Interpreter','latex')
    legend('$\sigma_1$','$\sigma_2$','SARS-CoV-2', 'B.1.1.7 - ALPHA VARIANT', 'P.1 - GAMMA VARIANT', ' B.1.617.2 - DELTA VARIANT', 'Interpreter','latex','Location','northwest')
    grid on
    xlim([SmoothVar.date(1), SmoothVar.date(end)])
    ylim([0, max(opti_coefficients.sigma1)*1.1])

    % Taus and Variants
    figure()
    yyaxis right;
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.SC2 zeros(1,57)], [0 0.4470 0.7410], 'LineStyle', '-','FaceAlpha',.1,'EdgeColor', [0 0.4470 0.7410])
    hold on
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.alpha zeros(1,57)], [0.8500 0.3250 0.0980], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.8500 0.3250 0.0980])
    hold on 
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.gamma zeros(1,57)], [0.9290 0.6940 0.1250], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.9290 0.6940 0.1250])
    hold on
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.delta zeros(1,57)], [0.4940 0.1840 0.5560], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.4940 0.1840 0.5560])
    ax2 = gca;
    ax2.YColor = 'k';
    yyaxis left;
    plot(new_time, opti_coefficients.tau1, 'k-', LineWidth=2)
    hold on
    plot(new_time, opti_coefficients.tau2, 'r-', LineWidth=2)
    ylabel('$\tau$ Coefficient Values','Interpreter','latex')
    ax1 = gca;
    ax1.YColor = 'k';
    title('\textbf{}','Interpreter','latex')
    title('\textbf{SARS-CoV-2 Variants - Coefficient $\tau$}','Interpreter','latex')
    legend('$\tau_1$','$\tau_2$','SARS-CoV-2', 'B.1.1.7 - ALPHA VARIANT', 'P.1 - GAMMA VARIANT', ' B.1.617.2 - DELTA VARIANT', 'Interpreter','latex','Location','northwest')
    grid on
    xlim([SmoothVar.date(1), SmoothVar.date(end)])
    ylim([0, max(opti_coefficients.tau2)*1.2])


    % Lambda and Variants
    figure()
    yyaxis right;
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.SC2 zeros(1,57)], [0 0.4470 0.7410], 'LineStyle', '-','FaceAlpha',.1,'EdgeColor', [0 0.4470 0.7410])
    hold on
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.alpha zeros(1,57)], [0.8500 0.3250 0.0980], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.8500 0.3250 0.0980])
    hold on 
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.gamma zeros(1,57)], [0.9290 0.6940 0.1250], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.9290 0.6940 0.1250])
    hold on
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.delta zeros(1,57)], [0.4940 0.1840 0.5560], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.4940 0.1840 0.5560])
    ax2 = gca;
    ax2.YColor = 'k';
    yyaxis left;
    plot(new_time, opti_coefficients.lambda, 'k','LineWidth', 2);
    ylabel('$\lambda$ Coefficient Values','Interpreter','latex')
    ax1 = gca;
    ax1.YColor = 'k';
    title('\textbf{}','Interpreter','latex')
    title('\textbf{SARS-CoV-2 Variants - Coefficient $\lambda$}','Interpreter','latex')
    legend('$\lambda$','SARS-CoV-2', 'B.1.1.7 - ALPHA VARIANT', 'P.1 - GAMMA VARIANT', ' B.1.617.2 - DELTA VARIANT', 'Interpreter','latex','Location','northwest')
    grid on
    xlim([SmoothVar.date(1), SmoothVar.date(end)])
    ylim([0, max(opti_coefficients.lambda)*1.2])


    % Epsi and Variants
    figure()
    yyaxis right;
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.SC2 zeros(1,57)], [0 0.4470 0.7410], 'LineStyle', '-','FaceAlpha',.1,'EdgeColor', [0 0.4470 0.7410])
    hold on
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.alpha zeros(1,57)], [0.8500 0.3250 0.0980], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.8500 0.3250 0.0980])
    hold on 
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.gamma zeros(1,57)], [0.9290 0.6940 0.1250], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.9290 0.6940 0.1250])
    hold on
    fill([SmoothVar.date; flip(SmoothVar.date)], [SmoothVar.delta zeros(1,57)], [0.4940 0.1840 0.5560], 'LineStyle', '-', 'FaceAlpha',.1,'EdgeColor', [0.4940 0.1840 0.5560])
    ax2 = gca;
    ax2.YColor = 'k';
    yyaxis left;
    plot(new_time, opti_coefficients.epsi, 'k','LineWidth', 2);
    ylabel('$\epsilon$ Coefficient Values','Interpreter','latex')
    ax1 = gca;
    ax1.YColor = 'k';
    title('\textbf{}','Interpreter','latex')
    title('\textbf{SARS-CoV-2 Variants - Coefficient $\epsilon$}','Interpreter','latex')
    legend('$\epsilon$','SARS-CoV-2', 'B.1.1.7 - ALPHA VARIANT', 'P.1 - GAMMA VARIANT', ' B.1.617.2 - DELTA VARIANT', 'Interpreter','latex','Location','northwest')
    grid on
    xlim([SmoothVar.date(1), SmoothVar.date(end)])


    % Coefficients phi
    figure(14)
    plot(new_time, opti_coefficients.phi, LineWidth=1.5)
    ylabel('Coefficients Values','Interpreter','latex')
    title('\textbf{Coefficient $\phi$}','Interpreter','latex')
    grid on
    legend('$\phi$', 'Interpreter','latex', 'Location','northeast')
    xlim([new_time(1), new_time(end)])
    ylim([0, max(opti_coefficients.phi)*1.2])
    set(gca, 'TickLabelInterpreter', 'Latex')
    
    % Coefficient kappa
    figure(15)
    plot(new_time, opti_coefficients.kappa, LineWidth=1.5)
    ylabel('Coefficients Values','Interpreter','latex')
    title('\textbf{Coefficient $\kappa$}','Interpreter','latex')
    grid on
    legend('$\kappa$', 'Interpreter','latex', 'Location','northeast')
    xlim([new_time(1), new_time(end)])
    ylim([0, max(opti_coefficients.kappa)*1.2])
    set(gca, 'TickLabelInterpreter', 'Latex')

end