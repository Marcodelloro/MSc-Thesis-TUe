clc; close all; clear;

opts = odeset('MaxStep',1e-3,'RelTol',1e-6);
[t, y] = ode45(@(t,y) -3*y, [0,5], 1, opts);

y = y + 1e-1*(rand(size(y))-0.5);
y = y(1:100:end);
t = t(1:100:end);

par1 = 0:1e-2:4;
par2 = 0:1e-2:4;
% obj = zeros(length(par1), length(par2));
% for ii = 1:length(par1)
%     for jj = 1:length(par2)
%         [~,x] = ode45(@(t,x) -(par1(ii)+par2(jj))*x, t, y(1), opts);
%         e = y-x;
%         obj(ii,jj) = sum(e.^2);
%     end
% end

load('toy_example_cooked.mat');

%% Data cooking
obj2 = obj;
for ii = 1:length(par1)
    for jj = 1:length(par2)
        if obj2(ii,jj) > 2
            obj2(ii,jj) = NaN;
        end
    end
end

m = min(obj, [], 'all');
[row, col] = find(obj == m);
par(1) = par1(row(1));
par(2) = par2(col(1));

%% Perturbed objective function
obj3 = obj;
for ii = 1:length(par1)
    for jj = 1:length(par2)
        obj3(ii,jj) = obj3(ii,jj) + (norm([par1(ii), par2(jj)] - par) - 1)^2;
        if obj3(ii,jj) > 2
            obj3(ii,jj) = NaN;
        end
    end
end

m2 = min(obj3(1:200,1:200), [], 'all');
[row, col] = find(obj3 == m2);
par_p1(1) = par1(row(1));
par_p1(2) = par2(col(1));
par_p2(1) = par1(row(2));
par_p2(2) = par2(col(2));

%% Plotting
fig1 = figure(1);
c = [0.8500 0.3250 0.0980];
fig1.Color = 'w';
scatter(t, y, 100, 'filled', 'MarkerFaceColor', c);
grid on
xlim([0 3])
xlabel('Time', 'Interpreter', 'Latex');
ylabel('Data', 'Interpreter', 'Latex');
% lgd.FontSize = 18;
box on
set(gca, 'TickLabelInterpreter', 'Latex')


%%
figure(2)
surf(par1, par2, obj2', 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
hold on;
surf(par1, par2, obj3', 'EdgeColor', 'none');
plot3(par(1), par(2), m, 'r*', 'LineWidth', 4, handleVisibility='off');
plot3(par_p1(1), par_p1(2), m2, 'r*', 'LineWidth', 4, HandleVisibility='off');
plot3(par_p2(1), par_p2(2), m2, 'r*', 'LineWidth', 4, HandleVisibility='off');
fplot3(@(v) par(1) + sin(v), @(v) par(2) + cos(v), @(v) m, 'r:', 'LineWidth', 2);
zlim([0, 3]);
xlabel('$\alpha$', 'Interpreter', 'Latex');
ylabel('$\beta$', 'Interpreter', 'Latex');
set(gca, 'TickLabelInterpreter', 'Latex');
legend('$V(\theta)$','$V^R_{pen}$','$R$ penalty radius','Interpreter', 'Latex');
% Set a different view for the second subplot
az2 = 20;
el2 = 40; 
view(az2, el2); % Apply the new view settings
box on 

figure(3)
% Second subplot with a different view
surf(par1, par2, obj2', 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
hold on;
surf(par1, par2, obj3', 'EdgeColor', 'none');
plot3(par(1), par(2), m, 'r*', 'LineWidth', 4);
plot3(par_p1(1), par_p1(2), m2, 'r*', 'LineWidth', 4);
plot3(par_p2(1), par_p2(2), m2, 'r*', 'LineWidth', 4);
fplot3(@(v) par(1) + sin(v), @(v) par(2) + cos(v), @(v) m, 'r:', 'LineWidth', 2);
zlim([0, 3]);
xlabel('$\alpha$', 'Interpreter', 'Latex');
ylabel('$\beta$', 'Interpreter', 'Latex');
set(gca, 'TickLabelInterpreter', 'Latex');
box on
% legend('$V^R_{pen}$', '$V(\theta)$','Interpreter', 'Latex')
% legend('Interpreter', 'Latex')

% Set a different view for the second subplot
az2 = -30; % New azimuth for the second view
el2 = 20; % New elevation for the second view
view(az2, el2); % Apply the new view settings

%%
figure(3)
surf(par1,par2,obj3','EdgeColor','none');
xlabel('$\alpha$','Interpreter', 'Latex');
ylabel('$\beta$','Interpreter', 'Latex');
set(gca, 'TickLabelInterpreter', 'Latex')
az = 25;
el = 35; 
view(az, el); 

%%
% Create meshgrid for plotting
[x1, x2] = meshgrid(par1, par2);

figure(10)
surf(x1, x2, obj2','EdgeColor','none');
xlabel('$\alpha$', 'Interpreter', 'Latex');
ylabel('$\beta$', 'Interpreter', 'Latex');
zlabel('Objective', 'Interpreter', 'Latex');
% colorbar('TickLabelInterpreter', 'Latex');
set(gca, 'TickLabelInterpreter', 'Latex')
az = 20;
el = 50; 
view(az, el); 
box on 

figure(5)
contourf(x1, x2, obj2');
xlabel('$\alpha$', 'Interpreter', 'Latex');
ylabel('$\beta$', 'Interpreter', 'Latex');
colorbar('TickLabelInterpreter', 'Latex');
% title('\textbf{(b)}', 'Interpreter', 'Latex');
set(gca, 'TickLabelInterpreter', 'Latex')
box on
