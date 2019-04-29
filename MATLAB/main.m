close all; clear all; clc;

%% Load Coactivation Matrix
load('./MatlabData/Coactivation_matrix.mat')
load('./MatlabData/Coactivation_matrix_coord.mat')

%% Plot Coordinates Coactivation Matrix
% figure;
% scatter3(Coactivation_matrix_coord(:,1),Coactivation_matrix_coord(:,2),Coactivation_matrix_coord(:,3))
% val = 120;
% xlim([-val val]);
% ylim([-val val]);
% zlim([-val val]);

%% Plot Network Coactivation Matrix with communities

figure;
A = Coactivation_matrix;
threshold = 0.4*max(max(A));
for i=1:length(A)
for j=1:length(A)
if A(i,j) < threshold
   A(i,j) = 0;
end
end
end

G_coact_un = graph(A);
G_coact_w = graph(Coactivation_matrix);

% Plot 
p = plot(G_coact_un, 'XData', Coactivation_matrix_coord(:,1), 'YData', Coactivation_matrix_coord(:,2), 'ZData', Coactivation_matrix_coord(:,3));
p.Marker = 'o';
p.NodeColor = 'r';
p.MarkerSize = 10;
title('Coactivation Network','FontSize',18)
grid on
% % Nodes with highest degree
% deg = sum(Coactivation_matrix,2);
% [sortdeg, indsortdeg] = sort(deg);
% highlight(p,indsortdeg(end-10:end),'NodeColor',[0, 0.5, 0])

[list1, list2, list3, list4, list5] = communitiesCoact();
highlight(p,list1,'NodeColor','r')
highlight(p,list2,'NodeColor','b')
highlight(p,list3,'NodeColor',[0, 0.5, 0])
highlight(p,list4,'NodeColor','c')
highlight(p,list5,'NodeColor','y')

%write_matrix_to_pajek(Coactivation_matrix,'./NetworksPajek/Coactivation.net','weighted',true,'directed',false,'coords',Coactivation_matrix_coord);

%% Load Resting-State Connection Matrix
load('./MatlabData/RestingState_matrix.mat')
load('./MatlabData/RestingState_matrix_coord.mat')

%% Plot Coordinates Resting-State Connection Matrix
% figure;
% scatter3(RestingState_matrix_coord(:,1),RestingState_matrix_coord(:,2),RestingState_matrix_coord(:,3))
% val = 120;
% xlim([-val val]);
% ylim([-val val]);
% zlim([-val val]);

%% Plot Network Resting-State Connection Matrix with communities
figure;
A = RestingState_matrix;
threshold = 0.6*max(max(A));
for i=1:length(A)
for j=1:length(A)
if A(i,j) < threshold
   A(i,j) = 0;
end
end
end
G_rest_un = graph(A);
G_rest_w = graph(RestingState_matrix);

% Plot
p = plot(G_rest_un, 'XData', RestingState_matrix_coord(:,1), 'YData', RestingState_matrix_coord(:,2), 'ZData', RestingState_matrix_coord(:,3));
p.Marker = 'o';
p.NodeColor = 'r';
p.MarkerSize = 8;
title('Resting State Network','FontSize',18)
grid on

[list1,list2,list3,list4,list5,list6,list7,list8,list9,list10] = communitiesRest();
highlight(p,list1,'NodeColor','r')
highlight(p,list2,'NodeColor','b')
highlight(p,list3,'NodeColor','g')
highlight(p,list4,'NodeColor','c')
highlight(p,list5,'NodeColor','y')

highlight(p,list6,'NodeColor',[0.6350, 0.0780, 0.1840])
highlight(p,list7,'NodeColor',[0.8500, 0.3250, 0.0980])
highlight(p,list8,'NodeColor',[0, 0.5, 0])
highlight(p,list9,'NodeColor',[0.75, 0.75, 0])
highlight(p,list10,'NodeColor',[0.25, 0.25, 0.25])

% % Nodes with highest degree
% deg = sum(RestingState_matrix,2);
% [sortdeg, indsortdeg] = sort(deg);
% highlight(p,indsortdeg(end-10:end),'NodeColor',[0, 0.5, 0])

%write_matrix_to_pajek(RestingState_matrix,'./NetworksPajek/RestingState.net','weighted',true,'directed',false,'coords',RestingState_matrix_coord);