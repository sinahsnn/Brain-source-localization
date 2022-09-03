clear, close all, clc ;
load('ElecPosXYZ') ;
%Forward Matrix
ModelParams.R = [8 8.5 9.2] ;
ModelParams.Sigma = [3.3e-3 8.25e-5 3.3e-3];
ModelParams.Lambda = [.5979 .2037 .0237];
ModelParams.Mu = [.6342 .9364 1.0362];
Resolution = 1 ;
[LocMat,GainMat] = ForwardModel_3shell(Resolution, ModelParams) ;
for i = 1:21
    ElecName{i} = ElecPos{1,i}.Name;
end
G = GainMat;
save("G.mat",'G')
save("LocMat.mat",'LocMat')
%% PART A
figure()
for i = 1:length(LocMat(1,:))
    plot3(LocMat(1,i),LocMat(2,i),LocMat(3,i),'*');
    hold on
end
title("possible dipoles with reslution" + num2str(Resolution));
xlabel('x');ylabel('y');zlabel('z');
%% Part B
figure()
for i = 1:length(LocMat(1,:))
    plot3(LocMat(1,i),LocMat(2,i),LocMat(3,i),'co');
    hold on
    if i<=21
        locations =ElecPos{1,i}.XYZ; 
        temp_pose = ModelParams.R(3)*locations;
        temp_name = ElecPos{1,i}.Name;
        plot3(temp_pose(1),temp_pose(2),temp_pose(3),'r*');
        text(temp_pose(1), temp_pose(2), temp_pose(3), num2str(temp_name))
        hold on
        xlabel('x');ylabel('y');zlabel('z');
        title("Dipoles & Electrodes (wirh names)");

    end
end
%% PART C 
selected = 1203;

loc_vec = [LocMat(1,selected) LocMat(2,selected) LocMat(3,selected)];
loc_vec = loc_vec/norm(loc_vec);
endpts = loc_vec +[LocMat(1,selected) LocMat(2,selected) LocMat(3,selected)];
figure()
for i = 1:length(LocMat(1,:))
    plot3(LocMat(1,i),LocMat(2,i),LocMat(3,i),'co');
    hold on
    if i<=21
        locations =ElecPos{1,i}.XYZ; 
        temp_pose = ModelParams.R(3)*locations;
        temp_name = ElecPos{1,i}.Name;
        plot3(temp_pose(1),temp_pose(2),temp_pose(3),'r*');
        text(temp_pose(1), temp_pose(2), temp_pose(3), num2str(temp_name))
        hold on
        xlabel('x');ylabel('y');zlabel('z');
        title("Dipoles & Electrodes (wirh names) and one selected dipole");
    end
end
plot3(LocMat(1,selected),LocMat(2,selected),LocMat(3,selected),'bd','linewidth',2);
hold on
plot3([0,1*LocMat(1,selected)],[0,1*LocMat(2,selected)],[0,1*LocMat(3,selected)],'--k','LineWidth',1);
plot3([1*LocMat(1,selected),endpts(1)],[1*LocMat(2,selected),endpts(2)],[1*LocMat(3,selected),endpts(3)],'r','LineWidth',1);
%% PART D 
load('Interictal.mat');
selected_idx = 1;
selected_interictal = Interictal(selected_idx,:);
Q = zeros(3951,10240);
Q((selected-1)*3+1:selected*3,:) = [selected_interictal*loc_vec(1);selected_interictal*loc_vec(2);selected_interictal*loc_vec(3)];
EEG = GainMat*Q;
disp_eeg(EEG,max(abs(EEG(:))),256,ElecName,'EEG signal');
%% PART E 
mean_ineric_m = zeros(21,5);
for i = 1:21
    mean_ineric_m(i,1) = mean(EEG(i,1761-3:1761+3));
    mean_ineric_m(i,2) = mean(EEG(i,5429-3:5429+3));
    mean_ineric_m(i,3) = mean(EEG(i,6302-3:6302+3));
    mean_ineric_m(i,4) = mean(EEG(i,8304-3:8304+3));
    mean_ineric_m(i,5) = mean(EEG(i,9417-3:9417+3));
end
mean_mean_ineric_m = abs(mean(mean_ineric_m,2));
M = mean(mean_ineric_m,2);
scale_means = rescale(mean_mean_ineric_m);

for i = 1:length(LocMat(1,:))
    plot3(LocMat(1,i),LocMat(2,i),LocMat(3,i),'co');
    hold on
    if i<=21
        locations =ElecPos{1,i}.XYZ; 
        temp_pose = ModelParams.R(3)*locations;
        temp_name = ElecPos{1,i}.Name;
        scatter3(temp_pose(1),temp_pose(2),temp_pose(3),50,scale_means(i),'filled');
        text(temp_pose(1), temp_pose(2), temp_pose(3), num2str(temp_name))
        hold on
        xlabel('x');ylabel('y');zlabel('z');
        title("Electrodes (wirh corresponding voltage color)");
    end
end
%% PART F
figure()
Display_Potential_3D(ModelParams.R(3),scale_means);
%% part G , H , I  shallow FOR MNE
%G
G = GainMat;
alpha = 0.5;
Q_MNE = G'*inv(G*G' + alpha*eye(21,21))*M;
% H
T = 1;
P = 1037;
moment_domain_MNE = zeros(P,T);
for i = 1:length(LocMat)
    for j = 1:T
        moment_domain_MNE(i,j) = norm(Q_MNE((i-1)*3+1:i*3,j));
    end
end
[max_MNE,idx_MNE] = max((moment_domain_MNE));
loc_mne_est = LocMat(:,idx_MNE);
%I
rmse_mne = norm([LocMat(1,selected) LocMat(2,selected) LocMat(3,selected)]-loc_mne_est');
disp(['MSE error in predicted location of MNE method is ' num2str(rmse_mne)]);
temp_mne_direct = Q_MNE(3*(idx_MNE-1)+1:3*idx_MNE,1);
temp_mne_direct = temp_mne_direct/norm(temp_mne_direct);
disp(['predicted direction of MNE method is ', num2str(temp_mne_direct')]);
disp(['predicted location of MNE method is ', num2str(loc_mne_est')]);
disp(['real location ', num2str([LocMat(1,selected) LocMat(2,selected) LocMat(3,selected)])]);
disp(['real direction ', num2str(loc_vec)]);
%%  part G , H , I  shallow  FOR WMNE
% G
p = 1317;
omega = zeros(p,p);
for i = 1:1317
    temp = 0;
    for j = 1:21
        temp = temp + G(j,(i-1)*3+1:i*3)*G(j,(i-1)*3+1:i*3)';
    end
    omega(i,i) = sqrt(temp);
end
W = kron(omega, eye(3));
Q_WMNE = inv(W'*W)*G'*inv(G*(inv(W'*W))*G'+alpha*eye(21,21))*M; %#ok
% H
T = 1;
P = 1037;
moment_domain_WMNE = zeros(P,T);
for i = 1:length(LocMat)
    for j = 1:T
        moment_domain_WMNE(i,j) = norm(Q_WMNE((i-1)*3+1:i*3,j));
    end
end
[max_WMNE , idx_WMNE] = max((moment_domain_WMNE));
loc_wmne_est = LocMat(:,idx_WMNE);
% I
rmse_wmne = norm([LocMat(1,selected) LocMat(2,selected) LocMat(3,selected)]-loc_wmne_est');
disp(['MSE error in predicted location of WMNE method is ' num2str(rmse_wmne)]);
temp_wmne_direct = Q_WMNE(3*(idx_WMNE-1)+1:3*idx_WMNE,1);
temp_wmne_direct = temp_wmne_direct/norm(temp_wmne_direct);

disp(['predicted direction of WMNE method is ', num2str(temp_wmne_direct')]);
disp(['predicted location of WMNE method is ', num2str(loc_wmne_est')]);
disp(['real location ', num2str([LocMat(1,selected) LocMat(2,selected) LocMat(3,selected)])]);
disp(['real direction ', num2str(loc_vec)]);

%% %%  part G , H , I  shallow FOR LORETA
%G

d = 1;
A = zeros(p,p);
for i = 1:p
    for j = 1:p
        dist = norm(LocMat(:,i) - LocMat(:,j));
        if dist == d
            A(i,j) = 1/6;
        end
    end
end

temp_A = A*ones(p);
A0 = inv(temp_A.*eye(size(temp_A)))*A; %#ok
B = (6/(d^2))*(kron(A0, eye(3)) - eye(3*p,3*p));
W = kron(omega, eye(3)) * B'*B*kron(omega, eye(3)); %#ok 
Q_LORETA = inv(W'*W)*G'*inv(G*(inv(W'*W))*G'+alpha*eye(21,21))*M; %#ok 
%H
T = 1;
P = 1037;
moment_domain_LORETA = zeros(P,T);
for i = 1:length(LocMat)
    for j = 1:T
        moment_domain_LORETA(i,j) = norm(Q_LORETA((i-1)*3+1:i*3,j));
    end
end
[~ , idx_LORETA] = max((moment_domain_LORETA));
loc_loreta_est = LocMat(:,idx_LORETA);
%I
rmse_loreta = norm([LocMat(1,selected) LocMat(2,selected) LocMat(3,selected)]-loc_loreta_est');
disp(['MSE error in predicted location of LoRETA method is ' num2str(rmse_loreta)]);

temp_loreta_direct = Q_LORETA(3*(idx_LORETA-1)+1:3*idx_LORETA,1);
temp_loreta_direct = temp_loreta_direct/norm(temp_loreta_direct);

%disp(['difference between angle of two vectors in loreta is equal to ' num2str(abs(angleVec(temp_loreta_direct,loc_vec)))]);
disp(['predicted direction of LORETA method is ', num2str(temp_loreta_direct')]);
disp(['predicted location of LORETA method is ', num2str(loc_loreta_est')]);
disp(['real location ', num2str([LocMat(1,selected) LocMat(2,selected) LocMat(3,selected)])]);
disp(['real direction ', num2str(loc_vec)]);
%% %%  part G , H , I  shallow FOR SLORETA
S_Q = G'*(G*G' + alpha*eye(21))^(-1)*G;
Q_SLORETA = zeros(3*p, 1);
for i = 1:p
    Q_SLORETA(3*i-2:3*i) = Q_MNE(3*i-2:3*i)'*S_Q(3*i-2:3*i, 3*i-2:3*i)^(-1)*Q_MNE(3*i-2:3*i);
end
%%%
T = 1;
P = 1037;
moment_domain_SLORETA = zeros(P,T);
for i = 1:length(LocMat)
    for j = 1:T
        moment_domain_SLORETA(i,j) = norm(Q_SLORETA((i-1)*3+1:i*3,j));
    end
end
[max_Sloreta , idx_SLORETA] = max((moment_domain_SLORETA));
loc_SLORETA_est = LocMat(:,idx_SLORETA);
% I
rmse_SLORETA = norm([LocMat(1,selected) LocMat(2,selected) LocMat(3,selected)]-loc_SLORETA_est');
disp(['MSE error in predicted location of SLORETA method is ' num2str(rmse_SLORETA)]);
temp_SLORETA_direct = Q_SLORETA(3*(idx_SLORETA-1)+1:3*idx_SLORETA,1);
temp_SLORETA_direct = temp_SLORETA_direct/norm(temp_SLORETA_direct);

%disp(['difference between angle of two vectors in SLORETA is equal to ' num2str(abs(angleVec(temp_SLORETA_direct,loc_vec)))]);
disp(['predicted direction of SLORETA method is ', num2str(temp_SLORETA_direct')]);
disp(['predicted location of SLORETA method is ', num2str(loc_SLORETA_est')]);
disp(['real location ', num2str([LocMat(1,selected) LocMat(2,selected) LocMat(3,selected)])]);
disp(['real direction ', num2str(loc_vec)]);
%% phase
disp(['difference between angle of two vectors in mne is equal to ' num2str(abs(angleVec(temp_mne_direct,loc_vec)))]);
disp(['difference between angle of two vectors in wmne is equal to ' num2str(abs(angleVec(temp_wmne_direct,loc_vec)))]);
disp(['difference between angle of two vectors in loreta is equal to ' num2str(abs(angleVec(temp_loreta_direct,loc_vec)))]);
disp(['difference between angle of two vectors in Sloretta is equal to ' num2str(abs(angleVec(temp_SLORETA_direct,loc_vec)))]);

%% DEEP DIPOLE 
%% PART C 
selected = 246;

loc_vec = [LocMat(1,selected) LocMat(2,selected) LocMat(3,selected)];
loc_vec = loc_vec/norm(loc_vec);
endpts = loc_vec +[LocMat(1,selected) LocMat(2,selected) LocMat(3,selected)];
figure()
for i = 1:length(LocMat(1,:))
    plot3(LocMat(1,i),LocMat(2,i),LocMat(3,i),'co');
    hold on
    if i<=21
        locations =ElecPos{1,i}.XYZ; 
        temp_pose = ModelParams.R(3)*locations;
        temp_name = ElecPos{1,i}.Name;
        plot3(temp_pose(1),temp_pose(2),temp_pose(3),'r*');
        text(temp_pose(1), temp_pose(2), temp_pose(3), num2str(temp_name))
        hold on
        xlabel('x');ylabel('y');zlabel('z');
        title("Dipoles & Electrodes (wirh names) and one selected dipole");
    end
end
plot3(LocMat(1,selected),LocMat(2,selected),LocMat(3,selected),'bd','linewidth',2);
hold on
plot3([0,1*LocMat(1,selected)],[0,1*LocMat(2,selected)],[0,1*LocMat(3,selected)],'--k','LineWidth',1);
plot3([1*LocMat(1,selected),endpts(1)],[1*LocMat(2,selected),endpts(2)],[1*LocMat(3,selected),endpts(3)],'r','LineWidth',1);
%% PART D 
load('Interictal.mat');
selected_idx = 1;
selected_interictal = Interictal(selected_idx,:);
Q = zeros(3951,10240);
Q((selected-1)*3+1:selected*3,:) = [selected_interictal*loc_vec(1);selected_interictal*loc_vec(2);selected_interictal*loc_vec(3)];
EEG = GainMat*Q;
disp_eeg(EEG,max(abs(EEG(:))),256,ElecName,'EEG signal');
%% PART E 
mean_ineric_m = zeros(21,5);
for i = 1:21
    mean_ineric_m(i,1) = mean(EEG(i,1762-3:1762+3));
    mean_ineric_m(i,2) = mean(EEG(i,5430-3:5430+3));
    mean_ineric_m(i,3) = mean(EEG(i,6302-3:6302+3));
    mean_ineric_m(i,4) = mean(EEG(i,8305-3:8304+5));
    mean_ineric_m(i,5) = mean(EEG(i,9416-3:9416+3));
end
mean_mean_ineric_m = abs(mean(mean_ineric_m,2));
M2 = mean(mean_ineric_m,2);
scale_means = rescale(mean_mean_ineric_m);

for i = 1:length(LocMat(1,:))
    plot3(LocMat(1,i),LocMat(2,i),LocMat(3,i),'co');
    hold on
    if i<=21
        locations =ElecPos{1,i}.XYZ; 
        temp_pose = ModelParams.R(3)*locations;
        temp_name = ElecPos{1,i}.Name;
        scatter3(temp_pose(1),temp_pose(2),temp_pose(3),50,scale_means(i),'filled');
        text(temp_pose(1), temp_pose(2), temp_pose(3), num2str(temp_name))
        hold on
        xlabel('x');ylabel('y');zlabel('z');
        title("Electrodes (wirh corresponding voltage color)");
    end
end
%% PART F
figure()
Display_Potential_3D(ModelParams.R(3),scale_means);
%% part G , H , I  DEEP FOR MNE
%G
G = GainMat;
alpha = 0.5;
Q_MNE = G'*inv(G*G' + alpha*eye(21,21))*M;
% H
T = 1;
P = 1037;
moment_domain_MNE = zeros(P,T);
for i = 1:length(LocMat)
    for j = 1:T
        moment_domain_MNE(i,j) = norm(Q_MNE((i-1)*3+1:i*3,j));
    end
end
[max_MNE,idx_MNE] = max((moment_domain_MNE));
loc_mne_est = LocMat(:,idx_MNE);
%I
rmse_mne = norm([LocMat(1,selected) LocMat(2,selected) LocMat(3,selected)]-loc_mne_est');
disp(['MSE error in estimated location of MNE method is ' num2str(rmse_mne)]);
temp_mne_direct = Q_MNE(3*(idx_MNE-1)+1:3*idx_MNE,1);
temp_mne_direct = temp_mne_direct/norm(temp_mne_direct);
%disp(['difference between angle of two vectors in mne is equal to ' num2str(abs(angleVec(temp_mne_direct,loc_vec)))]);
disp(['estimated direction of MNE method is ', num2str(temp_mne_direct')]);
disp(['estimated location of MNE method is ', num2str(loc_mne_est')]);
disp(['real location ', num2str([LocMat(1,selected) LocMat(2,selected) LocMat(3,selected)])]);
disp(['real direction ', num2str(loc_vec)]);
%%  part G , H , I  DEEP  FOR WMNE
% G
p = 1317;
omega = zeros(p,p);
for i = 1:1317
    temp = 0;
    for j = 1:21
        temp = temp + G(j,(i-1)*3+1:i*3)*G(j,(i-1)*3+1:i*3)';
    end
    omega(i,i) = sqrt(temp);
end
W = kron(omega, eye(3));
Q_WMNE = inv(W'*W)*G'*inv(G*(inv(W'*W))*G'+alpha*eye(21,21))*M; %#ok
% H
T = 1;
P = 1037;
moment_domain_WMNE = zeros(P,T);
for i = 1:length(LocMat)
    for j = 1:T
        moment_domain_WMNE(i,j) = norm(Q_WMNE((i-1)*3+1:i*3,j));
    end
end
[max_WMNE , idx_WMNE] = max((moment_domain_WMNE));
loc_wmne_est = LocMat(:,idx_WMNE);
% I
rmse_wmne = norm([LocMat(1,selected) LocMat(2,selected) LocMat(3,selected)]-loc_wmne_est');
disp(['MSE error in estimated location of WMNE method is ' num2str(rmse_wmne)]);
temp_wmne_direct = Q_WMNE(3*(idx_WMNE-1)+1:3*idx_WMNE,1);
temp_wmne_direct = temp_wmne_direct/norm(temp_wmne_direct);

%disp(['difference between angle of two vectors in wmne is equal to ' num2str(abs(angleVec(temp_wmne_direct,loc_vec)))]);
disp(['estimated direction of WMNE method is ', num2str(temp_wmne_direct')]);
disp(['estimated location of WMNE method is ', num2str(loc_wmne_est')]);
disp(['real location ', num2str([LocMat(1,selected) LocMat(2,selected) LocMat(3,selected)])]);
disp(['real direction ', num2str(loc_vec)]);
%%
disp(['difference between angle of two vectors in mne is equal to ' num2str(abs(angleVec(temp_mne_direct,loc_vec)))]);
disp(['difference between angle of two vectors in wmne is equal to ' num2str(abs(angleVec(temp_wmne_direct,loc_vec)))]);
disp(['difference between angle of two vectors in loreta is equal to ' num2str(abs(angleVec(temp_loreta_direct,loc_vec)))]);
disp(['difference between angle of two vectors in Sloretta is equal to ' num2str(abs(angleVec(temp_SLORETA_direct,loc_vec)))]);

%%   part G , H , I  DEEP FOR LORETA
%G
d = 1;
A = zeros(p,p);
for i = 1:p
    for j = 1:p
        dist = norm(LocMat(:,i) - LocMat(:,j));
        if dist == d
            A(i,j) = 1/6;
        end
    end
end

temp_A = A*ones(p);
A0 = inv(temp_A.*eye(size(temp_A)))*A; %#ok
B = (6/(d^2))*(kron(A0, eye(3)) - eye(3*p,3*p));
W = kron(omega, eye(3)) * B'*B*kron(omega, eye(3)); %#ok 
Q_LORETA = inv(W'*W)*G'*inv(G*(inv(W'*W))*G'+alpha*eye(21,21))*M; %#ok 
%H
T = 1;
P = 1037;
moment_domain_LORETA = zeros(P,T);
for i = 1:length(LocMat)
    for j = 1:T
        moment_domain_LORETA(i,j) = norm(Q_LORETA((i-1)*3+1:i*3,j));
    end
end
[~ , idx_LORETA] = max((moment_domain_LORETA));
loc_loreta_est = LocMat(:,idx_LORETA);
%I
rmse_loreta = norm([LocMat(1,selected) LocMat(2,selected) LocMat(3,selected)]-loc_loreta_est');
disp(['MSE error in estimated location of LoRETA method is ' num2str(rmse_loreta)]);

temp_loreta_direct = Q_LORETA(3*(idx_LORETA-1)+1:3*idx_LORETA,1);
temp_loreta_direct = temp_loreta_direct/norm(temp_loreta_direct);

%disp(['difference between angle of two vectors in loreta is equal to ' num2str(abs(angleVec(temp_loreta_direct,loc_vec)))]);
disp(['estimated direction of LORETA method is ', num2str(temp_loreta_direct')]);
disp(['estimated location of LORETA method is ', num2str(loc_loreta_est')]);
disp(['real location ', num2str([LocMat(1,selected) LocMat(2,selected) LocMat(3,selected)])]);
disp(['real direction ', num2str(loc_vec)]);
%% part G , H , I  DEEP FOR SLORETA
S_Q = G'*(G*G' + alpha*eye(21))^(-1)*G;
Q_SLORETA = zeros(3*p, 1);
for i = 1:p
    Q_SLORETA(3*i-2:3*i) = Q_MNE(3*i-2:3*i)'*S_Q(3*i-2:3*i, 3*i-2:3*i)^(-1)*Q_MNE(3*i-2:3*i);
end
%%%
T = 1;
P = 1037;
moment_domain_SLORETA = zeros(P,T);
for i = 1:length(LocMat)
    for j = 1:T
        moment_domain_SLORETA(i,j) = norm(Q_SLORETA((i-1)*3+1:i*3,j));
    end
end
[max_Sloreta , idx_SLORETA] = max((moment_domain_SLORETA));
loc_SLORETA_est = LocMat(:,idx_SLORETA);
% I
rmse_SLORETA = norm([LocMat(1,selected) LocMat(2,selected) LocMat(3,selected)]-loc_SLORETA_est');
disp(['MSE error in estimated location of SLORETA method is ' num2str(rmse_SLORETA)]);
temp_SLORETA_direct = Q_SLORETA(3*(idx_SLORETA-1)+1:3*idx_SLORETA,1);
temp_SLORETA_direct = temp_SLORETA_direct/norm(temp_SLORETA_direct);

%disp(['difference between angle of two vectors in SLORETA is equal to ' num2str(abs(angleVec(temp_SLORETA_direct,loc_vec)))]);
disp(['estimated direction of SLORETA method is ', num2str(temp_SLORETA_direct')]);
disp(['estimated location of SLORETA method is ', num2str(loc_SLORETA_est')]);
disp(['real location ', num2str([LocMat(1,selected) LocMat(2,selected) LocMat(3,selected)])]);
disp(['real direction ', num2str(loc_vec)]);

%% PART J : parametric model 
peredicted_source1 = 1098;
Q_predicted1 = 9.347;
%% not fixed 
peredicted_source2 = 1098;
Q_predicted2 = [-10.671,3.672,1.272]
%% part K 
Set_selec = [1143 1156 1130 1119 1107 1131 1179 1167 1155 1229 1210 1250 1219 1240 1220 1201 1260 1239];
Set_loc = LocMat(:, Set_selec);
figure()
for i = 1:length(LocMat(1,:))
    plot3(LocMat(1,i),LocMat(2,i),LocMat(3,i),'co');
    hold on
	if i <= 21
		temp_pose = ModelParams.R(3)*ElecPos{1,i}.XYZ;
		temp_name = ElecPos{1,i}.Name;
		plot3(temp_pose(1),temp_pose(2),temp_pose(3),'r*');
		text(temp_pose(1), temp_pose(2), temp_pose(3), num2str(temp_name))
		hold on
	end
	if i <= length(Set_selec)
		scatter3(Set_loc(1,i),Set_loc(2,i), Set_loc(3,i), 'bd',"linewidth",2);
		title('locations with selected dipole');
		xlabel('x');ylabel('y');zlabel('z');
		hold on
		plot3([0,1.2*Set_loc(1,i)],[0,1.2*Set_loc(2,i)],[0,1.2*Set_loc(3,i)],'y','LineWidth',2);
		hold on
	end
end
%% part L
load('Interictal.mat');
Q = zeros(3951,10240);
for selected_idx = 1 : 18
    selected_idx = 1;
    selected_interictal = Interictal(selected_idx,:);
    selected = Set_selec(selected_idx);
    loc_vec = [LocMat(1,selected) LocMat(2,selected) LocMat(3,selected);];
    loc_vec = loc_vec/norm(loc_vec);
    Q((selected-1)*3+1:selected*3,:) = [selected_interictal*loc_vec(1);selected_interictal*loc_vec(2);selected_interictal*loc_vec(3)];
end
EEG = GainMat*Q;
disp_eeg(EEG,max(abs(EEG(:))),256,ElecName,'EEG signal');

%% 

mean_ineric_m = zeros(21,4);
for i = 1:21
    mean_ineric_m(i,1) = mean(EEG(i,1746-3:1746+3));
    mean_ineric_m(i,2) = mean(EEG(i,5381-3:5381+3));
    mean_ineric_m(i,3) = mean(EEG(i,8306-3:8306+3));
    mean_ineric_m(i,4) = mean(EEG(i,9350-3:9350+3));
end
mean_mean_ineric_m = abs(mean(mean_ineric_m,2));
M = mean(mean_ineric_m,2);
scale_means = rescale(mean_mean_ineric_m);

for i = 1:length(LocMat(1,:))
    plot3(LocMat(1,i),LocMat(2,i),LocMat(3,i),'co');
    hold on
    if i<=21
        locations =ElecPos{1,i}.XYZ; 
        temp_pose = ModelParams.R(3)*locations;
        temp_name = ElecPos{1,i}.Name;
        scatter3(temp_pose(1),temp_pose(2),temp_pose(3),50,scale_means(i),'filled');
        text(temp_pose(1), temp_pose(2), temp_pose(3), num2str(temp_name))
        hold on
        xlabel('x');ylabel('y');zlabel('z');
        title("Electrodes (wirh corresponding voltage color)");
    end
end

%
figure()
Display_Potential_3D(ModelParams.R(3),scale_means);
%% PART M , N MNE algorithm
%M
G = GainMat;
alpha = 0.5;
Q_MNE = G'*inv(G*G' + alpha*eye(21,21))*M; %#ok
%
T = 1;
P = 1037;
moment_domain_MNE = zeros(P,T);
for i = 1:length(LocMat)
    for j = 1:T
        moment_domain_MNE(i,j) = norm(Q_MNE((i-1)*3+1:i*3,j));
    end
end
%N
label_set = zeros(length(LocMat),1);
label_set(Set_selec) = 1;
figure()
plotroc(label_set',moment_domain_MNE');
[tpr_mne,fpr_mne,thresholds_mne] = roc(label_set',moment_domain_MNE');

title('ROC for MNE');
figure()
subplot(2,1,1)
plot(thresholds_mne,tpr_mne);
title('TPR for mne based on threshould',"linewidth",1);
xlabel('Thresholds');
ylabel('TPR');
grid on
subplot(2,1,2)
plot(thresholds_mne,fpr_mne,"linewidth",1);
title('FPR for mne based on threshould');
xlabel('Thresholds');
ylabel('FPR');
grid on
%% PART M , N WMNE algorithm
%M
p = 1317;
omega = zeros(p,p);
for i = 1:1317
    temp = 0;
    for j = 1:21
        temp = temp + G(j,(i-1)*3+1:i*3)*G(j,(i-1)*3+1:i*3)';
    end
    omega(i,i) = sqrt(temp);
end
W = kron(omega, eye(3));
Q_WMNE = inv(W'*W)*G'*inv(G*(inv(W'*W))*G'+alpha*eye(21,21))*M; %#ok
%
T = 1;
P = 1037;
moment_domain_WMNE = zeros(P,T);
for i = 1:length(LocMat)
    for j = 1:T
        moment_domain_WMNE(i,j) = norm(Q_WMNE((i-1)*3+1:i*3,j));
    end
end
%N

label_set = zeros(length(LocMat),1);
label_set(Set_selec) = 1;
figure()
plotroc(label_set',moment_domain_WMNE');
title('ROC for WMNE');
[tpr_wmne,fpr_wmne,thresholds_wmne] = roc(label_set',moment_domain_WMNE');


figure()
subplot(2,1,1)
plot(thresholds_wmne,tpr_wmne,"linewidth",1);
title('TPR for wmne based on threshould');
xlabel('Thresholds');
ylabel('TPR');
grid on 
subplot(2,1,2)
plot(thresholds_wmne,fpr_wmne);
title('FPR for wmne based on threshould');
xlabel('Thresholds');
ylabel('FPR');
grid on
%% PART M , N LORETA
%M
d = 1;
A = zeros(p,p);
for i = 1:p
    for j = 1:p
        dist = norm(LocMat(:,i) - LocMat(:,j));
        if dist == d
            A(i,j) = 1/6;
        end
    end
end
temp_A = A*ones(p);
A0 = inv(temp_A.*eye(size(temp_A)))*A; %#ok

B = (6/(d^2))*(kron(A0, eye(3)) - eye(3*p,3*p));

W = kron(omega, eye(3)) * B'*B*kron(omega, eye(3)); %#ok 
Q_LORETA = inv(W'*W)*G'*inv(G*(inv(W'*W))*G'+alpha*eye(21,21))*M; %#ok 
%
T = 1;
P = 1037;
moment_domain_LORETA = zeros(P,T);
for i = 1:length(LocMat)
    for j = 1:T
        moment_domain_LORETA(i,j) = norm(Q_LORETA((i-1)*3+1:i*3,j));
    end
end
% N

label_set = zeros(length(LocMat),1);
label_set(Set_selec) = 1;
figure()
plotroc((~label_set)',moment_domain_LORETA');
title('ROC for LORETA');
[tpr_loreta,fpr_loreta,thresholds_loreta] = roc(label_set',-moment_domain_LORETA');

figure()
subplot(2,1,1)
plot(thresholds_loreta,tpr_loreta);
title('TPR for loreta based on threshould');
xlabel('Thresholds');
ylabel('TPR');
grid on
subplot(2,1,2)
plot(thresholds_loreta,fpr_loreta);
title('FPR for loreta based on threshould');
xlabel('Thresholds');
ylabel('FPR');
grid on

%% PARTM m,N SLORETA
%M
S_Q = G'*(G*G' + alpha*eye(21))^(-1)*G;
Q_SLORETA = zeros(3*p, 1);
for i = 1:p
    Q_SLORETA(3*i-2:3*i) = Q_MNE(3*i-2:3*i)'*S_Q(3*i-2:3*i, 3*i-2:3*i)^(-1)*Q_MNE(3*i-2:3*i);
end

%%%
T = 1;
P = 1037;
moment_domain_SLORETA = zeros(P,T);
for i = 1:length(LocMat)
    for j = 1:T
        moment_domain_SLORETA(i,j) = norm(Q_SLORETA((i-1)*3+1:i*3,j));
    end
end
label_set = zeros(length(LocMat),1);
label_set(Set_selec) = 1;
figure()
plotroc((~label_set)',-moment_domain_SLORETA');
title('ROC for SLORETA');
%N
[tpr_Sloreta,fpr_Sloreta,thresholds_Sloreta] = roc(label_set',-moment_domain_SLORETA');

figure()
subplot(2,1,1)
plot(thresholds_Sloreta,tpr_Sloreta);
title('TPR for sloreta based on threshould');
xlabel('Thresholds');
ylabel('TPR');
grid on
subplot(2,1,2)
plot(thresholds_Sloreta,fpr_Sloreta);
title('FPR for sloreta based on threshould');
xlabel('Thresholds');
ylabel('FPR');
grid on
function theta = angleVec(u,v)
CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
theta  = real(acosd(CosTheta));
end