% Load the simulated and experimental data
power = load('powersimulations');
data = load('experimentaldata');

% Prepare the tables
psi_table=zeros(3,6,4);

% There are three tables, one for each type of simulation
for agent_type=1:3 

% Each table row is PSi for a given theory, each column - treatment in the
% experiment
psi_table(agent_type,1,:)=(data.count_PO./data.num_obs - power.count_PO(agent_type)/power.max_iter)*100;
psi_table(agent_type,2,:)=(data.count_lott_PO./data.num_obs - power.count_lott_PO(agent_type)/power.max_iter)*100;

% Indices of efficient assinment [312] in the list of matchings in the two
% .mat files
[~, power312ind]=ismember([3 1 2], power.sorted_pairings,'rows');
[~, data312ind]=ismember([3 1 2], data.sorted_pairings,'rows');

% Remaining rows
psi_table(agent_type,3,:)=(data.sorted_pairings_cnt(data312ind,:)./data.num_obs - power.sorted_pairings_cnt(power312ind,agent_type)/power.max_iter)*100;
psi_table(agent_type,4,:)=(data.count_ce./data.num_obs - power.count_ce(agent_type)/power.max_iter)*100;
psi_table(agent_type,5,:)=(data.count_ne./data.num_obs - power.count_ne(agent_type)/power.max_iter)*100;
psi_table(agent_type,6,:)=(data.count_Ncore./data.num_obs - power.count_Ncore(agent_type)/power.max_iter)*100;

% Produce three heatmap tables within MatLab
subplot(1,3,agent_type);
hm = heatmap(squeeze(psi_table(agent_type,:,:)));
hm.YData={'PO' 'PO lott', '312', 'ce', 'ne', 'PN-core'};
hm.XData=data.treatment_names;
end

% Print the tables as LaTeX (Table F.8)
for agent_type=1:3
    sprd=max(max(max(psi_table(agent_type,:,:))))-min(min(min(psi_table(agent_type,:,:))));
    for r=1:6
        row_to_print=zeros(1,8);
        row_to_print(2:2:8)=psi_table(agent_type,r,:);
        row_to_print(1:2:7)=round((psi_table(agent_type,r,:)-min(min(min(psi_table(agent_type,:,:)))))/sprd/2*100);
        fprintf('{\\cellcolor{black!%d} %0.2f} & {\\cellcolor{black!%d}%0.2f} & {\\cellcolor{black!%d}%0.2f} & {\\cellcolor{black!%d}%0.2f} \n', row_to_print)
    end
end