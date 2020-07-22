%%%%%%%%%%%%%% Automatic Image Analysis - real-time version %%%%%%%%%%%%%
%{ 
data structure
    CurrFrame: [particle index, x( in unit of a), y]

    P_Type: partical category (<0 perpheral particles (unused), =-3 grain
    boundary, >0 crystalline particles differentiate by domains indice

    GB_INFO: [grain boundary orientation (-67,112), grain boundary center x,
    grain boundary center y]. GB orientation between +-90 from 22.5 field
    direction

    DOMAIN_INFO: [domain center x, domain center y, domain orientation 
    (-30 ~ +30), domain size] captured up to two domains

    G_PSI6: global psi_6
    G_C6: global C_6
    
    L_PSI6: local psi_6 of each particle
    L_C6: local C_6 of each particle
    L_theta6: local lattice orientation of each particle
%}
function [P_Type, G_PSI6, G_C6, L_PSI6, L_C6] = IMAGE_ANALYSIS(CurrFrame)
%% format the input data
CurrFrame(:,1) = (1:size(CurrFrame,1))';
    
%% differentiate particle types
[P_Type, neighboridx, L_C6, G_PSI6, G_C6,L_PSI6] = CATEGORIZATION(CurrFrame);
    
%% identify major grain boundary
if ~isempty(CurrFrame(P_Type == 0,:))
    GB_list = CONNECTIVITY2(CurrFrame(P_Type == 0,1), neighboridx);
%     for i = 1:length(GB_list)
%         if length(GB_list{i})>6
%             GB = cat(2,GB,GB_list{i});
%         end
%     end
    GB = GB_list{1};
    P_Type(GB) = -3;
end
end
