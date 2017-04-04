%Copyright (c) 2017 Birga Soetje.
% * 
% * This script is part of the 'spheroid polarity' source code/program.
% *
% * This program is free software: you can redistribute it and/or modify  
% * it under the terms of the GNU General Public License as   
% * published by the Free Software Foundation, version 3.
% *
% * This program is distributed in the hope that it will be useful, but 
% * WITHOUT ANY WARRANTY; without even the implied warranty of 
% * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
% * General Public License for more details.
% *
% * You should have received a copy of the GNU General Public License
% * along with this program. If not, see <http://www.gnu.org/licenses/>.

%Part IV (of IV) of 'spheroid polarity': MATLAB
%This macro is used to classify spheroids based on the results from Parts II and III 
%and a manually trained decision tree. Necessarity for reclassification is checked 
%and results including classification and reclassification are saved. 
 
close all;
clear all;
 %Input of the files: reduced results table and decision tree
 [Results,path_res] = uigetfile('*.mat','Select the Results mat file');
 [DecTree,path_tree] = uigetfile('*.mat','Select the Decision Tree mat file');
 cd(path_tree);
 DecisionTree=load(DecTree);
 name = fieldnames(DecisionTree);
 Tree_Str=strcat('Tree=DecisionTree.',name{1,1},';');
 eval(Tree_Str);
 
 cd(path_res);
 res=load(Results);
 name = fieldnames(res);
 Res_Str=strcat('res_cell=res.',name{1,1},';');
 eval(Res_Str);
  
% clearvars res DecisionTree path_tree

% At first the resultsing parameters of each picture are copied 
% into an array of numbers without containing the image name
% resulting array is <res>

 [m,n]=size(res_cell);
 %if res_cell does not contain 3 emty columns, change to 'res_red=zeros(m-1,n-1)'
 res_red=zeros(m-4,n-1); 
 [m,n]=size(res_red);
for i=1:m
     for j=1:n
         res_red(i,j)=res_cell{i+1,j+1};
     end
end

% The resulting array is reduced by eleminating the columns 2,3 and 4
% because the classification was more appropiate using only 
% the residual 12 columns/parameters
% This step is adaptable for another set of parameters for classification
res_red2=res_red;
res_red2(:,2:4)=[];
%%
% Prediction of the polarity group and saving of the results in an 
% additional column (n) in the variable <res_class>
yfit=predict(Tree,res_red2);
res_class=res_red;
[m,n]=size(res_red);
res_class(:,n+1)=yfit;
for (i=1:length(yfit))
res_cell{i+1,17}=yfit(i,1);
end
res_cell(1,17)={'First_Classification'};

%% From here the reclassification progress begins.
%If initial classification is unequivocal the reclassification step can be skipped by converting to a comment

cutoffs= [0.006181217 0.59885 1.502288487 0.324035 0.782];
res_class(:,17)=0;
[a,b]=size(res_class);
for (i=1:a)
    if res_red2(i,1)<cutoffs(1)
        res_class(i,17)=res_class(i,17)+1
    end
    if res_red2(i,4)<cutoffs(2)
        res_class(i,17)=res_class(i,17)+1
    end
    if res_red2(i,5)>cutoffs(3)
        res_class(i,17)=res_class(i,17)+1
    end
    if res_red2(i,6)>cutoffs(4)
        res_class(i,17)=res_class(i,17)+1
    end
    if res_red2(i,7)>cutoffs(5)
        res_class(i,17)=res_class(i,17)+1
    end
end
%%
potWrong=0;
res_class(:,18)=0;
for (i=1:a)
    if ((res_class(i,17)<2 & res_class(i,16)==3) | (res_class(i,17)>2 & res_class(i,16)<=2))
        res_class(i,18)=1;
        potWrong=potWrong+1;
    end
end

%% 
yes=0; 
promptMessage = strcat(num2str(potWrong),' of  ',num2str(a),' spheroids are potentially clasified wrong. Do you want to continue with reclassification?');
button = questdlg(promptMessage, 'Continue', 'Continue', 'Cancel', 'Continue');
if strcmpi(button, 'Cancel')
	yes=1;
end
 save('Results_Classified','res_cell');
%% 
if yes==0
    sort_gr1_question=uigetfile('*.mat','Select the DecisionTree mat file for questioning the reclassification of group1 spheroids');
    sort_gr1=uigetfile('*.mat','Select the DecisionTree mat file for reclassification of group1 spheropids');
    sort_gr2_question=uigetfile('*.mat','Select the DecisionTree mat file for questioning the reclassification of group2 spheroids');
    sort_gr2=uigetfile('*.mat','Select the DecisionTree mat file for reclassification of group2 spheropids');
    sort_gr3_question=uigetfile('*.mat','Select the DecisionTree mat file for questioning the reclassification of group3 spheroids');
    sort_gr3=uigetfile('*.mat','Select the DecisionTree mat file for reclassification of group3 spheropids');
%%
    gr1_d=load(sort_gr1);
    name = fieldnames(gr1_d);
    gr1_decision=strcat('gr1_dec=gr1_d.',name{1,1},';');
    eval(gr1_decision);
    
    gr1_q=load(sort_gr1_question);
    name = fieldnames(gr1_q);
    gr1_question=strcat('gr1_ques=gr1_q.',name{1,1},';');
    eval(gr1_question);
    
    gr2_d=load(sort_gr2);
    name = fieldnames(gr2_d);
    gr2_decision=strcat('gr2_dec=gr2_d.',name{1,1},';');
    eval(gr2_decision);
    
    gr2_q=load(sort_gr2_question);
    name = fieldnames(gr2_q);
    gr2_question=strcat('gr2_ques=gr2_q.',name{1,1},';');
    eval(gr2_question);
    
    gr3_d=load(sort_gr3);
    name = fieldnames(gr3_d);
    gr3_decision=strcat('gr3_dec=gr3_d.',name{1,1},';');
    eval(gr3_decision);
    
    gr3_q=load(sort_gr3_question);
    name = fieldnames(gr3_q);
    gr3_question=strcat('gr3_ques=gr3_q.',name{1,1},';');
    eval(gr3_question);
    
%%   
% Reduction of the <res_class> array to 8 parameters, 
% informative enough to decide if a reclassification is necessary <res_reclass>
% If a reclassification is necessary, in an additional column (y+1) of  <res_class> 
% the value 0 indicates no necessarity and 1 the necessarity of reclassification
% result of the classification is saved in one additional column (y+2) of  <res_class>

res_reclass=res_class;
[x,y]=size(res_class);
res_reclass(:,17:end)=[];
res_reclass(:,14:15)=[];
res_reclass(:,12)=[];
res_reclass(:,3:6)=[];
%% 
%variable <res_reclass_q> contains all parameters except the classification result

res_reclass_q=res_reclass(:,1:end-1);
[m,n]=size(res_reclass_q);
%%
% for determining the reclassification necessarity or performing the reclassification
% often just a subset of parameters appeared to be most effective

for i=1:m
    if res_reclass(i,n+1)==1;
        tmp=res_reclass_q(i,:);
        q_sort=predict(gr1_ques,tmp);
        res_class(i,y+1)=q_sort;
        if q_sort==1;
            tmp2=zeros(1,3); %array for subset of parameters
            tmp2(1,1)=tmp(1,4);
            tmp2(1,2)=tmp(1,5); 
            tmp2(1,3)=tmp(1,8);
            res_sort=predict(gr1_dec,tmp2);
            res_class(i,y+2)=res_sort;
            clearvars tmp tmp2 q_sort; 
        else res_class(i,y+2)=0;
        end
    
    elseif res_reclass(i,n+1)==2;
        tmp=res_reclass_q(i,:);
        q_sort=predict(gr2_ques,tmp);
        res_class(i,y+1)=q_sort;
        if q_sort==1;
            res_sort=predict(gr2_dec,tmp);
            res_class(i,y+2)=res_sort;
            clearvars tmp q_sort;
        else res_class(i,y+2)=0;
        end        
    else
        if res_reclass(i,n+1)==3;
            tmp=res_reclass_q(i,:);
            tmp2=zeros(1,3);
            tmp2(1,1)=tmp(1,4);
            tmp2(1,2)=tmp(1,5); 
            tmp2(1,3)=tmp(1,8);
            q_sort=predict(gr3_ques,tmp2);
            res_class(i,y+1)=q_sort;
            if q_sort==1;
                tmp3=zeros(1,2);
                tmp3(1,1)=tmp(1,5); 
                tmp3(1,2)=tmp(1,8);
                res_sort=predict(gr3_dec,tmp3);
                res_class(i,y+2)=res_sort;
                clearvars tmp tmp2 tmp3 q_sort;
            else res_class(i,y+2)=0;
            end
        end
    end
end
%%

res_cell(1,18)={'Score_possible_wrong_classified'};
res_cell(1,19)={'possible_wrong_classified'};
res_cell(1,20)={'Reclassification_from_tree_recommended'};
res_cell(1,21)={'Reclassification'};
res_cell(1,22)={'End_Classification'};
res_cell(1,23)={'Comments'};
for (i=1:a)
    res_cell{i+1,18}=res_class(i,17);
	res_cell{i+1,19}=res_class(i,18);
    res_cell{i+1,20}=res_class(i,19);
    res_cell{i+1,21}=res_class(i,20);
end

for (i=1:a)
   if res_class(i,20)==0
    res_cell{i+1,22}=res_class(i,16);
   elseif ((res_class(i,20)==res_class(i,16)) & (res_class(i,18)+res_class(i,19)==2))
        res_cell{i+1,23}={'Reclassification_approved_group'};   
        res_cell{i+1,22}=res_class(i,20);
   elseif ((res_class(i,20)~=res_class(i,16)) & (res_class(i,18)+res_class(i,19)==2))
        res_cell{i+1,23}={'Reclassified_to_another_group'};
        res_cell{i+1,22}=res_class(i,20);
   elseif ((res_class(i,20)~=res_class(i,16)) & (res_class(i,17)>=2) & (res_class(i,20)==1) )
        res_cell{i+1,23}={'Reclassified_to_group1_but_reevaluation_recommended'};
        res_cell{i+1,22}=res_class(i,20);
   else res_cell{i+1,22}=res_class(i,20);
   end
end

end
%%
save('Results_Classified','res_cell');
T=cell2table(res_cell);
T.Properties.VariableNames=res_cell(1,:);
writetable(T, 'EndResults.csv');







        