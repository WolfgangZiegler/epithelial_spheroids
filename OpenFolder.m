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

%Part II (of IV) of 'spheroid polarity': MATLAB
%This macro is used to open the folder and results table from Part I
%calls the function readstack.m (Part III), which processes images, 
%and saves results for analysis in Part IV.

close all;
clear all;
%Input of the folder containing all midplane images and the ImageJ macro results file
path = uigetdir;
liste = dir(path);
cd(path);
files = {liste.name};

[results,direct_results] = uigetfile('*.csv', 'Select Spheroid Parameter Results');

%Initialisation of a Cell variable containing all results data
Cell_results=cell([numel(files)+1 16]);
Cell_results(1,1)={'Image_name'};
Cell_results(1,2)={'Difference_gp58_gp135'};
Cell_results(1,3)={'Difference_actin_gp135'};
Cell_results(1,4)={'Position_gp58'};
Cell_results(1,5)={'Position_gp135'};
Cell_results(1,6)={'Position_actin'};
Cell_results(1,7)={'Position_nuclei'};
Cell_results(1,8)={'Circularity_Spheroid'};
Cell_results(1,9)={'Circularity_Nuclei'};
Cell_results(1,10)={'Slope_Nuclei'};
Cell_results(1,11)={'Actin_Area_Amount'};
Cell_results(1,12)={'Actin_Particle_Count'};
Cell_results(1,13)={'Actin_Particle_Count_after_WS'};
Cell_results(1,14)={'Circularity_Actin'};
Cell_results(1,15)={'Distance_CoM_nuclei_spheroid'};
Cell_results(1,16)={'Distance_CoM_nuclei_actin'};
cell_r=cell([1 16]);

number=1;
%loop over all files, filenames should be longer than 10 characters
%if file is a tif, call of the function 'readstack'
%results are saved in the variable 'cell_results'
for i=1:numel(files)
    name=cell2string(files(i));
    if (length(name)>=10)
     name=name(3:length(name)-3);
     if (strcmp(name(length(name)-3:length(name)),'.tif'))
     number=number+1;    

     
     [cell_r]=readstack(name, results);

     Cell_results(number,1)=cell_r(1,1); %Image name
     Cell_results(number,2)=cell_r(1,2); %Difference Markers (gp58-gp135)
     Cell_results(number,3)=cell_r(1,3); %Difference actin-gp135 
     Cell_results(number,4)=cell_r(1,4); %gp58 Position
     Cell_results(number,5)=cell_r(1,5); %gp135 Position
     Cell_results(number,6)=cell_r(1,6); %actin Position
     Cell_results(number,7)=cell_r(1,7); %nuclei Position
     Cell_results(number,8)=cell_r(1,8); %Circularity Spheroid
     Cell_results(number,9)=cell_r(1,9); %Circularity nuclei
     Cell_results(number,10)=cell_r(1,10);%Slope nuclei
     Cell_results(number,11)=cell_r(1,11);%Amount Actin Area
     Cell_results(number,12)=cell_r(1,12);%Actin Particle Count
     Cell_results(number,13)=cell_r(1,13);%Actin Particle Count after Watershed
     Cell_results(number,14)=cell_r(1,14);%Cirtularity actin
     Cell_results(number,15)=cell_r(1,15);%Distance CoM nuclei-spheroid
     Cell_results(number,16)=cell_r(1,16);%Distance CoM nuclei-actin
     end
     save('Results','Cell_results');
    end
    
end

%after finishing the loop/analysis of all images
%all variables except the image name are saved in a table
%and the table is saved as text document

T=cell2table(Cell_results);
T.Properties.VariableNames=Cell_results(1,:);
writetable(T, 'ResultsTable of Parameters');






