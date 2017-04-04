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

%Part III (of IV) of 'spheroid polarity': MATLAB
%This macro is called by part II, processes images and returns results. 


% stack='BBOpt3_025.tif';
% results='Spheroid_Parameters_Results.txt';
% cell=cell([1 16]);
function [cell]=readstack(stack, results)
%% 1 Initialisation

%Definition of scale
%PixelPerMicron=1; Necessary if scale is not in px.

%Get image infos
info1 = imfinfo (stack);
r = max([info1.Height]);
c = max([info1.Width]);
z = length([info1.FileSize]);
for i = 1:z; 
I(:,:,i) = imread(stack,i);
end 
%Order of the channels, adoptable for other stainings/channels order
gp58=imread(stack,1);
gp135=imread(stack,2);
aktin=imread(stack,3);
nuclei=imread(stack,4);


%% Import results table

delimiter = ',';
startRow = 2;

%% Format string for each line of text:
%   Input order: Image_name X_nuclei Y_nuclei Area_nuclei Width_nuclei Height_nuclei Circularity_nuclei 
%	X_spheroid Y_spheroid Width_spheroid Height_spheroid Area_spheroid Circularity_spheroid MaximumRadius_spheroid 
%	Area_actin X_Actin Y_Actin Count_Actin_Particles Count_Actin_Particles_afterWatershed Circularity_Actin_Particles
%   column1: text (%s) Image_name
%	column2: double (%f) X_nuclei 
%   column3: double (%f) Y_nuclei 
%	column4: double (%f) Area_nuclei 
%   column5: double (%f) Width_nuclei 
%	column6: double (%f) Height_nuclei 
%   column7: double (%f) Circularity_nuclei 
%	column8: double (%f) X_spheroid 
%   column9: double (%f) Y_spheroid 
%	column10: double (%f) Width_spheroid 
%   column11: double (%f) Height_spheroid 
%	column12: double (%f) Area_spheroid 
%   column13: double (%f) Circularity_spheroid 
%	column14: double (%f) MaximumRadius_spheroid 
%	column15: double (%f) Area_actin
%	column16: double (%f) X_Actin 
%	column17: double (%f) Y_Actin
%	column18: double (%f) Count_Actin_Particles
%	column19: double (%f) Count_Actin_Particles_afterWatershed
%	column20: double (%f) Circularity_Actin_Particles
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(results,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Create output variable
CenterCoordinates = table(dataArray{1:end-1}, 'VariableNames', {'Image','X_nuclei','Y_nuclei','Area_nuclei','Width_nuclei','Height_nuclei','Circularity_nuclei','X_spheroid','Y_spheroid','Width_spheroid','Height_spheroid','Area_spheroid','Circularity_spheroid','MaximumRadius_spheroid','Area_actin', 'X_actin', 'Y_actin','Count_Actin_Particles','Count_Actin_Particles_after_WS', 'Circularity_Actin_Particles'});

%% Clear temporary variables
clearvars delimiter startRow formatSpec fileID dataArray ans;

%% Read Center Coordinates for opened stack
[row,~]=size(CenterCoordinates);
%Imagename -4 characters for extraction of image-name without '.tif'
stack_search=stack(1:length(stack)-4); 
row_id=0;
%Searching for position of image in table 'CenterCoordinates'
for(i=1:row) 
  Coor_search=CenterCoordinates.Image{i} ; 
  Coor_search=Coor_search(1:length(Coor_search)-4);
  if(strcmp(Coor_search,stack_search))
         row_id=i;
  end
end

%Save Imagename in cell (1,1)
cell(1,1)={stack_search};

%% Clear temporary variables
clearvars row col Coor_search stack_search;

%%Parameters from results table 'CenterCoordinates' are saved as variables

X_Coor_nu=CenterCoordinates.X_nuclei(row_id);
Y_Coor_nu=CenterCoordinates.Y_nuclei(row_id);
Area_nu=CenterCoordinates.Area_nuclei(row_id);
Width_nu=CenterCoordinates.Width_nuclei(row_id);
Heigth_nu=CenterCoordinates.Height_nuclei(row_id);
Circ_nu=CenterCoordinates.Circularity_nuclei(row_id);

X_Coor_s=CenterCoordinates.X_spheroid(row_id); 
Y_Coor_s=CenterCoordinates.Y_spheroid(row_id);
Area_s=CenterCoordinates.Area_spheroid(row_id);
Width_s=CenterCoordinates.Width_spheroid(row_id);
Heigth_s=CenterCoordinates.Height_spheroid(row_id);
Circ_s=CenterCoordinates.Circularity_spheroid(row_id);

Area_act=CenterCoordinates.Area_actin(row_id);

%Save Circularity values for the spheroid, nuclei and actin in cells (1,8)/(1,9) and (1,14)
cell(1,8)={Circ_s};
cell(1,9)={Circ_nu};
Circularity_Actin=CenterCoordinates.Circularity_Actin_Particles(row_id);
cell(1,14)={Circularity_Actin};

%Save amount of actin area on whole spheroid area in cell (1,11)
actin_percentage=Area_act/Area_s;
cell(1,11)={actin_percentage};

%Get and save actin particle counts before and after watershed segmentation in cells (1,12) and (1,13) 
Actin_Particle_Count=CenterCoordinates.Count_Actin_Particles(row_id);
cell(1,12)={Actin_Particle_Count};
Actin_Particle_Count_after_WS=CenterCoordinates.Count_Actin_Particles_after_WS(row_id);
cell(1,13)={Actin_Particle_Count_after_WS};

%Get Center of mass (CoM) for the actin signal and calculate distances for 
%Centres of mass for the different signals and save in cells (1,15) and (1,16)
X_Coor_a=CenterCoordinates.X_actin(row_id); 
Y_Coor_a=CenterCoordinates.Y_actin(row_id); 
Distance_CoM_nu_sph=round(abs(sqrt(((Y_Coor_nu-Y_Coor_s)^2)+((X_Coor_nu-X_Coor_s)^2))));
cell(1,15)={Distance_CoM_nu_sph};
Distance_CoM_nu_act=round(abs(sqrt(((Y_Coor_nu-Y_Coor_a)^2)+((X_Coor_nu-X_Coor_a)^2))));
cell(1,16)={Distance_CoM_nu_act};

%Get maximum radius of the spheroid, needed for calculation of the relativ radius
MaxRad_s=CenterCoordinates.MaximumRadius_spheroid(row_id); 


%% 2a - Mean Intensity on radius  
%Radius array and Intensity array for each Channel are compiled 
%Sum of Intensities on each radius in each channel

radius_a=zeros(r,c);
m=max(r,c);
actin_int_mean=int64(zeros(m,4));
gp135_int_mean=int64(zeros(m,4));
gp58_int_mean=int64(zeros(m,4));
nuclei_int_mean=int64(zeros(m,4));

%radius: distance to centre of mass of the spheroid to each image pixel coordinate
%loop over each image pixel and calculation of the radius value
for (i=1:r)
    for(j=1:c)
    radius_a(i,j)=round(abs(sqrt(((i-Y_Coor_s)^2)+((j-X_Coor_s)^2))));
    end
end

%%Loop over all pixels of the image, summation of intensities and counting of summed pixels
for (i=1:r)
    for(j=1:c)
    radius=radius_a(i,j)+1;   
    actin_int_mean(radius,1)=int64(radius); %column 1: radius
    gp135_int_mean(radius,1)=int64(radius);
    gp58_int_mean(radius,1)=int64(radius);
    nuclei_int_mean(radius,1)=int64(radius);    
    
    actin_int_mean(radius,2)=actin_int_mean(radius,2)+int64(aktin(i,j)); %column 2: summed intensity on the radius
    gp135_int_mean(radius,2)=gp135_int_mean(radius,2)+int64(gp135(i,j));
    gp58_int_mean(radius,2)=gp58_int_mean(radius,2)+int64(gp58(i,j));
    nuclei_int_mean(radius,2)=nuclei_int_mean(radius,2)+int64(nuclei(i,j));
    
    actin_int_mean(radius,3)=actin_int_mean(radius,3)+1; %column 3: number of pixels on the radius
    gp135_int_mean(radius,3)=gp135_int_mean(radius,3)+1;
    gp58_int_mean(radius,3)=gp58_int_mean(radius,3)+1;
    nuclei_int_mean(radius,3)=nuclei_int_mean(radius,3)+1;
    end
end 

%% 2b - Mean Intensity/px on radius

%The variables 'channel'_int_mean contain the uncalculated summed intensities and pixelcounts
%Via division of the summed intensity (column 2) by the pixelcount on the radius (column 3)
%the mean intensity per pixel is calculated and saved in 'channel'_mean_int
%Division of the radius (px) by maximum radius MaxRad (MaxRad=1)

actin_mean_int=zeros(c,3);
gp135_mean_int=zeros(c,3);
gp58_mean_int=zeros(c,3);
nuclei_mean_int=zeros(c,3);

for (i=1:c)
    actin_mean_int(i,1)=i; %column 1: radius (px)
    gp135_mean_int(i,1)=i;
    gp58_mean_int(i,1)=i;
    nuclei_mean_int(i,1)=i;
    actin_mean_int(i,2)=actin_int_mean(i,2)/actin_int_mean(i,3); %column 2: intensity/px    
    gp135_mean_int(i,2)=gp135_int_mean(i,2)/gp135_int_mean(i,3);
    gp58_mean_int(i,2)=gp58_int_mean(i,2)/gp58_int_mean(i,3);
    nuclei_mean_int(i,2)=nuclei_int_mean(i,2)/nuclei_int_mean(i,3);
    actin_mean_int(i,3)=i/MaxRad_s; %column 3: relative radius
    gp135_mean_int(i,3)=i/MaxRad_s;
    gp58_mean_int(i,3)=i/MaxRad_s;
    nuclei_mean_int(i,3)=i/MaxRad_s;
end




%% 2c - Normalisation mean intensity / whole intensity 

%Division of the mean intensity/px 'channel'_mean_int by the summed intensity per channel
%calculates the normalized intensity. The sum of the intensity is set to 1 for each channel.
%Saved as 'channel'_mean_int_corr_norm (column 1: radius(px), column 2: normalized intensity, column 3: relative radius)

actin_mean_int_corr_norm=actin_mean_int;
gp135_mean_int_corr_norm=gp135_mean_int;
gp58_mean_int_corr_norm=gp58_mean_int;
nuclei_mean_int_corr_norm=nuclei_mean_int;

for (i=1:c)
    actin_mean_int_corr_norm(i,2)=actin_mean_int(i,2)/sum(actin_mean_int(:,2));  %normalisation by dividing through whole Intensity
    gp135_mean_int_corr_norm(i,2)=gp135_mean_int(i,2)/sum(gp135_mean_int(:,2));
    gp58_mean_int_corr_norm(i,2)=gp58_mean_int(i,2)/sum(gp58_mean_int(:,2));
    nuclei_mean_int_corr_norm(i,2)=nuclei_mean_int(i,2)/sum(nuclei_mean_int(:,2));
end

%% Clear temporary variables

clearvars actin_int_mean gp135_int_mean gp58_int_mean nuclei_int_mean i j actin_mean_int gp135_mean_int gp58_mean_int nuclei_mean_int;


%% 3 Cumulative Plotting

%The cumulative intensity is calculated by summation of column 2 
%of the variables'channel'_mean_int_corr_norm and saved in column 2 
%of the variable 'channel'_mean_int_cumul (column 1: radius(px), column 2: cumulative intensity, column 3: relative radius)

actin_mean_int_cumul=actin_mean_int_corr_norm;
gp135_mean_int_cumul=gp135_mean_int_corr_norm;
gp58_mean_int_cumul=gp58_mean_int_corr_norm;
nuclei_mean_int_cumul=nuclei_mean_int_corr_norm;

for (i=1:c)
     actin_mean_int_cumul(i,2)=sum(actin_mean_int_corr_norm(1:i,2));
     gp135_mean_int_cumul(i,2)=sum(gp135_mean_int_corr_norm(1:i,2));
     gp58_mean_int_cumul(i,2)=sum(gp58_mean_int_corr_norm(1:i,2));
     nuclei_mean_int_cumul(i,2)=sum(nuclei_mean_int_corr_norm(1:i,2));
end

%Plot of cumulative intensity (y-axis) vs. relative radius (x-axis)
%Plot is saved as pdf file

figure(1)
P=plot(actin_mean_int_cumul(:,3),actin_mean_int_cumul(:,2),'m')
hold on
P=plot(gp135_mean_int_cumul(:,3),gp135_mean_int_cumul(:,2),'r')
hold on
P=plot(gp58_mean_int_cumul(:,3),gp58_mean_int_cumul(:,2),'g')
hold on
P=plot(nuclei_mean_int_cumul(:,3),nuclei_mean_int_cumul(:,2),'b')

set(gca,'XLim',[0 1.4],'YLim',[0 1.2],'XTick',[0:0.2:1.4],'YTick',[0:0.1:1.2]);
filename=cat(2,stack(1:length(stack)-4),'.pdf');

saveas(gca,filename)

close all;

%% Clear temporary variables

clearvars actin_mean_int_corr_norm gp135_mean_int_corr_norm gp58_mean_int_corr_norm nuclei_mean_int_corr_norm;


%% 4 - Save cumulative signals as txt
%Order: rel. radius, actin, gp135, gp58, nuclei

cumul_allMarkers=zeros(c,5); 

cumul_allMarkers(:,1)=actin_mean_int_cumul(:,3); %relative radius
cumul_allMarkers(:,2)=actin_mean_int_cumul(:,2); %actin
cumul_allMarkers(:,3)=gp135_mean_int_cumul(:,2); %gp135
cumul_allMarkers(:,4)=gp58_mean_int_cumul(:,2);	 %gp58
cumul_allMarkers(:,5)=nuclei_mean_int_cumul(:,2);%nuclei

filename=cat(2,stack(1:length(stack)-4),'.txt');
save(filename,'cumul_allMarkers','-ascii','-tabs');

%% 5 - Calculation radius positions at 60% intensity 

n_actin=int8(0);
n_gp135=int8(0);
n_gp58=int8(0);
n_nuclei=int8(0);

%Determination of the row/radius where intensity is >60% the first time 
for(i=1:c)
    if (n_actin==0)
        if(actin_mean_int_cumul(i,2)>=0.6)
           n_actin=i;
        end
    end
    if (n_gp135==0)
        if(gp135_mean_int_cumul(i,2)>=0.6)
           n_gp135=i;
        end
    end
    if (n_gp58==0)
        if(gp58_mean_int_cumul(i,2)>=0.6)
           n_gp58=i;
        end
    end
    if (n_nuclei==0)
        if(nuclei_mean_int_cumul(i,2)>=0.6)
           n_nuclei=i;
        end
    end  
end

%Linear interpolation of the exact radius at 60%
y=0.6;

x2=actin_mean_int_cumul(n_actin,3);
y2=actin_mean_int_cumul(n_actin,2);
x1=actin_mean_int_cumul(n_actin-1,3);
y1=actin_mean_int_cumul(n_actin-1,2);

r_actin=x1+((x2-x1)*(y-y1))/(y2-y1);

clearvars x2 x1 y2 y1;

x2=gp135_mean_int_cumul(n_gp135,3);
y2=gp135_mean_int_cumul(n_gp135,2);
x1=gp135_mean_int_cumul(n_gp135-1,3);
y1=gp135_mean_int_cumul(n_gp135-1,2);

r_gp135=x1+((x2-x1)*(y-y1))/(y2-y1);

clearvars x2 x1 y2 y1;

x2=gp58_mean_int_cumul(n_gp58,3);
y2=gp58_mean_int_cumul(n_gp58,2);
x1=gp58_mean_int_cumul(n_gp58-1,3);
y1=gp58_mean_int_cumul(n_gp58-1,2);

r_gp58=x1+((x2-x1)*(y-y1))/(y2-y1);

clearvars x2 x1 y2 y1;

x2=nuclei_mean_int_cumul(n_nuclei,3);
y2=nuclei_mean_int_cumul(n_nuclei,2);
x1=nuclei_mean_int_cumul(n_nuclei-1,3);
y1=nuclei_mean_int_cumul(n_nuclei-1,2);

r_nuclei=x1+((x2-x1)*(y-y1))/(y2-y1);

clearvars x2 x1 y2 y1;

%Save radius of 60% intensity for each channel
cell(1,4)={r_actin};
cell(1,5)={r_gp135};
cell(1,6)={r_gp58};
cell(1,7)={r_nuclei};

%% Clear temporary variables

clearvars n_actin n_gp135 n_gp58 n_nuclei;

%% 6 - Slope Nuclei - Heigth of rel nuclei signal at 30% of spheroid radius

x_rad=int8(0);
y_nuclei=int8(0);

%Determination of the nuclei intensity at radius >30% the first time 
for(i=1:c)
    if (x_rad==0)
        if(nuclei_mean_int_cumul(i,3)>=0.3)
           x_rad=i;
        end
    end
end
%Linear Interpolation and calculation of the initial nuclei slope

x=0.3;

x2=nuclei_mean_int_cumul(x_rad,3);
y2=nuclei_mean_int_cumul(x_rad,2);
x1=nuclei_mean_int_cumul(x_rad-1,3);
y1=nuclei_mean_int_cumul(x_rad-1,2);

y_nuclei=y1+((y2-y1)*(x-x1))/(x2-x1);
slope=y_nuclei/0.3;
cell(1,10)={slope};

%% Calculation of the differences between gp58 and gp135/ actin and gp135
% 
diff_Markers=r_gp58-r_gp135;
cell(1,2)={diff_Markers};
diff_actin_gp135=abs(r_actin-r_gp135);
cell(1,3)={diff_actin_gp135};


end