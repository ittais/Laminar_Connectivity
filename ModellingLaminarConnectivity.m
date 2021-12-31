
function [ADJ] = ModellingLaminarConnectivity(connectivity,SG2,G2,IG2,granularity,labels)
%This function models cortical laminar connections based on cortical
%connectivity, cortical laminar composition, and cortical granularity
%indices. 
%Model implementation necessitates all of three components per
%cortical region in order to properly compute.
%
%Written by: Ittai Shamir (Mail: ittaisha@mail.tau.ac.il)
%
%Reference: 
%Shamir I, Assaf Y (2021a). An MRI-based, data-driven model of cortical laminar connectivity. Neuroinformatics, vol. 19, 205–218. 
%(doi: 10.1007/s12021-020-09491-7)
%
%
%The defualt atlas used in this example is the BRAINNETOME atlas.
%Reference:
%Fan L, Li H, Zhuo J, Zhang Y, Wang J, Chen L, Yang Z, Chu C, Xie S, 
%Laird AR, Fox PT, Eickhoff SB, Yu C, Jiang T (2016) The Human Brainnetome 
%Atlas: A New Brain Atlas Based on Connectional Architecture. 
%Cerebral Cortex, 26-8, 3508–3526. 
%doi: 10.1093/cercor/bhw157
%
%Brainnetome region order:
%1:123: left hemisphere regions (LH)
    %1:105: cortical regions in left hemisphere  
    %106:123: subcortical regions of left hemisphere
%124:246: right hemisphere regions (RH)
    %124:228: cortical regions in right hemisphere
    %229:246: subcortical regions in right hemisphere
%
%IMPORTANT: FOR OTHER ATLASES, CHANGE INDICES OF START/END OF
%CORTICAL/SUBCORTICAL REGIONS (ACCORDING TO FITTING ORDER OF REGIONS).
%
%Abbreviations:  
%General to brain structure:
    %LH- left hemisphere
    %RH- right hemisphere
    %C- cortical
    %SC- subcortex  
    %S- connection strength
%
%Specific to cortical laminar structure:
    %SG- supragranular layers 
    %G- granular layer
    %IG- infragranular layers
%
%Structure of all bihemispheric variables: 1st half: left hemisphere, 
%2nd half: right hemisphere***
%
%Input arguments: 
    %1. connectivity- cortical connectivity matrix (LxL, L-number of regions)
    %2. SG2- supragranular regional content (Lx1)
    %3. G2- granular regional content (Lx1)
    %4. IG2- infragrnular regional content (Lx1)
        %(sum of SG2,G2, and IG2 per cortical region should equal 1) 
    %5. granularity- regional granularity indices, (Lx1) 
    %where: 
        %6- granular, 
        %5- increasing granular presence, 
        %4- increasing granular presence, 
        %3- increasing granular presence, 
        %2- slightly granular, 
        %1- agranular, 
        %0- irrelevant (SC) 
    %6. labels- region labels (Lx1)

%
%Output arguments:
    %ADJ22- adjacency matrix representing laminar connectivity (3Lx3L)
%Adjacency matrix structure:
%       IG:     G:     SG:
% IG: (IG,IG) (IG,G) (IG,SG)
% G:  (G,IG)   (G,G)  (G,SG)
% SG: (SG,IG) (SG,G) (SG,SG)
%

%Initialize matrices:
len=length(connectivity); %length of matrix 
ADJ3=zeros(len*3,len*3); %adjacency matrix

%Individual components for adjacency matrix:
IG2IG2=zeros(len,len); %infragranular to infragranular connections
IG2G2=zeros(len,len); %infragranular to granular connections
IG2SG2=zeros(len,len); %infragranular to supragranular connections

G2IG2=zeros(len,len); %granular to infragranular connections
G2G2=zeros(len,len); %granular to granular connections
G2SG2=zeros(len,len); %granular to supragranular connections

SG2IG2=zeros(len,len); %supragranular to infragranular connections
SG2G2=zeros(len,len); %supragranular to granular connections
SG2SG2=zeros(len,len); %supragranular to supragranular connections

%LAMINAR CONNECTIVITY:
%
%Indices of start/end of cortical/subcortical regions (change according to atlas): 
LHendC=105; %end of cortical regions LH
LHstartSC=106; %start of subcortical regions LH
LHendSC=123; %end of subcortical regions LH

RHstartC=124; %start of cortical regions RH
RHendC=228; %end of cortical regions RH
RHstartSC=229; %start of subcortical regions RH

%TRACTOGRAPHY-BASED (LONG-RANGE) LAMINAR CONNECTIVITY:
for i=1:size(connectivity,1)
    for j=1:size(connectivity,2)
        %HORIZONTAL CONNECTIONS:
        if i>=j %to phalf, matrix symmetric
        if connectivity(i,j)~=0 && i~=j %if tract i-j exists
            %intrahemipshere (single hemi)
            if (i<=LHendC&& j<=LHendC) || (i>=RHstartC && i<=RHendC) %LH / RH
               S=connectivity(i,j);
                %
                %RULE A:
                if granularity(i)==1 && granularity(j)==6
                    IG2SG2(i,j)=S; 
                %
                %RULE B (symmetric to A):
                elseif granularity(i)==6 && granularity(j)==1
                    SG2IG2(i,j)=S; 
                %
                %RULE C:
                elseif granularity(i)<granularity(j) %i- low, j-high
                    
                    IG2IG2(i,j)=S*(1/3)*(2/3)*IG2(j)/(IG2(j)+G2(j));
                    IG2G2(i,j)=S*(1/3)+(2/3)*G2(j)/(IG2(j)+G2(j));
                    IG2SG2(i,j)=S*(2/3)^2;
                    SG2IG2(i,j)=S*((1/3)^2)*IG2(j)/(IG2(j)+G2(j));
                    SG2G2(i,j)=S*((1/3)^2)*G2(j)/(IG2(j)+G2(j));
                    SG2SG2(i,j)=S*(1/3)*(2/3);
                %    
                %RULE D (symmetric to C):
                elseif granularity(i)>granularity(j) %i-high, j-low
                                       
                    IG2IG2(j,i)=S*(1/3)*(2/3)*IG2(i)/(IG2(i)+G2(i));
                    IG2G2(j,i)=S*(1/3)+(2/3)*G2(i)/(IG2(i)+G2(i));
                    IG2SG2(j,i)=S*(2/3)^2;
                    SG2IG2(j,i)=S*((1/3)^2)*IG2(i)/(IG2(i)+G2(i));
                    SG2G2(j,i)=S*((1/3)^2)*G2(i)/(IG2(i)+G2(i));
                    SG2SG2(j,i)=S*(1/3)*(2/3);
                    
                %
                %RULE E:
                elseif granularity(i)==granularity(j) && granularity(i)~=0 && granularity(j)~=0
                    if granularity(i)==5 || granularity(i)==6
                        SG2SG2(i,j)=S; %S*SG2(i)/(SG2(i)+SG2(j));
                    elseif granularity(i)==3 || granularity(i)==4
                        G2G2(i,j)=S;
                    elseif granularity(i)==1 || granularity(i)==2
                        IG2IG2(i,j)=S; 
                    end
                end
            end  
            %Interhemispheric connections (between two hemispheres):
            if (i<=LHendC && (j>=RHstartC && j<=RHendC)) || ((i>=RHstartC && i<=RHendC) && j<=LHendC) %LH+RH / RH+LH
                SG2SG2(i,j)=S; %S*SG2(i)/(SG2(i)+SG2(j));
            end
        %VERTICAL CONNECTIONS (C-SC):

            if (j>=LHstartSC && j<=LHendSC) || (i>=RHstartSC && j<=246) %CORTEX-SC
               G2G2(i,j)=S; %?????????????? SC2G
            end         
         end
        end
    end
end
%
%ASSUMED LAMINAR (SHORT-RANGE) CONNECTIVITY (NOT TRACTOGRAPHY-BASED):
for i=1:size(connectivity,1)
    for j=1:size(connectivity,2)
        if i==j %intrinsic connection (within a single region)
           S2=1; %intrinsic connection weight (assumed)
           G2SG2(i,j)=S2; %
           SG2G2(i,j)=S2; %S2*G2(i);
           SG2SG2(i,j)=S2;%*SG2(i);
           SG2IG2(i,j)=S2; %S2*SG2(i);
           IG2SG2(i,j)=S2; %S2*SG2(i);
        end
    end
end

%Replace any NaNs:
IG2IG2(isnan(IG2IG2))=0;
IG2G2(isnan(IG2G2))=0;
IG2SG2(isnan(IG2SG2))=0; 

G2IG2(isnan(G2IG2))=0;
G2G2(isnan(G2G2))=0;
G2SG2(isnan(G2SG2))=0;

SG2IG2(isnan(SG2IG2))=0; 
SG2IG2(isnan(SG2G2))=0; 
SG2IG2(isnan(SG2SG2))=0; 

%Rebuild adjacency matrix:
ADJ3(1:len,1:len)=IG2IG2;
ADJ3((len+1):len*2,1:len)=G2IG2;
ADJ3((len*2+1):end,1:len)=SG2IG2;

ADJ3(1:len,(len+1):len*2)=IG2G2;
ADJ3((len+1):len*2,(len+1):len*2)=G2G2;
ADJ3((len*2+1):end,(len+1):len*2)=SG2G2;

ADJ3(1:len,(len*2+1):end)=IG2SG2;
ADJ3((len+1):len*2,(len*2+1):end)=G2SG2;
ADJ3((len*2+1):end,(len*2+1):end)=SG2SG2;

ADJ3(isnan(ADJ3))=0;

ADJ=(ADJ3+ADJ3')/2; %make matrix symmetric
 

%Plot resulting adjacency matrix:
myLabel=cat(1,labels,labels,labels); 
labs={'IG','G','SG'}; %laminar components

figure; 
imagesc(log(ADJ)); 
set(gca, 'XTick', (len/2):len:(len/2)+len*2, 'XTickLabel', labs,'Fontsize',10);
set(gca, 'YTick', (len/2):len:(len/2)+len*2, 'YTickLabel', labs,'Fontsize',10); grid on;
set(gca, 'XTick',1:length(ADJ),'XTickLabel', myLabel,'Fontsize',2); xtickangle(90);
set(gca, 'YTick',1:length(ADJ),'YTickLabel', myLabel,'Fontsize',2); 

hold on; line([len,len], [3*len,0], 'Color', 'w','LineStyle','--');
hold on; line([2*len,2*len], [3*len,0], 'Color', 'w','LineStyle','--');
hold on; line([3*len,0], [len,len], 'Color', 'w','LineStyle','--');
hold on; line([3*len,0], [2*len,2*len], 'Color', 'w','LineStyle','--');
set(gca,'xaxisLocation','top');
axis square;

colorbar(gca,'FontSize',8);
set(gcf,'Color','w'); %set background to white
caxis([0 max(max(log(ADJ)))]); %set colorbar min and max

grid on; 
colormap(hot); %colormap
%Additional colormap option: colorbrewer
%For colorbrewer uncomment the following line and replace the line 
%'colormap(hot)':
%colormap(flipud(cell2mat(colorbrewer.seq.Blues(9))./255))%seq.YlOrRd(9))./255))

title('Log adjacency matrix','Fontsize',14); %title

end

