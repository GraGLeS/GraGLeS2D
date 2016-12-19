% ----------------------------
% Install MTEX first
% ----------------------------


% ----------------------------
% Reminder: You have to set the number of grains generated from EBSD
% manualy in your parameter.xml file!!!!
% Check the smGrains object in the end to find the number.
% ----------------------------

%load ebsd data
[ebsd,~,~] = loadEBSD('V19800C1hSubset.ctf',...
    'convertEuler2SpatialReferenceFrame');

%Start MTEX
startup_mtex


% determine the biggest square if it is a rectangular:
% decide which coordinate is biggest
mmax = max(get(ebsd,'y'))
ebsd1 = ebsd(get(ebsd,'x')< mmax)

mmax = max(get(ebsd,'y'))
ebsd2 = ebsd(get(ebsd,'x')< mmax)

% calculate grain objects
ebsd = ebsd1;
grains  = calcGrains(ebsd,'angle',3*degree);

% compute TWIN boundraies and estimate their share
%TwinBoundary = perimeter(grains,CSL(3));
%BoundaryLength = sum(perimeter(grains));
%frac = TwinBoundary/BoundaryLength;

grain=grains;
mymean=mean(area(grain));
mymax=max(area(grain));
mymax/mymean

opts.stepsize = get(ebsd(20),'x')/19;
opts.segAngle = 3;
opts.nFilter = 15;
opts.iter = 5;
opts.sm = 0;

% smooth and correct grain boundaries
[ebsd_cor,grains_cor] =  correctEBSD(ebsd,opts);
smGrains = smooth(grains_cor,10);



%figure
%plot(grains) %,'colorcoding','ipdfHSV','r',zvector,'antipodal');
plotx2east;
figure
plot(smGrains)

%get(ebsd,'quaternion');
% get(grains,'V');
%get maximum coordinates:

%grains2=merge(smGrains,CSL(3));
%grains3=merge(grains2,CSL(9));

file = fopen('TWIP.txt','wt');
vertices=get(smGrains,'V');
[q]=quaternion(get(smGrains,'orientation'));
[phi]=Euler(get(smGrains,'orientation'))

vertices(1,:)=0;
xymax = max(vertices);
xymin = min(vertices);
nn=xymax -xymin
normalized_V =zeros(length(vertices),2);

%normiert beide auf ein square
normalized_V(:,1) = (vertices(:,1)-xymin(1))/(nn(1));
normalized_V(:,2) = (vertices(:,2)-xymin(2))/(nn(2));

% normiert die längere achse auf 1
if(nn(1) > nn(2))
    nominator = nn(1)
else
    nominator = nn(2)
end



normalized_V(:,1) = (vertices(:,1)-xymin(1))/(nominator);
normalized_V(:,2) = (vertices(:,2)-xymin(2))/(nominator);


for i = 1 : numel(smGrains)
    
    %get the orientation and write
    [positions]=get(smGrains(i),'boundaryEdgeorder');
    
    if( iscell(positions{1}))
        i
        %kill entire grainboundaries like wholes
        [data]=normalized_V(positions{1}{1},:);
        
    else
        [data]=normalized_V(positions{1},:);
        
    end
    
    NumberFaces=length(data); 
    fprintf(file, '%d %d %6.6f %6.6f %6.6f %6.6f \n',i, NumberFaces, get(q(i), 'a'), get(q(i), 'b'), get(q(i), 'c'), get(q(i), 'd')); 
    %fprintf(file, '%d %d %6.6f %6.6f %6.6f  \n',i, NumberFaces, phi(i,1), phi(i,2), phi(i,3)); 
    for j = 1 : NumberFaces
       fprintf(file,'%6.6f %6.6f \n',data(j,1), data(j,2));
    end
    fprintf(file,'\n');
end


fclose(file); 
