% Copyright 2024 by Georgy Moshkin <georgy.moshkin@outlook.com>
% More openEMS examples: https://www.udemy.com/antenna-arrays
%
% Notes:
% 1. openEMS functions are defined in C:\openEMS\matlab (.m files)
% 2. We can use Octave or Matlab
% 3. You can learn from examples C:\openEMS\matlab\examples
%    and C:\openEMS\matlab\Tutorials

% 1. Init "CSX" 3D CAD data structure
CSX = InitCSX();

% 2. Define materials
CSX = AddMetal(CSX, 'my_line');         % Microstrip line metal
CSX = AddMetal(CSX, 'my_ground');       % Ground plane metal
CSX = AddMaterial(CSX, 'my_substrate'); % Substrate material
CSX = SetMaterialProperty(CSX, 'my_substrate', 'Epsilon', 4.8); % Er = 4.8 for substrate

% 3. Define geometry
% AddBox(CSX, name, priority, [x1 y1 z1], [x2 y2 z2])
CSX = AddBox(CSX, 'my_substrate', 0, [-10 -20 -1], [10 20 0]);  % z1=-1, substrate thickness is 1 mm
CSX = AddBox(CSX, 'my_ground', 0, [-10 -20 -1], [10 20 -1]);    % z1=z2, zero thickness ideal ground plane
CSX = AddBox(CSX, 'my_line', 0, [-1.8/2 -20 0], [1.8/2 20 0]);  % z1=z2, because it's zero thickness microstrip line
CSX = AddBox(CSX, 'my_line', 0, [1.8/2 -1.8/2 0], [6.7 1.8/2 0]);  % quarterwave stub

% Field dump for electromagnetic field visualization
CSX = AddDump(CSX,'E_field','FileType',1);
CSX = AddBox(CSX,'E_field',10,[-10 -20 -0.5], [10 20 -0.5]); % -0.5mm is inside substrate


% 4. Add two lumped ports
[CSX port{1}] = AddLumpedPort(CSX, 1, 1, 50, [-1.8/2 -20 -1], [1.8/2 -20 0], [0 0 1], true);
[CSX port{2}] = AddLumpedPort(CSX, 1, 2, 50, [-1.8/2 20 -1], [1.8/2 20 0], [0 0 1], false);

% 5. Slice 3D space using 2D plaes through shape edges added with AddBox() function
mesh = DetectEdges(CSX);

% 6. Append mesh with surrounding empty space boundary coordinates
% appending principle:
% Assume that mesh.x contains array [1 2 3]
% Then mesh.x = [mesh.x 4 5] creates a new array [1 2 3 4 5] with two new elements [4 5]
mesh.x = [mesh.x -25 25]; % two YZ planes at X=-25 and X=25
mesh.y = [mesh.y -25 25]; % two XZ planes at Y=-25 and Y=25
mesh.z = [mesh.z -15 15]; % two XY planes at Z=-15 and Z=15

% 7. Smooth rough mesh
% Basically subdivide the space by more 2D planes in a such way
% that distance between adjacent planes is not larger than "max_res"
% and distance ratio between adjacent planes doesn't exceed "ratio"
mesh = SmoothMesh(mesh, 0.5, 1.25); % (mesh, max_res, ratio)


% 8. Set FDTD parameters
F0 = 5.8*10^9;  % Center frequency 5.8 GHz
FC = 3*10^9;    % Corner frequency +-3 GHz
FDTD = InitFDTD('EndCriteria', 10^-4);
FDTD = SetGaussExcite(FDTD, F0, FC);
FDTD = SetBoundaryCond(FDTD, {'MUR','MUR','MUR','MUR','MUR','MUR'});

% 9. Define rectangular grid
% Basically "mesh" was filled with coordinates in millimeters
% We use deltaUnit=1/1000 to convert this to meters, because
% FDTD calculations use meters.
CSX = DefineRectGrid(CSX, 1/1000, mesh);

% 10. Save data file that can be used by openEMS.exe and AppCSXCAD.exe
mkdir('temp');
WriteOpenEMS('temp/test.xml', FDTD, CSX);

% 11. Display 3D model using CSX CAD
CSXGeomPlot('temp/test.xml');

% 12. Run FDTD simulation
RunOpenEMS('temp', 'test.xml');

% 13. Process and display FDTD output

set(0, "defaulttextfontsize", 24)  % title
set(0, "defaultaxesfontsize", 16)  % axes labels
set(0, "defaultlinelinewidth", 1)

close all % close existing windows if any
freq = linspace(F0-FC, F0+FC, 201);
port = calcPort(port, 'temp', freq);

s11 = port{1}.uf.ref./ port{1}.uf.inc;
s21 = port{2}.uf.ref./ port{1}.uf.inc;

figure
hold on;
plot( freq/1e6, 20*log10(abs(s11)), 'r-');
plot( freq/1e6, 20*log10(abs(s21)), 'b-');
grid on
title({'Reflection Coefficients {\color{red}|S_{11}|} and {\color{blue}|S_{21}|}','Copyright 2024 by Georgy Moshkin'});
xlabel('frequency f / MHz');
ylabel('Magnitude, dB');
ylim([-50 5]);


% draw electromagnetic field distribution
[myField myMesh] = ReadHDF5Dump(['temp' '/E_field.h5']);
myField2 = GetField_TD2FD(myField, F0);
sx=size(myField2.FD.values{1})(1);
sy=size(myField2.FD.values{1})(2);


figure;
hold on;

colormap('jet');
[xx,yy]=meshgrid(myMesh.lines{1},myMesh.lines{2});
cc=zeros(sy,sx);

for ii = 1:sx
  for kk = 1:sy
    fz = myField2.FD.values{1}(ii,kk,1,3);
    amp=abs(fz);
    cc(kk,ii)=sin(arg(fz))*abs(fz);
  endfor
endfor

ss=pcolor(xx,yy,cc);

set(ss,'FaceColor','interp','EdgeColor','none'); % replace 'none' with 'black' to view mesh

SCALE=1/1000; % to meters

% draw layout (metal layers)
metalN=size(CSX.Properties.Metal)(2);
for nn = 1:metalN
  primitives=CSX.Properties.Metal{1,nn}.Primitives.Box;
  primitivesN=size(primitives)(2);
  for tt= 1:primitivesN
    X1=primitives{1,tt}.P1.ATTRIBUTE.X;
    Y1=primitives{1,tt}.P1.ATTRIBUTE.Y;
    X2=primitives{1,tt}.P2.ATTRIBUTE.X;
    Y2=primitives{1,tt}.P2.ATTRIBUTE.Y;
    SX=X2-X1;
    SY=Y2-Y1;
    rectangle("Position", [X1*SCALE, Y1*SCALE, SX*SCALE, SY*SCALE], "EdgeColor", "black", "FaceColor", "none");
  endfor
endfor

title(sprintf("Ez field distribution @ %.2f GHz, \nCopyright 2024 by Georgy Moshkin",F0/1e9));

X1=abs(CSX.Properties.Material{1,1}.Primitives.Box{1,1}.P1.ATTRIBUTE.X);
Y1=abs(CSX.Properties.Material{1,1}.Primitives.Box{1,1}.P1.ATTRIBUTE.Y);
X2=abs(CSX.Properties.Material{1,1}.Primitives.Box{1,1}.P2.ATTRIBUTE.X);
Y2=abs(CSX.Properties.Material{1,1}.Primitives.Box{1,1}.P2.ATTRIBUTE.Y);

% scale
DIM1=X1*SCALE;
DIM2=X2*SCALE;
DIM3=Y1*SCALE;
DIM4=Y2*SCALE;
DIM=max([DIM1,DIM2,DIM3,DIM4])*1.25; % leave 25% empty from the sides

disp(DIM)
axis ([-DIM, DIM, -DIM, DIM], "square");


