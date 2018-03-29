function [T,R, lambdas, P1] = rekonstruktion(T1,T2,R1,R2, Korrespondenzen, K)
% Funktion zur Bestimmung der korrekten euklidischen Transformation, der
% Tiefeninformation und der 3D Punkte der Merkmalspunkte in Bild 1

R3 = R1;
R4 = R2;
T3 = T2;
T4 = T1;
num = length(Korrespondenzen);
M1 = zeros(3*num,num+1);
M2 = zeros(3*num,num+1);
M3 = zeros(3*num,num+1);
M4 = zeros(3*num,num+1);
M5 = zeros(3*num,num+1);
M6 = zeros(3*num,num+1);
M7 = zeros(3*num,num+1);
M8 = zeros(3*num,num+1);
x1 = K\[Korrespondenzen(1:2,:);ones(1,num)];
x2 = K\[Korrespondenzen(3:4,:);ones(1,num)];
hat = @(x) [0,-x(3),x(2);x(3),0,-x(1);-x(2),x(1),0];
X1 = zeros(3*num,3);
X2 = zeros(3*num,3);
E_Bewegung = cell(2,4);
for i = 1:num
    X1(3*i-2:3*i,:) = hat(x1(:,i));
    X2(3*i-2:3*i,:) = hat(x2(:,i));
end

Aa1 = X2*R1*x1;
Ab1 = X2*T1;
Aa2 = X2*R2*x1;
Ab2 = X2*T2;
Aa3 = X2*R3*x1;
Ab3 = X2*T3;
Aa4 = X2*R4*x1;
Ab4 = X2*T4;
Ba1 = X1*(R1\x2);
Bb1 = X1*R1*T1;
Ba2 = X1*(R2\x2);
Bb2 = X1*R2*T2;
Ba3 = X1*(R3\x2);
Bb3 = X1*R3*T3;
Ba4 = X1*(R4\x2);
Bb4 = X1*R4*T4;

for i = 1:num
    M1(3*i-2:3*i,i) = Aa1(3*i-2:3*i,i);
    M2(3*i-2:3*i,i) = Aa2(3*i-2:3*i,i);
    M3(3*i-2:3*i,i) = Aa3(3*i-2:3*i,i);
    M4(3*i-2:3*i,i) = Aa4(3*i-2:3*i,i);
    M5(3*i-2:3*i,i) = Ba1(3*i-2:3*i,i);
    M6(3*i-2:3*i,i) = Ba2(3*i-2:3*i,i);
    M7(3*i-2:3*i,i) = Ba3(3*i-2:3*i,i);
    M8(3*i-2:3*i,i) = Ba4(3*i-2:3*i,i);
end

M1(:,num+1) = Ab1; % lamdas1_R1,T1
M2(:,num+1) = Ab2; % lamdas1_R2,T2
M3(:,num+1) = Ab3; % lamdas1_R3,T3
M4(:,num+1) = Ab4; % lamdas1_R4,T4
M5(:,num+1) = Bb1; % lamdas2_R1,T1
M6(:,num+1) = Bb2; % lamdas2_R2,T2
M7(:,num+1) = Bb3;% lamdas2_R3,T3
M8(:,num+1) = Bb4; %lamdas2_R4,T4

E_Bewegung{1,1} = R1;
E_Bewegung{1,2} = R2;
E_Bewegung{1,3} = R3;
E_Bewegung{1,4} = R4;
E_Bewegung{2,1} = T1;
E_Bewegung{2,2} = T2;
E_Bewegung{2,3} = T3;
E_Bewegung{2,4} = T4;

[~,~,V1] = svd(M1);
V1 = V1(:,end);
[~,~,V2] = svd(M2);
V2 = V2(:,end);
[~,~,V3] = svd(M3);
V3 = V3(:,end);
[~,~,V4] = svd(M4);
V4 = V4(:,end);
V = [V1,V2,V3,V4];
% V = [V1,V2];

[~,~,L1] = svd(M5);
L1 = L1(:,end);
[~,~,L2] = svd(M6);
L2 = L2(:,end);
[~,~,L3] = svd(M7);
L3 = L3(:,end);
[~,~,L4] = svd(M8);
L4 = L4(:,end);
L = [L1,L2,L3,L4];
% L = [L1,L2];
WT1 = V(num+1,:);
WT2 = L(num+1,:);
% WNummer1 = find((WT1 > 0) == 1);
if isempty(WT1 < 0) == 0
    WNummer = find((WT1 < 0) == 1);
    V(:,WNummer) = -V(:,WNummer);
end

if isempty(WT2 < 0) == 0
    WNummer = find((WT2 < 0) == 1);
    L(:,WNummer) = -L(:,WNummer);
end

Vergleich = sum( [V;L] > 0,1 );
[~,Nummer] = max(Vergleich);
lambdas = zeros(num+1,2);
lambdas(:,1) = V(:,Nummer);
lambdas(:,2) = L(:,Nummer);
R = E_Bewegung{1,Nummer};
T = lambdas(num+1,1)*E_Bewegung{2,Nummer};

P1 = zeros(3,num);
for i = 1:num
    P1(:,i) = lambdas(i,1)*x1(:,i);
end
P2 = bsxfun(@plus,R*P1,T);
p0 = [0 0 0];
pp0 = R*p0'+T;
px = [0.05 0 0];
py = [0 0.05 0];
pz = [0 0 0.05];
p0 = [p0;p0;p0];
pxyz = [px;py;pz];
ppx = R*px'+T-pp0;
ppy = R*py'+T-pp0;
ppz = R*pz'+T-pp0;
pp0 = [pp0';pp0';pp0'];
ppxyz = [ppx';ppy';ppz'];

% hold on;
% figure('name','P1 und P2');
% line([p0(1),px(1)],[p0(2),px(2)],[p0(3),px(3)]);
% hold on;
% line([p0(1),py(1)],[p0(2),py(2)],[p0(3),py(3)]);
% hold on;
% line([p0(1),pz(1)],[p0(2),pz(2)],[p0(3),pz(3)]);
% hold on;
% scatter3(P1(1,:),P1(2,:),P1(3,:),'r*');
% hold on;
% line([pp0(1),ppx(1)],[pp0(2),ppx(2)],[pp0(3),ppx(3)]);
% hold on;
% line([pp0(1),ppy(1)],[pp0(2),ppy(2)],[pp0(3),ppy(3)]);
% hold on;
% line([pp0(1),ppz(1)],[pp0(2),ppz(2)],[pp0(3),ppz(3)]);
% hold off;
figure('name','P1 mit beiden Kameras');
quiver3(p0(:,1),p0(:,2),p0(:,3),pxyz(:,1),pxyz(:,2),pxyz(:,3));
hold on;
quiver3(pp0(:,1),pp0(:,2),pp0(:,3),ppxyz(:,1),ppxyz(:,2),ppxyz(:,3));
hold on;
scatter3(P1(1,:),P1(2,:),P1(3,:),'r*');
hold on;
scatter3(P2(1,:),P2(2,:),P2(3,:),'g*');
hold off;






end