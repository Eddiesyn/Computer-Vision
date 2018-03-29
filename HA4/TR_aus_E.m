function [T1,R1,T2,R2] = TR_aus_E(E)
% In dieser Funktion sollen die moeglichen euklidischen Transformationen
% aus der Essentiellen Matrix extrahiert werden

[U,S,V] = svd(E);
R = [0 -1 0;1 0 0;0 0 1];
if det(U)==-1 || det(V)==-1
    [U,S,V] = svd(-E);
end
R1 = U*R'*V';
Td1 = U*R*S*U';
T1 = [Td1(3,2),Td1(1,3),Td1(2,1)]';
R2 = U*R*V';
Td2 = U*R'*S*U';
T2 = [Td2(3,2),Td2(1,3),Td2(2,1)]';








end