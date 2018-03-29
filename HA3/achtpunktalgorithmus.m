function [EF] = achtpunktalgorithmus(Korrespondenzen,K)
% Diese Funktion berechnet die Essentielle Matrix oder Fundamentalmatrix
% mittels 8-Punkt-Algorithmus, je nachdem, ob die Kalibrierungsmatrix 'K'
% vorliegt oder nicht
P = inputParser;

P.addOptional('K',eye(3),@(x)isnumeric(x) && size(x,1) == 3 && size(x,2) == 3)
P.parse(K);
K = P.Results.K;

if size(Korrespondenzen,2) < 8
    disp('kein genug Korrespondenzpunkten')
% elseif size(Korrespondenzen,2) > 8
%     label = randperm(length(Korrespondenzen),8);
%     I = Korrespondenzen(:,label);
else
    I = Korrespondenzen;
end
num = length(I);
I = [I(1:2,:);ones(1,num);I(3:4,:);ones(1,num)];
for i = 1:num-1
    A = vertcat( kron(I(1:3,i),I(4:6,i))', kron(I(1:3,i+1),I(4:6,i+1))' );
end

[~,~,V] = svd(A);
G = vec2mat(V(:,9),3)';
[U,S,V] = svd(G);
S(3,3) = 0;
EF = K'*(U*S*V')*K;

end