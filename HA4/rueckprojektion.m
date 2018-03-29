function repro_error = rueckprojektion(Korrespondenzen, P1, I2, T, R, K)
% Diese Funktion berechnet die projizierten Punkte in Kamera 2 und den
% mittleren Rueckprojektionsfehler

P2 = bsxfun(@plus,R*P1,T);
X2 = zeros(3,length(P2));

P2 = [P2;ones(1,length(P2))];
PI = [1 0 0 0;0 1 0 0;0 0 1 0];
P2 = K*PI*P2;
% for i = 1:3
%     X2(i,:) = P2(i,:)./lambdas(1:end-1,2)';
% end
for i = 1:length(P2)
    X2(:,i) = P2(:,i)/P2(3,i);
end
X2(3,:) = [];
x2 = Korrespondenzen(3:4,:);
repro_error = sum(sqrt(sum((X2-x2).^2,1)))/length(P2);
disp(repro_error)

figure('name','rueckprojektion');
imshow(uint8(I2));
hold on
plot(X2(1,:),X2(2,:),'r*')
hold on
plot(x2(1,:),x2(2,:),'g*')

for i = 1:length(P2)
    l1 = [X2(1,i),x2(1,i)];
    l2 = [X2(2,i),x2(2,i)];
    line(l1,l2);
end
hold off


end