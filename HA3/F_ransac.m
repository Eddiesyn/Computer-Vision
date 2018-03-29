function [Korrespondenzen_robust] = F_ransac(Korrespondenzen,varargin)
% Diese Funktion implementiert den RANSAC-Algorithmus zur Bestimmung von
% robusten Korrespondenzpunktpaaren

P = inputParser;

P.addOptional('p',0.95,@(x)isnumeric(x) && x > 0 && x < 1);
P.addOptional('epsilon',0.5,@(x)isnumeric(x) && x > 0 && x < 1);
P.addOptional('tolerance',10,@isnumeric);

P.parse(varargin{:});

p              = P.Results.p;
epsilon        = P.Results.epsilon;
tolerance      = P.Results.tolerance;

s = (log(1-p))/log(1-(1-epsilon)^8);
n1 = length(Korrespondenzen);
I = [Korrespondenzen(1:2,:);ones(1,n1);Korrespondenzen(3:4,:);ones(1,n1)];

e = [0 -1 0;1 0 0;0 0 0];
ranpick = @(x) x(:,randperm(length(x),8));
sampson_dis = @(I,F,e) (diag(I(4:6,:)'*F*I(1:3,:)).^2)'./(sum((e*F*I(1:3,:)).^2,1)+sum((I(4:6,:)'*F*e).^2,2)');

Itest = ranpick(Korrespondenzen);
F = achtpunktalgorithmus(Itest,eye(3));
d1 = sampson_dis(I,F,e);
m1 = sum(d1 < tolerance);
Korrespondenzen_robust = I(:,d1 < tolerance);
for i = 1:round(s)
    Itest = ranpick(Korrespondenzen);
    F = achtpunktalgorithmus(Itest,eye(3));
    d2 = sampson_dis(I,F,e);
    m2 = sum(d2 < tolerance);
    if m2 > m1
        m1 = m2;
        Korrespondenzen_robust = I(:,d2 < tolerance);
        d1 = d2;
    elseif m2 == m1
           d1s = sum(d1,2);
           d2s = sum(d2,2);
           if d2s < d1s
               m1 = m2;
               Korrespondenzen_robust = I(:,d2 < tolerance);
               d1 = d2;
           end
    end
end
Korrespondenzen_robust = [Korrespondenzen_robust(1:2,:);Korrespondenzen_robust(4:5,:)];

hold on
figure(1);
% plot(Korrespondenzen(1,:),Korrespondenzen(2,:),'.g');
% plot(Korrespondenzen(3,:),Korrespondenzen(4,:),'.r');
plot(Korrespondenzen_robust(1,:),Korrespondenzen_robust(2,:),'^m');
plot(Korrespondenzen_robust(3,:),Korrespondenzen_robust(4,:),'^b');
hold off
% If = ranpick(Korrespondenzen_robust);
% F = achtpunktalgorithmus(If);




end