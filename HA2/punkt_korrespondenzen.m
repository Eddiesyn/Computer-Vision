function [Korrespondenzen] = punkt_korrespondenzen(I1,I2,Mpt1,Mpt2,varargin)
% In dieser Funktion sollen die extrahierten Merkmalspunkte aus einer
% Stereo-Aufnahme mittels NCC verglichen werden um Korrespondenzpunktpaare
% zu ermitteln.
tic

window_length = varargin{1};
min_corr = varargin{2};
It1 = double(I1);
It2 = double(I2);
n1 = size(Mpt1,2);
n2 = size(Mpt2,2);
N = ones(window_length);
label = zeros(n1,n2);

for i = 1:n1
    test1 = Mpt1(:,i);
    if test1(2)-floor(window_length/2)>=1 && test1(2)+floor(window_length/2)<=2000 && test1(1)-floor(window_length/2)>=1 && test1(1)+floor(window_length/2)<=3000
        W = It1(test1(2)-floor(window_length/2):test1(2)+floor(window_length/2),test1(1)-floor(window_length/2):test1(1)+floor(window_length/2));
    else
        continue;
    end
    Wa = (N*W*N)/(window_length^2);
    sigma1 = sqrt(1/(window_length^2-1)*trace((W-Wa)'*(W-Wa)));
    Wn = (W-Wa)/sigma1;
    for j = 1:n2
        test2 = Mpt2(:,j);
        if test2(2)-floor(window_length/2)>=1 && test2(2)+floor(window_length/2)<=2000 && test2(1)-floor(window_length/2)>=1 && test2(1)+floor(window_length/2)<=3000
            V = It2(test2(2)-floor(window_length/2):test2(2)+floor(window_length/2),test2(1)-floor(window_length/2):test2(1)+floor(window_length/2));
        else
            continue;
        end
        Va = (N*V*N)/(window_length^2);
        sigma2 = sqrt(1/(window_length^2-1)*trace((V-Va)'*(V-Va)));
        Vn = (V-Va)/sigma2;
        Ncc = 1/(window_length^2-1)*trace(Wn'*Vn);
        if Ncc > min_corr
            label(i,j) = Ncc;
        end
    end
    label(i,label(i,:) < max(label(i,:))) = 0;
end

[~,label1] = max(label);
label2 = find(label1 > 1);
label1(label1 == 1)=[];
Inew = [I1,I2];
num = size(label1,2);

Korrespondenzen = zeros(4,num);
for i = 1:num
    Korrespondenzen(1,i) = Mpt1(1,label1(i));
    Korrespondenzen(2,i) = Mpt1(2,label1(i));
    Korrespondenzen(3,i) = Mpt2(1,label2(i));
    Korrespondenzen(4,i) = Mpt2(2,label2(i));
end

imshow(Inew);
hold on
plot(Mpt1(1,label1),Mpt1(2,label1),'*r');
hold on
plot(Mpt2(1,label2)+3000,Mpt2(2,label2),'*g');
hold on
for i = 1:num
    plot([Mpt1(1,label1(i)),Mpt2(1,label2(i))+3000],[Mpt1(2,label1(i)),Mpt2(2,label2(i))])
    hold on
end
hold off

toc
end


