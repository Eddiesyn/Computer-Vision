function  Merkmale = harris_detektor(Image,varargin) 
% In dieser Funktion soll der Harris-Detektor implementiert werden, der
% Merkmalspunkte aus dem Bild extrahiert
tic
Gray_image=rgb_to_gray(Image);
% Gray_image=imresize(Gray_image,0.3);

%{ 
Definieren wir die Groesse des Bildsegments, Variable k, sigma von Gaussian 
Operator, minimaler Abstand, Breite und Hoehe der Kachelgroesse, maximale Anzahl an Merkmalen innerhalb einer Kachel
%} 
%%
segment_length=varargin{1}; 
k=varargin{2};              
sigma=varargin{3};
min_dist=varargin{4};
tile_sizex=varargin{5};
tile_sizey=varargin{6};
N=varargin{7};
% tau=varargin{8};
if N > tile_sizex*tile_sizey
    fprintf('Wrong input parameter!!!\n\n\n'); % festzustellen, ob die Eingabevariable N rational ist.
end
%%
[Fx,Fy] = sobel_xy(Gray_image);
F=sqrt(Fx.^2+Fy.^2);
figure(1);
subplot(2,2,1);
imshow(Gray_image);
subplot(2,2,2);
imshow(Fx);
subplot(2,2,3);
imshow(Fy);
subplot(2,2,4);
imshow(F);
% Angabe 3,1 faengt an.
% [Fx,Fy] = sobel_xy(F);
Fx2=Fx.^2;
Fy2=Fy.^2;
Fxy=Fx.*Fy;
h=fspecial('gaussian',[segment_length,segment_length],sigma); % Setze ich den Gaussian Fenster Groesse und sigma
FX2=conv2(Fx2,h,'same'); 
FY2=conv2(Fy2,h,'same');
FXY=conv2(Fxy,h,'same');

% R=FX2.*FY2-FXY.^2-k*(FX2+FY2).^2; % die Matrix wird aus dem Grauwertbild Image Harris-Merkmale ueber das Kriterium extrahiert.
% R=conv2(R,h,'same');
[ii,jj]=size(Gray_image);
R=zeros(ii,jj);
sita=3; % Setze ich den Schwellenwert fuer das Filtern bei Singulaerwertzerlegung
for i=1:ii
    for j=1:jj
        M=[FX2(i,j) FXY(i,j);FXY(i,j) FY2(i,j)];
        [~,S,V]=svd(M);
        T=eig(S);
        if abs(T(1)/T(2))<sita && abs(T(2)/T(1))<sita % Wenn damit zufrieden ist, werden sie nicht falsch detektierende Ecken.
            R(i,j)=det(M)-k*(trace(M))^2;
        else
            continue;
        end
    end
end
Rm=max(max(R)); % finden wir den groessten Wert in Matrix R
Merkmale=zeros(ii,jj);
for i=1:ii
    for j=1:jj
        if R(i,j)>0.01*Rm
            Merkmale(i,j)=R(i,j); % detektieren wir die Ecken ( benutzen wir 0.01*Rm als Schwellenwert 0.01*Rm )
        end
    end
end

[a,b]=find(Merkmale~=0); % plotten wir das Resultat
num1=length(a)
figure(2);
imshow(Gray_image);
hold on;
plot(b,a,'go');
% Angabe 3,1 beendet
%%
% Angabe 3,2 faengt an.
clear a;
clear b;
for i=1:ii-min_dist+1
    for j=1:jj-min_dist+1
        if sum(sum(Merkmale(i:i+min_dist-1,j:j+min_dist-1)))~=0 % wenn es hier Ecken gibt
            test1=Merkmale(i:i+min_dist-1,j:j+min_dist-1);
            Test1=zeros(min_dist);
            P= test1==max(max(test1));   % finden wir den groessten Wert in jedem kleinen Fenster
            Test1(P)=test1(P);           % halten den groessenten Wert und beseitigen den anderen
            Merkmale(i:i+min_dist-1,j:j+min_dist-1)=Test1;
        else
            continue;
        end
    end
end
[e,f]=find(Merkmale~=0);
num2=length(e)
figure(3);
imshow(Gray_image);
hold on;
plot(f,e,'go');
% Angabe 3,2 beendet.
%%
% Angabe 3,3 faengt an.
tile_sizexx=rem(ii,tile_sizex); % Berechnen wir den moeglicherweise existierenden Rest in Zeilen 
tile_sizeyy=rem(jj,tile_sizey); % Berechnen wir den moeglicherweise existierenden Rest in Spalten
% Durch zwei Kreislaeufe filtern wir das ganze Bild, damit in jeder Kachel
% nur N maximale Merkmale existieren. 
for i=1:tile_sizex:ii-tile_sizexx 
    for j=1:tile_sizey:jj-tile_sizeyy
        test2=Merkmale(i:i+tile_sizex-1,j:j+tile_sizey-1);
        if max(size(find(test2>0)))>N
            [~,bbb]=sort(test2(:)); % anorden die Punkte nach ihrer Groesse
            test2(bbb(1:end-N))=0; % beseitigen die andere Punkte und nur die maximale N Punkte festhalten
            Merkmale(i:i+tile_sizex-1,j:j+tile_sizey-1)=test2;
        else
            continue;
        end
    end
end
clear test1;
clear test2;
clear bbb;
% Die entsprechende Behandlung von Kaechelgroesse, die nicht ganzzahlig in
% das Bild passen. 
if tile_sizexx~=0 || tile_sizeyy~=0
    test1=Merkmale(ii-tile_sizexx+1:ii,1:jj);
    test2=Merkmale(1:ii,jj-tile_sizeyy+1:jj);
    if max(size(find(test1>0)))>N
        [~,bbb]=sort(test1(:));
        test1(bbb(1:end-N))=0;
        Merkmale(ii-tile_sizexx+1:ii,1:jj)=test1;
    end
    clear bbb;
    if max(size(find(test2>0)))>N
        [~,bbb]=sort(test2(:));
        test2(bbb(1:end-N))=0;
        Merkmale(1:ii,jj-tile_sizeyy+1:jj)=test2;
    end
end
[x,y]=find(Merkmale~=0);
num3=length(x)
figure(4);
imshow(Gray_image);
hold on;
for i=1:num3
    plot(y(i),x(i),'go');
end
% Angabe 3,3 beendet.
%%
toc
end