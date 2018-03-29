function [Fx,Fy] = sobel_xy(Image)
% In dieser Funktion soll das Sobel-Filter implementiert werden, welches
% ein Graustufenbild einliest und den Bildgradienten in x- sowie in
% y-Richtung zurueckgibt.
Sabelx=[1 0 -1;2 0 -2;1 0 -1]; % Definieren wir Sabel-Operator in x-Richtung.
Sabely=[1 2 1;0 0 0;-1 -2 -1]; % Definieren wir Sabel-Operator in y-Richtung.
% Sabelx=[5 0 -5;8 0 -8;5 0 -5];
% Sabely=[5 0 -5;8 0 -8;5 0 -5];
Fx=conv2(Image,Sabelx,'same'); % Berechnen wir die erste Ableitung in x-Richtung mit Sabelfilter-Operator.
Fy=conv2(Image,Sabely,'same'); % Berechnen wir die erste Ableitung in y-Richtung mit Sabelfilter-Operator.
% figure(1)
% imshow(Fx);
% figure(2)
% imshow(Fy);
% figure(3)
% imshow(sqrt(Fx.^2+Fy.^2));
end

