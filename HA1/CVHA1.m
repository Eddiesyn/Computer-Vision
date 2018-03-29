%  Gruppennummer: M23
%  Gruppenmitglieder: Shi Yinan; Liu Siyuan; Niu Xianrui; Dong Chao; Hua Jia

%% Hausaufgabe 1
%  Einlesen und Konvertieren von Bildern sowie Bestimmung von 
%  Merkmalen mittels Harris-Detektor. 

%  Für die letztendliche Abgabe bitte die Kommentare in den folgenden Zeilen
%  enfernen und sicherstellen, dass alle optionalen Parameter über den
%  entsprechenden Funktionsaufruf fun('var',value) modifiziert werden können.


%% Bild laden
Image = imread('szene.jpg');
Gray_image = rgb_to_gray(Image);




%% Harris-Merkmale berechnen
Merkmale = harris_detektor(Gray_image,7,0.05,1.5,25,325,325,1);
