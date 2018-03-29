function [Gray_image] = rgb_to_gray(Image)
% Diese Funktion soll ein RGB-Bild in ein Graustufenbild umwandeln. Falls
% das Bild bereits in Graustufen vorliegt, soll es direkt zur?ckgegeben werden.
a=size(size(Image),2); % festzustellen, ob das uebergebene Bild ein Grauwertbild ist. 
if a==2                % wenn ja, dieses Bild wird sofort als Resultat zurueckgegeben.
    Gray_image=Image;
else                   % wenn nicht, dieses Bild wird in ein Grauwertbild konvertiert. 
    R=Image(:,:,1);
    G=Image(:,:,2);
    B=Image(:,:,3);
    Gray_image=0.299*R+0.587*G+0.114*B;
end

Gray_image=im2double(Gray_image);

end
