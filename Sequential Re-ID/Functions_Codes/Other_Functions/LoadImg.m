function I = LoadImg(Image_address,k)

file = dir(Image_address);
filename = strcat([Image_address,'\'],file(k+2).name);
I = imread(filename);

end