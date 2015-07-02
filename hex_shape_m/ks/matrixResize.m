function result = matrixResize(image, new_size)
if size(image, 1) > size(image, 2)
   my_size = size(image, 2);
else
   my_size = size(image, 1);
end
image = image(1:my_size, 1:my_size);
coef = new_size / my_size;
result = imresize(image, coef);
