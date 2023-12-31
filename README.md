# What is Hexelate?

The name "Hexelate" is a wordplay combining the words "hexagon" and "pixelate". This is a small set of code (currently only in Matlab) that takes a digital image, samples its intensities in a hexagonal grid, and then draws colored dots of various sizes at the grid points to create a "downsampled(pixelated)" representation of the image. Alternatively, an image drawn with hexagons of different intensity/color can be created. 


# Why did I wrote Hexelate and who may want to use it?

I wrote this code in 2022 to aid create a physical art work for my wife's birthday. I used a family portrait as the input and stamped dots of various sizes on a canvas using items like pens, a permanent marker and a coin. Thus the section of "instruction" which color-codes the different dot size.

Hexelate is useful for anybody who wants to generate a similar kind of hexelated images or artwork with their own picture.

# Instructions

## Installation
This code is written as a Matlab script as a single .m file. All custom functions are included as local functions at the bottom of file. To run this you will need Matlab installed on your computer. Simply download the file "HexImageConversions.m" to a directory of your choice and open it in Matlab. 

If you are not familiar with the idea of a script you can read more [here](https://www.mathworks.com/help/matlab/learn_matlab/scripts.html)

## Generate a dot image with canvas dimensions
Open the file "HexImageConversions.m" in Matlab. Run the cells one by one by the "Run and Advance." button in the editor, or by pressing Ctrl+Enter (Windows) or Replace the path to the example image with the path to your own image. Like this:

``` Matlab
ImgMat = imread('\path\to\your\own\image.JPG'); hexWidthPx = 30;
```

The variable ```hexWidthPx``` determines how big the hexagons will be (in pixels of the original image). This is the distance between the center of neighboring hexagons, or equivalently the distance between opposite sides of the hexagons.

In the third cell calculate output size, you can specify the size of the canvas this image is supposed to occupy by adjusting ```canvasWidth``` and ```canvasHeight```. 

``` Matlab
%% calculate output size
% specify canvas size
canvasWidth = 50; % cm
canvasHeight = 50; % cm
```


# Credit

The part of the code that deals with hexagonal coordinate systems was adapted from codes created by [Red Blob Games](https://www.redblobgames.com/grids/hexagons), who released their codes under a CC0 license. 
The artistic style of dot image was inspired by two artists as of 2022 residing in [Ateliers Ackersdijkstraat 20](http://ackersdijkstraat20.nl/), Rotterdam, Netherlands: [Michael Bom](http://michaelbom.com/) and Boris Pas. 