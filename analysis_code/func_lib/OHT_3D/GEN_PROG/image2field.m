function [Kfield] = image2field(image,blackvalue, whitevalue)

%image2field: Reads a 24-bit bitmap image file and converts the results to
%a  field with minimum and maximum specified by the user. Dimensions of the
%field will match the dimensions of the original image. This program
%combines all RGB color channels, effectively resulting in the brightness
%(black to white) value being used for mapping each pixel
%
% © 2012-2013 Michael Cardiff, Warren Barrash, and Peter K. Kitanidis, All
% Rights Reserved.
%
%   [Kfield] = image2field(image,blackvalue, whitevalue)
%   where:
%       -image is a text string with the name of the image (must be in the
%       current directory)
%       -blackvalue is the value to be assigned to the cells which are
%       black
%       -whitevalue is the value to be assigned to the cells that are
%       white. All values between whitevalue and blackvalue are
%       interpolated linearly.

numcolors = 256;

%Load.
Kfield = imread(image);

%Values are loaded as integers. This converts everything to normal real
%numbers
   
%If three channels or more, assume first three channels are RGB
if size(Kfield,3) >= 3 %color case
    Kfield = str2num(int2str(sum(Kfield(:,:,(1:3)),3)))./3;
    warning('More than 3 channels, assumed first 3 represent RGB');
elseif size(Kfield,3) == 1 %B&W case
    Kfield = str2num(int2str(Kfield));
else
    error('No protocol for dealing with only 2 channels in image')
end
    
Kfield = Kfield./(numcolors-1).*(whitevalue - blackvalue) + blackvalue;
Kfield = flipdim(Kfield,1);