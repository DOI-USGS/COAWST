function crop(filename,append)
%CROP gets rid of whitespace around an image
%   CROP('filename.ext') crops the image in the file and saves it using
%   the original filename (overwrites the old image). The extension (ext)
%   can be anything IMREAD supports.
%
%   CROP('filename.ext',1) saves the cropped image as filenamecropped.ext
%   in the same directory as the original.
%
%   Changes since version 1:
%       1. Now accepts directories as input.
%       2. Just prior to saving version 1, I added the 'filenamecropped' option.
%
%   Example:
%       crop('C:\MATLAB7\toolbox\matlab\demos\html\cruller_01.png',1)
%
%   See also: IMREAD, IMWRITE

%   Requirements: your FIND must allow the 'last' option (version 7+?)
%   Copyright Andy Bliss Sept 8th, 2006

%set the margin width in pixels
margin=15;

%Get image names
%   if the input is a directory, get all the image files from the directory
if isdir(filename)
    currentdir=pwd;
    cd (filename)
    files=[dir('*.png'); dir('*.gif'); dir('*.bmp'); dir('*.jpg')];
else %if it is a single file:
    files.name=filename;
end

%loop over all the files
for n=1:length(files)
    filename=files(n).name;

    %get the image
    T=imread(filename);

    %sum the RGB values of the image
    xsum=sum(sum(T,3));
    ysum=sum(sum(T,3),2);
    %xsum will be equal to max(xsum) wherever there is a blank column in 
    %   the image (rgb white is [255,255,255]). The left edge for the 
    %   cropped image is found by looking for the first column in which 
    %   xsum is less than max(xsum) and then subtracting the margin
    xleftedge=492;
%     xleftedge=find(xsum<max(xsum),1,'first')-margin;
%     if xleftedge<1
%         xleftedge=1;
%     end
    xrightedge=1602;
%     xrightedge=find(xsum<max(xsum),1,'last')+margin;
%     if xrightedge>length(xsum)
%         xrightedge=length(xsum);
%     end
ytopedge=128;
%     ytopedge=find(ysum<max(ysum),1,'first')-margin;
%     if ytopedge<1
%         ytopedge=1;
%     end
ybottomedge=961;
%     ybottomedge=find(ysum<max(ysum),1,'last')+margin;
%     if ybottomedge>length(ysum)
%         ybottomedge=length(ysum);
%     end

    %resave the image
    if nargin==2 && append
        imwrite(T(ytopedge:ybottomedge,xleftedge:xrightedge,:),[filename(1:end-4) 'cropped' filename(end-3:end)])
    else
        imwrite(T(ytopedge:ybottomedge,xleftedge:xrightedge,:),filename)
    end
end
if exist('currentdir','var')
    cd(currentdir)
end
