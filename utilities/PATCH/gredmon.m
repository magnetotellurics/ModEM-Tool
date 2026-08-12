temp = computer;
if (temp(1:4) =='PCWI')
   %  assume notebook ...
%   size of window
width =  750;
height = 550; 
% location of lower left corner
loc = [30,30];
%  bottom/side margin for control buttons (normal coordinates)
marg_bot = .20;
%   side margins
marg = .05; ax = 1 - 2*marg;
%  font size for GUI buttons
guiButtonFontSize = 9;
undo_button = 'alt'
else
%   size of window
width =  1050;
height = 840; 
% location of lower left corner
loc = [10,10];
%  bottom/side margin for control buttons (normal coordinates)
marg_bot = .20;
%   side margins
marg = .025; ax = 1 - 2*marg;
%  font size for GUI buttons
guiButtonFontSize = 12;
undo_button = 'extend';
end
