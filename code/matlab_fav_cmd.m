%% runline windows
currentEditor = matlab.desktop.editor.getActive;
originalSelection = currentEditor.Selection;
assert(originalSelection(1)==originalSelection(3));%Check that multiple lines are not selected
currentEditor.Selection = [originalSelection(1) 1 originalSelection(1) Inf];%Select the whole line
disp(currentEditor.SelectedText)
eval(currentEditor.SelectedText);%Run the whole line
% currentEditor.Selection = originalSelection + 1;%Reset selection to original state + 1 line (go to next line)
currentEditor.Selection = [originalSelection(1) + 1, 1, originalSelection(1) + 1, 1];%Reset selection to original state + 1 line (go to next line)
clear currentEditor originalselection
disp(' ')

%% startup 
addpath(genpath('D:\Lan_temp'))
disp('path added')

currentdir = 'C:\Users\lan\Documents\repos';
% currentdir = 'D:\Lan_temp';
svnroot = [currentdir, '\ImagingCode-Glickfeld-Hull'];
ijroot = 'C:\ProgramFiles\ImageJ';
coreInitJavaPath(svnroot,ijroot);
coreInitMatlabPath(svnroot,ijroot);
disp('imageJ added')
disp(' ')