%% figure tips
figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf, 'Position', get(0, 'Screensize'));
cd C:\Users\lan\Documents\repos\inter\code
saveas(gcf, ['ori_tuning_fit_', num2str(icell)], 'jpg')

line([0-5, 180+5], [0, 0], 'Color', 'g', 'LineWidth', 1);
line([1, 1], [0, 300], 'Color', 'g', 'LineWidth',1, 'LineStyle','--');
yl = ylim; % query [ymin ymax]
line([ori_pref, ori_pref], [yl(1), (b_hat + R1_hat)], 'Color', 'r', 'LineWidth', 1);

errorbar(ori_pref_binned_list, osi_shift_avg(:, igap), osi_shift_ste(:, igap),...
        'color', color_list{igap} , 'LineStyle','none')
legend('isi 750', 'isi 250', 'Location','southeast')

for itext = 1 : length(dis_list)
    text(dis_list(itext), ...
       dfof_dis_norm_avg(itext, 1) + dfof_dis_norm_ste(itext, 1) + 0.02, ...
        ['n=', num2str(ntrial_dis(itext))], 'HorizontalAlignment', 'center')
end

%% preference setting
% evaluate selection: Ctrl + Shift + E
% command window wrap lines
% general - confirmation dialogs - confirm before matlab exit 

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

%% openvar
currentEditor = matlab.desktop.editor.getActive;
originalSelection = currentEditor.Selection;
disp(currentEditor.SelectedText)
openvar(currentEditor.SelectedText);
% currentEditor.Selection = [originalSelection(1) + 1, 1, originalSelection(1) + 1, 1];%Reset selection to original state + 1 line (go to next line)
clear currentEditor originalselection
disp(' ')

%% edit function for detailed help
currentEditor = matlab.desktop.editor.getActive;
originalSelection = currentEditor.Selection;

prefix = 'edit ';
function_name = currentEditor.SelectedText;
disp([prefix function_name])

call_help = [prefix function_name];
eval(call_help);

%% startup 
addpath(genpath('D:\Lan_temp'))
addpath(genpath('C:\Users\lan\Documents\repos'))
disp('path added')

currentdir = 'C:\Users\lan\Documents\repos';
% currentdir = 'D:\Lan_temp';
svnroot = [currentdir, '\ImagingCode-Glickfeld-Hull'];
ijroot = 'C:\ProgramFiles\ImageJ';
coreInitJavaPath(svnroot,ijroot);
coreInitMatlabPath(svnroot,ijroot);
disp('imageJ added')

cd C:\Users\lan\Documents\repos\inter\code\
disp('directory set home')
disp(' ')
clc

%% saving var
save('data_reg.mat', 'data_reg', '-v7.3') % force save >2GB .mat

%% search field name in struct

tt = fieldnames(input);
index = cellfun(@(x) any(contains(x, 'Contr')),tt); sum(index)
id = find(index > 0);
for i = 1 : length(id)
    fprintf(['input.', tt{id(i)}])
    fprintf('\n')
end
