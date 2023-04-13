%% preference setting
% keyboard shortcut for `evaluate selection` (not `evaluate CURRENT selection`!): Ctrl + Shift + E
% command window - wrap lines
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
restoredefaultpath;matlabrc
repo_dir = 'C:\Users\ll357\Documents\'
addpath(genpath(repo_dir))
% addpath(genpath('C:\Users\GlickfeldLab\Documents\test\inter'))
% addpath(genpath('C:\Users\GlickfeldLab\Documents\test\ImagingCode-Glickfeld-Hull'))
% addpath(genpath('C:\Users\GlickfeldLab\Documents\test\BehaviorCode-Glickfeld-Hull'))
% addpath(genpath('C:\Users\GlickfeldLab\Documents\test\Scanbox'))

% svnroot = [repo_dir, 'ImagingCode-Glickfeld-Hull'];
% ijroot = 'C:\Program Files\ImageJ';
% coreInitJavaPath(svnroot,ijroot);
% coreInitMatlabPath(svnroot,ijroot);

cd(repo_dir)
clc

%% saving var
save('data_reg.mat', 'data_reg', '-v7.3') % force save >2GB .mat
save('fit_bootstrap.mat', 'well_fit_cell', '-append') % append var to saved mat

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

%% search field name in struct

tt = fieldnames(input_behav);
index = cellfun(@(x) any(contains(x, 'iti')),tt); sum(index)
id = find(index > 0);
for i = 1 : length(id)
    fprintf(['input_behav.', tt{id(i)}])
    fprintf('\n')
end