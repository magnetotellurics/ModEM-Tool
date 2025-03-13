val = get(gcbo,'Value');
str = get(gcbo,'String');
current_patch = char(str(val));
patch_num = val;
mask_color = mask_nums(patch_num);
set(findobj('Tag','patch_name'),'String',current_patch);

set(findobj('Tag','MC'),'BackGroundColor',map(mask_color,:));
