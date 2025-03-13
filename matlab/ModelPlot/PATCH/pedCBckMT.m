editTask = get(gcbo,'Tag');

switch editTask
    case 'Read'
        %   get pathname for starting patch file to edit
        [cfile,cpath] = uigetfile('*pat.mat','Input Patch File To Edit');
        cfile = [cpath cfile];
        eval(['load ' cfile ]);
        patch_names = patchStruct.names;
        npatch = patchStruct.nPatches;
        mask = patchStruct.mask;
        mask_nums = patchStruct.mask_nums;
        patch_num = 1;
        set(findobj('Tag','ListBox'),'String',...
            patch_names,'Value',patch_num)
        mask_color = mask_nums(patch_num);
        set(findobj('Tag','MC'),'BackGroundColor',map(mask_color,:));
        initPltMT
        next_num = max(mask_nums);
    case 'Mask'
        mask_on = 1;
        set(findobj('Tag','Prev_Grid'),'Enable','off','Value',0);
        set(findobj('Tag','Mask_Blk'),'Enable','off','Value',0);
        set(findobj('Tag','Unmsk_Blk'),'Enable','off','Value',0);
        set(findobj('Tag','Sub_grid'),'Enable','off','Value',0);
        set(findobj('Tag','Full_Grid'),'Enable','off','Value',0);
        x0 = lims_c(1);
        y0 = lims_c(3);
        xs = (n2_c-n1_c)/(lims_c(2) - lims_c(1));
        ys = (m2_c-m1_c)/(lims_c(4) - lims_c(3));
        findobj('Tag','GridPlotAxis');
        axes(findobj('Tag','GridPlotAxis'));
        hold on
        while mask_on
            waitforbuttonpress
            point = get(gca,'CurrentPoint');
            x = point(1,1);
            y = point(1,2);
            button = get(gcf,'SelectionType');
            %        [x,y,button] = ginput(1)
            if button(1:3) == 'nor'
                ii = n1_c + round((x-x0)*xs);
                jj = m1_c + round((y-y0)*ys);
                if ( mask(ii,jj) == 0 )
                    x1 = (ii-n1_c-.5)/xs + x0;
                    x2 = x1 + 1./xs;
                    y1 = (jj-m1_c-.5)/ys + y0;
                    y2 = y1 + 1./ys;
                    px = [x1,x2,x2,x1,x1];
                    py = [y1,y1,y2,y2,y1];
                    patch(px,py,mask_color,'EdgeColor','none','CDataMapping','Direct')
                    mask(ii,jj) = mask_color;
                end
            elseif button(1:3) == undo_button(1:3)
                ii = n1_c + round((x-x0)*xs);
                jj = m1_c + round((y-y0)*ys);
                if(mask(ii,jj) == mask_color )
                    x1 = (ii-n1_c-.5)/xs + x0;
                    x2 = x1 + 1./xs;
                    y1 = (jj-m1_c-.5)/ys + y0;
                    y2 = y1 + 1./ys;
                    px = [x1,x2,x2,x1,x1];
                    py = [y1,y1,y2,y2,y1];
                    cind=ceil((hz(ii,jj)-hLims(1))*nbath/(hLims(2)-hLims(1)));
                    cind = min(cind,32);
                    patch(px,py,map(cind,:),'EdgeColor','none','CDataMapping','Direct');
                    mask(ii,jj) = 0;
                end
            else
                mask_on = 0;
            end
        end
        set(gcbo,'Value',0);
        set(findobj('Tag','Prev_Grid'),'Enable','on');
        set(findobj('Tag','Mask_Blk'),'Enable','on');
        set(findobj('Tag','Unmsk_Blk'),'Enable','on');
        set(findobj('Tag','Sub_grid'),'Enable','on');
        set(findobj('Tag','Full_Grid'),'Enable','on');
        set(findobj('Tag','Open_BC'),'Enable','on');
        
    case 'Mask_Blk'
        box_ind
        mask(n1:n2,m1:m2) = mask_color*(mask(n1:n2,m1:m2)==0) + ...
            mask(n1:n2,m1:m2);
        pltmaskMT(hz(n1_c:n2_c,m1_c:m2_c),lims_c,gcf,mask(n1_c:n2_c,m1_c:m2_c),hLims)
        hold on
        plot(info{1}.lon,info{1}.lat,'k*','Markersize',10,'linewidth',2)
        set(gcbo,'Value',0);
    case 'Unmsk_Blk'
        box_ind
        mask(n1:n2,m1:m2) = ...
            (mask(n1:n2,m1:m2) ~= mask_color).*mask(n1:n2,m1:m2);
        pltmaskMT(hz(n1_c:n2_c,m1_c:m2_c),lims_c,gcf,mask(n1_c:n2_c,m1_c:m2_c),hLims)
        hold on
        plot(info{1}.lon,info{1}.lat,'k*','Markersize',10,'linewidth',2)
        set(gcbo,'Value',0);
    case 'Full_Grid'
        %    restore to initial full grid
        if(m1_c ~= 1 & m2_c ~= m)
            initPltMT
        end
        set(gcbo,'Value',0);
    case 'Prev_Grid'
        %    restore to previous subgrid
        if (m1_c ~= m1_old) | (n1_c ~= n1_old) | ...
                (m2_c ~= m2_old) | (n2_c ~= n2_old)
            m1_c = m1_old;
            n1_c = n1_old;
            m2_c = m2_old;
            n2_c = n2_old;
            pltmaskMT(hz(n1_c:n2_c,m1_c:m2_c),lims_c,gcf,mask(n1_c:n2_c,m1_c:m2_c),hLims);
            hold on
            plot(info{1}.lon,info{1}.lat','k*','Markersize',10,'linewidth',2)
        end
        set(gcbo,'Value',0);
        hold off
    case 'Save'
        [file,path] = uiputfile('*.mat');
        outFile = [ path file];
        nPatches = length(patch_names);
        patchStruct = struct('nPatches',nPatches,'ll_lims',ll_lims,...
            'mask',mask,'names',{patch_names},'mask_nums',mask_nums)
        eval(['save ' outFile ' patchStruct']);
    case 'Add'
        %  if this looks weird, it's because I don't understand cell arrays ye
        npatch = npatch+1;
        next_num = next_num+1;
        patch_num = npatch;
        ctemp = cell(npatch,1);
        mask_color = next_num;
        mask_nums =  [ mask_nums mask_color ];
        if(npatch == 1)
            ctemp(1) = {current_patch};
        else
            ctemp(1:npatch-1) = patch_names(1:npatch-1);
            ctemp(npatch) = {current_patch};
        end
        patch_names = ctemp;
        set(findobj('Tag','ListBox'),'String',...
            patch_names,'Value',patch_num)
        set(findobj('Tag','MC'),'BackGroundColor',map(mask_color,:));
   
    case 'Del'
        %  delete an entry from the list of patch names
        mask = mask.*(mask ~= mask_nums(patch_num));
        ctemp = cell(npatch-1,1);
        if(patch_num == 1)
            ctemp(1:npatch-1) = patch_names(1:npatch-1);
            mask_nums(1:npatch-1) = mask_nums(1:npatch-1);
            patch_num = 1;
        else
            ctemp(1:patch_num-1) = patch_names(1:patch_num-1);
            ctemp(patch_num:npatch-1) = patch_names(patch_num+1:npatch);
            mask_nums(1:patch_num-1) = mask_nums(1:patch_num-1);
            mask_nums(patch_num:npatch-1) = mask_nums(patch_num+1:npatch);
            patch_num = patch_num - 1;
        end
        npatch = npatch-1;
        mask_nums = mask_nums(1:npatch);
        patch_names = ctemp;
        set(findobj('Tag','ListBox'),'String',...
            patch_names,'Value',patch_num)
        mask_color = mask_nums(patch_num);
        set(findobj('Tag','MC'),'BackGroundColor',map(mask_color,:));
        if(npatch > 0 )
            current_patch = char(patch_names(patch_num));
        end
        set(findobj('Tag','patch_name'),'String',current_patch);
        pltmaskMT(hz(n1_c:n2_c,m1_c:m2_c),lims_c,gcf,...
            mask(n1_c:n2_c,m1_c:m2_c),hLims)
        hold on
        plot(info{1}.lon,info{1}.lat','k*','Markersize',10,'linewidth',2)
    case 'depth'
        patch_names = ctemp;
        set(findobj('Tag','ListBox'),'String',...
            patch_names,'Value',patch_num)
        set(findobj('Tag','MC'),'BackGroundColor',map(mask_color,:));
   
        initPltMT
    otherwise
        ['Not coded']
end
