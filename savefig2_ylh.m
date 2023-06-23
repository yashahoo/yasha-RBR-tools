function []= savefig_ylh(save_matfig,print_fig,use_export_fig,img_type,fname) 
% export figs with different requirements using either export_fig or
% inbuilt matlab print.
% [USAGE]
% savefig_ylh(save_matfig,print_fig,use_export_fig,img_type,fname)
% save_matfig = 'N'; % save .fig ? 'Y' or 'N'
% print_fig ='Y';    % export plot to file ? 'Y' or 'N'
% use_export_fig ='Y' % attemts to use export_fig that requires ghostscriptinstalled)
% img_type = {'pdf'}; % format to save figure if print_fig = 'Y' Options
% are pdf, png, eps. Can be multiple, e.g.  {'pdf','png','eps'}. pdf can
% cause problems, so default saves to eps unless change code
% fname = complete filename excluding file extension but inculding path. If no path it will save in working directory.
switch save_matfig
    case 'Y'
        savefig(fname);
end

switch print_fig
    case 'Y'
        try
            for k=1:length(img_type)
                if exist('export_fig','file') & strcmpi(use_export_fig,'Y')
                    if strcmpi(img_type{k},'pdf')
                        export_fig(fname,'-pdf','-painters');
                    elseif  strcmpi(img_type{k},'png')
                        export_fig(fname,'-png','-r300'); %, '-transparent'
                    elseif  strcmpi(img_type{k},'eps')
                        export_fig(fname,'-eps');
                    end
                elseif  ~exist('export_fig','file') | strcmpi(use_export_fig,'N')
                    disp('NOT using export_fig')
                    if strcmpi(img_type{k},'pdf')
                        print(fname,'-depsc','-painters');
                        disp('saving to eps instead of pdf to avoid problems...');
                        %                     print(fname,'-dpdf','-painters','-bestfit');
                        
                    elseif  strcmpi(img_type{k},'eps')
                        print(fname,'-depsc','-painters');
                    elseif  strcmpi(img_type{k},'png')
                        print(fname,'-dpng','-r300')
                    end
                end
            end
            disp([img_type{:} ' saved to ' fname '.' img_type{:}])
        catch
            disp('ERROR saving fig to file...')
        end
    case 'N'
        disp('Not printing fig to file...')
        
end

end