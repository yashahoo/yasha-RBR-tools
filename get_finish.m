function [finish]=get_finish(ts,y)

% [begin finish]=get_startend(ts,y)
% Funtion to manually select begin and finish times 
% 
% Yasha Hetzel jan 2017



% rename to fit
ts1=ts;
res=y;

% set up so can loop through points
badi=nan; clear p

% plot initial data
fig=figure;
plot(ts1,res,'b')
hold on
% plot(ts1,res,'k.')
ylabel('height (m)','fontsize',14)
%  xlabel('Date (2016)','fontsize',14)
set(gca, 'fontsize',14,'ylim',[min(res) max(res)])
grid on
box on
% title('Data to clean')
% set(gca, 'fontsize',14,'xticklabel',[],'xlim',[min(ts1) max(ts1)],'xtick',XX,'ylim',[min(res)-1 max(res)+.5])
set(gca, 'fontsize',14,'xlim',[min(ts1) max(ts1)],'ylim',[min(res) max(res)])
datetick('x','keeplimits')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Find the start of good data
% title('ZOOM in to START of good data, then press Return')
% zoom on;
% % disp('ZOOM in...Then click any key to continue...')
% 
% pause % you can zoom with your mouse and when your image is okay, you press any key
% datetick('x','yyyy-mm-dd HH:MM','keeplimits')
% doc_datacursormode(gcf);
% dcm_obj = datacursormode(gcf);
% set(dcm_obj,'DisplayStyle','datatip',...
%     'SnapToDataVertex','off','Enable','on')
% title('Click on START of good data, then press Return')
% disp('Click on START of good data, then press Return')
% % Wait while the user does this.
% pause 
% c_info = getCursorInfo(dcm_obj)
% p1= c_info.DataIndex;       
%         hold on;
%         plot(ts1(p1),res(p1),'g*')
% begin=c_info.Position(1);
% zoom out;
% zoom off; % to escape the zoom mode
% datetick('x','yyyy-mm-dd HH:MM','keeplimits')

 hold on;
        plot(ts1(1),res(1),'g*')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the end of good data
title('ZOOM in to END of good data, then press Return')
zoom on;
% disp('ZOOM in...Then click any key to continue...')

pause % you can zoom with your mouse and when your image is okay, you press any key
datetick('x','yyyy-mm-dd HH:MM','keeplimits')
doc_datacursormode(gcf);
dcm_obj = datacursormode(gcf);
set(dcm_obj,'DisplayStyle','datatip',...
    'SnapToDataVertex','off','Enable','on')
title('Click on END of good data, then press Return')
disp('Click on END of good data, then press Return')
% Wait while the user does this.
pause 

c_info = getCursorInfo(dcm_obj)
p2= c_info.DataIndex;       
        hold on;
        plot(ts1(p2),res(p2),'g*')

finish=c_info.Position(1);

zoom out;
zoom off; % to escape the zoom mode


% plot(ts1(1:p1),res(1:p1),'r')
plot(ts1(p2:end),res(p2:end),'r')
% plot(ts1(p1),res(p1),'g*')
plot(ts1(p2),res(p2),'g*')
title('Done selecting finish!')
delete(findall(gcf,'Type','hggroup'))
datetick('x','yyyy-mm-dd HH:MM','keeplimits')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







