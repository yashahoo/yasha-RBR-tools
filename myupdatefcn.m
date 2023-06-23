function txt = myupdatefcn(empt,event_obj)
% Customizes text of data tips

pos = get(event_obj,'Position');
txt = {['Time: ',datestr(pos(1),'yyyy-mm-dd HH:MM:SS')],...
	      ['Amplitude: ',num2str(pos(2))]};