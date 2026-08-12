action = get(gcbo,'Tag');
%  Handle for main figure which generated this subwindow
hWindow = get(gcbo,'Parent');
h0 = get(hWindow,'UserData');
P = get(h0,'UserData');
switch action
case 'No Interp'
   P.nointerp = get(gcbo,'Value');
   set(h0,'UserData',P);
case 'Flip Color Scale'
   P.flipud = get(gcbo,'Value');
   P.cmap = flipud(P.cmap);
   set(h0,'UserData',P);
case 'X limits 1'
   P.i1 = str2num(get(gcbo,'String'));
   set(h0,'UserData',P);
case 'X limits 2'
   P.i2 = str2num(get(gcbo,'String'))+1;
   set(h0,'UserData',P);
case 'Y limits 1'
   P.j1 = str2num(get(gcbo,'String'));
   set(h0,'UserData',P);
case 'Y limits 2'
   P.j2 = str2num(get(gcbo,'String'))+1;
   set(h0,'UserData',P);
case 'Z limits 1'
   P.k1 = str2num(get(gcbo,'String'));
   set(h0,'UserData',P);
case 'Z limits 2'
   P.k2 = str2num(get(gcbo,'String'))+1;
   set(h0,'UserData',P);
case 'Sigma limits 1'
   P.Clims(1) = str2num(get(gcbo,'String'));
   set(h0,'UserData',P);
case 'Sigma limits 2'
   P.Clims(2) = str2num(get(gcbo,'String'));
   set(h0,'UserData',P);
case 'Plot Lat Lon'
   P.PlotLatLon = get(gcbo,'Value');
   set(h0,'UserData',P);
case 'Lat'
   P.lat0 = str2num(get(gcbo,'String'));
   set(h0,'UserData',P);
case 'Lon'
   P.lon0 = str2num(get(gcbo,'String'));
   set(h0,'UserData',P);
case 'Limits OK'
   CondReplot
   close(hWindow)
case 'Cancel'
   close(hWindow)
end
