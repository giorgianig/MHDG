% function readyforprint([W H],fntsize,fgc,bgc,lwidth)
%
% this script resizes and reformats images for inclusion in documents,
% presentations or web pages.
%
% What this function modifies;
% - paper size (size of output document when output usin 'print')
% - font sizes
% - colours of axes, legend boxes and colour bar boxes
% - widths of lines in the plot, including the axes, legend and colour bar
% boxes and
%
% What this function does not modify;
% - marker sizes
%
% Inputs are;
% [W H];    width and height of the image (inches)
% fntsize;	Font size (pts). All text will be set to this size.
% fgc;      foreground colour, either as 'r' or [0.5 0.5 0.5] style. This
%			colour is for lines in the foreground (data)
% bgc;      background colour, either as 'r' or [0.5 0.5 0.5] style. This
%			colour is for axes, legend boxes, etc.
% lwidth;	Width of lines in figures. use 0.5 for journal plots, 2 for
%			overhead presentations, etc.
%
% example;
% readyforprint([6 3],10,'k','w',0.5)
%
% Any or all inputs can be left empty, in which case a default set is used.
% This is equivalent to
% > readyforprint([6 3],10,'k','w',0.5)
% These defaults are approximately suitable for a full-width journal image.
%
% Suggested image sizes are;
% [7 4];        for a picture spanning a whole A4 or letter page
% [3.33 3.33];  for a picture half the width of a page (subfigure)
%
% Note that typically used paper sizes are;
% [8.5 11];     letter paper
% [8.27 11.69]; A4 paper
% [5.83 8.27];  A5 paper
%
% and in fact the paper size should be fitted to the text width, not the
% paper width. This is the user's concern, and may be specified by a
% Journal or University, etc, rather than the user having any say in the
% matter... A consistent use of the same widths or heights should help any
% document look much smarter.
%
% follow this command with something like;
% print('-depsc2','-loose','GLR_fig2a')     % for a colour eps
% print('-djpeg','GRL_fig2a')               % for a colour jpg
%
% This version written by Andy, January 2010 for submission to MatlabCentral

function readyforprintnew(wh,fs,fgcin,bgcin,lw,fh,ms,name)



  
%% input check

% first deal with empty ones
if isempty(wh)
    wh = [6 3];
elseif(numel(wh)) == 1
    disp('')
    disp('Please set the required width and height (in inches)!')
    disp('')
    error('Width and height not defined in inputs to size_for_pub!')
    return
end
if isempty(fs)
    fs = 10;
end
if isempty(fgcin)
    fgc = 'k';
else
    fgc = fgcin;
end
if isempty(bgcin)
    bgc = 'w';
else
    bgc = bgcin;
end
if isempty(ms)
    ms=10;
end
if exist('fh','var')
    if isempty(fh)
        % Run on current figure
        fh = gcf;
    end
else
    fh = gcf;
end
set(0,'CurrentFigure',fh)

%specify some figure properties
defaults.figure_props.InvertHardcopy = 'off';

% specify the paper properties
defaults.paper_props.Color = bgc;
defaults.paper_props.PaperUnits ='inches';
defaults.paper_props.PaperOrientation ='Portrait';
defaults.paper_props.PaperPosition = [0 0 wh(1) wh(2)];
defaults.paper_props.PaperPositionMode = 'manual';
defaults.paper_props.PaperSize = [wh(1) wh(2)];

% specify the axis properties
defaults.axes_props.FontSize = max(floor(0.8*fs),6);
defaults.axes_props.Units = 'normalized';
defaults.axes_props.Color = 'none';
defaults.axes_props.Xcolor = fgc;					% foreground colours
defaults.axes_props.Ycolor = fgc;
defaults.axes_props.Zcolor = fgc;

% specify the line properties
%defaults.line_props.Color = fgc;

% specify line width
if ~isempty(lw)
	defaults.axes_props.LineWidth = lw*0.5;
	defaults.line_props.LineWidth = lw;
end

% specify the text properties
% defaults.text_props.Color = fgc;					% text colour
defaults.text_props.FontSize = fs;

% specify the legend and colorbar properties
defaults.legend_props.FontSize = max(floor(0.8*fs),6);
defaults.legend_props.Xcolor = fgc;
defaults.legend_props.Ycolor = fgc;
defaults.legend_props.TextColor = fgc;
defaults.legend_props.Color = bgc;

defaults.colorbar_props.FontSize = max(floor(0.8*fs),6);
defaults.colorbar_props.Xcolor = fgc;
defaults.colorbar_props.Ycolor = fgc;
defaults.colorbar_props.Color = bgc;

% specify the title properties
defaults.title_props.Color = fgc;
defaults.title_props.FontSize = max(ceil(1.2*fs),6);
defaults.title_props.FontName = 'Times New Roman';
defaults.title_props.FontWeight = 'normal';

% specify the label properties
defaults.label_props.Color =  fgc;
defaults.label_props.FontSize = max(fs,6);
defaults.label_props.BackgroundColor = 'none';
defaults.label_props.FontName = 'Times New Roman';
defaults.label_props.FontWeight = 'normal';


%%%%%%%%%%%%%%% END OF INPUTS %%%%%%%%%%%%%%%

% now need to concatenate these inputs into a useful format; prefered
% option for setting properties is a cell array of strings for the
% properties, and a cell array of properties

%%%%%%%%%%%%%%% START MODIFICATION %%%%%%%%%%%%%%%

%% PAPER PROPERTIES %%%%%%%%%%%%%%%
set(fh,fieldnames(defaults.figure_props)',struct2cell(defaults.figure_props)')
set(fh,fieldnames(defaults.paper_props)',struct2cell(defaults.paper_props)')

%% ORDER OF AXES
myaxesraw = findobj(fh,'Type','Axes');
% figure out what these axes are, and sort so that legends and colorbars
% are dealt with at the end
mytypes = get(findobj(fh,'Type','Axes'),'Tag');

naxes = numel(myaxesraw);
legends = strcmp('legend', mytypes);
Colorbars = strcmp('Colorbar', mytypes);
nothers = sum(legends)+sum(Colorbars);

myaxes = zeros(naxes,1);
myaxes(end-nothers+1:end) = vertcat(myaxesraw(legends),...
    myaxesraw(Colorbars));

myaxes(1:end-nothers) = myaxesraw(~(Colorbars | legends));
    


  
%% CORRECT AXES
for myaxe = 1:length(myaxes)
    % set the focus onto the current axis
    Ah = myaxes(myaxe);
    if ishandle(Ah)
        
        % copy default axes properties into a new structure
        axes_props = defaults.axes_props;
        % work out what controls the position of the axes
        axes_props.ActivePositionProperty = get(Ah,'ActivePositionProperty');
        
        
        % get the current foreground and background colours
        old_fgc = get(Ah,'Xcolor');
        old_bgc = get(Ah,'Color');        
        
        %%%%%%%%%%%%%%% LINE PROPERTIES %%%%%%%%%%%%%%%
        
        % find all line items on this axis
        mylines = findobj(Ah,'Type','line');
        % update these so they correspond to the values passed in
		if isfield(defaults,'line_props')
			set(mylines,fieldnames(defaults.line_props)',struct2cell(defaults.line_props)')
		end
        set(findobj(mylines,'MarkerFaceColor','w','-and',...
            'MarkerEdgeColor','k'),...
            'MarkerEdgeColor',fgc,...
            'MarkerFaceColor',bgc,...
            'MarkerSize',ms)
        
        %%%%%%%%%%%%%%% TEXT PROPERTIES %%%%%%%%%%%%%%%
        % Find text items on this axis; note that we use the axes handle
        % as this explicitly ignores the legend. The legend is set and
        % controlled seperately
        mytexts = findobj(Ah,'Type','text','-and','BackgroundColor','none');
        set(mytexts,fieldnames(defaults.text_props)',struct2cell(defaults.text_props)');
        mytexts = findobj(Ah,'Type','text','-and','BackgroundColor',old_bgc);
        set(mytexts,fieldnames(defaults.text_props)',struct2cell(defaults.text_props)');
        set(mytexts,'BackgroundColor',bgc);                
        
        
        %%%%%%%%%%%%%%% FONTS %%%%%%%%%%%%%%%%
        % change the font on numbers and labels to something sans-serif
%         set(findobj(fh,'-regexp','FontName','[^'']'),'FontName','Times New Roman','Fontweight','Bold');
        set(findobj(fh,'-regexp','FontName','[^'']'),'FontName','Times New Roman');
        
        %%%%%%%%%%%%%%% TITLE PROPERTIES %%%%%%%%%%%%%%
        try
            if ~isempty(get(get(Ah,'Title'),'String'))
                set(get(Ah,'Title'),...
                    fieldnames(defaults.title_props)',struct2cell(defaults.title_props)')
            end
        end
        
        %%%%%%%%%%%%%% LABELS %%%%%%%%%%%%%%%%%%
        try
            if ~isempty(get(get(Ah,'XLabel'),'String'))
                set(get(Ah,'XLabel'),...
                    fieldnames(defaults.label_props)',struct2cell(defaults.label_props)')
            end
        end
        try
            if ~isempty(get(get(Ah,'YLabel'),'String'))
                set(get(Ah,'YLabel'),...
                    fieldnames(defaults.label_props)',struct2cell(defaults.label_props)')
            end
        end
        try
            if ~isempty(get(get(Ah,'ZLabel'),'String'))
                set(get(Ah,'ZLabel'),...
                    fieldnames(defaults.label_props)',struct2cell(defaults.label_props)')
            end
		end
		%%%%%%%%%%%%%%%% CHECK AXES PROPERTIES %%%%%%%%%%%%%%%
		% get the current positions
        %axes_props.Position = get(Ah,'Position');
        %axes_props.Position(axes_props.Position <= 0) = 0.001;
        %axes_props.Position(axes_props.Position > 1) = 1;
        %axes_props.OuterPosition = get(Ah,'OuterPosition');
        %axes_props.OuterPosition(axes_props.OuterPosition <= 0) = 0.001;
        %axes_props.OuterPosition(axes_props.OuterPosition >1 ) = 1;
        
        %%%%%%%%%%%%%%% UPDATE AXES PROPERTIES %%%%%%%%%%%%%%%
        switch lower(get(Ah,'Tag'))
            case 'legend'
                legend_props = defaults.legend_props;
                legend_props.Location = get(Ah,'Location');
                set(Ah,fieldnames(legend_props)',struct2cell(legend_props)');
            case 'colorbar'   
                colorbar_props = defaults.colorbar_props;
                set(Ah,fieldnames(colorbar_props)',struct2cell(colorbar_props)');
            otherwise
                set(Ah,fieldnames(axes_props)',struct2cell(axes_props)');                
        end
    end
end



  
%% PLAY NICELY WITH BOX PLOTS
if ischar(bgc)
    if strcmp(bgc,'none')
        boxfacecolor = 'w';
    else
        boxfacecolor = bgc;
    end
else
    boxfacecolor = bgc;
end
set(findobj(fh,'Tag','Box'),...
    'LineWidth',lw,...
    'Color',fgc,...
    'MarkerFaceColor',boxfacecolor)

% set(findobj(fh,'Type','line'),'Markersize',ms);
% set(findobj(fh,'Type','text'),'Fontsize',fs-2);
% set(findobj(fh,'Type','axis'),'Fontname','times new roman');
print('-dpdf','-r300',['/home/giorgio/Dropbox/PostDoc_Marseille/Latex/Paper_HDG/figures/figure_2107/' name '.pdf'])
% print('-deps','-r300',['/home/giorgio/Dropbox/PostDoc_Marseille/Latex/ReportHDGMHD/' name '.eps']),system(['epstopdf /home/giorgio/Dropbox/PostDoc_Marseille/Latex/ReportHDGMHD/' name '.eps'])

system(['pdfcrop /home/giorgio/Dropbox/PostDoc_Marseille/Latex/Paper_HDG/figures/figure_2107/'  name '.pdf /home/giorgio/Dropbox/PostDoc_Marseille/Latex/Paper_HDG/figures/figure_2107/' name '.pdf'])
% system(['pdfcrop /home/giorgio/Dropbox/PostDoc_Marseille/Latex/Report/'  name '.pdf  --margins "-20 -10 -35 0 " /home/giorgio/Dropbox/PostDoc_Marseille/Latex/Report/' name '.pdf'])

%% code for testing; uncomment before use!
% figure
% plot(rand(10),rand(10),'ko','MarkerFaceColor','r')
% legend('test')
% set(legend,'location','best')
% saveas(gcf,'this_is_my_raw_figure','fig')
% readyforprint([4 3],10,[1 1 1],[0 0.1 0.8],3)
% print('-djpeg','this_is_my_printed_jpg')  % for a jpg
% print('-depsc2','-loose','this_is_my_printed_eps') % for an eps