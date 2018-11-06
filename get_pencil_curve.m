function P = get_pencil_curve(f)


% Get the input figure or get current one (creates new one if none exist)
if nargin == 0 || isempty(f)
    f = gcf;
end
figure(f);

% get axes of current figure (creates on if doesn't exist)
a = gca;

% set equal axis
axis equal;
% freeze axis
axis manual;
% set view to XY plane
view(2);

set(gcf,'windowbuttondownfcn',@ondown);
% living variables
P = [];
p = [];

% loop until mouse up or ESC is pressed
done = false;
while(~done)
    drawnow;
end

% We've been also gathering Z coordinate, drop it
P = P(:,1:2);

% Callback for mouse press
    function ondown(src,ev)
        % Tell window that we'll handle drag and up events
        set(gcf,'windowbuttonmotionfcn', @ondrag);
        set(gcf,'windowbuttonupfcn',     @onup);
        append_current_point();
    end

% Callback for mouse drag
    function ondrag(src,ev)
        append_current_point();
    end

% Callback for mouse release
    function onup(src,ev)
        finish();
    end


    function append_current_point()
        % get current mouse position
        cp = get(gca,'currentpoint');
        % append to running list
        P = [P;cp(1,:)];
        if isempty(p)
            % init plot
            hold on;
            p = plot(P(:,1),P(:,2));
            hold off;
        else
            % update plot
            set(p,'Xdata',P(:,1),'Ydata',P(:,2));
        end
    end

    function finish()
        done = true;
        set(gcf,'windowbuttonmotionfcn','');
        set(gcf,'windowbuttonupfcn','');
        set(gcf,'windowbuttondownfcn','');
        set(gcf,'keypressfcn','');
    end

end