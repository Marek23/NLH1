function [Ps] = getPencils(f)

SSL = 0;

% Get the input figure or get current one (creates new one if none exist)
if nargin == 0 || isempty(f)
    f = gcf;
end
figure(f);


set(gcf,'windowbuttondownfcn',@ondown);
set(gcf,'keypressfcn',        @onkeypress);
% living variables

% loop until mouse up or ESC is pressed
done = false;
while(~done)
    drawnow;
end

% We've been also gathering Z coordinate, drop it

% Callback for mouse press
    function ondown(src,ev)
        SSL = SSL+1;
        Ps{SSL} = get_pencil_curve(f);
        set(gcf,'windowbuttondownfcn',@ondown);
        set(gcf,'keypressfcn',        @onkeypress);
    end



    function onkeypress(src,ev)
        % escape character id
        ESC = char(27);
        switch ev.Character
            case ESC
                finish();
            otherwise
                error(['Unknown key: ' ev.Character]);
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

