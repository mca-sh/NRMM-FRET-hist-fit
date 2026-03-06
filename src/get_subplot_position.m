function pos = get_subplot_position(varargin)
    if nargin==3
        R = varargin{1};
        C = varargin{2};
        n = varargin{3};
    else
        R = varargin{1}(1);
        C = varargin{1}(2);
        n = varargin{1}(3);
    end
    r = ceil(n/C);
    c = n-(r-1)*C;

    w = 1/C;
    h = 1/R;

    x = w*(c-1);
    y = 1-r*h;

    pos = [x,y,w,h];
end