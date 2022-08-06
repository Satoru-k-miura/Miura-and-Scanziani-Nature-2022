function varargout = runlen(varargin)

if length(varargin)==1
    data=varargin{1};
    breaks = [true logical(diff(data))];
    runStart = find([breaks true]);
    elem = data(breaks);
    len = diff(runStart);
    varargout = {elem, len};
elseif length(varargin)==2
    elem = varargin{1};
    len = varargin{2};
    elem=elem(:).';
    len=len(:);
    lc=cumsum(len);
    lx=zeros(1,lc(end));
    lx([1;lc(1:end-1)+1])=1;
    lc=cumsum(lx);
    varargout{1}=elem(lc);
end