function x = load_texture(filename,varargin)

if numel(varargin)
    ext = varargin{1};
    x = double(imread([filename,'.',ext]))/255;
else
    x = double(imread(filename))/255;
end