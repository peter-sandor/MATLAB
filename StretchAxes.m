function StretchAxes(varargin)
N=length(varargin);
tweak=[0 0 0 0];
for ind1=1:N
    if strcmp(varargin{ind1},'top2up')
      tweak=tweak+[0 0 0 0.01];
    end
    if strcmp(varargin{ind1},'top2down')
      tweak=tweak+[0 0 0 -0.01];
    end
    if strcmp(varargin{ind1},'bottom2down')
      tweak=tweak+[0 -0.01 0 0.01];
    end
    if strcmp(varargin{ind1},'bottom2up')
      tweak=tweak+[0 +0.01 0 -0.01];
    end
    if strcmp(varargin{ind1},'left2left')
      tweak=tweak+[-0.01 0 0.01 0];
    end
    if strcmp(varargin{ind1},'left2right')
      tweak=tweak+[+0.01 0 -0.01 0];
    end
    if strcmp(varargin{ind1},'right2right')
      tweak=tweak+[0 0 0.01 0];
    end
    if strcmp(varargin{ind1},'right2left')
      tweak=tweak+[0 0 -0.01 0];
    end
    if strcmp(varargin{ind1},'up')
      tweak=tweak+[0 0.01 0 0];
    end
    if strcmp(varargin{ind1},'down')
      tweak=tweak+[0 -0.01 0 0];
    end
    if strcmp(varargin{ind1},'left')
      tweak=tweak+[-0.01 0 0 0];
    end
    if strcmp(varargin{ind1},'right')
      tweak=tweak+[0.01 0 0 0];
    end
    if strcmp(varargin{ind1},'shrink')
      tweak=tweak+[0.01 0.01 -0.02 -0.02];
    end
    if strcmp(varargin{ind1},'expand')
      tweak=tweak+[-0.01 -0.01 0.02 0.02];
    end
end
axpos=get(gca,'position');
set(gca,'position',axpos+tweak);
end