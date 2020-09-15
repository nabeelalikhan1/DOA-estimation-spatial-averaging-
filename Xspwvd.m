%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 Boualem Boashash
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
% Authors: Prof. Boualem Boashash        (boualem.boashash@gmail.com)
%          Dr. Abdeldjalil Aissa-El-Bey  (abdeldjalil.aissaelbey@telecom-bretagne.eu)
%          RA: Md.F.A
%
% The following references should be cited whenever this script is used:
% [1] B. Boashash, A. Aissa-El-Bey, Multisensor Time-Frequency Signal Processing:
%     A tutorial review with illustrations in selected application areas, Digital
%     Signal Processing, In Press.
% [2] B. Boashash, A. Aissa-El-Bey, M. F. Al-Sa'd, Multisensor time-frequency
%     signal processing software Matlab package: An analysis tool for multichannel
%     non-stationary data , SoftwareX, In Press.
%
% Last Modification: 21-04-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Cross Smoothed Pseudo Wigner-Ville Distribution (XSPWVD)
%
% Syntax : TFR = Xspwvd(s1, s2, varargin);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <INPUTs>
% s1, s2   : Input signals. These can be real or analytic.
% varargin : This is a dynamic input argument that has to be supplied in
%            the following sequence:
%            varargin{1} : The Lag window name. See help of the Matlab function "window" to
%                          list all available windows (default : 'hamming').
%            varargin{2} : Length of the Lag window. This has to be odd (default : N-1 or N).
%            varargin{3} : Smoothing parameter for the following Lag windows:
%                          1) 'gausswin' : Reciprocal of the standard deviation (default : 1).
%                          2) 'kaiser'   : Beta parameter default : 1).
%                          3) 'tukeywin' : specifies the ratio of the length of taper section to the
%                                          total length of the window (default : 1).
%         *** NOTE THAT IF OTHER LAG WINDOWS ARE USED, TAKE THIS VARIABLE OUT OF THE SEQUENCE ***
%
%            varargin{4} : The Time window name. See help of the Matlab function "window" to
%                          list all available windows (default : 'hamming').
%            varargin{5} : Length of the Time window. This has to be odd (default : N-1 or N).
%            varargin{6} : Smoothing parameter for the following Time windows:
%                          1) 'gausswin' : Reciprocal of the standard deviation (default : 1).
%                          2) 'kaiser'   : Beta parameter default : 1).
%                          3) 'tukeywin' : specifies the ratio of the length of taper section to the
%                                          total length of the window (default : 1).
%         *** NOTE THAT IF OTHER TIME WINDOWS ARE USED, TAKE THIS VARIABLE OUT OF THE SEQUENCE ***
%
%            varargin{7} : Number of Frequency bins. This has to be larger than
%                          the Lag window length (default : M = 2^nextpow2(N)).
% <OUTPUTs>
% TFR  : Reduced Interference Cross-TFD using Cross-Smoothed Pseudo Wigner-Ville Distribution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <EXAMPLE>
% s1 = chirp(0:255,0.4,255,0.1);
% s2 = chirp(0:255,0.1,255,0.4);
% s = s1 + s2;
% TFD1 = Xspwvd(s, s, 'hann', 51, 'hann', 21, 512);
% TFD2 = Xspwvd(s, s, 'gausswin', 51, 1, 'hann', 21, 512);
% TFD3 = Xspwvd(s, s, 'hann', 51, 'gausswin', 21, 1, 512);
% TFD4 = Xspwvd(s, s, 'gausswin', 51, 1, 'gausswin', 21, 1, 512);
% figure; imagesc(0:1/1023:1/2,0:255,abs(TFD1')); axis xy
% figure; imagesc(0:1/1023:1/2,0:255,abs(TFD2')); axis xy
% figure; imagesc(0:1/1023:1/2,0:255,abs(TFD3')); axis xy
% figure; imagesc(0:1/1023:1/2,0:255,abs(TFD4')); axis xy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TFR = Xspwvd(s1, s2, varargin)
TFR = [];
ww  = {'bartlett','barthannwin','blackman','blackmanharris','bohmanwin','chebwin','flattopwin','gausswin',...
    'hamming','hann','kaiser','nuttallwin','parzenwin','rectwin','taylorwin','tukeywin','triang'}';
ww1 = {'bartlett','barthannwin','blackman','blackmanharris','bohmanwin','chebwin','flattopwin','hamming',...
    'hann','nuttallwin','parzenwin','rectwin','taylorwin','triang'}';
ww2 = {'gausswin','kaiser','tukeywin'};

%% Main Inputs Checkup
if(nargin == 0), error_msg(1); return;
elseif(nargin == 1), error_msg(2); return;
elseif(isempty(s1) || isempty(s2)), error_msg(2); return;
elseif(~isa(s1,'double') || ~isa(s2,'double')), error_msg(3); return;
end
if(iscolumn(s1)), s1 = s1'; elseif(iscolumn(s2)), s2 = s2'; end
[row1, col1] = size(s1);
[row2, col2] = size(s2);
if(row1 > 1 || row2 > 1), error_msg(4); return;
elseif(col1 ~= col2), error_msg(5); return; end

%% Auxiliary Inputs Checkup
N = length(s1);
if(nargin==2)
    win1 = 'hamming'; win2 = 'hamming';
    if(~mod(N,2)), L1 = N-1; L2 = N-1; else L1 = N; L2 = N; end; M = 2^nextpow2(N);
elseif(nargin==3)
    if(isempty(varargin{1})), win1 = 'hamming';
    elseif(~sum(strcmp(varargin{1},ww))); error_msg(6); return;
    elseif(sum(strcmp(varargin{1},ww2))), OPT1 = 1; win1 = varargin{1};
    else win1 = varargin{1};
    end
    win2 = 'hamming'; if(~mod(N,2)), L1 = N-1; L2 = N-1; else L1 = N; L2 = N; end; M = 2^nextpow2(N);
elseif(nargin==4)
    if(isempty(varargin{1})), win1 = 'hamming';
    elseif(~sum(strcmp(varargin{1},ww))); error_msg(6); return;
    elseif(sum(strcmp(varargin{1},ww2))), OPT1 = 1; win1 = varargin{1};
    else win1 = varargin{1};
    end
    if(isempty(varargin{2}))
        if(~mod(N,2)), L1 = N-1; else L1 = N; end;
    elseif(length(varargin{2})>1), error_msg(7); return;
    elseif(varargin{2} <= 0), error_msg(8); return;
    elseif(varargin{2} > N)
        if(~mod(N,2)), warning_msg(1,N-1); L1 = N-1;
        else warning_msg(1,N); L1 = N; end
    elseif(~mod(varargin{2},2))
        warning_msg(2,varargin{2}-1); L1 = varargin{2}-1;
    else L1 = varargin{2};
    end
    win2 = 'hamming'; if(~mod(N,2)), L2 = N-1; else L2 = N; end; M = 2^nextpow2(N);
elseif(nargin==5)
    if(isempty(varargin{1})), win1 = 'hamming';
    elseif(~sum(strcmp(varargin{1},ww))); error_msg(6); return;
    elseif(sum(strcmp(varargin{1},ww2))), OPT1 = 1; win1 = varargin{1};
    else win1 = varargin{1};
    end
    if(isempty(varargin{2}))
        if(~mod(N,2)), L1 = N-1; else L1 = N; end;
    elseif(length(varargin{2})>1), error_msg(7); return;
    elseif(varargin{2} <= 0), error_msg(8); return;
    elseif(varargin{2} > N)
        if(~mod(N,2)), warning_msg(1,N-1); L1 = N-1;
        else warning_msg(1,N); L1 = N; end
    elseif(~mod(varargin{2},2))
        warning_msg(2,varargin{2}-1); L1 = varargin{2}-1;
    else L1 = varargin{2};
    end
    if(sum(strcmp(win1,ww1)))
        if(isempty(varargin{3})), win2 = 'hamming';
        elseif(~sum(strcmp(varargin{3},ww))); error_msg(6); return;
        elseif(sum(strcmp(varargin{3},ww2))), OPT2 = 1; win2 = varargin{3};
        else win2 = varargin{3};
        end
        if(~mod(N,2)), L2 = N-1; else L2 = N; end; M = 2^nextpow2(N);
    elseif(sum(strcmp(win1,ww2)))
        if(isempty(varargin{3})), OPT1 = 1;
        elseif(length(varargin{3})>1), error_msg(9); return;
        else OPT1 = varargin{3};
        end
        win2 = 'hamming'; if(~mod(N,2)), L2 = N-1; else L2 = N; end; M = 2^nextpow2(N);
    end
elseif(nargin==6)
    if(isempty(varargin{1})), win1 = 'hamming';
    elseif(~sum(strcmp(varargin{1},ww))); error_msg(6); return;
    elseif(sum(strcmp(varargin{1},ww2))), OPT1 = 1; win1 = varargin{1};
    else win1 = varargin{1};
    end
    if(isempty(varargin{2}))
        if(~mod(N,2)), L1 = N-1; else L1 = N; end;
    elseif(length(varargin{2})>1), error_msg(7); return;
    elseif(varargin{2} <= 0), error_msg(8); return;
    elseif(varargin{2} > N)
        if(~mod(N,2)), warning_msg(1,N-1); L1 = N-1;
        else warning_msg(1,N); L1 = N; end
    elseif(~mod(varargin{2},2))
        warning_msg(2,varargin{2}-1); L1 = varargin{2}-1;
    else L1 = varargin{2};
    end
    if(sum(strcmp(win1,ww1)))
        if(isempty(varargin{3})), win2 = 'hamming';
        elseif(~sum(strcmp(varargin{3},ww))); error_msg(6); return;
        elseif(sum(strcmp(varargin{3},ww2))), OPT2 = 1; win2 = varargin{3};
        else win2 = varargin{3};
        end
        if(~mod(N,2)), L2 = N-1; else L2 = N; end; M = 2^nextpow2(N);
    elseif(sum(strcmp(win1,ww2)))
        if(isempty(varargin{3})), OPT1 = 1;
        elseif(length(varargin{3})>1), error_msg(9); return;
        else OPT1 = varargin{3};
        end
        win2 = 'hamming'; if(~mod(N,2)), L2 = N-1; else L2 = N; end; M = 2^nextpow2(N);
    end
    if(sum(strcmp(win1,ww2)))
        if(isempty(varargin{4})), win2 = 'hamming';
        elseif(~sum(strcmp(varargin{4},ww))); error_msg(6); return;
        elseif(sum(strcmp(varargin{4},ww2))), OPT2 = 1; win2 = varargin{4};
        else win2 = varargin{4};
        end
        if(~mod(N,2)), L2 = N-1; else L2 = N; end; M = 2^nextpow2(N);
    elseif(sum(strcmp(win1,ww1)))
        if(isempty(varargin{4}))
            if(~mod(N,2)), L2 = N-1; else L2 = N; end;
        elseif(length(varargin{4})>1), error_msg(7); return;
        elseif(varargin{4} <= 0), error_msg(8); return;
        elseif(varargin{4} > N)
            if(~mod(N,2)), warning_msg(1,N-1); L2 = N-1;
            else warning_msg(1,N); L2 = N; end
        elseif(~mod(varargin{4},2))
            warning_msg(2,varargin{4}-1); L2 = varargin{4}-1;
        else L2 = varargin{4};
        end
        M = 2^nextpow2(N);
    end
elseif(nargin==7)
    if(isempty(varargin{1})), win1 = 'hamming';
    elseif(~sum(strcmp(varargin{1},ww))); error_msg(6); return;
    elseif(sum(strcmp(varargin{1},ww2))), OPT1 = 1; win1 = varargin{1};
    else win1 = varargin{1};
    end
    if(isempty(varargin{2}))
        if(~mod(N,2)), L1 = N-1; else L1 = N; end;
    elseif(length(varargin{2})>1), error_msg(7); return;
    elseif(varargin{2} <= 0), error_msg(8); return;
    elseif(varargin{2} > N)
        if(~mod(N,2)), warning_msg(1,N-1); L1 = N-1;
        else warning_msg(1,N); L1 = N; end
    elseif(~mod(varargin{2},2))
        warning_msg(2,varargin{2}-1); L1 = varargin{2}-1;
    else L1 = varargin{2};
    end
    if(sum(strcmp(win1,ww1)))
        if(isempty(varargin{3})), win2 = 'hamming';
        elseif(~sum(strcmp(varargin{3},ww))); error_msg(6); return;
        elseif(sum(strcmp(varargin{3},ww2))), OPT2 = 1; win2 = varargin{3};
        else win2 = varargin{3};
        end
        if(~mod(N,2)), L2 = N-1; else L2 = N; end; M = 2^nextpow2(N);
    elseif(sum(strcmp(win1,ww2)))
        if(isempty(varargin{3})), OPT1 = 1;
        elseif(length(varargin{3})>1), error_msg(9); return;
        else OPT1 = varargin{3};
        end
        win2 = 'hamming'; if(~mod(N,2)), L2 = N-1; else L2 = N; end; M = 2^nextpow2(N);
    end 
    if(sum(strcmp(win1,ww2)))
        if(isempty(varargin{4})), win2 = 'hamming';
        elseif(~sum(strcmp(varargin{4},ww))); error_msg(6); return;
        elseif(sum(strcmp(varargin{4},ww2))), OPT2 = 1; win2 = varargin{4};
        else win2 = varargin{4};
        end
        if(~mod(N,2)), L2 = N-1; else L2 = N; end; M = 2^nextpow2(N);
    elseif(sum(strcmp(win1,ww1)))
        if(isempty(varargin{4}))
            if(~mod(N,2)), L2 = N-1; else L2 = N; end;
        elseif(length(varargin{4})>1), error_msg(7); return;
        elseif(varargin{4} <= 0), error_msg(8); return;
        elseif(varargin{4} > N)
            if(~mod(N,2)), warning_msg(1,N-1); L2 = N-1;
            else warning_msg(1,N); L2 = N; end
        elseif(~mod(varargin{4},2))
            warning_msg(2,varargin{4}-1); L2 = varargin{4}-1;
        else L2 = varargin{4};
        end
        M = 2^nextpow2(N);
    end
    if(sum(strcmp(win1,ww2)))
        if(isempty(varargin{5}))
            if(~mod(N,2)), L2 = N-1; else L2 = N; end;
        elseif(length(varargin{5})>1), error_msg(7); return;
        elseif(varargin{5} <= 0), error_msg(8); return;
        elseif(varargin{5} > N)
            if(~mod(N,2)), warning_msg(1,N-1); L2 = N-1;
            else warning_msg(1,N); L2 = N; end
        elseif(~mod(varargin{5},2))
            warning_msg(2,varargin{5}-1); L2 = varargin{5}-1;
        else L2 = varargin{5};
        end
        M = 2^nextpow2(N);
    elseif(sum(strcmp(win1,ww1)))
        if(sum(strcmp(win2,ww2)))
            if(isempty(varargin{5})), OPT2 = 1;
            elseif(length(varargin{5})>1), error_msg(9); return;
            else OPT2 = varargin{5};
            end
            M = 2^nextpow2(N);
        elseif(sum(strcmp(win2,ww1)))
            if(isempty(varargin{5})), M = 2^nextpow2(N);
            elseif(length(varargin{5})>1), error_msg(10); return;
            elseif(varargin{5} <= 0), error_msg(11); return;
            elseif(varargin{5} < L1), warning_msg(3,2^nextpow2(L1)); M = 2^nextpow2(L1);
            elseif(mod(varargin{5},2)), warning_msg(3,2^nextpow2(varargin{5})); M = 2^nextpow2(varargin{5});
            else M = varargin{5};
            end
        end
    end
elseif(nargin==8)
    if(isempty(varargin{1})), win1 = 'hamming';
    elseif(~sum(strcmp(varargin{1},ww))); error_msg(6); return;
    elseif(sum(strcmp(varargin{1},ww2))), OPT1 = 1; win1 = varargin{1};
    else win1 = varargin{1};
    end
    if(isempty(varargin{2}))
        if(~mod(N,2)), L1 = N-1; else L1 = N; end;
    elseif(length(varargin{2})>1), error_msg(7); return;
    elseif(varargin{2} <= 0), error_msg(8); return;
    elseif(varargin{2} > N)
        if(~mod(N,2)), warning_msg(1,N-1); L1 = N-1;
        else warning_msg(1,N); L1 = N; end
    elseif(~mod(varargin{2},2))
        warning_msg(2,varargin{2}-1); L1 = varargin{2}-1;
    else L1 = varargin{2};
    end
    if(sum(strcmp(win1,ww1)))
        if(isempty(varargin{3})), win2 = 'hamming';
        elseif(~sum(strcmp(varargin{3},ww))); error_msg(6); return;
        elseif(sum(strcmp(varargin{3},ww2))), OPT2 = 1; win2 = varargin{3};
        else win2 = varargin{3};
        end
        if(~mod(N,2)), L2 = N-1; else L2 = N; end; M = 2^nextpow2(N);
    elseif(sum(strcmp(win1,ww2)))
        if(isempty(varargin{3})), OPT1 = 1;
        elseif(length(varargin{3})>1), error_msg(9); return;
        else OPT1 = varargin{3};
        end
        win2 = 'hamming'; if(~mod(N,2)), L2 = N-1; else L2 = N; end; M = 2^nextpow2(N);
    end 
    if(sum(strcmp(win1,ww2)))
        if(isempty(varargin{4})), win2 = 'hamming';
        elseif(~sum(strcmp(varargin{4},ww))); error_msg(6); return;
        elseif(sum(strcmp(varargin{4},ww2))), OPT2 = 1; win2 = varargin{4};
        else win2 = varargin{4};
        end
        if(~mod(N,2)), L2 = N-1; else L2 = N; end; M = 2^nextpow2(N);
    elseif(sum(strcmp(win1,ww1)))
        if(isempty(varargin{4}))
            if(~mod(N,2)), L2 = N-1; else L2 = N; end;
        elseif(length(varargin{4})>1), error_msg(7); return;
        elseif(varargin{4} <= 0), error_msg(8); return;
        elseif(varargin{4} > N)
            if(~mod(N,2)), warning_msg(1,N-1); L2 = N-1;
            else warning_msg(1,N); L2 = N; end
        elseif(~mod(varargin{4},2))
            warning_msg(2,varargin{4}-1); L2 = varargin{4}-1;
        else L2 = varargin{4};
        end
        M = 2^nextpow2(N);
    end
    if(sum(strcmp(win1,ww2)))
        if(isempty(varargin{5}))
            if(~mod(N,2)), L2 = N-1; else L2 = N; end;
        elseif(length(varargin{5})>1), error_msg(7); return;
        elseif(varargin{5} <= 0), error_msg(8); return;
        elseif(varargin{5} > N)
            if(~mod(N,2)), warning_msg(1,N-1); L2 = N-1;
            else warning_msg(1,N); L2 = N; end
        elseif(~mod(varargin{5},2))
            warning_msg(2,varargin{5}-1); L2 = varargin{5}-1;
        else L2 = varargin{5};
        end
        M = 2^nextpow2(N);
    elseif(sum(strcmp(win1,ww1)))
        if(sum(strcmp(win2,ww2)))
            if(isempty(varargin{5})), OPT2 = 1;
            elseif(length(varargin{5})>1), error_msg(9); return;
            else OPT2 = varargin{5};
            end
            M = 2^nextpow2(N);
        elseif(sum(strcmp(win2,ww1)))
            if(isempty(varargin{5})), M = 2^nextpow2(N);
            elseif(length(varargin{5})>1), error_msg(10); return;
            elseif(varargin{5} <= 0), error_msg(11); return;
            elseif(varargin{5} < L1), warning_msg(3,2^nextpow2(L1)); M = 2^nextpow2(L1);
            elseif(mod(varargin{5},2)), warning_msg(3,2^nextpow2(varargin{5})); M = 2^nextpow2(varargin{5});
            else M = varargin{5};
            end
        end
    end
    if(sum(strcmp(win1,ww2)) && sum(strcmp(win2,ww2)))
        if(isempty(varargin{6})), OPT2 = 1;
        elseif(length(varargin{6})>1), error_msg(9); return;
        else OPT2 = varargin{6};
        end
    elseif(sum(strcmp(win1,ww2)) && sum(strcmp(win2,ww1)))
        if(isempty(varargin{6})), M = 2^nextpow2(N);
        elseif(length(varargin{6})>1), error_msg(10); return;
        elseif(varargin{6} <= 0), error_msg(11); return;
        elseif(varargin{6} < L1), warning_msg(3,2^nextpow2(L1)); M = 2^nextpow2(L1);
        elseif(mod(varargin{6},2)), warning_msg(3,2^nextpow2(varargin{6})); M = 2^nextpow2(varargin{6});
        else M = varargin{6};
        end
    elseif(sum(strcmp(win1,ww1)) && sum(strcmp(win2,ww2)))
        if(isempty(varargin{6})), M = 2^nextpow2(N);
        elseif(length(varargin{6})>1), error_msg(10); return;
        elseif(varargin{6} <= 0), error_msg(11); return;
        elseif(varargin{6} < L1), warning_msg(3,2^nextpow2(L1)); M = 2^nextpow2(L1);
        elseif(mod(varargin{6},2)), warning_msg(3,2^nextpow2(varargin{6})); M = 2^nextpow2(varargin{6});
        else M = varargin{6};
        end
    end
elseif(nargin==9)
    if(isempty(varargin{1})), win1 = 'hamming';
    elseif(~sum(strcmp(varargin{1},ww))); error_msg(6); return;
    elseif(sum(strcmp(varargin{1},ww2))), OPT1 = 1; win1 = varargin{1};
    else win1 = varargin{1};
    end
    if(isempty(varargin{2}))
        if(~mod(N,2)), L1 = N-1; else L1 = N; end;
    elseif(length(varargin{2})>1), error_msg(7); return;
    elseif(varargin{2} <= 0), error_msg(8); return;
    elseif(varargin{2} > N)
        if(~mod(N,2)), warning_msg(1,N-1); L1 = N-1;
        else warning_msg(1,N); L1 = N; end
    elseif(~mod(varargin{2},2))
        warning_msg(2,varargin{2}-1); L1 = varargin{2}-1;
    else L1 = varargin{2};
    end
    if(sum(strcmp(win1,ww1)))
        if(isempty(varargin{3})), win2 = 'hamming';
        elseif(~sum(strcmp(varargin{3},ww))); error_msg(6); return;
        elseif(sum(strcmp(varargin{3},ww2))), OPT2 = 1; win2 = varargin{3};
        else win2 = varargin{3};
        end
        if(~mod(N,2)), L2 = N-1; else L2 = N; end; M = 2^nextpow2(N);
    elseif(sum(strcmp(win1,ww2)))
        if(isempty(varargin{3})), OPT1 = 1;
        elseif(length(varargin{3})>1), error_msg(9); return;
        else OPT1 = varargin{3};
        end
        win2 = 'hamming'; if(~mod(N,2)), L2 = N-1; else L2 = N; end; M = 2^nextpow2(N);
    end 
    if(sum(strcmp(win1,ww2)))
        if(isempty(varargin{4})), win2 = 'hamming';
        elseif(~sum(strcmp(varargin{4},ww))); error_msg(6); return;
        elseif(sum(strcmp(varargin{4},ww2))), OPT2 = 1; win2 = varargin{4};
        else win2 = varargin{4};
        end
        if(~mod(N,2)), L2 = N-1; else L2 = N; end; M = 2^nextpow2(N);
    elseif(sum(strcmp(win1,ww1)))
        if(isempty(varargin{4}))
            if(~mod(N,2)), L2 = N-1; else L2 = N; end;
        elseif(length(varargin{4})>1), error_msg(7); return;
        elseif(varargin{4} <= 0), error_msg(8); return;
        elseif(varargin{4} > N)
            if(~mod(N,2)), warning_msg(1,N-1); L2 = N-1;
            else warning_msg(1,N); L2 = N; end
        elseif(~mod(varargin{4},2))
            warning_msg(2,varargin{4}-1); L2 = varargin{4}-1;
        else L2 = varargin{4};
        end
        M = 2^nextpow2(N);
    end
    if(sum(strcmp(win1,ww2)))
        if(isempty(varargin{5}))
            if(~mod(N,2)), L2 = N-1; else L2 = N; end;
        elseif(length(varargin{5})>1), error_msg(7); return;
        elseif(varargin{5} <= 0), error_msg(8); return;
        elseif(varargin{5} > N)
            if(~mod(N,2)), warning_msg(1,N-1); L2 = N-1;
            else warning_msg(1,N); L2 = N; end
        elseif(~mod(varargin{5},2))
            warning_msg(2,varargin{5}-1); L2 = varargin{5}-1;
        else L2 = varargin{5};
        end
        M = 2^nextpow2(N);
    elseif(sum(strcmp(win1,ww1)))
        if(sum(strcmp(win2,ww2)))
            if(isempty(varargin{5})), OPT2 = 1;
            elseif(length(varargin{5})>1), error_msg(9); return;
            else OPT2 = varargin{5};
            end
            M = 2^nextpow2(N);
        elseif(sum(strcmp(win2,ww1)))
            if(isempty(varargin{5})), M = 2^nextpow2(N);
            elseif(length(varargin{5})>1), error_msg(10); return;
            elseif(varargin{5} <= 0), error_msg(11); return;
            elseif(varargin{5} < L1), warning_msg(3,2^nextpow2(L1)); M = 2^nextpow2(L1);
            elseif(mod(varargin{5},2)), warning_msg(3,2^nextpow2(varargin{5})); M = 2^nextpow2(varargin{5});
            else M = varargin{5};
            end
        end
    end
    if(sum(strcmp(win1,ww2)) && sum(strcmp(win2,ww2)))
        if(isempty(varargin{6})), OPT2 = 1;
        elseif(length(varargin{6})>1), error_msg(9); return;
        else OPT2 = varargin{6};
        end
    elseif(sum(strcmp(win1,ww2)) && sum(strcmp(win2,ww1)))
        if(isempty(varargin{6})), M = 2^nextpow2(N);
        elseif(length(varargin{6})>1), error_msg(10); return;
        elseif(varargin{6} <= 0), error_msg(11); return;
        elseif(varargin{6} < L1), warning_msg(3,2^nextpow2(L1)); M = 2^nextpow2(L1);
        elseif(mod(varargin{6},2)), warning_msg(3,2^nextpow2(varargin{6})); M = 2^nextpow2(varargin{6});
        else M = varargin{6};
        end
    elseif(sum(strcmp(win1,ww1)) && sum(strcmp(win2,ww2)))
        if(isempty(varargin{6})), M = 2^nextpow2(N);
        elseif(length(varargin{6})>1), error_msg(10); return;
        elseif(varargin{6} <= 0), error_msg(11); return;
        elseif(varargin{6} < L1), warning_msg(3,2^nextpow2(L1)); M = 2^nextpow2(L1);
        elseif(mod(varargin{6},2)), warning_msg(3,2^nextpow2(varargin{6})); M = 2^nextpow2(varargin{6});
        else M = varargin{6};
        end
    end 
    if(sum(strcmp(win1,ww2)) && sum(strcmp(win2,ww2)))
        if(isempty(varargin{7})), M = 2^nextpow2(N);
        elseif(length(varargin{7})>1), error_msg(10); return;
        elseif(varargin{7} <= 0), error_msg(11); return;
        elseif(varargin{7} < L1), warning_msg(3,2^nextpow2(L1)); M = 2^nextpow2(L1);
        elseif(mod(varargin{7},2)), warning_msg(3,2^nextpow2(varargin{7})); M = 2^nextpow2(varargin{7});
        else M = varargin{7};
        end
    end
else error_msg(12); return;
end

%% Computing the Analytic Associate of the Input Signals S1 and S2
P = nextpow2(N); K = 2*(2^P);
z1 = fft(s1,K);           % Transforming S1 into the frequency domain
z1(1) = z1(1);            % Holding the zero frequency amplitude/energy
z1(2:K/2) = 2*z1(2:K/2);  % Accepting and doubling half of the frequencies
z1(K/2+1:K) = 0;          % Rejecting half of the frequencies
z1 = ifft(z1);            % Transforming Z into the time domain
z2 = fft(s2,K);           % Transforming S2 into the frequency domain
z2(1) = z2(1);            % Holding the zero frequency amplitude/energy
z2(2:K/2) = 2*z2(2:K/2);  % Accepting and doubling half of the frequencies
z2(K/2+1:K) = 0;          % Rejecting half of the frequencies
z2 = ifft(z2);            % Transforming Z into the time domain

%% Computing the Signal Kernel and Smoothing in the Time-Lag Domain
N = length(s1);
K_TL = zeros(K, N);
L_half = fix(L1/2);
if(sum(strcmp(win1,ww1)))
    g = window(str2func(win1),L1);
elseif(sum(strcmp(win1,ww2)))
    g = window(str2func(win1),L1,OPT1);
end
for n = 1:K
    for tau = -L_half:L_half
        G  = g(1 + tau + L_half);
        Z1 = z1(1 + rem(2*K + n - 1 + tau, K));
        Z2 = z2(1 + rem(2*K + n - 1 - tau, K));
        mm = 1 + rem(2*K + tau, M);
        K_TL(mm,n) = G*Z1*conj(Z2); % Smoothing by the Lag window G
    end
end

%% Smoothing by the Doppler window H
if(sum(strcmp(win2,ww1)))
    h = window(str2func(win2),L2);
elseif(sum(strcmp(win2,ww2)))
    h = window(str2func(win2),L2,OPT2);
end
for tau = -L_half:L_half
    mm = 1 + rem(2*K + tau, M);
    K_TL(mm,:) = conv(K_TL(mm,:), h, 'same');
end

%% Computing the Cross Smoothed Psuedo Wigner-Ville Distribution
TFR = zeros(M, N);
temp = fft(K_TL, M);
TFR(:,2:N) = temp(:,1:N-1);

%% Supplementary Functions
    function error_msg(n)
        switch n,
            case 1,  fprintf(2,'ERROR (No input signals)\n');
            case 2,  fprintf(2,'ERROR (Two input signals are required)\n');
            case 3,  fprintf(2,'ERROR (Input signals class must be double)\n');
            case 4,  fprintf(2,'ERROR (Input signals must be 1D)\n');
            case 5,  fprintf(2,'ERROR (Input signals must have the same size)\n');
            case 6,  fprintf(2,'ERROR (Window type is not valid. Type ''help window'')\n');
            case 7,  fprintf(2,'ERROR (Window length should be a 1x1 positive odd number)\n');
            case 8,  fprintf(2,'ERROR (Window length should be a positive odd number)\n');
            case 9,  fprintf(2,'ERROR (Window options should be a 1x1 number)\n');
            case 10, fprintf(2,'ERROR (FFT length should be a 1x1 power of 2)\n');
            case 11, fprintf(2,'ERROR (FFT length should be positive power of 2)\n'); 
            case 12, fprintf(2,'ERROR (Extra inputs, please check the function help)\n');
        end
    end
    function warning_msg(n, x)
        switch n,
            case 1,  fprintf('WARNING (Window length is truncated to %d)\n',x);
            case 2,  fprintf('WARNING (Window length must be odd, thus it is truncated to %d)\n',x);
            case 3,  fprintf('WARNING (FFT length is truncated to %d)\n',x);
        end
    end
end