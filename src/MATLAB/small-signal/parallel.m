function equival = parallel(x)
%
% Syntax   : equival = parallel(x)
%
% Purpose  : Finds the value of parallel combination of component values
%
% Input    : x - any vector containing the data
%
% Output   : equival - value of the equivalent combination
%
% See also :
%
% Calls    :
%
% Call by  :
%

% (c) Copyright 1991 Joe H. Chow - All Rights Reserved
%
% History ( in reverse chronologocal order )
%
% Version  : 1.0
% Author   : Pierre N. Accari
% Date     : 8 October, 1991

equival= 1 / sum( ones(x)./x);
