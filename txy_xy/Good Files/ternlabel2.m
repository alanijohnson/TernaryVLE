% TERNLABEL label ternary phase diagram
%   TERNLABEL('ALABEL', 'BLABEL', 'CLABEL') labels a ternary phase diagram created using TERNPLOT
%   
%   H = TERNLABEL('ALABEL', 'BLABEL', 'CLABEL') returns handles to the text objects created.
%   with the labels provided.  TeX escape codes are accepted.
%
%   See also TERNPLOT

% Author: Carl Sandrock 20020827

% To Do

% Modifications

% Modifiers

function h = ternlabel2(A, B, C)
r(1) = text(0.5, -0.45, A, 'horizontalalignment', 'center');
r(2) = text(1.4, 0.1, B, 'horizontalalignment', 'center');
r(3) = text(0.2*sin(deg2rad(30)), 0.5, C, 'rotation', 105, 'horizontalalignment', 'center');

if nargout > 0
    h = r;
end;