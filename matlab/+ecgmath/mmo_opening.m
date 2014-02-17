function Result = mmo_opening(Signal, varargin)
% Aplica a operaçao de abertura no sinal de acordo com um elemento
% estruturador B, ou com um par (B1,B2)

if nargin == 2
    B = varargin{1};
    Result = mmo_dilation(mmo_erosion(Signal, B), B);
else
    B1 = varargin{1};
    B2 = varargin{2};
    Result = mmo_dilation(mmo_erosion(Signal, B1), B2);
end
