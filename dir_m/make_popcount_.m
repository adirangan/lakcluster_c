function popcount_ = make_popcount_();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
%
% function output = make_popcount_();
%
% This function calculates the lookup table popcount_(uint8_input). ;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

popcount_ = uint8(sum(dec2bin(0:255)-'0',2));
