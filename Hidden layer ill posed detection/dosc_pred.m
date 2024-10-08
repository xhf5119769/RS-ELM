function [Znew,Tnew] = dosc_pred(Xnew,W,P);
%
% [Znew,Tnew] = dosc_pred(Xnew,W,P);
%
% dosc_pred is used to remove DOSC components from new data
%
% Input
%    Xnew    New data (usually spectra)            (Inew x J)
%    W       Weights to calculate DOSC component   (J x nocomp)
%    P       Loadings to remove DOSC component     (J x nocomp)
%
% Output
%    Znew    DOSC corrected new data               (Inew x J)
%    Tnew    value of DOSC component for new data  (Inew x nocomp)
% 
% See Reference 
% Westerhuis JA, de Jong S and Smilde AK, Direct orthogonal signal correction, 
% Chemometrics and Intelligent Laboratory Systems, 56, (2001), 13-25.

% Johan Westerhuis
%   ==========================================================================
%   Copyright 2005 Biosystems Data Analysis Group ; Universiteit van Amsterdam
%   ==========================================================================

Tnew = Xnew * W;
Znew = Xnew - Tnew * P';