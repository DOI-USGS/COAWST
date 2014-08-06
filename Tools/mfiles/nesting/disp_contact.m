function disp_contact (S);
  
%
% DISP_CONTACT:  Displays Nested Grids Contact Points unique values
%
% disp_contact(S)
%
% This function displays the Nested Grids Contact Points unique values
% for each Contact Region. The Contact Points Structure is rich and this
% simple function facilitates summarizing the Contact Points values.
%
% On Input:
%
%    S          Contact Points structure (struct array)
%                 See "contact.m" for full list of fields
%

% svn $Id: disp_contact.m 711 2014-01-23 20:36:13Z arango $
%=========================================================================%
%  Copyright (c) 2002-2014 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

for cr = 1:S.Ncontact,
  disp(' ');
  disp(['Contact Region = ', num2str(cr), blanks(4),                    ...
	'Donor Grid = ', num2str(S.contact(cr).donor_grid), blanks(4),  ...
	'Receiver Grid = ', num2str(S.contact(cr).receiver_grid),       ...
        blanks(4), '(only unique values are displayed)']);

  disp(' ');
  disp(['Irg_rho:', blanks(4), '[',                                     ...
	sprintf('%4i', unique(S.contact(cr).point.Irg_rho)), '  ]']);
  disp(['Idg_rho:', blanks(4), '[',                                     ...
	sprintf('%4i', unique(S.contact(cr).point.Idg_rho)), '  ]']);

  disp(' ');
  disp(['Jrg_rho:', blanks(4), '[',                                     ...
	sprintf('%4i', unique(S.contact(cr).point.Jrg_rho)), '  ]']);
  disp(['Jdg_rho:', blanks(4), '[',                                     ...
	sprintf('%4i', unique(S.contact(cr).point.Jdg_rho)), '  ]']);

  disp(' ');
  disp(['Irg_u:  ', blanks(4), '[',                                     ...
	sprintf('%4i', unique(S.contact(cr).point.Irg_u)), '  ]']);
  disp(['Idg_u:  ', blanks(4), '[',                                     ...
	sprintf('%4i', unique(S.contact(cr).point.Idg_u)), '  ]']);

  disp(' ');
  disp(['Jrg_u:  ', blanks(4), '[',                                     ...
	sprintf('%4i', unique(S.contact(cr).point.Jrg_u)), '  ]']);
  disp(['Jdg_u:  ', blanks(4), '[',                                     ...
	sprintf('%4i', unique(S.contact(cr).point.Jdg_u)), '  ]']);

  disp(' ');
  disp(['Irg_v:  ', blanks(4), '[',                                     ...
	sprintf('%4i', unique(S.contact(cr).point.Irg_v)), '  ]']);
  disp(['Idg_v:  ', blanks(4), '[',                                     ...
	sprintf('%4i', unique(S.contact(cr).point.Idg_v)), '  ]']);

  disp(' ');
  disp(['Jrg_v:  ', blanks(4), '[',                                     ...
	sprintf('%4i', unique(S.contact(cr).point.Jrg_v)), '  ]']);
  disp(['Jdg_v:  ', blanks(4), '[',                                     ...
	sprintf('%4i', unique(S.contact(cr).point.Jdg_v)), '  ]']);

  disp(' ');
end

return