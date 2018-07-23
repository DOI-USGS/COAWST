% SW_COPY    Copyright and licence information on SEAWATER library.
% =================================================================
% SW_COPY  $Id: sw_copy.m 330 2009-03-10 05:57:42Z arango $
%          Copyright (C) Phil Morgan 1993
%
%      SOFTWARE LICENCE AGREEMENT
%
%      1.0 Grant of Licence
%
%      1.1  The CSIRO Division of  Oceanography (herein referred to as
%           "CSIRO") hereby grants you (hereinafter  referred to  as
%           the  "Licencee"),  subject  to  the  Licencee agreeing  to
%           comply with the terms and  conditions of this Agreement, a
%           non-transferable,   non-exclusive  licence   to  use   the
%           computer  programs described in this document (hereinafter
%           referred to  as the   "Software")  for  the  purpose   of
%           the  Licencee's computing activity.
%
%      1.2  CSIRO hereby grants the Licencee the right to  make copies
%           of  the  Software  for   the  purpose  of  the  Licencee's
%           computing activity only.
%
%      1.3  The benefit of  the rights granted to the Licencee  by the
%           Licence and this Agreement  generally shall be personal to
%           the Licencee and the  Licencee shall not mortgage, charge,
%           assign,  rent,  lease,  sell or  otherwise  dispose  of or
%           transfer the same or any part to any third party.
%
%      1.4  Unless  otherwise agreed  in  writing or  provided  for in
%           this  Agreement, CSIRO  shall  be under  no  obligation or
%           responsibility  to provide the Licencee with any training,
%           maintenance  services,  enhancements  or  updates  of  the
%           Software or any services whatsoever.
%
%      2.0 Acknowledgment by the Licencee
%
%      2.1  The Licencee acknowledges and agrees that it shall not:
%
%           (i)   sell,  let for  hire or  by way  of  trade, offer  or
%                 exhibit  or  expose  for sale  or  hire  or otherwise
%             distribute the Software for  the purposes of trade or
%                 any other purpose;
%
%           (ii)  authorise  or assist any third person  to do any
%                 of the acts set out in (i)  above;
%
%           (iii) modify  the  Software source  code  without  advising
%                 CSIRO.
%
%      2.2 The Licencee agrees that:
%
%           (a)  CSIRO  is  the  owner  of  all  copyright  and  other
%                Intellectual  Property   Rights  subsisting  in   the
%                Software;
%
%           (b)  this   document  must  be   properly   cited  in  any
%                publication  reporting  results  derived   from  this
%                document or obtained from application and use of this
%                software. Any  of   the  Licencee's  documentation
%                describing  results  generated  by  the  Licencee's
%                use  of  the Software will contain an acknowledgement
%                of CSIRO's  ownership of the Software;
%
%           (c)  CSIRO reserves all rights  in the Software other than
%                the  rights   granted  to   the   Licencee  by   this
%                Agreement;
%
%           (d)  each  item  of  the Software  will  display  a  banner
%                summarising   the  terms   of   this   Agreement  and
%                acknowledging  the source  of  the Software,  and the
%                contents of  a banner  will not be  modified and  its
%                display  will  not  be  inactivated  by  the Licencee
%                without the approval of CSIRO.
%
%      3.0 Indemnity
%
%      3.1  To the full  extent permitted by  law, CSIRO  excludes any
%           and  all liability  in  respect  of any  loss  or  damage,
%           whether  personal  (includes  death  or  illness)  or  of
%           property  and  whether direct,  consequential  or  special
%           (including  consequential  financial  loss or  damage) of
%           the Licencee, its  officers, agents and  employees or  any
%           third party  howsoever caused,  which may  be suffered  or
%           incurred  or which  may  arise directly  or  indirectly in
%           respect of  or  arising  out  of  the  Licencee's  use  or
%           inability to use  the Software or the failure or  omission
%           on the  part of CSIRO  to comply  with the conditions  and
%           warranties under  this  Licence  Agreement.    Insofar  as
%           liability for  loss or damages  under or  pursuant to such
%           legislation cannot  be  excluded,  CSIRO's  liability  for
%           loss or  damages shall  be limited  to the  amount of  One
%           Dollar ($1.00).
%
%      3.2  CSIRO  make  no  warranties,  expressed  or  implied,  and
%           excludes all  other warranties  representations, terms  or
%           conditions, whether  express or implied,  oral or written,
%           statutory  or  otherwise,  relating  in  any  way  to  the
%           Software, or  to  this  Agreement, including  any  implied
%           warranty of merchantability  or of fitness for  particular
%           purpose.   To the full extent permitted  by the law of the
%           Commonwealth of  Australia  or the  laws of  any State  or
%           Territory  of  Australia,  any  conditions  or  warranties
%           imposed by such  legislation are hereby  excluded.   In so
%           far as  liability under  or pursuant  to such  legislation
%           may not  be excluded,  CSIRO's liability  to the  Licencee
%           pursuant to this Agreement shall be limited as  set out in
%           clause 3.1 hereof.
%
%      3.3  The  Licencee acknowledges  and agrees  that  the Software
%           was developed  for CSIRO  research purposes  and may  have
%           inherent  defects, errors or deficiencies, and  that it is
%           the  responsibility  of  the   Licencee  to  make  its  own
%           assessment  of the  suitability  of the  Software  for the
%           purpose  of  the   Licencee's  computing  activity.    The
%           Licencee will use  the Software, and  advice, opinions  or
%           information supplied by CSIRO,  its officers, employees or
%           agents  concerning  the  Software  at the  Licencee's  own
%           risk.
%
%      3.4  The  Licencee hereby  releases and  indemnifies and  shall
%           continue to  release and  indemnify  CSIRO, its  officers,
%           employees  and  agents  from  and  against  all   actions,
%           claims, proceedings  or demands  (including those  brought
%           by third parties) which may be bought against  it or them,
%           whether  on their  own or  jointly  with the  Licencee and
%           whether at common law, in equity or pursuant to statute or
%           otherwise, in respect of  any loss, death, injury, illness
%           or  damage  (whether  personal  or  property,  and whether
%           direct   or    consequential,   including    consequential
%           financial  loss)   and  any  infringement  of   copyright,
%           patents,   trade  marks,  designs  or  other  Intellectual
%           Property Rights, howsoever  arising out of the  Licencee's
%           exercise of its  rights under this  Agreement and from and
%           against  all  damages,  costs  and  expenses  incurred  in
%           defending  or  settling  any  such  claim,  proceeding  or
%           demand.
%
%      3.5  The  Licencee's  obligation to  indemnify  CSIRO  and  its
%           officers,  employees  and  agents set  out  in  clause 3.4
%           hereof   is   a   continuing   obligation   separate  from
%           and independent of the Licencee's other obligations  under
%           this  Agreement,  and  shall  survive  all  expiration  or
%           termination of this Agreement.
%
%      4.0 Termination
%
%      4.1  The Licence shall terminate  immediately upon the Licencee
%           breaching any term or  condition of this Agreement whether
%           or  not CSIRO is aware of  the occurrence of the breach at
%           the time that it happens.
%
%      4.2  CSIRO  may terminate the Licence on  reasonable grounds by
%           notice in  writing  to the  Licencee, and  such notice  of
%           termination shall  be effective  immediately upon  receipt
%           by  the  Licencee.
%
%

more on
help sw_copy
more off
return
