1) Extensions to unstructured mesh implementation:
   *) Alleviation of the garden-sprinkler effect,
   *) Phase-decoupled diffraction.
2) An alternative to the JONSWAP formulation for bottom friction
   is included (friction coefficient depends on the frequency-
   dependent directional spreading).
3) The SORDUP scheme is made more simpler and more consistent.
4) Bug fixes:
   *) correction of outputting the wave force on unstructured meshes
      in the case of spherical coordinates,
   *) remove small inconsistency in the concatenation program for
      hotfiles,
   *) small correction in collecting data for MPI parallel runs,
   *) no use of the Hersbach and Janssen limiter in the case of
      stationary runs,
   *) the user-value of water density will not be overwritten
      with the default value 1025 kg/m3.

--- CV.eps	1970-01-01 01:00:00.000000000 +0100
+++ CV.eps	2009-03-03 17:50:34.000000000 +0100
@@ -0,0 +1,650 @@
+%!PS-Adobe-2.0 EPSF-2.0
+%%Title: CV.fig
+%%Creator: fig2dev Version 3.2 Patchlevel 4
+%%CreationDate: Thu Feb 26 16:02:19 2009
+%%For: mzijlema@TUD11321 (Marcel Zijlema)
+%%BoundingBox: 0 0 379 319
+%%Magnification: 1.0000
+%%EndComments
+/MyAppDict 100 dict dup begin def
+/$F2psDict 200 dict def
+$F2psDict begin
+$F2psDict /mtrx matrix put
+/col-1 {0 setgray} bind def
+/col0 {0.000 0.000 0.000 srgb} bind def
+/col1 {0.000 0.000 1.000 srgb} bind def
+/col2 {0.000 1.000 0.000 srgb} bind def
+/col3 {0.000 1.000 1.000 srgb} bind def
+/col4 {1.000 0.000 0.000 srgb} bind def
+/col5 {1.000 0.000 1.000 srgb} bind def
+/col6 {1.000 1.000 0.000 srgb} bind def
+/col7 {1.000 1.000 1.000 srgb} bind def
+/col8 {0.000 0.000 0.560 srgb} bind def
+/col9 {0.000 0.000 0.690 srgb} bind def
+/col10 {0.000 0.000 0.820 srgb} bind def
+/col11 {0.530 0.810 1.000 srgb} bind def
+/col12 {0.000 0.560 0.000 srgb} bind def
+/col13 {0.000 0.690 0.000 srgb} bind def
+/col14 {0.000 0.820 0.000 srgb} bind def
+/col15 {0.000 0.560 0.560 srgb} bind def
+/col16 {0.000 0.690 0.690 srgb} bind def
+/col17 {0.000 0.820 0.820 srgb} bind def
+/col18 {0.560 0.000 0.000 srgb} bind def
+/col19 {0.690 0.000 0.000 srgb} bind def
+/col20 {0.820 0.000 0.000 srgb} bind def
+/col21 {0.560 0.000 0.560 srgb} bind def
+/col22 {0.690 0.000 0.690 srgb} bind def
+/col23 {0.820 0.000 0.820 srgb} bind def
+/col24 {0.500 0.190 0.000 srgb} bind def
+/col25 {0.630 0.250 0.000 srgb} bind def
+/col26 {0.750 0.380 0.000 srgb} bind def
+/col27 {1.000 0.500 0.500 srgb} bind def
+/col28 {1.000 0.630 0.630 srgb} bind def
+/col29 {1.000 0.750 0.750 srgb} bind def
+/col30 {1.000 0.880 0.880 srgb} bind def
+/col31 {1.000 0.840 0.000 srgb} bind def
+
+end
+save
+newpath 0 319 moveto 0 0 lineto 379 0 lineto 379 319 lineto closepath clip newpath
+-298.8 450.7 translate
+1 -1 scale
+
+% This junk string is used by the show operators
+/PATsstr 1 string def
+/PATawidthshow { 	% cx cy cchar rx ry string
+  % Loop over each character in the string
+  {  % cx cy cchar rx ry char
+    % Show the character
+    dup				% cx cy cchar rx ry char char
+    PATsstr dup 0 4 -1 roll put	% cx cy cchar rx ry char (char)
+    false charpath		% cx cy cchar rx ry char
+    /clip load PATdraw
+    % Move past the character (charpath modified the
+    % current point)
+    currentpoint			% cx cy cchar rx ry char x y
+    newpath
+    moveto			% cx cy cchar rx ry char
+    % Reposition by cx,cy if the character in the string is cchar
+    3 index eq {			% cx cy cchar rx ry
+      4 index 4 index rmoveto
+    } if
+    % Reposition all characters by rx ry
+    2 copy rmoveto		% cx cy cchar rx ry
+  } forall
+  pop pop pop pop pop		% -
+  currentpoint
+  newpath
+  moveto
+} bind def
+/PATcg {
+  7 dict dup begin
+    /lw currentlinewidth def
+    /lc currentlinecap def
+    /lj currentlinejoin def
+    /ml currentmiterlimit def
+    /ds [ currentdash ] def
+    /cc [ currentrgbcolor ] def
+    /cm matrix currentmatrix def
+  end
+} bind def
+% PATdraw - calculates the boundaries of the object and
+% fills it with the current pattern
+/PATdraw {			% proc
+  save exch
+    PATpcalc			% proc nw nh px py
+    5 -1 roll exec		% nw nh px py
+    newpath
+    PATfill			% -
+  restore
+} bind def
+% PATfill - performs the tiling for the shape
+/PATfill { % nw nh px py PATfill -
+  PATDict /CurrentPattern get dup begin
+    setfont
+    % Set the coordinate system to Pattern Space
+    PatternGState PATsg
+    % Set the color for uncolored pattezns
+    PaintType 2 eq { PATDict /PColor get PATsc } if
+    % Create the string for showing
+    3 index string		% nw nh px py str
+    % Loop for each of the pattern sources
+    0 1 Multi 1 sub {		% nw nh px py str source
+	% Move to the starting location
+	3 index 3 index		% nw nh px py str source px py
+	moveto			% nw nh px py str source
+	% For multiple sources, set the appropriate color
+	Multi 1 ne { dup PC exch get PATsc } if
+	% Set the appropriate string for the source
+	0 1 7 index 1 sub { 2 index exch 2 index put } for pop
+	% Loop over the number of vertical cells
+	3 index 		% nw nh px py str nh
+	{			% nw nh px py str
+	  currentpoint		% nw nh px py str cx cy
+	  2 index oldshow	% nw nh px py str cx cy
+	  YStep add moveto	% nw nh px py str
+	} repeat		% nw nh px py str
+    } for
+    5 { pop } repeat
+  end
+} bind def
+
+% PATkshow - kshow with the current pattezn
+/PATkshow {			% proc string
+  exch bind			% string proc
+  1 index 0 get			% string proc char
+  % Loop over all but the last character in the string
+  0 1 4 index length 2 sub {
+				% string proc char idx
+    % Find the n+1th character in the string
+    3 index exch 1 add get	% string proc char char+1
+    exch 2 copy			% strinq proc char+1 char char+1 char
+    % Now show the nth character
+    PATsstr dup 0 4 -1 roll put	% string proc chr+1 chr chr+1 (chr)
+    false charpath		% string proc char+1 char char+1
+    /clip load PATdraw
+    % Move past the character (charpath modified the current point)
+    currentpoint newpath moveto
+    % Execute the user proc (should consume char and char+1)
+    mark 3 1 roll		% string proc char+1 mark char char+1
+    4 index exec		% string proc char+1 mark...
+    cleartomark			% string proc char+1
+  } for
+  % Now display the last character
+  PATsstr dup 0 4 -1 roll put	% string proc (char+1)
+  false charpath		% string proc
+  /clip load PATdraw
+  neewath
+  pop pop			% -
+} bind def
+% PATmp - the makepattern equivalent
+/PATmp {			% patdict patmtx PATmp patinstance
+  exch dup length 7 add		% We will add 6 new entries plus 1 FID
+  dict copy			% Create a new dictionary
+  begin
+    % Matrix to install when painting the pattern
+    TilingType PATtcalc
+    /PatternGState PATcg def
+    PatternGState /cm 3 -1 roll put
+    % Check for multi pattern sources (Level 1 fast color patterns)
+    currentdict /Multi known not { /Multi 1 def } if
+    % Font dictionary definitions
+    /FontType 3 def
+    % Create a dummy encoding vector
+    /Encoding 256 array def
+    3 string 0 1 255 {
+      Encoding exch dup 3 index cvs cvn put } for pop
+    /FontMatrix matrix def
+    /FontBBox BBox def
+    /BuildChar {
+	mark 3 1 roll		% mark dict char
+	exch begin
+	Multi 1 ne {PaintData exch get}{pop} ifelse  % mark [paintdata]
+	  PaintType 2 eq Multi 1 ne or
+	  { XStep 0 FontBBox aload pop setcachedevice }
+	  { XStep 0 setcharwidth } ifelse
+	  currentdict		% mark [paintdata] dict
+	  /PaintProc load	% mark [paintdata] dict paintproc
+	end
+	gsave
+	  false PATredef exec true PATredef
+	grestore
+	cleartomark		% -
+    } bind def
+    currentdict
+  end				% newdict
+  /foo exch			% /foo newlict
+  definefont			% newfont
+} bind def
+% PATpcalc - calculates the starting point and width/height
+% of the tile fill for the shape
+/PATpcalc {	% - PATpcalc nw nh px py
+  PATDict /CurrentPattern get begin
+    gsave
+	% Set up the coordinate system to Pattern Space
+	% and lock down pattern
+	PatternGState /cm get setmatrix
+	BBox aload pop pop pop translate
+	% Determine the bounding box of the shape
+	pathbbox			% llx lly urx ury
+    grestore
+    % Determine (nw, nh) the # of cells to paint width and height
+    PatHeight div ceiling		% llx lly urx qh
+    4 1 roll				% qh llx lly urx
+    PatWidth div ceiling		% qh llx lly qw
+    4 1 roll				% qw qh llx lly
+    PatHeight div floor			% qw qh llx ph
+    4 1 roll				% ph qw qh llx
+    PatWidth div floor			% ph qw qh pw
+    4 1 roll				% pw ph qw qh
+    2 index sub cvi abs			% pw ph qs qh-ph
+    exch 3 index sub cvi abs exch	% pw ph nw=qw-pw nh=qh-ph
+    % Determine the starting point of the pattern fill
+    %(px, py)
+    4 2 roll				% nw nh pw ph
+    PatHeight mul			% nw nh pw py
+    exch				% nw nh py pw
+    PatWidth mul exch			% nw nh px py
+  end
+} bind def
+
+% Save the original routines so that we can use them later on
+/oldfill	/fill load def
+/oldeofill	/eofill load def
+/oldstroke	/stroke load def
+/oldshow	/show load def
+/oldashow	/ashow load def
+/oldwidthshow	/widthshow load def
+/oldawidthshow	/awidthshow load def
+/oldkshow	/kshow load def
+
+% These defs are necessary so that subsequent procs don't bind in
+% the originals
+/fill	   { oldfill } bind def
+/eofill	   { oldeofill } bind def
+/stroke	   { oldstroke } bind def
+/show	   { oldshow } bind def
+/ashow	   { oldashow } bind def
+/widthshow { oldwidthshow } bind def
+/awidthshow { oldawidthshow } bind def
+/kshow 	   { oldkshow } bind def
+/PATredef {
+  MyAppDict begin
+    {
+    /fill { /clip load PATdraw newpath } bind def
+    /eofill { /eoclip load PATdraw newpath } bind def
+    /stroke { PATstroke } bind def
+    /show { 0 0 null 0 0 6 -1 roll PATawidthshow } bind def
+    /ashow { 0 0 null 6 3 roll PATawidthshow }
+    bind def
+    /widthshow { 0 0 3 -1 roll PATawidthshow }
+    bind def
+    /awidthshow { PATawidthshow } bind def
+    /kshow { PATkshow } bind def
+  } {
+    /fill   { oldfill } bind def
+    /eofill { oldeofill } bind def
+    /stroke { oldstroke } bind def
+    /show   { oldshow } bind def
+    /ashow  { oldashow } bind def
+    /widthshow { oldwidthshow } bind def
+    /awidthshow { oldawidthshow } bind def
+    /kshow  { oldkshow } bind def
+    } ifelse
+  end
+} bind def
+false PATredef
+% Conditionally define setcmykcolor if not available
+/setcmykcolor where { pop } {
+  /setcmykcolor {
+    1 sub 4 1 roll
+    3 {
+	3 index add neg dup 0 lt { pop 0 } if 3 1 roll
+    } repeat
+    setrgbcolor - pop
+  } bind def
+} ifelse
+/PATsc {		% colorarray
+  aload length		% c1 ... cn length
+    dup 1 eq { pop setgray } { 3 eq { setrgbcolor } { setcmykcolor
+  } ifelse } ifelse
+} bind def
+/PATsg {		% dict
+  begin
+    lw setlinewidth
+    lc setlinecap
+    lj setlinejoin
+    ml setmiterlimit
+    ds aload pop setdash
+    cc aload pop setrgbcolor
+    cm setmatrix
+  end
+} bind def
+
+/PATDict 3 dict def
+/PATsp {
+  true PATredef
+  PATDict begin
+    /CurrentPattern exch def
+    % If it's an uncolored pattern, save the color
+    CurrentPattern /PaintType get 2 eq {
+      /PColor exch def
+    } if
+    /CColor [ currentrgbcolor ] def
+  end
+} bind def
+% PATstroke - stroke with the current pattern
+/PATstroke {
+  countdictstack
+  save
+  mark
+  {
+    currentpoint strokepath moveto
+    PATpcalc				% proc nw nh px py
+    clip newpath PATfill
+    } stopped {
+	(*** PATstroke Warning: Path is too complex, stroking
+	  with gray) =
+    cleartomark
+    restore
+    countdictstack exch sub dup 0 gt
+	{ { end } repeat } { pop } ifelse
+    gsave 0.5 setgray oldstroke grestore
+  } { pop restore pop } ifelse
+  newpath
+} bind def
+/PATtcalc {		% modmtx tilingtype PATtcalc tilematrix
+  % Note: tiling types 2 and 3 are not supported
+  gsave
+    exch concat					% tilingtype
+    matrix currentmatrix exch			% cmtx tilingtype
+    % Tiling type 1 and 3: constant spacing
+    2 ne {
+	% Distort the pattern so that it occupies
+	% an integral number of device pixels
+	dup 4 get exch dup 5 get exch		% tx ty cmtx
+	XStep 0 dtransform
+	round exch round exch			% tx ty cmtx dx.x dx.y
+	XStep div exch XStep div exch		% tx ty cmtx a b
+	0 YStep dtransform
+	round exch round exch			% tx ty cmtx a b dy.x dy.y
+	YStep div exch YStep div exch		% tx ty cmtx a b c d
+	7 -3 roll astore			% { a b c d tx ty }
+    } if
+  grestore
+} bind def
+/PATusp {
+  false PATredef
+  PATDict begin
+    CColor PATsc
+  end
+} bind def
+
+% this is the pattern fill program from the Second edition Reference Manual
+% with changes to call the above pattern fill
+% left30
+11 dict begin
+/PaintType 1 def
+/PatternType 1 def
+/TilingType 1 def
+/BBox [0 0 1 1] def
+/XStep 1 def
+/YStep 1 def
+/PatWidth 1 def
+/PatHeight 1 def
+/Multi 2 def
+/PaintData [
+  { clippath } bind
+  { 32 16 true [ 32 0 0 -16 0 16 ]
+	{<c000c000300030000c000c000300030000c000c000300030
+	000c000c00030003c000c000300030000c000c0003000300
+	00c000c000300030000c000c00030003>}
+     imagemask } bind
+] def
+/PaintProc {
+	pop
+	exec fill
+} def
+currentdict
+end
+/P1 exch def
+
+/cp {closepath} bind def
+/ef {eofill} bind def
+/gr {grestore} bind def
+/gs {gsave} bind def
+/sa {save} bind def
+/rs {restore} bind def
+/l {lineto} bind def
+/m {moveto} bind def
+/rm {rmoveto} bind def
+/n {newpath} bind def
+/s {stroke} bind def
+/sh {show} bind def
+/slc {setlinecap} bind def
+/slj {setlinejoin} bind def
+/slw {setlinewidth} bind def
+/srgb {setrgbcolor} bind def
+/rot {rotate} bind def
+/sc {scale} bind def
+/sd {setdash} bind def
+/ff {findfont} bind def
+/sf {setfont} bind def
+/scf {scalefont} bind def
+/sw {stringwidth} bind def
+/tr {translate} bind def
+/tnt {dup dup currentrgbcolor
+  4 -2 roll dup 1 exch sub 3 -1 roll mul add
+  4 -2 roll dup 1 exch sub 3 -1 roll mul add
+  4 -2 roll dup 1 exch sub 3 -1 roll mul add srgb}
+  bind def
+/shd {dup dup currentrgbcolor 4 -2 roll mul 4 -2 roll mul
+  4 -2 roll mul srgb} bind def
+ /DrawEllipse {
+	/endangle exch def
+	/startangle exch def
+	/yrad exch def
+	/xrad exch def
+	/y exch def
+	/x exch def
+	/savematrix mtrx currentmatrix def
+	x y tr xrad yrad sc 0 0 1 startangle endangle arc
+	closepath
+	savematrix setmatrix
+	} def
+
+/$F2psBegin {$F2psDict begin /$F2psEnteredState save def} def
+/$F2psEnd {$F2psEnteredState restore end} def
+
+$F2psBegin
+10 setmiterlimit
+0 slj 0 slc
+ 0.06299 0.06299 sc
+%
+% Fig objects follow
+%
+%
+% here starts figure with depth 50
+% Ellipse
+7.500 slw
+n 7635 3360 67 67 0 360 DrawEllipse gs col0 s gr
+
+% Ellipse
+n 6645 5415 67 67 0 360 DrawEllipse gs col0 s gr
+
+% Ellipse
+n 9120 4800 67 67 0 360 DrawEllipse gs col0 s gr
+
+% Ellipse
+n 8955 3570 67 67 0 360 DrawEllipse gs col0 s gr
+
+% Ellipse
+n 8070 5715 67 67 0 360 DrawEllipse gs col0 s gr
+
+% Ellipse
+n 6426 4232 67 67 0 360 DrawEllipse gs col0 s gr
+
+% Polyline
+n 8835 2280 m 6105 3120 l 5430 4935 l 6705 6885 l 9465 5880 l 10140 3990 l
+
+ 8790 2280 l gs col0 s gr
+% Polyline
+n 7631 3361 m 6431 4261 l 6641 5431 l 8051 5731 l 9116 4816 l 8951 3571 l
+
+ 7631 3361 l  cp gs /PC [[1.00 1.00 1.00] [0.00 0.00 0.00]] def
+15.00 15.00 sc P1 [16 0 0 -8 428.73 224.07] PATmp PATsp ef gr PATusp gs col0 s gr
+% Polyline
+n 7813 4532 m
+ 9433 5867 l gs col0 s gr
+% Polyline
+n 7816 4504 m
+ 8791 2299 l gs col0 s gr
+% Polyline
+n 7846 4505 m
+ 10156 3980 l gs col0 s gr
+% Polyline
+n 7810 4534 m
+ 5455 4939 l gs col0 s gr
+% Polyline
+15.000 slw
+gs  clippath
+8635 6906 m 8746 6861 l 8637 6595 l 8673 6840 l 8526 6640 l cp
+eoclip
+n 7950 5085 m
+ 8685 6870 l gs col0 s gr gr
+
+% arrowhead
+n 8526 6640 m 8673 6840 l 8637 6595 l 8526 6640 l  cp gs 0.00 setgray ef gr  col0 s
+% Polyline
+7.500 slw
+ [60] 0 sd
+n 6645 5445 m 5820 5850 l
+ 5850 5835 l gs col0 s gr  [] 0 sd
+% Polyline
+ [60] 0 sd
+n 6420 4230 m
+ 5475 3900 l gs col0 s gr  [] 0 sd
+% Polyline
+ [60] 0 sd
+n 7650 3345 m
+ 7395 2460 l gs col0 s gr  [] 0 sd
+% Polyline
+ [60] 0 sd
+n 8970 3570 m
+ 9630 2910 l gs col0 s gr  [] 0 sd
+% Polyline
+ [60] 0 sd
+n 9120 4800 m
+ 9915 4905 l gs col0 s gr  [] 0 sd
+% Polyline
+ [60] 0 sd
+n 8070 5745 m
+ 7785 6645 l gs col0 s gr  [] 0 sd
+% Polyline
+ [60] 0 sd
+n 5535 4140 m 4755 4140 l 4755 3705 l 5535 3705 l
+ cp gs col0 s gr  [] 0 sd
+% Polyline
+ [60] 0 sd
+n 7710 2490 m 6900 2490 l 6900 2115 l 7710 2115 l
+ cp gs col0 s gr  [] 0 sd
+% Polyline
+ [60] 0 sd
+n 10365 3075 m 9555 3075 l 9555 2670 l 10365 2670 l
+ cp gs col0 s gr  [] 0 sd
+% Polyline
+ [60] 0 sd
+n 10740 5100 m 9930 5100 l 9930 4710 l 10740 4710 l
+ cp gs col0 s gr  [] 0 sd
+% Polyline
+ [60] 0 sd
+n 8085 7080 m 7380 7080 l 7380 6675 l 8085 6675 l
+ cp gs col0 s gr  [] 0 sd
+% Polyline
+ [60] 0 sd
+n 5910 6240 m 5220 6240 l 5220 5820 l 5910 5820 l
+ cp gs col0 s gr  [] 0 sd
+% Polyline
+n 7784 4528 m
+ 6134 3133 l gs col0 s gr
+% Polyline
+n 7815 4500 m
+ 6720 6855 l gs col0 s gr
+% Polyline
+2 slj
+n 7815 4425 m 7814 4425 l 7810 4426 l 7801 4427 l 7789 4430 l 7776 4433 l
+ 7763 4438 l 7751 4445 l 7740 4455 l 7733 4466 l 7730 4474 l
+ 7730 4479 l 7731 4481 l 7733 4484 l 7735 4489 l 7737 4500 l
+ 7740 4515 l 7744 4531 l 7747 4544 l 7749 4552 l 7750 4558 l
+ 7751 4563 l 7754 4570 l 7760 4579 l 7770 4590 l 7779 4596 l
+ 7787 4599 l 7792 4601 l 7795 4601 l 7798 4600 l 7800 4599 l
+ 7804 4598 l 7810 4596 l 7819 4594 l 7830 4590 l 7841 4585 l
+ 7851 4581 l 7857 4578 l 7862 4576 l 7865 4575 l 7868 4574 l
+ 7872 4571 l 7877 4566 l 7883 4557 l 7890 4545 l 7896 4528 l
+ 7898 4515 l 7899 4506 l 7899 4500 l 7898 4494 l 7897 4485 l
+ 7894 4472 l 7890 4455 l 7883 4440 l 7876 4429 l 7868 4422 l
+ 7860 4417 l 7853 4413 l 7848 4411 l
+ 7845 4410 l gs col7 0.00 shd ef gr gs col0 s gr
+/Times-Roman ff 375.00 scf sf
+8715 7050 m
+gs 1 -1 sc (A) col0 sh gr
+/Times-Roman ff 300.00 scf sf
+9000 7155 m
+gs 1 -1 sc (cv) col0 sh gr
+/Times-Roman ff 375.00 scf sf
+5205 3975 m
+gs 1 -1 sc (y) col0 sh gr
+/Times-Roman ff 225.00 scf sf
+5370 4095 m
+gs 1 -1 sc (0) col0 sh gr
+/Times-Roman ff 225.00 scf sf
+5745 6195 m
+gs 1 -1 sc (1) col0 sh gr
+/Times-Roman ff 375.00 scf sf
+5625 6075 m
+gs 1 -1 sc (y) col0 sh gr
+/Times-Roman ff 225.00 scf sf
+7890 7050 m
+gs 1 -1 sc (2) col0 sh gr
+/Times-Roman ff 375.00 scf sf
+7740 6900 m
+gs 1 -1 sc (y) col0 sh gr
+/Times-Roman ff 225.00 scf sf
+7605 7050 m
+gs 1 -1 sc (2) col0 sh gr
+/Times-Roman ff 375.00 scf sf
+9990 4980 m
+gs 1 -1 sc (x) col0 sh gr
+/Times-Roman ff 225.00 scf sf
+10155 5085 m
+gs 1 -1 sc (3) col0 sh gr
+/Times-Roman ff 375.00 scf sf
+10380 4995 m
+gs 1 -1 sc (y) col0 sh gr
+/Times-Roman ff 225.00 scf sf
+10575 5070 m
+gs 1 -1 sc (3) col0 sh gr
+/Times-Roman ff 375.00 scf sf
+9630 2955 m
+gs 1 -1 sc (x) col0 sh gr
+/Times-Roman ff 225.00 scf sf
+9810 3045 m
+gs 1 -1 sc (4) col0 sh gr
+/Times-Roman ff 375.00 scf sf
+9990 2955 m
+gs 1 -1 sc (y) col0 sh gr
+/Times-Roman ff 225.00 scf sf
+10185 3030 m
+gs 1 -1 sc (4) col0 sh gr
+/Times-Roman ff 225.00 scf sf
+7185 2475 m
+gs 1 -1 sc (5) col0 sh gr
+/Times-Roman ff 375.00 scf sf
+7365 2355 m
+gs 1 -1 sc (y) col0 sh gr
+/Times-Roman ff 225.00 scf sf
+7545 2460 m
+gs 1 -1 sc (5) col0 sh gr
+/Times-Roman ff 225.00 scf sf
+5460 6195 m
+gs 1 -1 sc (1) col0 sh gr
+/Times-Roman ff 225.00 scf sf
+5025 4095 m
+gs 1 -1 sc (0) col0 sh gr
+/Times-Roman ff 375.00 scf sf
+4815 3975 m
+gs 1 -1 sc (x) col0 sh gr
+/Times-Roman ff 375.00 scf sf
+5295 6090 m
+gs 1 -1 sc (x) col0 sh gr
+/Times-Roman ff 375.00 scf sf
+7005 2370 m
+gs 1 -1 sc (x) col0 sh gr
+/Times-Roman ff 375.00 scf sf
+7440 6915 m
+gs 1 -1 sc (x) col0 sh gr
+% here ends figure;
+$F2psEnd
+rs
+end
+showpage
--- Makefile	2009-03-03 17:49:37.000000000 +0100
+++ Makefile	2009-03-03 17:50:34.000000000 +0100
@@ -69,6 +69,8 @@
 SwanPropvelS.$(EXTO) \
 SwanTranspAc.$(EXTO) \
 SwanTranspX.$(EXTO) \
+SwanDiffPar.$(EXTO) \
+SwanGSECorr.$(EXTO) \
 SwanInterpolatePoint.$(EXTO) \
 SwanInterpolateAc.$(EXTO) \
 SwanInterpolateOutput.$(EXTO) \
--- m_fileio.ftn90	2008-05-25 14:59:38.000000000 +0200
+++ m_fileio.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -83,7 +83,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -371,7 +371,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -427,7 +427,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- mod_xnl4v5.ftn90	2008-05-25 14:59:39.000000000 +0200
+++ mod_xnl4v5.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -469,7 +469,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -853,7 +853,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1075,7 +1075,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1204,7 +1204,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1461,7 +1461,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1676,7 +1676,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1804,7 +1804,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1911,7 +1911,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2219,7 +2219,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2705,7 +2705,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2858,7 +2858,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -3059,7 +3059,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -3542,7 +3542,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -3855,7 +3855,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -4254,7 +4254,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -4675,7 +4675,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -5145,7 +5145,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -5302,7 +5302,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -5579,7 +5579,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -5882,7 +5882,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -6158,7 +6158,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -6289,7 +6289,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -6526,7 +6526,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -6627,7 +6627,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -6992,7 +6992,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -7219,7 +7219,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -7369,7 +7369,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -7806,7 +7806,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -7908,7 +7908,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -8025,7 +8025,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -8137,7 +8137,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -8256,7 +8256,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -8360,7 +8360,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -8465,7 +8465,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -8904,7 +8904,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- ocpcre.ftn	2008-05-25 14:59:50.000000000 +0200
+++ ocpcre.ftn	2009-03-03 17:50:34.000000000 +0100
@@ -44,7 +44,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -132,7 +132,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -242,7 +242,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -391,7 +391,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -516,7 +516,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -712,7 +712,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -874,7 +874,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1070,7 +1070,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1264,7 +1264,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1409,7 +1409,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1792,7 +1792,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1901,7 +1901,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2002,7 +2002,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2110,7 +2110,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2222,7 +2222,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2336,7 +2336,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2425,7 +2425,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- ocpids.ftn	2008-05-25 14:59:50.000000000 +0200
+++ ocpids.ftn	2009-03-03 17:50:34.000000000 +0100
@@ -27,7 +27,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -339,7 +339,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -438,7 +438,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -562,7 +562,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- ocpmix.ftn	2009-03-03 17:49:37.000000000 +0100
+++ ocpmix.ftn	2009-03-03 17:50:34.000000000 +0100
@@ -41,7 +41,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -181,7 +181,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -364,7 +364,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -455,7 +455,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -552,7 +552,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -783,7 +783,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1153,7 +1153,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1265,7 +1265,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1411,7 +1411,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1509,7 +1509,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1608,7 +1608,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1949,7 +1949,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2076,7 +2076,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2268,7 +2268,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2364,7 +2364,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- platform.pl	2009-03-03 17:49:37.000000000 +0100
+++ platform.pl	2009-03-03 17:50:34.000000000 +0100
@@ -296,7 +296,7 @@
     print OUTFILE "F90_OMP = \n";
     print OUTFILE "F90_MPI = mpif90\n";
     print OUTFILE "FLAGS_OPT = -O\n";
-    print OUTFILE "FLAGS_MSC = -w\n";
+    print OUTFILE "FLAGS_MSC = -w -ffree-line-length-none\n";
     print OUTFILE "FLAGS_SER =\n";
     print OUTFILE "FLAGS_OMP =\n";
     print OUTFILE "FLAGS_MPI =\n";
--- serv_xnl4v5.ftn90	2008-05-25 14:59:40.000000000 +0200
+++ serv_xnl4v5.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -54,7 +54,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -130,7 +130,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -341,7 +341,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -451,7 +451,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -520,7 +520,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -738,7 +738,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -808,7 +808,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanBpntlist.ftn90	2009-03-03 17:49:37.000000000 +0100
+++ SwanBpntlist.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanCheckGrid.ftn90	2009-03-03 17:49:37.000000000 +0100
+++ SwanCheckGrid.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- swancom1.ftn	2009-03-03 17:49:37.000000000 +0100
+++ swancom1.ftn	2009-03-03 17:50:34.000000000 +0100
@@ -76,7 +76,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -434,7 +434,7 @@
 !                          2        (variable breaking coefficient
 !                                    according to Nelson (1994))
 !     PSURF(1)  = ALFA     1.0      (Battjes Janssen)
-!     PSURF(2)  = GAMMA    0.8      (Breaking criterium)
+!     PSURF(2)  = GAMMA    0.73     (breaking criterium)
 !
 !     PWCAP(1)  = ALFAWC   2.36e-5  (Empirical coefficient)
 !     PWCAP(2)  = ALFAPM   3.02E-3  (Alpha of Pierson Moskowitz frequency)
@@ -2070,7 +2070,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2371,10 +2371,11 @@
 !       Call SOURCE  to compute the source terms for each bin which fall
 !                    within a sweep:
 !                   1. Dissipation by wave-bottom effects
-!                   2. Dissipation due to whitecapping
-!                   3. Generation of wave energy by wind effects
-!                   4. Nonlinear wave-wave interactions (quadruplets)
-!                   5. Nonlinear wave-wave interactions (triads)
+!                   2. Dissipation due to surf breaking
+!                   3. Dissipation due to whitecapping
+!                   4. Generation of wave energy by wind effects
+!                   5. Nonlinear wave-wave interactions (quadruplets)
+!                   6. Nonlinear wave-wave interactions (triads)
 !       -----------------------------------------------------------------
 !       Call ACTION  calculate the derivatives in x,y,s,d space and store
 !                    the results in the corresponding arrays
@@ -3116,7 +3117,7 @@
 !
 !TIMG        CALL SWTSTA(122)                                                  40.23
         IF (PNUMS(20).LT.100.) THEN
-           IF (IWIND.NE.4) THEN                                           40.61
+           IF (IWIND.NE.4 .OR. NSTATC.NE.1) THEN                          40.96 40.61
 !             default limiter
               CALL PHILIM (AC2, SWMATR(1,1,JAOLD),                        40.31 40.16
      &                     CGO, KWAVE,                                    40.00
@@ -3189,7 +3190,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -3376,14 +3377,19 @@
          WRITE (PRINTF, *) 'Triads is off'
       ENDIF                                                               40.41
       IF (IBOT.EQ.2) THEN
-         WRITE(PRINTF,7005) PBOT(2),PBOT(1)
+         WRITE(PRINTF,7005) PBOT(2), PBOT(1)
  7005    FORMAT(' Collins (`72)        : CFW    ',E12.4,' CFC   ',E12.4)
       ELSEIF (IBOT.EQ.3) THEN
-         WRITE(PRINTF,7335) PBOT(4),PBOT(5)
+         WRITE(PRINTF,7335) PBOT(4), PBOT(5)
  7335    FORMAT(' Madsen et al. (`84)  : MF     ',E12.4,' KN    ',E12.4)
       ELSEIF (IBOT.EQ.1) THEN
          WRITE(PRINTF,7325) PBOT(3)
  7325    FORMAT(' JONSWAP (`73)        : GAMMA  ',E12.4)
+      ELSEIF (IBOT.EQ.4) THEN                                             41.04
+         WRITE(PRINTF,7326) PBOT(6), PBOT(7)
+ 7326    FORMAT(' JONSWAP (`73)        : GAMMA1 ',E12.4,' GAMMA2',E12.4)
+         WRITE(PRINTF,7327) PBOT(8), PBOT(9)
+ 7327    FORMAT('                      : DSPR1  ',E12.4,' DSPR2 ',E12.4)
       ELSE
          WRITE (PRINTF, *) 'Bottom friction is off'
       ENDIF
@@ -3427,9 +3433,29 @@
  7012    FORMAT(' Battjes&Janssen (`78): ALPHA  ',E12.4,
      &          ' GAMMA ',E12.4)
       ELSEIF (ISURF.EQ.2) THEN                                            970219
-         WRITE(PRINTF,7212) PSURF(1), PSURF(4), PSURF(5)
- 7212    FORMAT(' Nelson (`94): ALPHA  ',E12.4,
-     &          ' GAMmin ',E12.4, ' GAMmax ',E12.4)
+         WRITE(PRINTF,7212) PSURF(1), PSURF(4)
+ 7212    FORMAT(' Nelson (`94)         : ALPHA  ',E12.4,
+     &          ' GAMmin',E12.4)
+         WRITE(PRINTF,7213) PSURF(5)
+ 7213    FORMAT('                        GAMmax ',E12.4)
+      ELSEIF (ISURF.EQ.3) THEN                                            41.03
+         WRITE(PRINTF,7214) PSURF(1), PSURF(4)
+ 7214    FORMAT(' Ruessink et al (2003): ALPHA  ',E12.4,
+     &          ' A     ',E12.4)
+         WRITE(PRINTF,7215) PSURF(5)
+ 7215    FORMAT('                        B      ',E12.4)
+      ELSEIF (ISURF.EQ.4) THEN                                            41.03
+         WRITE(PRINTF,7216) PSURF(1), PSURF(2)
+ 7216    FORMAT(' Thornton&Guza (`83)  : ALPHA  ',E12.4,
+     &          ' GAMMA ',E12.4)
+         WRITE(PRINTF,7217) PSURF(4)
+ 7217    FORMAT('                        N      ',E12.4)
+      ELSEIF (ISURF.EQ.5) THEN                                            41.03
+         WRITE(PRINTF,7218) PSURF(1), PSURF(4)
+ 7218    FORMAT(' Surf biphase (2009)  : ALPHA  ',E12.4,
+     &          ' N     ',E12.4)
+         WRITE(PRINTF,7219) PSURF(5)
+ 7219    FORMAT('                        BREF   ',E12.4)
       ELSE
          WRITE (PRINTF, *) 'Surf breaking is off'
       ENDIF
@@ -3520,7 +3546,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -3988,7 +4014,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -4208,7 +4234,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -4602,7 +4628,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -4753,6 +4779,7 @@
 !
 !     AB2                : Sum of E_DSHKD2_DS_DD
 !     ACTOT_DSIG         : Integration term for calculating ACTOT
+!     BIPH               : parameterized biphase of the spectrum
 !     EDRKTOT_SWELL      : Swell-part of EDRKTOT
 !     EMAX               : Maximum energy according calculated HM
 !     ETOT_DRK_DSIG      : Integration term for calculating EDRKTOT
@@ -4766,8 +4793,10 @@
 !     FRINT_X_DDIR       : FRINTF * DDIR
 !     SINH_K_X_DEP_2     : SINH(KWAVE*DEP2)**2
 !     UB2                : Sum of ETOT_SIG2_DSHKD2_DSIG
+!     WH                 : fraction of breaking waves
 !
       REAL              :: AB2, EMAX, UB2
+      REAL              :: BIPH, WH
       REAL              :: FRINTF_X_DDIR
       REAL              :: BRCOEF   ! variable breaking coefficient (calc. in BRKPAR)  40.22
 !
@@ -4834,12 +4863,12 @@
 !
       ETOT = ETOT + ETOT_DSIG(MSC) * PWTAIL(6) / FRINTF
 !
-      IF (ISURF .EQ. 1) THEN
+      IF (ISURF .EQ. 1 .OR. ISURF.EQ.4 .OR. ISURF.EQ.5) THEN              41.03
         HM   = PSURF(2) * DEP2(KCGRD(1))                                  40.22
-      ELSEIF (ISURF .GE. 2) THEN
+      ELSEIF (ISURF .EQ. 2 .OR. ISURF.EQ.3) THEN                          41.03
 !        Calulate the correct breaking coefficient BRCOEF
          CALL BRKPAR (BRCOEF  ,SPCDIR(1,2), SPCDIR(1,3), AC2     ,
-     &                SIGPOW(:,1), DEP2, RDX, RDY            )            40.22
+     &                SIGPOW(:,1), DEP2, RDX, RDY, KWAVE         )        41.03 40.22
         HM   = BRCOEF * DEP2(KCGRD(1))                                    40.22
       ELSE
 !       breaking disabled, assign very high value to HM
@@ -4915,9 +4944,18 @@
           HS       = 4. * SQRT (ETOT)
         END IF
 !
-!     Calculate QB, when breaking is activated
+!       --- calculate Qb when Battjes/Janssen breaking is activated
+!
+        IF ( ISURF.GT.0 .AND. ISURF.LT.4 ) THEN
+           CALL FRABRE (HM, ETOT, QB(KCGRD(1)))
 !
-        IF (ISURF.GT.0) CALL FRABRE (HM, ETOT, QB(KCGRD(1)))
+!       --- calculate Qb when Thornton/Guza breaking is activated         41.03
+!
+        ELSEIF (ISURF.EQ.4) THEN
+           WH = (2.*SQRT(2.*ETOT)/HM)**PSURF(4)
+           WH = MIN(1.,WH)
+           QB(KCGRD(1)) = WH
+        ENDIF
 !
 !     Calculate the orbital velocity UBOT, orbital excursion ABRBOT and
 !     bottom wave period TMBOT                                            40.51
@@ -4945,6 +4983,19 @@
         URSELL(KCGRD(1)) = 0.                                             40.41
       END IF                                                              40.41
 !
+!     --- calculate Qb when biphase scaling is activated                  41.03
+!
+      IF (ISURF.EQ.5) THEN
+         BIPH = 0.5*PI*(TANH(PTRIAD(4)/URSELL(KCGRD(1)))-1.)
+         IF ( BIPH.LT.0. ) THEN
+            WH = (BIPH/PSURF(5))**PSURF(4)
+            WH = MIN(1.,WH)
+            QB(KCGRD(1)) = WH
+         ENDIF
+      ENDIF
+!
+      QB_LOC = QB(KCGRD(1))
+!
 !     *** test output ***
 !
       IF (TESTFL .AND. ITEST.GE.60) THEN
@@ -4952,13 +5003,11 @@
  901     FORMAT (' SINTGRL: ETOT Hs Sigma K Aorb', 5(1X, E11.4))
       END IF
 !
-!
 !     Set variables used outside the whitecapping scope
 !
       AC2TOT = ACTOT
       KMESPC = KM_WAM
       SMEBRK = SIGM01
-      QB_LOC = QB(KCGRD(1))
 !
       RETURN
 !
@@ -4997,7 +5046,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -5217,7 +5266,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -5463,7 +5512,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -5819,7 +5868,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -5990,10 +6039,8 @@
 !     PBOT(4)   = MF      -0.08     (Madsen et al. equation)
 !     PBOT(5)   = KN       0.05     (Madsen et al. bottom roughness)
 !
-!     PNUMS(*)  =
-!
 !     PSURF(1)  = ALFA     1.0      (Battjes & Janssen, 1978)
-!     PSURF(2)  = GAMMA    0.8      (Breaking criterium)
+!     PSURF(2)  = GAMMA    0.73     (breaking criterium)
 !
 !     PWCAP(1)  = ALFAWC   2.36e-5  (Emperical coefficient)
 !     PWCAP(2)  = ALFAPM   3.02E-3  (Alpha of Pierson Moskowitz frequency)
@@ -6020,7 +6067,7 @@
 !     -----------
 !     PWIND(14) 1D    alfa (which is tuned at 0.01)
 !     PWIND(15) 1D    Kappa ( 0.41)
-!     PWIND(16) 1D    Rho air (1.32)
+!     PWIND(16) 1D    Rho air (1.28)
 !     PWIND(17) 1D    Rho water (1025)
 !
 !  6. Local variables
@@ -6049,7 +6096,6 @@
 !   ------------------------------------------------------------
 !   If SSURF is on (ISURF > 0 ) then,
 !     Call SSURF to compute the source term due to wave breaking
-!       according to Battjes and Janssen (1978)
 !   ------------------------------------------------------------
 !   IF IWIND =1 OR IWIND =2 THEN
 !     Call WNDPAR (first or second generation mode of source terms
@@ -6185,7 +6231,7 @@
 !
 !       *** wave-bottom interactions ***
 !
-        CALL SBOT (ABRBOT   ,DEP2     ,ECOS     ,ESIN     ,
+        CALL SBOT (ABRBOT   ,DEP2     ,ECOS     ,ESIN     ,AC2      ,     41.04
      &             IMATDA   ,KWAVE    ,SPCSIG   ,UBOT     ,UX2      ,     30.72
      &             UY2      ,IDCMIN   ,IDCMAX   ,
      &             PLBTFR   ,ISSTOP   ,DISSC1   ,VARFR    ,FRCOEF   )     40.67
@@ -6195,10 +6241,7 @@
 !TIMG      CALL SWTSTA(131)                                                    40.23
       IF (ISURF .GE. 1) THEN
 !
-!       *** wave breaking with Kirby type formulation (f/fm)^2 ***
-!
-!
-!         *** wave breaking according to Battjes and Janssen (1978) ***
+!         *** calculate surf breaking source term (5 formulations) ***    41.03
 !
           CALL SSURF (ETOT    ,HM      ,
      &                QBLOC   ,SMEBRK  ,AC2     ,IMATRA  ,                30.81
@@ -6631,7 +6674,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -6782,7 +6825,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -6926,7 +6969,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -7115,7 +7158,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -7726,7 +7769,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -8236,7 +8279,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -8354,7 +8397,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -8716,7 +8759,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -9271,7 +9314,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- swancom2.ftn	2009-03-03 17:49:37.000000000 +0100
+++ swancom2.ftn	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !     SBOT    (Bottom friction)
 !     FRABRE  (Fraction of breaking waves)                                30.77
-!     SSURF   (Wave breaking: Battjes and Janssen (1978))
+!     SSURF   (Wave breaking: five formulations)
 !     SWCAP   (White capping: seven formulations)                         40.53
 !     BRKPAR  (wave breaking criterion according to Nelson (1987))
 !     CNTAIL  (contributions to the spectrum of the high frequency tail)
@@ -19,10 +19,10 @@
 !
 !****************************************************************
 !
-      SUBROUTINE SBOT (ABRBOT  ,DEP2    ,ECOS    ,ESIN    ,
+      SUBROUTINE SBOT (ABRBOT  ,DEP2    ,ECOS    ,ESIN    ,AC2     ,      41.04
      &                 IMATDA  ,KWAVE   ,SPCSIG  ,UBOT    ,UX2     ,      30.72
      &                 UY2     ,IDCMIN  ,IDCMAX  ,
-     &                 PLBTFR  ,ISSTOP  ,DISSC1  ,VARFR   , FRCOEF )      40.67
+     &                 PLBTFR  ,ISSTOP  ,DISSC1  ,VARFR   ,FRCOEF  )      40.67
 !
 !****************************************************************
 !
@@ -42,7 +42,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -67,6 +67,7 @@
 !     40.41: Marcel Zijlema
 !     40.61: Marcel Zijlema
 !     40.67: Nico Booij
+!     41.04: Marcel Zijlema
 !
 !  1. Updates
 !
@@ -76,6 +77,7 @@
 !     40.41, Oct. 04: common blocks replaced by modules, include files removed
 !     40.61, Sep. 06: introduce DISBOT variable for output purposes
 !     40.67, Jun. 07: more accurate computation fo dissipation terms
+!     41.04, Mar. 09: frequency-dependent JONSWAP formulation
 !
 !  2. Purpose
 !
@@ -196,8 +198,6 @@
 !     ---------
 !
 !     ABRBOT      Near bottom excursion amplitude
-!     DD          Spectral direction band width
-!     DS          Width of the frequency band
 !     FACB        an auxiliary factor contributing to bottom friction
 !     FW          Friction factor
 !     GRAV        Gravitational acceleration
@@ -210,6 +210,7 @@
 !     one and more dimensional arrays:
 !     ---------------------------------
 !
+!     AC2       2D    Action density
 !     DEP2      2D    Depth
 !     ESIN      1D    Sin per spectral direction (id)
 !     ECOS      1D    Cos per spectral direction (id)
@@ -255,13 +256,16 @@
 !
       INTEGER  ID     ,IS     ,ISSTOP
 !
-      REAL     XDUM   ,KD     ,SBOTEO ,CFBOT  ,FACB
-     &         CFW    ,FW     ,CURR   ,UC     ,ABRBOT ,
+      REAL     XDUM   ,KD     ,SBOTEO ,FACB  ,
+     &         CFW    ,FW     ,CURR   ,UC    ,ABRBOT,
      &         ADUM   ,CDUM   ,DDUM
+      REAL     CFBOT(MSC)
+      REAL     DSP    ,ETOT   ,EEX    ,EEY   ,EAD
 !
       LOGICAL  VARFR
 !
-      REAL     DEP2(MCGRD)               ,
+      REAL     AC2(MDC,MSC,MCGRD)        ,                                41.04
+     &         DEP2(MCGRD)               ,
      &         ECOS(MDC)                 ,
      &         ESIN(MDC)                 ,
      &         IMATDA(MDC,MSC)           ,
@@ -339,19 +343,51 @@
             FW = 0.3
           ENDIF
           CFBOT =  UBOT(KCGRD(1)) * FW / (SQRT(2.) * GRAV)
+        ELSE IF ( IBOT.EQ.4 ) THEN
+!
+!            *** Jonswap model with variable friction coefficient  ***
+!            *** as function of frequency-dependent directional    ***
+!                spreading (varies linearly between 0.038 - 0.067) ***
+!
+          DO IS = 1, MSC
+             ETOT = 0.
+             EEX  = 0.
+             EEY  = 0.
+             DO ID = 1, MDC
+                EAD  = SPCSIG(IS)*AC2(ID,IS,KCGRD(1))
+                ETOT = ETOT + EAD
+                EEX  = EEX  + EAD * ECOS(ID)
+                EEY  = EEY  + EAD * ESIN(ID)
+             ENDDO
+             IF ( ETOT.GT.0. ) THEN
+                XDUM = 1.-MIN(1.,SQRT(EEX*EEX+EEY*EEY)/ETOT)
+                DSP  = SQRT(2.*XDUM) *180./PI
+             ELSE
+                DSP  = 0.
+             ENDIF
+             IF ( DSP.LT.PBOT(8) ) THEN
+                CFBOT(IS) = PBOT(6)
+             ELSEIF ( DSP.GT.PBOT(9) ) THEN
+                CFBOT(IS) = PBOT(7)
+             ELSE
+                CFBOT(IS) = PBOT(6) + (PBOT(7)-PBOT(6))*(DSP-PBOT(8))/
+     &                                (PBOT(9)-PBOT(8))
+             ENDIF
+             CFBOT(IS) = CFBOT(IS) / GRAV**2
+          ENDDO
         ENDIF
 !
 !       *** test output ***
 !
         IF (TESTFL .AND. ITEST.GE.60) THEN
-          WRITE (PRTEST, 910) IBOT, KCGRD(1), DEP2(KCGRD(1)), CFBOT
+          WRITE (PRTEST, 910) IBOT, KCGRD(1), DEP2(KCGRD(1)), CFBOT(1)
  910      FORMAT (' SBOT :IBOT INDX DEP CFBOT:', 2I5, 2E12.4)
         END IF
 !
         DO 700 IS = 1, ISSTOP
-          KD      = KWAVE(IS,1) * DEP2(KCGRD(1))
+          KD = KWAVE(IS,1) * DEP2(KCGRD(1))
           IF ( KD .LT. 10. ) THEN
-            FACB = CFBOT * (SPCSIG(IS) / SINH(KD)) **2                    40.55 30.72
+            FACB = CFBOT(IS) * (SPCSIG(IS) / SINH(KD)) **2                41.04 40.55 30.72
 !
             DO 690 IDDUM = IDCMIN(IS) , IDCMAX(IS)
               ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
@@ -405,7 +441,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -600,7 +636,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -627,6 +663,7 @@
 !     40.41: Marcel Zijlema
 !     40.61: Marcel Zijlema
 !     40.67: Nico Booij
+!     41.03: Andre van der Westhuysen
 !
 !  1. Updates
 !
@@ -638,30 +675,36 @@
 !     40.41, Oct. 04: common blocks replaced by modules, include files removed
 !     40.61, Sep. 06: introduce DISSRF variable for output purposes
 !     40.67, Jun. 07: more accurate computation of dissipation terms
+!     41.03, Feb. 09: extension to alternative surf breaking formula's
 !
 !  2. Purpose
 !
-!     Computation of the source term due to wave breaking.
-!     White capping is not taken into account
+!     Computation of the source term due to wave breaking with one of three
+!     formulation:
+!                 1) Battjes and Janssen (1978)
+!                 2) Thornton and Guza (1983)
+!                 3) biphase scaling of Van der Westhuysen (2009)
+!
+!     Note: white capping is not taken into account
 !
 !  3. Method
 !
-!     The source term for surf breaking is implemented following
+!     Basically, the source term for surf breaking is implemented following
 !     the approach of Battjes/Janssen (1978) for the energy dissipation:
 !
 !             Alpha      -     2                  -   SMEBRK
 !     Dtot =  ----  Qb * f * Hm              with f = ------
 !              4                                      2 * Pi
 !
-!     Now the source term is:
+!     Now, the source term is:
 !
 !                      SIGMA * AC2(ID,IS,IX,IY)
-!     Sbr =  - Dtot *  ------------------------  =
+!     Sbr =    Dtot *  ------------------------  =
 !                              Etot
 !
 !
 !              Alpha * SMEBRK * Qb * Hm * Hm    SIGMA * AC2(ID,IS,IX,IY)
-!         =  - ------------------------------ * -------------------------
+!         =    ------------------------------ * -------------------------
 !                       8 * Pi                            Etot
 !
 !
@@ -704,16 +747,16 @@
 !     SbrD =  -----  =  ---- * N + WS
 !             d N       d N
 !
-!     Since BB and SMEBRK * N are proportional, we have
+!     Since BB and N are proportional, we have
 !
 !     d Sbr     d WS                   SMEBRK  (d Qb/ d BB) *BB - Qb
-!     -----  =  ---- * BB + WS = Alpha ------  --------------------- * BB + WS =
+!     -----  =  ---- * BB + WS = Alpha ------  --------------------- * BB + WS
 !     d N       d BB                     Pi           sqr(BB)
 !
 !
-!           SMEBRK d Qb
-!     Alpha ------ ----
-!            Pi    d BB
+!                    SMEBRK d Qb
+!            = Alpha ------ ----
+!                     Pi    d BB
 !
 !     With:
 !
@@ -730,6 +773,55 @@
 !          = ------------                  ;
 !            BB (BB - Qb)
 !
+!     Hence,
+!
+!     d Sbr       1 - Qb
+!     ----- = WS  -------
+!     d N         BB - Qb
+!
+!
+!     Alternatively, the source term for surf breaking is implemented following
+!     the approach of Thornton and Guza (1983) for energy dissipation:
+!
+!               3  -
+!              B * f                3                          -   SMEBRK
+!     Dtot =  ------- * INT(0,inf){H * W(H) * p(H)}dH     with f = ------
+!              4 * d                                               2 * Pi
+!
+!                    3
+!               3 * B * SMEBRK             3
+!          =  ------------------ * W * Hrms
+!              32 * sqrt(Pi) * d
+!
+!                           3                             3
+!          with INT(0,inf){H * p(H)}dH = 3/4*sqrt(Pi)*Hrms
+!
+!
+!          after Thornton and Guza (1983), Eqs. (24),(25)
+!
+!     For the weighting function W (or fraction of breaking waves), two choices can be made.
+!     The first one is after Thornton and Guza (1983):
+!
+!           Hrms  n
+!     W = (------)
+!           Hmax
+!
+!     and the second one is based on the biphase scaling of Van der Westhuysen (2009):
+!
+!               Beta    n
+!     W =   (----------)
+!             Beta_ref
+!
+!     where
+!                                                                                      1/2
+!              Pi   Pi         0.2               g                   sqr(Tm01) * 4*Etot
+!     Beta = - -- + -- * tanh (---)   and   Ur = ----------------- * ---------------------
+!               2    2          Ur               8*sqrt(2)*sqr(Pi)         sqr(d)
+!
+!
+!     For implementation details, see the Scientific/Technical documentation.
+!
+!
 !  4. Argument variables
 !
 !     AC2     input :   Action density array
@@ -767,6 +859,7 @@
 !  6. Local variables
 !     BB      Rate between the total energy and the energy
 !             according to the maximum wave height HM
+!     BIPH    parameterized biphase of the spectrum
 !     DIS0    Dummy variable
 !     ID      Counter for directional steps
 !     IDDUM   Counter
@@ -775,11 +868,12 @@
 !     SURFA0  Coefficient for old source term in matrix equation
 !             (i.e. SURFA0 * Nold = right hand side of matrix equation)
 !     SURFA1  Coefficient for new source term in matrix equation
+!     WH      fraction of breaking waves
 !     WS      Wavebreaking source term coefficient = DTOT/ETOT
 !     SbrD    Derivative of source term for surf breaking (Sbr) to action density
 !
-      INTEGER ID,      IDDUM,   IENT,   IS
-!
+      INTEGER          ID,      IDDUM,   IENT,   IS
+      REAL             BIPH,    DEPLOC,  URSLOC, WH
       DOUBLE PRECISION BB,      DIS0,    SbrD,
      &                 SURFA0,  SURFA1,  WS
 !
@@ -829,25 +923,86 @@
 !
 !     ALFA = PSURF(1)   <default = 1.0>
 !
-      BB = 8. * DBLE(ETOT) / ( DBLE(HM)**2 )                              30.82
-      SURFA0 = 0.
-      SURFA1 = 0.
-      IF (REAL(BB) .GT. 0. .AND.                                          30.82
-     &    REAL(ABS(BB - DBLE(QB))) .GT. 0.) THEN                          30.82
-        IF ( BB .LT. 1. ) THEN
-          WS  = ( DBLE(PSURF(1)) / DBLE(PI)) *                            30.82
-     &            DBLE(QB) * DBLE(SMEBRK) / BB                            30.82
-          SbrD = WS * (1. - DBLE(QB)) / (BB - DBLE(QB))                   30.82 40.00
-        ELSE
-          WS  = ( DBLE(PSURF(1)) / DBLE(PI)) * DBLE(SMEBRK)               30.82
-          SbrD = 0.
-        END IF
-        SURFA0 = SbrD
-        SURFA1 = WS + SbrD
-      ELSE
-        SURFA0 = 0.
-        SURFA1 = 0.
-      END IF
+      BB = 8D0 * DBLE(ETOT) / ( DBLE(HM)**2 )                             41.03 30.82
+      SURFA0 = 0D0
+      SURFA1 = 0D0
+!
+      IF ( ISURF.LE.3 ) THEN                                              41.03
+!
+!        --- Battjes and Janssen (1978)
+!
+         IF (REAL(BB) .GT. 0. .AND.                                       30.82
+     &       REAL(ABS(BB - DBLE(QB))) .GT. 0.) THEN                       30.82
+            IF ( BB .LT. 1D0 ) THEN                                       41.03
+               WS  = ( DBLE(PSURF(1)) / DBLE(PI)) *                       30.82
+     &                 DBLE(QB) * DBLE(SMEBRK) / BB                       30.82
+               SbrD = WS * (1D0 - DBLE(QB)) / (BB - DBLE(QB))             41.03 30.82 40.00
+            ELSE
+               WS  = ( DBLE(PSURF(1)) / DBLE(PI)) * DBLE(SMEBRK)          30.82
+               SbrD = 0D0
+            END IF
+            SURFA0 = SbrD
+            SURFA1 = WS + SbrD
+         ELSE
+            SURFA0 = 0D0
+            SURFA1 = 0D0
+         ENDIF
+!
+      ELSEIF ( ISURF.EQ.4 ) THEN                                          41.03
+!
+!        --- Thornton and Guza (1983)
+!
+         IF ( BB.GT.0D0 ) THEN
+            IF ( BB.LT.1D0 ) THEN
+               WS = 75D-2*DBLE(PSURF(2))*DBLE(PSURF(1))**3*DBLE(SMEBRK)*
+     &              BB**(0.5*(PSURF(4)+1))/DBLE(SQRT(PI))
+               SbrD = 5D-1*DBLE(3.+PSURF(4))*WS
+            ELSE
+               WS = 75D-2*DBLE(PSURF(2))*DBLE(PSURF(1))**3*DBLE(SMEBRK)/
+     &              DBLE(SQRT(PI))
+               SbrD = WS
+            ENDIF
+            SURFA0 = SbrD - WS
+            SURFA1 = SbrD
+         ELSE
+            SURFA0 = 0D0
+            SURFA1 = 0D0
+         ENDIF
+!
+      ELSEIF ( ISURF.EQ.5 ) THEN                                          41.03
+!
+!        --- biphase scaling of Van der Westhuysen (2009)
+!
+         DEPLOC = HM/PSURF(2)
+         URSLOC = (GRAV*2.*SQRT(ETOT))/(SQRT(2.)*SMEBRK**2*DEPLOC**2)
+         BIPH = 0.5*PI*(TANH(PTRIAD(4)/URSLOC)-1.)
+         IF ( BIPH.LT.0. ) THEN
+            WH = (BIPH/PSURF(5))**PSURF(4)
+            WH = MIN(1.,WH)
+!
+            IF ( BB.GT.0D0 ) THEN
+               IF ( BB.LT.1D0 ) THEN
+                  WS = 75D-2*DBLE(WH)*DBLE(PSURF(2))*DBLE(PSURF(1))**3*
+     &                 DBLE(SMEBRK)*SQRT(BB)/DBLE(SQRT(PI))
+                  SbrD = (15D-1-25D-2*DBLE(PSURF(4))*
+     &                                          (1D0+5D0/DBLE(BIPH)))*WS
+               ELSE
+                  WS = 75D-2*DBLE(PSURF(2))*DBLE(PSURF(1))**3*
+     &                 DBLE(SMEBRK)/DBLE(SQRT(PI))
+                  SbrD = WS
+               ENDIF
+               SURFA0 = SbrD - WS
+               SURFA1 = SbrD
+            ELSE
+               SURFA0 = 0D0
+               SURFA1 = 0D0
+            ENDIF
+         ELSE
+            SURFA0 = 0D0
+            SURFA1 = 0D0
+         ENDIF
+!
+      ENDIF
 !
 !     *** store the results for surf wave breaking  ***
 !     *** in the matrices IMATDA and IMATRA         ***
@@ -908,7 +1063,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1461,7 +1616,8 @@
 !****************************************************************
 !
       SUBROUTINE BRKPAR (BRCOEF  ,ECOS    ,ESIN    ,AC2     ,             40.22
-     &                   SPCSIG  ,DEP2    ,RDX     ,RDY     )             30.72
+     &                   SPCSIG  ,DEP2    ,RDX     ,RDY     ,             41.03 30.72
+     &                   KWAVE                              )             41.03
 !
 !****************************************************************
 !
@@ -1483,7 +1639,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1507,6 +1663,7 @@
 !     40.02: IJsbrand Haagsma
 !     40.22: Nico Booij
 !     40.41: Marcel Zijlema
+!     41.03: Andre van der Westhuysen
 !
 !  1. Updates
 !
@@ -1517,6 +1674,8 @@
 !     40.08, Mar. 03: Dimensioning of RDX, RDX changed to be consistent
 !                     with other subroutines
 !     40.41, Oct. 04: common blocks replaced by modules, include files removed
+!     41.03, Feb. 09: extension with spatially varying breaker parameter
+!                     according to Ruessink et al (2003)
 !
 !  2. Purpose
 !
@@ -1541,9 +1700,18 @@
 !     slope of 0.01) is imposed on the model of Nelson.
 !
 !     For negative bottom slopes (not considered by Nelson) a value
-!     op 0.73 is imposed (which is the average value in Table 2 of
+!     of 0.73 is imposed (which is the average value in Table 2 of
 !     Battjes and Janssen (1978).
 !
+!
+!     Alternatively, the breaker index can be computed according to
+!     Ruessink et al. (2003):
+!
+!     Hm / d = 0.76 kp * d + 0.29
+!
+!     with kp the peak wave number
+!
+!
 !  4. Argument variables
 !
       REAL, INTENT(OUT) :: BRCOEF    ! variable breaker coefficient       40.22
@@ -1556,6 +1724,7 @@
       REAL, INTENT(IN)  :: DEP2(MCGRD)         ! depths at grid points    40.22
 !     RDX, RDY:  coefficients to obtain spatial derivatives               40.22
       REAL, INTENT(IN)  :: RDX(10), RDY(10)                               40.08
+      REAL, INTENT(IN)  :: KWAVE(MSC,MICMAX)                              41.03
 !
 !        INTEGERS :
 !        ----------
@@ -1622,62 +1791,97 @@
 !************************************************************************
 !
       INTEGER :: ID    ,IS      ! counters                                40.22
+      INTEGER :: ISIGM                                                    41.03
 !
 !
       REAL  :: ETOTS ,EEX   ,EEY   ,
      &         EAD   ,SIGMA1,COSDIR,SINDIR,DDDX  ,                        40.22
      &         DDDY  ,DDDS  ,DETOT                                        40.22
+      REAL  :: EMAX, ETD, KP, KPD                                         41.03
 !
       INTEGER, SAVE :: IENT=0
       IF (LTRACE) CALL STRACE (IENT,'BRKPAR')
 !
-!     *** determine the average wave direction ***
+      IF ( ISURF.EQ.2 ) THEN                                              41.03
 !
-      EEX   = 0.
-      EEY   = 0.
-      ETOTS = 0.
-      DO ID = 1, MDC
-        EAD = 0.
-        DO IS = 1, MSC
-          SIGMA1 = SPCSIG(IS)                                             30.72
-          DETOT  = SIGMA1**2 * AC2(ID,IS,KCGRD(1))
-          EAD    = EAD + DETOT
-        ENDDO
-        ETOTS = ETOTS + EAD
-        EEX   = EEX + EAD * ECOS(ID)
-        EEY   = EEY + EAD * ESIN(ID)
-      ENDDO
-!
-      IF ( ETOTS .GT. 0.) THEN
-        COSDIR = EEX / ETOTS
-        SINDIR = EEY / ETOTS
-      ELSE
-        COSDIR = 1.
-        SINDIR = 0.
-      ENDIF
+!        *** determine the average wave direction ***
 !
-!     *** Determine bottom slope in average wave propagation direction ***
+         EEX   = 0.
+         EEY   = 0.
+         ETOTS = 0.
+         DO ID = 1, MDC
+            EAD = 0.
+            DO IS = 1, MSC
+               SIGMA1 = SPCSIG(IS)                                        30.72
+               DETOT  = SIGMA1**2 * AC2(ID,IS,KCGRD(1))
+               EAD    = EAD + DETOT
+            ENDDO
+            ETOTS = ETOTS + EAD
+            EEX   = EEX + EAD * ECOS(ID)
+            EEY   = EEY + EAD * ESIN(ID)
+         ENDDO
+!
+         IF ( ETOTS .GT. 0.) THEN
+            COSDIR = EEX / ETOTS
+            SINDIR = EEY / ETOTS
+         ELSE
+            COSDIR = 1.
+            SINDIR = 0.
+         ENDIF
+!
+!        *** Determine bottom slope in average wave propagation direction ***
+!
+         DDDX =  RDX(1) * (DEP2(KCGRD(1)) - DEP2(KCGRD(2)))
+     &         + RDX(2) * (DEP2(KCGRD(1)) - DEP2(KCGRD(3)))
+         DDDY =  RDY(1) * (DEP2(KCGRD(1)) - DEP2(KCGRD(2)))
+     &         + RDY(2) * (DEP2(KCGRD(1)) - DEP2(KCGRD(3)))
+!
+         DDDS = -1. * ( DDDX * COSDIR + DDDY * SINDIR )
+!
+!        *** calculate breaking coefficient according to Nelson (1987) ***
+!
+         IF ( DDDS .GE. 0. ) THEN
+            DDDS   = MAX ( 1.E-6 , DDDS)
+            BRCOEF = PSURF(4) + PSURF(7) * EXP ( -PSURF(8) / DDDS )       40.22
+            BRCOEF = MIN ( PSURF(5) , BRCOEF )                            40.22
+         ELSE
+            BRCOEF = PSURF(6)                                             40.22
+         ENDIF
+!
+      ELSE IF ( ISURF.EQ.3 ) THEN                                         41.03
 !
-      DDDX =  RDX(1) * (DEP2(KCGRD(1)) - DEP2(KCGRD(2)))
-     &      + RDX(2) * (DEP2(KCGRD(1)) - DEP2(KCGRD(3)))
-      DDDY =  RDY(1) * (DEP2(KCGRD(1)) - DEP2(KCGRD(2)))
-     &      + RDY(2) * (DEP2(KCGRD(1)) - DEP2(KCGRD(3)))
+!       calculate breaker index according to Ruessink et al (2003)
+!
+        EMAX = 0.
+        ISIGM = -1
+        DO IS = 1, MSC
+           ETD = 0.
+           DO ID = 1, MDC
+              ETD = ETD + SPCSIG(IS)*AC2(ID,IS,KCGRD(1))*DDIR
+           ENDDO
+           IF (ETD.GT.EMAX) THEN
+              EMAX  = ETD
+              ISIGM = IS
+           ENDIF
+        ENDDO
+        IF (ISIGM.GT.0) THEN
+           KP = KWAVE(ISIGM,1)
+        ELSE
+           KP = 0.
+        ENDIF
 !
-      DDDS = -1. * ( DDDX * COSDIR + DDDY * SINDIR )
+        KPD = KP*DEP2(KCGRD(1))
 !
-!     *** calculate breaking coefficient according to Nelson (1987) ***
+        IF ( KPD.LT.0.) THEN
+           BRCOEF = 0.73
+        ELSE
+           BRCOEF = PSURF(4)*KPD + PSURF(5)
+           BRCOEF = MIN( 1.2, BRCOEF)
+           BRCOEF = MAX( 0.3, BRCOEF)
+        ENDIF
 !
-      IF ( DDDS .GE. 0. ) THEN
-        DDDS   = MAX ( 1.E-6 , DDDS)
-        BRCOEF = PSURF(4) + PSURF(7) * EXP ( -PSURF(8) / DDDS )           40.22
-        BRCOEF = MIN ( PSURF(5) , BRCOEF )                                40.22
-      ELSE
-        BRCOEF = PSURF(6)                                                 40.22
       ENDIF
 !
-!      PSURF(2) = BRKVAR                              deleted             40.22
-!      PSURF(2) is no longer used to transmit br. coefficient             40.22
-!
 !     *** test output ***
 !
       IF ( TESTFL .AND. ITEST .GE. 40 ) THEN
@@ -1720,7 +1924,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- swancom3.ftn	2009-03-03 17:49:37.000000000 +0100
+++ swancom3.ftn	2009-03-03 17:50:34.000000000 +0100
@@ -53,7 +53,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -562,7 +562,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -927,7 +927,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1160,7 +1160,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1337,7 +1337,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1576,7 +1576,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1809,7 +1809,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2354,7 +2354,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- swancom4.ftn	2009-03-03 17:49:37.000000000 +0100
+++ swancom4.ftn	2009-03-03 17:50:34.000000000 +0100
@@ -70,7 +70,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -512,7 +512,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -705,7 +705,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1203,7 +1203,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1605,7 +1605,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1999,7 +1999,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2405,7 +2405,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2693,7 +2693,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2820,7 +2820,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -3555,7 +3555,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -3764,7 +3764,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -3896,7 +3896,7 @@
 !
 !     AUX1  :     auxiliary real
 !     AUX2  :     auxiliary real
-!     BIPH  :     parameterized bi-phase of the spectrum
+!     BIPH  :     parameterized biphase of the spectrum
 !     C0    :     phase velocity at central bin
 !     CM    :     phase velocity at interacting bin
 !     DEP   :     water depth
--- swancom5.ftn	2009-03-03 17:49:37.000000000 +0100
+++ swancom5.ftn	2009-03-03 17:50:34.000000000 +0100
@@ -55,7 +55,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -76,11 +76,13 @@
 !  0. Authors
 !
 !     40.41: Marcel Zijlema
+!     40.98: Marcel Zijlema
 !
 !  1. Updates
 !
 !     40.41, Sep. 04: New subroutine (taken from routine SWPSEL)
 !     40.41, Oct. 04: common blocks replaced by modules, include files removed
+!     40.98, Feb. 09: SORDUP scheme is made consistent
 !
 !  2. Purpose
 !
@@ -109,25 +111,17 @@
 !  6. Local variables
 !
 !     DET   :     determinant or volume of cell
-!     DET2  :     another determinant (for SORDUP scheme)
 !     DX1   :     first component of covariant base vector a_(1)
 !     DX2   :     second component of covariant base vector a_(1)
-!     DX3   :     =DX1 (to be needed in SORDUP scheme)
-!     DX4   :     =DX2 (to be needed in SORDUP scheme)
 !     DY1   :     first component of covariant base vector a_(2)
 !     DY2   :     second component of covariant base vector a_(2)
-!     DY3   :     =DY1 (to be needed in SORDUP scheme)
-!     DY4   :     =DY2 (to be needed in SORDUP scheme)
 !     IC    :     counter
 !     IENT  :     number of entries
 !     IXY   :     counter
 !     VIRT  :     indicates virtual point for 1D mode
 !
       INTEGER IC, IENT, IXY
-      REAL    VIRT,
-     &        DET, DET2,
-     &        DX1, DX2, DX3, DX4,
-     &        DY1, DY2, DY3, DY4
+      REAL    VIRT, DET, DX1, DX2, DY1, DY2
 !
 !  7. Common blocks used
 !
@@ -151,26 +145,19 @@
         IF (SWPDIR.EQ.1 .OR. SWPDIR.EQ.4) THEN
           DX1 = DX * COSPC                                                33.09
           DY1 = DX * SINPC                                                33.09
-          DX3 = DX1 ! 40.08 NOTE 3 and 4 won't be used unless PROPSL.eq.2
-          DY3 = DY1 ! 40.08 NOTE 3 and 4 won't be used unless PROPSL.eq.2
         ELSE
           DX1 = -DX * COSPC                                               33.09
           DY1 = -DX * SINPC                                               33.09
-          DX3 = DX1 ! 40.08 NOTE 3 and 4 won't be used unless PROPSL.eq.2
-          DY3 = DY1 ! 40.08 NOTE 3 and 4 won't be used unless PROPSL.eq.2
         ENDIF
-        IF (ONED) THEN                                                    33.09
+        IF ( ONED ) THEN                                                  33.09
+!         *** Inclusion of virtual point ***                              33.09
           VIRT = 1.E6                                                     33.09
           IF ( SWPDIR .EQ. 1 .OR. SWPDIR .EQ. 3 ) THEN                    33.09
             DX2 = -VIRT * DY1                                             33.09
             DY2 =  VIRT * DX1                                             33.09
-            DX4 = DX2 ! 40.08 NOTE 3 and 4 won't be used unless PROPSL.eq.2
-            DY4 = DY2 ! 40.08 NOTE 3 and 4 won't be used unless PROPSL.eq.2
           ELSE                                                            33.09
             DX2 =  VIRT * DY1                                             33.09
             DY2 = -VIRT * DX1                                             33.09
-            DX4 = DX2 ! 40.08 NOTE 3 and 4 won't be used unless PROPSL.eq.2
-            DY4 = DY2 ! 40.08 NOTE 3 and 4 won't be used unless PROPSL.eq.2
           ENDIF                                                           33.09
         ELSE                                                              33.09
           IF (SWPDIR.LE.2) THEN                                           40.13
@@ -180,44 +167,23 @@
             DX2 =   DY * SINPC                                            40.13
             DY2 = - DY * COSPC                                            40.13
           ENDIF                                                           40.13
-          DX4 = DX2 ! 40.08 NOTE 3 and 4 won't be used unless PROPSL.eq.2
-          DY4 = DY2 ! 40.08 NOTE 3 and 4 won't be used unless PROPSL.eq.2
         ENDIF                                                             33.09
       ELSE                                                                33.09
-
         DX1 = XCGRID(IXCGRD(1),IYCGRD(1)) - XCGRID(IXCGRD(2),IYCGRD(2))
         DY1 = YCGRID(IXCGRD(1),IYCGRD(1)) - YCGRID(IXCGRD(2),IYCGRD(2))
-        IF(PROPSL.EQ.2)THEN                                               40.08
-          DX3 = (XCGRID(IXCGRD(1),IYCGRD(1))                              40.08
-     &         - XCGRID(IXCGRD(4),IYCGRD(4)))/2.0                         40.08
-          DY3 = (YCGRID(IXCGRD(1),IYCGRD(1))                              40.08
-     &         - YCGRID(IXCGRD(4),IYCGRD(4)))/2.0                         40.08
-        ENDIF                                                             40.08
         IF  ( ONED ) THEN                                                 32.02
-
 !         *** Inclusion of virtual point ***                              32.02
-
           VIRT = 1.E6                                                     32.02
           IF ( SWPDIR .EQ. 1 .OR. SWPDIR .EQ. 3 ) THEN                    32.02
             DX2 = -VIRT * DY1                                             32.02
             DY2 =  VIRT * DX1                                             32.02
-            DX4 = DX2 ! 40.08 NOTE 3 and 4 won't be used unless PROPSL.eq.2
-            DY4 = DY2 ! 40.08 NOTE 3 and 4 won't be used unless PROPSL.eq.2
           ELSE IF ( SWPDIR .EQ. 2 .OR. SWPDIR .EQ. 4 ) THEN               32.02
             DX2 =  VIRT * DY1                                             32.02
             DY2 = -VIRT * DX1                                             32.02
-            DX4 = DX2 ! 40.08 NOTE 3 and 4 won't be used unless PROPSL.eq.2
-            DY4 = DY2 ! 40.08 NOTE 3 and 4 won't be used unless PROPSL.eq.2
           ENDIF                                                           32.02
         ELSE                                                              32.02
           DX2 = XCGRID(IXCGRD(1),IYCGRD(1))-XCGRID(IXCGRD(3),IYCGRD(3))
           DY2 = YCGRID(IXCGRD(1),IYCGRD(1))-YCGRID(IXCGRD(3),IYCGRD(3))
-          IF(PROPSL.EQ.2)THEN                                             40.08
-             DX4 = (XCGRID(IXCGRD(1),IYCGRD(1))                           40.08
-     &            - XCGRID(IXCGRD(5),IYCGRD(5)))/2.0                      40.08
-             DY4 = (YCGRID(IXCGRD(1),IYCGRD(1))                           40.08
-     &            - YCGRID(IXCGRD(5),IYCGRD(5)))/2.0                      40.08
-          ENDIF                                                           40.08
         ENDIF                                                             32.02
       ENDIF                                                               33.09
 !
@@ -227,19 +193,6 @@
       RDX(2) = -DY1/DET
       RDY(2) =  DX1/DET
 !
-      IF(PROPSL.EQ.2)THEN                                                 40.08
-         DET2   =  DY4*DX3 - DY3*DX4                                      40.08
-         RDX(3) =  DY4/DET2                                               40.08
-         RDY(3) = -DX4/DET2                                               40.08
-         RDX(4) = -DY3/DET2                                               40.08
-         RDY(4) =  DX3/DET2                                               40.08
-      ELSE
-         RDX(3) = 0.                                                      40.41
-         RDY(3) = 0.                                                      40.41
-         RDX(4) = 0.                                                      40.41
-         RDY(4) = 0.                                                      40.41
-      ENDIF                                                               40.08
-!
 !     in case of spherical coordinates determine cos of latitude
 !     note: latitude is in degrees
 !
@@ -248,8 +201,7 @@
           COSLAT(IC) =
      &    COS(DEGRAD*(YCGRID(IXCGRD(IC),IYCGRD(IC))+YOFFS))               33.09
         ENDDO
-!       --- "2" changed to "4". 3 and 4 won't be used unless PROPSL.eq.2  40.08
-        DO IXY = 1, 4                                                     40.08
+        DO IXY = 1, 2                                                     40.08
           RDY(IXY) = RDY(IXY) / LENDEG                                    33.09
           RDX(IXY) = RDX(IXY) / (COSLAT(1) * LENDEG)                      33.09
         ENDDO
@@ -309,7 +261,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -925,7 +877,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1212,7 +1164,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1393,7 +1345,6 @@
      &         KCG1  ,KCG2  ,KCG3  ,ISWEEP                                30.80
       INTEGER  IX1, IY1, IX2, IY2, IX3, IY3                               40.03
       INTEGER :: ISWP                                                     40.02
-      INTEGER IC, KCGI                                                    40.21
 !
 !     logical local variables
 !
@@ -1647,7 +1598,7 @@
 !
           IF (VALSWP) THEN
             KCG2 = KGRPNT(IX2,IY2)
-            IF (KCG2.LE.1 .OR. DEP2(KCG2).LE.DEPMIN) VALSWP = .FALSE.         40.14
+            IF (KCG2.LE.1 .OR. DEP2(KCG2).LE.DEPMIN) VALSWP = .FALSE.     40.14
             IF (KREPTX.GT.0) THEN                                         33.09
               DX2 = DX * COSPC                                            33.09
               DY2 = -DX * SINPC                                           33.09
@@ -1953,7 +1904,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2127,7 +2078,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2434,7 +2385,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2459,6 +2410,7 @@
 !        40.08: Erick Rogers
 !        40.41: Marcel Zijlema
 !        40.85: Marcel Zijlema
+!        40.98: Marcel Zijlema
 !
 !     1. UPDATE
 !
@@ -2469,16 +2421,17 @@
 !                          be a noticeable improvement. (I have not seen an example of a
 !                          case where the original SORDUP does poorly relative to BSBT,
 !                          so this is a speculative improvement).
-!                          Remove option for controllable 1st order diffusion ("MUXY",
+!                          Remove option for controllable 1st order diffusion ("XYMU",
 !                          "THETAK", etc.)
 !        40.41, Oct. 04: common blocks replaced by modules, include files removed
 !        40.85, Aug. 08: store xy-propagation for output purposes
+!        40.98, Feb. 09: SORDUP scheme is made consistent
 !
 !     2. PURPOSE
 !
 !        Purpose is to compute the space derivative of action transport using
 !        the SORDUP scheme.
-!        This is for stationary calcs only (no time derivative).
+!        This is for stationary runs only (no time derivative).
 !        The scheme is 2nd order accurate.                                40.08
 !        The scheme reduces to the "best" approximation of                40.08
 !             d/dx which can be determined using Taylor Series for the    40.08
@@ -2500,7 +2453,7 @@
 !        KCGRD(4) :   IX-2,IY
 !        KCGRD(5) :   IX  ,IY-2
 !
-!        The scheme is:                                                   40.08
+!        The scheme is:
 !
 !        @[CAX AC2]
 !        --------- =
@@ -2578,8 +2531,8 @@
 !     Changed ICMAX to MICMAX, since MICMAX doesn't vary over gridpoint   40.22
       REAL  :: CAX(MDC,MSC,MICMAX) ,CAY(MDC,MSC,MICMAX)                   40.22
       REAL  :: IMATRA(MDC,MSC)    ,IMATDA(MDC,MSC)            ,
-     &         RDX(10)            ,RDY(10)                                40.08
-     &         ,XMU(6),YMU(6)                                             40.08 33.10
+     &         RDX(10)            ,RDY(10)                    ,           40.08
+     &         XMU(5)             ,YMU(5)                                 40.08 33.10
       REAL  :: TRAC0(MDC,MSC,MTRNP)                                       40.85
       REAL  :: TRAC1(MDC,MSC,MTRNP)                                       40.85
 
@@ -2616,18 +2569,10 @@
 !         find Courant number values: XMU, YMU                            40.08
 !         depending on relative size of XMU and YMU, XNUM is true or      40.08
 !         false because of RDX and RDY, XMU and YMU are always positive   40.08
-          DO IXY=1,3                                                      40.08
-          XMU(IXY) = RDX(1)*CAX(ID,IS,IXY) + RDY(1)*CAY(ID,IS,IXY)        40.08
-          YMU(IXY) = RDX(2)*CAX(ID,IS,IXY) + RDY(2)*CAY(ID,IS,IXY)        40.08
+          DO IXY=1,5                                                      40.08
+             XMU(IXY) = RDX(1)*CAX(ID,IS,IXY) + RDY(1)*CAY(ID,IS,IXY)     40.08
+             YMU(IXY) = RDX(2)*CAX(ID,IS,IXY) + RDY(2)*CAY(ID,IS,IXY)     40.08
           END DO                                                          40.08
-          DO IXY=4,5                                                      40.08
-          XMU(IXY) = RDX(3)*CAX(ID,IS,IXY) + RDY(3)*CAY(ID,IS,IXY)        40.08
-          YMU(IXY) = RDX(4)*CAX(ID,IS,IXY) + RDY(4)*CAY(ID,IS,IXY)        40.08
-          END DO                                                          40.08
-!         XMU(6) is an alternate form of XMU(1). It uses the dx over 2    40.08
-!         grid spaces instead of 1. Similar for YMU(6)                    40.08
-          XMU(6) = RDX(3)*CAX(ID,IS,1) + RDY(3)*CAY(ID,IS,1) ! 40.08
-          YMU(6) = RDX(4)*CAX(ID,IS,1) + RDY(4)*CAY(ID,IS,1) ! 40.08
           IF(YMU(1).GT.XMU(1))THEN                                        33.10
 !            propagation mainly from grid point 3
              XNUM=.TRUE.                                                  33.10
@@ -2641,17 +2586,16 @@
           IF(XNUM)THEN                                                    33.10
 
 !           diagonal  FXY1
-            FXY1 = 2.0*XMU(1)-0.5*XMU(6)                                  40.08
-     &           + 2.0*YMU(1)-0.5*YMU(6)                                  40.08
+            FXY1 = 1.5*XMU(1) + 1.5*YMU(1)                                40.98 40.08
 
             IF (KSPHER.EQ.0) THEN                                         33.08
 !             Cartesian coordinates
 
 !             the known, rhs part FXY2                                    33.08
-              FXY2 = AC2(ID,IS,IND2) * XMU(2)*2.0                         33.10,40.08
-     &              -AC2(ID,IS,IND4) * 0.5*XMU(4)                         33.10,40.08
-     &              +AC2(ID,IS,IND3) * YMU(3)*2.0                         33.10,40.08
-     &              -AC2(ID,IS,IND5) * 0.5*YMU(5)                         33.10,40.08
+              FXY2 = AC2(ID,IS,IND2) * 2.0*XMU(2)                         40.08 33.10
+     &              -AC2(ID,IS,IND4) * 0.5*XMU(4)                         40.08 33.10
+     &              +AC2(ID,IS,IND3) * 2.0*YMU(3)                         40.08 33.10
+     &              -AC2(ID,IS,IND5) * 0.5*YMU(5)                         40.08 33.10
 
             ELSE                                                          33.10
 !             Spherical coordinates
@@ -2659,45 +2603,44 @@
 !             the known, rhs part FXY2                                    33.08
 
               FXY2 =
-     &        AC2(ID,IS,IND2) * CAX(ID,IS,2) * RDX(1) * 2.0 !33.10,40.08
-     &       -AC2(ID,IS,IND4) * CAX(ID,IS,4) * RDX(3) * 0.5 !33.10,40.08
-     &       +AC2(ID,IS,IND3) * CAX(ID,IS,3) * RDX(2) * 2.0 !33.10,40.08
-     &       -AC2(ID,IS,IND5) * CAX(ID,IS,5) * RDX(4) * 0.5 !33.10,40.08
-     &      +(AC2(ID,IS,IND2) * CAY(ID,IS,2) * RDY(1) * COSLAT(2) * 2.0   33.10,40.08
-     &       -AC2(ID,IS,IND4) * CAY(ID,IS,4) * RDY(3) * COSLAT(4) * 0.5   33.10,40.08
-     &       +AC2(ID,IS,IND3) * CAY(ID,IS,3) * RDY(2) * COSLAT(3) * 2.0   33.10,40.08
-     &       -AC2(ID,IS,IND5) * CAY(ID,IS,5) * RDY(4) * COSLAT(5) * 0.5   33.10,40.08
+     &        AC2(ID,IS,IND2) * CAX(ID,IS,2) * RDX(1) * 2.0               40.08 33.10
+     &       -AC2(ID,IS,IND4) * CAX(ID,IS,4) * RDX(1) * 0.5               40.98 40.08 33.10
+     &       +AC2(ID,IS,IND3) * CAX(ID,IS,3) * RDX(2) * 2.0               40.08 33.10
+     &       -AC2(ID,IS,IND5) * CAX(ID,IS,5) * RDX(2) * 0.5               40.98 40.08 33.10
+     &      +(AC2(ID,IS,IND2) * CAY(ID,IS,2) * RDY(1) * COSLAT(2) * 2.0   40.08 33.10
+     &       -AC2(ID,IS,IND4) * CAY(ID,IS,4) * RDY(1) * COSLAT(4) * 0.5   40.98 40.08 33.10
+     &       +AC2(ID,IS,IND3) * CAY(ID,IS,3) * RDY(2) * COSLAT(3) * 2.0   40.08 33.10
+     &       -AC2(ID,IS,IND5) * CAY(ID,IS,5) * RDY(2) * COSLAT(5) * 0.5   40.98 40.08 33.10
      &        ) / COSLAT(1) !33.10
             ENDIF
 
-          ELSE      ! switch 2<==>3, 4<==>5, and YMU<==>XMU  RDX(2)<==>RDX(3) RDX(4)<==>RDX(5) 33.10
+          ELSE      ! switch 2<==>3, 4<==>5 and YMU<==>XMU                33.10
 
 ! The diag part FXY1
-            FXY1 = 2.0*YMU(1)-0.5*YMU(6)                                  33.10,40.08
-     &           + 2.0*XMU(1)-0.5*XMU(6)                                  33.10,40.08
+            FXY1 = 1.5*YMU(1)+ 1.5*XMU(1)                                 40.98 40.08 33.10
 
             IF (KSPHER.EQ.0) THEN                                         33.08
 !             Cartesian coordinates
 
 !             the known, rhs part  FXY2                                   33.08
-              FXY2 = AC2(ID,IS,IND3) * YMU(3)*2.0                         33.10,40.08
-     &              -AC2(ID,IS,IND5) * 0.5*YMU(5)                         33.10,40.08
-     &              +AC2(ID,IS,IND2) * XMU(2)*2.0                         33.10,40.08
-     &              -AC2(ID,IS,IND4) * 0.5*XMU(4)                         33.10,40.08
+              FXY2 = AC2(ID,IS,IND3) * 2.0*YMU(3)                         40.08 33.10
+     &              -AC2(ID,IS,IND5) * 0.5*YMU(5)                         40.08 33.10
+     &              +AC2(ID,IS,IND2) * 2.0*XMU(2)                         40.08 33.10
+     &              -AC2(ID,IS,IND4) * 0.5*XMU(4)                         40.08 33.10
 
             ELSE
 !             Spherical coordinates
 
 !             the known, rhs part  FXY2                                   33.08
               FXY2 =
-     &        AC2(ID,IS,IND2) * CAX(ID,IS,2) * RDX(1) * 2.0               33.10,40.08
-     &       -AC2(ID,IS,IND4) * CAX(ID,IS,4) * RDX(3) * 0.5               33.10,40.08
-     &       +AC2(ID,IS,IND3) * CAX(ID,IS,3) * RDX(2) * 2.0               33.10,40.08
-     &       -AC2(ID,IS,IND5) * CAX(ID,IS,5) * RDX(4) * 0.5               33.10,40.08
-     &      +(AC2(ID,IS,IND2) * CAY(ID,IS,2) * RDY(1) * COSLAT(2) * 2.0   33.10,40.08
-     &       -AC2(ID,IS,IND4) * CAY(ID,IS,4) * RDY(3) * COSLAT(4) * 0.5   33.10,40.08
-     &       +AC2(ID,IS,IND3) * CAY(ID,IS,3) * RDY(2) * COSLAT(3) * 2.0   33.10,40.08
-     &       -AC2(ID,IS,IND5) * CAY(ID,IS,5) * RDY(4) * COSLAT(5) * 0.5   33.10,40.08
+     &        AC2(ID,IS,IND2) * CAX(ID,IS,2) * RDX(1) * 2.0               40.08 33.10
+     &       -AC2(ID,IS,IND4) * CAX(ID,IS,4) * RDX(1) * 0.5               40.98 40.08 33.10
+     &       +AC2(ID,IS,IND3) * CAX(ID,IS,3) * RDX(2) * 2.0               40.08 33.10
+     &       -AC2(ID,IS,IND5) * CAX(ID,IS,5) * RDX(2) * 0.5               40.98 40.08 33.10
+     &      +(AC2(ID,IS,IND2) * CAY(ID,IS,2) * RDY(1) * COSLAT(2) * 2.0   40.08 33.10
+     &       -AC2(ID,IS,IND4) * CAY(ID,IS,4) * RDY(1) * COSLAT(4) * 0.5   40.98 40.08 33.10
+     &       +AC2(ID,IS,IND3) * CAY(ID,IS,3) * RDY(2) * COSLAT(3) * 2.0   40.08 33.10
+     &       -AC2(ID,IS,IND5) * CAY(ID,IS,5) * RDY(2) * COSLAT(5) * 0.5   40.98 40.08 33.10
      &        )/ COSLAT(1)
             ENDIF
           END IF
@@ -2710,22 +2653,17 @@
 !         *** the term FXY2 is known, store in IMATRA ***
 !         *** the term FXY1 is unknown, store in IMATDA ***
 !
-!         --- Using an if statement like this might be helpful            40.08
-!            (see comments in STRSXY above).                              40.08
-!            I'm leaving it as-is for now.                                40.08
-!          IF((FXY1.GE.0).AND.(FXY2.GE.0))THEN
-!            IMATRA(ID,IS) = IMATRA(ID,IS) + FXY2
-!            IMATDA(ID,IS) = IMATDA(ID,IS) + FXY1
-!          ENDIF
           IMATRA(ID,IS) = IMATRA(ID,IS) + FXY2
           IMATDA(ID,IS) = IMATDA(ID,IS) + FXY1
           TRAC0(ID,IS,1) = TRAC0(ID,IS,1) - FXY2                          40.85
           TRAC1(ID,IS,1) = TRAC1(ID,IS,1) + FXY1                          40.85
+!
 !         *** test output ***
           IF ( ITEST .GE. 150 .AND. TESTFL ) THEN
               WRITE(PRINTF,6022) ID, FXY1, FXY2
  6022         FORMAT (' - ID FXY1 FXY2:', I4, 2(1X,E12.4))
           ENDIF
+!
  100    CONTINUE
  200  CONTINUE
 !
@@ -2775,7 +2713,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2807,7 +2745,7 @@
 !        33.08, July 98: SANDL: New subroutine using a Stelling and Leenderste
 !                        SANDL: scheme (Qo=0,Q1=1/6) is created.
 !        33.09, Aug. 99: extension with spherical coordinates
-!        40.02, Aug. 00: Avoide more than 19 continuatino lines
+!        40.02, Aug. 00: Avoid more than 19 continuation lines
 !        40.08, Feb. 03: Check for exceedence of soft CFL criterion
 !        40.41, Oct. 04: common blocks replaced by modules, include files removed
 !        40.85, Aug. 08: store xy-propagation for output purposes
@@ -3229,7 +3167,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -3547,7 +3485,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -3915,7 +3853,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -4198,7 +4136,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -4438,7 +4376,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -4648,7 +4586,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -4918,7 +4856,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -5228,7 +5166,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanCompdata.ftn90	2009-03-03 17:49:37.000000000 +0100
+++ SwanCompdata.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanCompUnstruc.ftn90	2009-03-03 17:49:37.000000000 +0100
+++ SwanCompUnstruc.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -32,11 +32,14 @@
 !   Authors
 !
 !   40.80: Marcel Zijlema
+!   40.85: Marcel Zijlema
+!   41.02: Marcel Zijlema
 !
 !   Updates
 !
-!   40.80, July 2007: New subroutine
-!   40.85, August 2008: add propagation, generation and redistribution terms for output purposes
+!   40.80,     July 2007: New subroutine
+!   40.85,   August 2008: add propagation, generation and redistribution terms for output purposes
+!   41.02, February 2009: implementation of diffraction
 !
 !   Purpose
 !
@@ -448,11 +451,11 @@
              !
              ! first guess settings
              !
-             IWIND     = 2            ! if first guess should be based on 1st generation mode, set IWIND = 1
+             IWIND     = 2               ! if first guess should be based on 1st generation mode, set IWIND = 1
              IWCAP     = 0
              IQUAD     = 0
-             PNUMS(20) = 1.E22        ! no limiter
-             PNUMS(30) = 0.           ! no under-relaxation
+             PNUMS(20) = 1.E22           ! no limiter
+             PNUMS(30) = 0.              ! no under-relaxation
              !
              write (PRINTF,101)
              !
@@ -480,6 +483,10 @@
           !
        endif
        !
+       ! calculate diffraction parameter and its derivatives
+       !
+       if ( IDIFFR /= 0 ) call SwanDiffPar ( ac2, compda(1,JDP2), spcsig )
+       !
        ! all vertices are set untagged except non-active ones and those where boundary conditions are given
        !
        vert(:)%fullupdated = .false.
@@ -766,10 +773,11 @@
                      !
 !TIMG                     call SWTSTA(118)
                      call SwanTranspAc ( amat  , rhs   , leakcf, ac2   , ac1   , &
-                                         cax   , cay   , cad   , cas   , anybin, &
-                                         rdx   , rdy   , obredf, spcsig, idcmin, &
-                                         idcmax, iscmin, iscmax, iddlow, iddtop, &
-                                         isslow, isstop, trac0 , trac1 )
+                                         cgo   , cax   , cay   , cad   , cas   , &
+                                         anybin, rdx   , rdy   , spcsig, spcdir, &
+                                         obredf, idcmin, idcmax, iscmin, iscmax, &
+                                         iddlow, iddtop, isslow, isstop, trac0 , &
+                                         trac1 )
 !TIMG                     call SWTSTO(118)
                      !
                      ! compute the source part of the action balance equation
--- SwanComputeForce.ftn90	2008-05-25 14:59:42.000000000 +0200
+++ SwanComputeForce.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -50,7 +50,9 @@
 !   Modules used
 !
     use ocpcomm4
+    use swcomm2
     use swcomm3
+    use swcomm4
     use SwanGriddata
     use SwanGridobjects
 !
@@ -82,7 +84,8 @@
     !
     integer, dimension(3)                 :: v        ! vertices in present cell
     !
-    real                                  :: area     ! twices the area of centroid dual around present vertex
+    double precision                      :: area     ! twices the area of centroid dual around present vertex
+    real                                  :: cslat    ! cosine of latitude
     real                                  :: deploc   ! local depth
     real                                  :: dsxxdx   ! x-gradient of sxx
     real                                  :: dsxydx   ! x-gradient of sxy
@@ -100,10 +103,10 @@
     real                                  :: sxxsums  ! cumulated sxx over frequency space
     real                                  :: sxysums  ! cumulated sxy over frequency space
     real                                  :: syysums  ! cumulated syy over frequency space
-    real                                  :: x0       ! x-coordinate of the centroid of present cell
-    real                                  :: x1       ! x-coordinate of the centroid of next cell
-    real                                  :: y0       ! y-coordinate of the centroid of present cell
-    real                                  :: y1       ! y-coordinate of the centroid of next cell
+    double precision                      :: x0       ! x-coordinate of the centroid of present cell
+    double precision                      :: x1       ! x-coordinate of the centroid of next cell
+    double precision                      :: y0       ! y-coordinate of the centroid of present cell
+    double precision                      :: y1       ! y-coordinate of the centroid of next cell
     !
     real, dimension(1)                    :: cg       ! group velocity
     real, dimension(1)                    :: k        ! wave number
@@ -205,7 +208,7 @@
        !
        ! first, compute gradients of the radiation stresses in vertices
        !
-       area   = 0.
+       area   = 0d0
        dsxxdx = 0.
        dsxydx = 0.
        dsxydy = 0.
@@ -261,28 +264,41 @@
           !
           ! compute contribution to x-gradient of radiation stresses sxx and sxy
           !
-          dsxxdx = dsxxdx + ( sxx0 + sxx1 ) * ( y1 - y0 )
-          dsxydx = dsxydx + ( sxy0 + sxy1 ) * ( y1 - y0 )
+          dsxxdx = dsxxdx + ( sxx0 + sxx1 ) * real( y1 - y0 )
+          dsxydx = dsxydx + ( sxy0 + sxy1 ) * real( y1 - y0 )
           !
           ! compute contribution to y-gradient of radiation stresses sxy and syy
           !
-          dsxydy = dsxydy + ( sxy0 + sxy1 ) * ( x1 - x0 )
-          dsyydy = dsyydy + ( syy0 + syy1 ) * ( x1 - x0 )
+          dsxydy = dsxydy + ( sxy0 + sxy1 ) * real( x1 - x0 )
+          dsyydy = dsyydy + ( syy0 + syy1 ) * real( x1 - x0 )
           !
        enddo
        !
        ! if area is non-positive, give error and go to next vertex
        !
-       if ( area <= 0. ) then
+       if ( area <= 0d0 ) then
           write (msgstr, '(a,i5)') ' Area of centroid dual is negative or zero in vertex ', ivert
           call msgerr( 2, trim(msgstr) )
           cycle vertexloop
        endif
        !
-       dsxxdx =  dsxxdx/area
-       dsxydx =  dsxydx/area
-       dsxydy = -dsxydy/area
-       dsyydy = -dsyydy/area
+       dsxxdx =  dsxxdx/real(area)
+       dsxydx =  dsxydx/real(area)
+       dsxydy = -dsxydy/real(area)
+       dsyydy = -dsyydy/real(area)
+       !
+       ! in case of spherical coordinates, transform back to Cartesian coordinates
+       !
+       if ( KSPHER > 0 ) then
+          !
+          cslat = cos(DEGRAD*(vert(ivert)%attr(VERTY) + YOFFS))
+          !
+          dsxxdx = dsxxdx/(cslat * LENDEG)
+          dsxydx = dsxydx/(cslat * LENDEG)
+          dsxydy = dsxydy/LENDEG
+          dsyydy = dsyydy/LENDEG
+          !
+       endif
        !
        ! finally, compute wave-induced force
        !
--- SwanConvAccur.ftn90	2008-05-25 14:59:42.000000000 +0200
+++ SwanConvAccur.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanConvStopc.ftn90	2009-03-03 17:49:37.000000000 +0100
+++ SwanConvStopc.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanCreateEdges.ftn90	2009-03-03 17:49:37.000000000 +0100
+++ SwanCreateEdges.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanCrossObstacle.ftn90	2008-05-25 14:59:43.000000000 +0200
+++ SwanCrossObstacle.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanDiffPar.ftn90	1970-01-01 01:00:00.000000000 +0100
+++ SwanDiffPar.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -0,0 +1,462 @@
+subroutine SwanDiffPar ( ac2, dep2, spcsig )
+!
+!   --|-----------------------------------------------------------|--
+!     | Delft University of Technology                            |
+!     | Faculty of Civil Engineering and Geosciences              |
+!     | Environmental Fluid Mechanics Section                     |
+!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
+!     |                                                           |
+!     | Programmer: Marcel Zijlema                                |
+!   --|-----------------------------------------------------------|--
+!
+!
+!     SWAN (Simulating WAves Nearshore); a third generation wave model
+!     Copyright (C) 2009  Delft University of Technology
+!
+!     This program is free software; you can redistribute it and/or
+!     modify it under the terms of the GNU General Public License as
+!     published by the Free Software Foundation; either version 2 of
+!     the License, or (at your option) any later version.
+!
+!     This program is distributed in the hope that it will be useful,
+!     but WITHOUT ANY WARRANTY; without even the implied warranty of
+!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
+!     GNU General Public License for more details.
+!
+!     A copy of the GNU General Public License is available at
+!     http://www.gnu.org/copyleft/gpl.html#SEC3
+!     or by writing to the Free Software Foundation, Inc.,
+!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
+!
+!
+!   Authors
+!
+!   41.02: Marcel Zijlema
+!
+!   Updates
+!
+!   41.02, February 2009: New subroutine
+!
+!   Purpose
+!
+!   Computes diffraction parameter and its derivatives in vertices
+!
+!   Method
+!
+!   Diffraction is approximated using the eikonal equation which
+!   relates the wavenumber K to the separation factor k. Several
+!   expressions of the eikonal equation have been presented in
+!   the literature:
+!                                 DH
+!   Battjes (1968):   K^2 = k^2 + --
+!                                 H
+!
+!   where H is the wave height and D is the Laplacian operator.
+!
+!                                 D.(pDH)
+!   Berkhoff (1972):  K^2 = k^2 + -------
+!                                   pH
+!
+!   where p = cc_g and D is the gradient operator in this case.
+!
+!   In both cases, the eikonal equation may be written as follows:
+!
+!   K = k (1+delta)^0.5
+!
+!   with
+!
+!           D.(pDH)
+!   delta = -------
+!           k^2 pH
+!
+!   From implementation point of view, the Battjes' eikonal
+!   equation can be obtained if
+!                                c_g = k
+!
+!   Modules used
+!
+    use ocpcomm4
+    use swcomm2
+    use swcomm3
+    use swcomm4
+    use m_diffr
+    use SwanGriddata
+    use SwanGridobjects
+!
+    implicit none
+!
+!   Argument variables
+!
+    real, dimension(MDC,MSC,nverts), intent(in) :: ac2      ! action density at current time level
+    real, dimension(nverts), intent(in)         :: dep2     ! water depth at current time level
+    real, dimension(MSC), intent(in)            :: spcsig   ! relative frequency bins
+!
+!   Local variables
+!
+    integer                               :: icell          ! index of present cell
+    integer, save                         :: ient = 0       ! number of entries in this subroutine
+    integer                               :: ivert          ! loop counter over vertices
+    integer                               :: jc             ! loop counter
+    integer                               :: jcell          ! index of next cell
+    integer, parameter                    :: jeiko=1        ! choice parameter:
+                                                            ! 0 = eikonal equation according to Battjes (1968)
+                                                            ! 1 = eikonal equation according to Berkhoff (1972)
+    !
+    integer, dimension(3)                 :: v              ! vertices in present cell
+    !
+    real                                  :: area0          ! area of present cell
+    real                                  :: area1          ! area of next cell
+    double precision                      :: carea          ! twices the area of centroid dual around present vertex
+    real                                  :: cg0            ! mean group velocity in centroid of present cell
+    real                                  :: cg1            ! mean group velocity in centroid of next cell
+    real                                  :: cslat          ! cosine of latitude
+    real                                  :: ctot           ! zeroth moment of energy times group velocity
+    real                                  :: delta          ! local diffraction parameter
+    real                                  :: denom          ! a denominator
+    real                                  :: deploc         ! local depth
+    real                                  :: dgx0           ! x-component of diffusion gradient inside present cell
+    real                                  :: dgx1           ! x-component of diffusion gradient inside next cell
+    real                                  :: dgxdx          ! x-gradient of x-diffusion gradient component
+    real                                  :: dgy0           ! y-component of diffusion gradient inside present cell
+    real                                  :: dgy1           ! y-component of diffusion gradient inside next cell
+    real                                  :: dgydy          ! y-gradient of y-diffusion gradient component
+    real                                  :: dhsdx          ! x-gradient of wave height
+    real                                  :: dhsdy          ! y-gradient of wave height
+    real                                  :: difp0          ! diffraction parameter in centroid of present cell
+    real                                  :: difp1          ! diffraction parameter in centroid of next cell
+    real                                  :: etot           ! zeroth moment of the variance spectrum
+    real                                  :: fmax           ! upper bound of frequency space for integration
+    real                                  :: fmin           ! lower bound of frequency space for integration
+    real                                  :: k0             ! mean wave number in centroid of present cell
+    real                                  :: k1             ! mean wave number in centroid of next cell
+    real                                  :: ktot           ! zeroth moment of energy times wave number
+    double precision                      :: x0             ! x-coordinate of the centroid of present cell
+    double precision                      :: x1             ! x-coordinate of the centroid of next cell
+    double precision                      :: y0             ! y-coordinate of the centroid of present cell
+    double precision                      :: y1             ! y-coordinate of the centroid of next cell
+    !
+    real, dimension(MSC)                  :: cgloc          ! group velocity
+    real, dimension(MSC)                  :: kloc           ! wave number
+    real, dimension(MSC)                  :: n              ! ratio of group and phase velocity
+    real, dimension(MSC)                  :: nd             ! derivative of n with respect to depth
+    real, dimension(MDC)                  :: ecs            ! help array containing (co)sine of spectral directions
+    !
+    real, dimension(:), allocatable       :: cg             ! mean group velocity
+    real, dimension(:), allocatable       :: hs             ! wave height
+    real, dimension(:), allocatable       :: k              ! mean wave number
+    !
+    real                                  :: SwanIntgratSpc ! integration of variance over a part of frequency space
+    !
+    type(celltype), dimension(:), pointer :: cell           ! datastructure for cells with their attributes
+    type(verttype), dimension(:), pointer :: vert           ! datastructure for vertices with their attributes
+!
+!   Structure
+!
+!   Description of the pseudo code
+!
+!   Source text
+!
+    if (ltrace) call strace (ient,'SwanDiffPar')
+    !
+    ! point to vertex and cell objects
+    !
+    vert => gridobject%vert_grid
+    cell => gridobject%cell_grid
+    !
+    ! allocation and initialization of wave height and mean wave parameters
+    !
+    allocate(hs(nverts))
+    allocate( k(nverts))
+    allocate(cg(nverts))
+    !
+    hs = 0.
+    k  = 10.
+    cg = 0.
+    !
+    ! compute total energy, mean wave number and mean group velocity in vertices
+    !
+    do ivert = 1, nverts
+       !
+       deploc = dep2(ivert)
+       !
+       if ( deploc <= DEPMIN ) cycle
+       !
+       ! compute group velocity and wave number for all frequencies
+       !
+       call KSCIP1 (MSC,spcsig,deploc,kloc,cgloc,n,nd)
+       !
+       ! integration over f in [0,infty]
+       !
+       fmin = 0.
+       fmax = 1000.
+       ecs  = 1.
+       !
+       etot = SwanIntgratSpc(0.   , fmin, fmax, spcsig, ecs,            &
+                             kloc , ecs , 0.  , 0.    , ac2(:,:,ivert), &
+                             1    )
+       !
+       ktot = SwanIntgratSpc(1.   , fmin, fmax, spcsig, ecs,            &
+                             kloc , ecs , 0.  , 0.    , ac2(:,:,ivert), &
+                             3    )
+       !
+       ctot = SwanIntgratSpc(1.   , fmin, fmax, spcsig, ecs,            &
+                             cgloc, ecs , 0.  , 0.    , ac2(:,:,ivert), &
+                             4    )
+       !
+       if ( etot > 0. ) then
+          hs(ivert) = 4.*sqrt(etot)
+          k (ivert) = ktot/etot
+          cg(ivert) = ctot/etot
+       endif
+       !
+    enddo
+    !
+    if ( jeiko == 0 ) cg = k
+    !
+    ! compute diffraction parameter in vertices
+    !
+    DIFPARAM = 1.
+    !
+    vertexloop : do ivert = 1, nverts
+       !
+       if ( vert(ivert)%atti(VMARKER) == 1 ) cycle vertexloop    ! boundary vertex
+       !
+       cslat = cos(DEGRAD*(vert(ivert)%attr(VERTY) + YOFFS))
+       !
+       ! compute contributions to the Laplacian in present vertex
+       !
+       carea = 0d0
+       dgxdx = 0.
+       dgydy = 0.
+       !
+       ! loop over cells around considered vertex
+       !
+       do jc = 1, vert(ivert)%noc
+          !
+          ! get present cell and its vertices
+          !
+          icell = vert(ivert)%cell(jc)%atti(CELLID)
+          !
+          v(1) = cell(icell)%atti(CELLV1)
+          v(2) = cell(icell)%atti(CELLV2)
+          v(3) = cell(icell)%atti(CELLV3)
+          !
+          if ( dep2(v(1)) <= DEPMIN .or. dep2(v(2)) <= DEPMIN .or. dep2(v(3)) <= DEPMIN ) cycle vertexloop
+          !
+          ! determine centroid of present cell
+          !
+          x0    = cell(icell)%attr(CELLCX)
+          y0    = cell(icell)%attr(CELLCY)
+          area0 = cell(icell)%attr(CELLAREA)
+          !
+          ! determine mean wave number and group velocity in centroid in present cell
+          !
+          cg0 = ( cg(v(1)) + cg(v(2)) + cg(v(3)) )/ 3.
+          k0  = ( k (v(1)) + k (v(2)) + k (v(3)) )/ 3.
+          !
+          ! determine derivatives of wave height inside present cell
+          !
+          dhsdx = 0.5*( hs(v(1))*(ycugrd(v(2))-ycugrd(v(3))) + &
+                        hs(v(2))*(ycugrd(v(3))-ycugrd(v(1))) + &
+                        hs(v(3))*(ycugrd(v(1))-ycugrd(v(2))) )/area0
+          !
+          dhsdy = 0.5*( hs(v(1))*(xcugrd(v(3))-xcugrd(v(2))) + &
+                        hs(v(2))*(xcugrd(v(1))-xcugrd(v(3))) + &
+                        hs(v(3))*(xcugrd(v(2))-xcugrd(v(1))) )/area0
+          !
+          ! in case of spherical coordinates, transform back to Cartesian coordinates
+          !
+          if ( KSPHER > 0 ) then
+             !
+             dhsdx = dhsdx/(cslat * LENDEG)
+             dhsdy = dhsdy/LENDEG
+             !
+          endif
+          !
+          ! determine diffusion gradients in centroid of present cell
+          !
+          dgx0 = cg0*dhsdx/k0
+          dgy0 = cg0*dhsdy/k0
+          !
+          ! get next cell in counterclockwise direction
+          !
+          jcell = vert(ivert)%cell(jc)%atti(NEXTCELL)
+          !
+          v(1) = cell(jcell)%atti(CELLV1)
+          v(2) = cell(jcell)%atti(CELLV2)
+          v(3) = cell(jcell)%atti(CELLV3)
+          !
+          ! determine centroid of next cell
+          !
+          x1    = cell(jcell)%attr(CELLCX)
+          y1    = cell(jcell)%attr(CELLCY)
+          area1 = cell(jcell)%attr(CELLAREA)
+          !
+          ! determine mean wave number and group velocity in centroid in next cell
+          !
+          cg1 = ( cg(v(1)) + cg(v(2)) + cg(v(3)) )/ 3.
+          k1  = ( k (v(1)) + k (v(2)) + k (v(3)) )/ 3.
+          !
+          ! determine derivatives of wave height inside next cell
+          !
+          dhsdx = 0.5*( hs(v(1))*(ycugrd(v(2))-ycugrd(v(3))) + &
+                        hs(v(2))*(ycugrd(v(3))-ycugrd(v(1))) + &
+                        hs(v(3))*(ycugrd(v(1))-ycugrd(v(2))) )/area1
+          !
+          dhsdy = 0.5*( hs(v(1))*(xcugrd(v(3))-xcugrd(v(2))) + &
+                        hs(v(2))*(xcugrd(v(1))-xcugrd(v(3))) + &
+                        hs(v(3))*(xcugrd(v(2))-xcugrd(v(1))) )/area1
+          !
+          ! in case of spherical coordinates, transform back to Cartesian coordinates
+          !
+          if ( KSPHER > 0 ) then
+             !
+             dhsdx = dhsdx/(cslat * LENDEG)
+             dhsdy = dhsdy/LENDEG
+             !
+          endif
+          !
+          ! determine diffusion gradients in centroid of next cell
+          !
+          dgx1 = cg1*dhsdx/k1
+          dgy1 = cg1*dhsdy/k1
+          !
+          ! compute contribution to area of centroid dual
+          !
+          carea = carea + x0*y1 - x1*y0
+          !
+          ! compute contribution to x-gradient of x-diffusion gradient in centroid dual
+          !
+          dgxdx = dgxdx + ( dgx0 + dgx1 ) * real( y1 - y0 )
+          !
+          ! compute contribution to y-gradient of y-diffusion gradient in centroid dual
+          !
+          dgydy = dgydy + ( dgy0 + dgy1 ) * real( x0 - x1 )
+          !
+       enddo
+       !
+       if ( carea > 0d0 ) then
+          !
+          dgxdx = dgxdx/real(carea)
+          dgydy = dgydy/real(carea)
+          !
+          ! in case of spherical coordinates, transform back to Cartesian coordinates
+          !
+          if ( KSPHER > 0 ) then
+             !
+             dgxdx = dgxdx/(cslat * LENDEG)
+             dgydy = dgydy/LENDEG
+             !
+          endif
+          !
+          denom = k(ivert)*cg(ivert)*hs(ivert)
+          !
+          if ( denom > 0. ) then
+             delta = (dgxdx + dgydy)/denom
+          else
+             delta = 0.
+          endif
+          !
+          if ( delta > -1. ) DIFPARAM(ivert) = sqrt(1.+delta)
+          !
+       endif
+       !
+    enddo vertexloop
+    !
+    ! deallocation of wave parameters
+    !
+    deallocate(cg)
+    deallocate(hs)
+    deallocate(k )
+    !
+    ! compute derivatives of diffraction parameter in vertices
+    !
+    DIFPARDX = 0.
+    DIFPARDY = 0.
+    !
+    vertexloop2 : do ivert = 1, nverts
+       !
+       if ( vert(ivert)%atti(VMARKER) == 1 ) cycle vertexloop2   ! boundary vertex
+       !
+       cslat = cos(DEGRAD*(vert(ivert)%attr(VERTY) + YOFFS))
+       !
+       carea = 0d0
+       dgxdx = 0.
+       dgydy = 0.
+       !
+       ! loop over cells around considered vertex
+       !
+       do jc = 1, vert(ivert)%noc
+          !
+          ! get present cell and its vertices
+          !
+          icell = vert(ivert)%cell(jc)%atti(CELLID)
+          !
+          v(1) = cell(icell)%atti(CELLV1)
+          v(2) = cell(icell)%atti(CELLV2)
+          v(3) = cell(icell)%atti(CELLV3)
+          !
+          if ( dep2(v(1)) <= DEPMIN .or. dep2(v(2)) <= DEPMIN .or. dep2(v(3)) <= DEPMIN ) cycle vertexloop2
+          !
+          ! determine centroid of present cell
+          !
+          x0 = cell(icell)%attr(CELLCX)
+          y0 = cell(icell)%attr(CELLCY)
+          !
+          ! determine diffraction parameter in centroid in present cell
+          !
+          difp0 = ( DIFPARAM(v(1)) + DIFPARAM(v(2)) + DIFPARAM(v(3)) )/ 3.
+          !
+          ! get next cell in counterclockwise direction
+          !
+          jcell = vert(ivert)%cell(jc)%atti(NEXTCELL)
+          !
+          v(1) = cell(jcell)%atti(CELLV1)
+          v(2) = cell(jcell)%atti(CELLV2)
+          v(3) = cell(jcell)%atti(CELLV3)
+          !
+          ! determine centroid of next cell
+          !
+          x1 = cell(jcell)%attr(CELLCX)
+          y1 = cell(jcell)%attr(CELLCY)
+          !
+          ! determine diffraction parameter in centroid of next cell
+          !
+          difp1 = ( DIFPARAM(v(1)) + DIFPARAM(v(2)) + DIFPARAM(v(3)) )/ 3.
+          !
+          ! compute contribution to area of centroid dual
+          !
+          carea = carea + x0*y1 - x1*y0
+          !
+          ! compute x-gradient of diffraction parameter
+          !
+          dgxdx = dgxdx + ( difp0 + difp1 ) * real( y1 - y0 )
+          !
+          ! compute y-gradient of diffraction parameter
+          !
+          dgydy = dgydy + ( difp0 + difp1 ) * real( x0 - x1 )
+          !
+       enddo
+       !
+       if ( carea > 0d0 ) then
+          !
+          dgxdx = dgxdx/real(carea)
+          dgydy = dgydy/real(carea)
+          !
+          ! in case of spherical coordinates, transform back to Cartesian coordinates
+          !
+          if ( KSPHER > 0 ) then
+             !
+             dgxdx = dgxdx/(cslat * LENDEG)
+             dgydy = dgydy/LENDEG
+             !
+          endif
+          !
+          DIFPARDX(ivert) = dgxdx
+          DIFPARDY(ivert) = dgydy
+          !
+       endif
+       !
+    enddo vertexloop2
+    !
+end subroutine SwanDiffPar
--- SwanDispParm.ftn90	2008-05-25 14:59:43.000000000 +0200
+++ SwanDispParm.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- swan.edt	2009-03-03 17:49:37.000000000 +0100
+++ swan.edt	2009-03-03 17:50:34.000000000 +0100
@@ -90,7 +90,7 @@
 !         | WESTHuysen                  |
 !
 !
-!         | ->KOMen    [cds2] [stpm] [powst] [delta] [powk]  |
+!         | -> KOMen   [cds2] [stpm] [powst] [delta] [powk]  |
 !         |                                                  |
 !         |   JANSsen  [cds1]  [delta] [pwtail]              |
 !         |                                                  |
@@ -110,15 +110,29 @@
 !   MDIA LAMbda < [lambda] > <                                 >
 !                             | CNL4_12 < [cnl4_1] [cnl4_2] > |
 !
-!   BREaking CONSTANT [alpha] [gamma]
+!        | -> CON [alpha] [gamma]                                      |
+!        |                                                             |
+!        |    VAR [alpha] [gammin] [gammax] [gamneg] [coeff1] [coeff2] |
+!        |                                                             |
+!   BRE <     RUE [alpha] [a] [b]                                       >
+!        |                                                             |
+!        |    TG  [alpha] [gamma] [pown]                               |
+!        |                                                             |
+!        !    BIP [alpha] [pown] [bref]                                |
+!
+!              |             | -> CONstant [cfjon]
+!              | -> JONswap <
+!              |             |    VARiable [cfj1] [cfj2] [dsp1] [dsp2]
+!   FRICtion  <
+!              |    COLLins [cfw] [cfc]
+!              |
+!              |    MADsen  [kn]
 !
-!   FRICTION  JONSWAP [cfjon] | COLLINS [cfw] | MADSEN [kn]
-!
-!   TRIAD   [trfac] [cutfr] [urcrit] [urslim]
+!   TRIad   [trfac] [cutfr] [urcrit] [urslim]
 !
 !   LIMiter [ursell] [qb]
 !
-!               |  TRANSm [trcoef]                          |
+!               | -> TRANSm [trcoef]                        |
 !   OBSTacle   <                                            |
 !               |       | -> GODA [hgt] [alpha] [beta]       >              &
 !               |  DAM <                                    |
--- SwanFindObstacles.ftn90	2008-05-25 14:59:43.000000000 +0200
+++ SwanFindObstacles.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanFindPoint.ftn90	2008-05-25 14:59:43.000000000 +0200
+++ SwanFindPoint.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanGridCell.ftn90	2009-03-03 17:49:37.000000000 +0100
+++ SwanGridCell.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanGriddata.ftn90	2009-03-03 17:49:37.000000000 +0100
+++ SwanGriddata.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanGridFace.ftn90	2009-03-03 17:49:37.000000000 +0100
+++ SwanGridFace.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanGridobjects.ftn90	2009-03-03 17:49:37.000000000 +0100
+++ SwanGridobjects.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanGridTopology.ftn90	2008-05-25 14:59:44.000000000 +0200
+++ SwanGridTopology.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanGridVert.ftn90	2008-05-25 14:59:44.000000000 +0200
+++ SwanGridVert.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanGSECorr.ftn90	1970-01-01 01:00:00.000000000 +0100
+++ SwanGSECorr.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -0,0 +1,282 @@
+subroutine SwanGSECorr ( rhs, ac2, cgo, spcdir, idcmin, idcmax, isslow, isstop, trac0 )
+!
+!   --|-----------------------------------------------------------|--
+!     | Delft University of Technology                            |
+!     | Faculty of Civil Engineering and Geosciences              |
+!     | Environmental Fluid Mechanics Section                     |
+!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
+!     |                                                           |
+!     | Programmer: Marcel Zijlema                                |
+!   --|-----------------------------------------------------------|--
+!
+!
+!     SWAN (Simulating WAves Nearshore); a third generation wave model
+!     Copyright (C) 2009  Delft University of Technology
+!
+!     This program is free software; you can redistribute it and/or
+!     modify it under the terms of the GNU General Public License as
+!     published by the Free Software Foundation; either version 2 of
+!     the License, or (at your option) any later version.
+!
+!     This program is distributed in the hope that it will be useful,
+!     but WITHOUT ANY WARRANTY; without even the implied warranty of
+!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
+!     GNU General Public License for more details.
+!
+!     A copy of the GNU General Public License is available at
+!     http://www.gnu.org/copyleft/gpl.html#SEC3
+!     or by writing to the Free Software Foundation, Inc.,
+!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
+!
+!
+!   Authors
+!
+!   41.00: Marcel Zijlema
+!
+!   Updates
+!
+!   41.00, February 2009: New subroutine
+!
+!   Purpose
+!
+!   Computes waveage-dependent diffusion terms in x-y space to counteract the garden-sprinkler effect
+!
+!   Modules used
+!
+    use ocpcomm4
+    use swcomm2
+    use swcomm3
+    use swcomm4
+    use SwanGriddata
+    use SwanGridobjects
+    use SwanCompdata
+!
+    implicit none
+!
+!   Argument variables
+!
+    integer, intent(in)                         :: isslow ! minimum frequency that is propagated within a sweep
+    integer, intent(in)                         :: isstop ! maximum frequency that is propagated within a sweep
+    !
+    integer, dimension(MSC), intent(in)         :: idcmax ! maximum frequency-dependent counter in directional space
+    integer, dimension(MSC), intent(in)         :: idcmin ! minimum frequency-dependent counter in directional space
+    !
+    real, dimension(MDC,MSC,nverts), intent(in) :: ac2    ! action density at current time level
+    real, dimension(MSC,ICMAX), intent(in)      :: cgo    ! group velocity
+    real, dimension(MDC,MSC), intent(inout)     :: rhs    ! right-hand side of system of equations in (sigma,theta) space
+    real, dimension(MDC,6), intent(in)          :: spcdir ! (*,1): spectral direction bins (radians)
+                                                          ! (*,2): cosine of spectral directions
+                                                          ! (*,3): sine of spectral directions
+                                                          ! (*,4): cosine^2 of spectral directions
+                                                          ! (*,5): cosine*sine of spectral directions
+                                                          ! (*,6): sine^2 of spectral directions
+    real, dimension(MDC,MSC,MTRNP), intent(out) :: trac0  ! explicit part of propagation in present vertex for output purposes
+!
+!   Local variables
+!
+    integer                               :: icell    ! index of present cell
+    integer                               :: id       ! loop counter over direction bins
+    integer                               :: iddum    ! counter in directional space for considered sweep
+    integer, save                         :: ient = 0 ! number of entries in this subroutine
+    integer                               :: is       ! loop counter over frequency bins
+    integer                               :: ivert    ! index of present vertex
+    integer                               :: jc       ! loop counter
+    integer                               :: jcell    ! index of next cell
+    !
+    integer, dimension(3)                 :: v        ! vertices in present cell
+    !
+    real                                  :: area0    ! area of present cell
+    real                                  :: area1    ! area of next cell
+    double precision                      :: carea    ! twices the area of centroid dual around present vertex
+    real                                  :: cslat    ! cosine of latitude
+    real                                  :: dac2dx   ! x-derivative of action density
+    real                                  :: dac2dy   ! y-derivative of action density
+    real                                  :: dcg      ! group velocity difference across frequency bin
+    real                                  :: dgx0     ! x-component of diffusion gradient inside present cell
+    real                                  :: dgx1     ! x-component of diffusion gradient inside next cell
+    real                                  :: dgxdx    ! x-gradient of x-diffusion gradient component
+    real                                  :: dgy0     ! y-component of diffusion gradient inside present cell
+    real                                  :: dgy1     ! y-component of diffusion gradient inside next cell
+    real                                  :: dgydy    ! y-gradient of y-diffusion gradient component
+    real                                  :: dnn      ! waveage-dependent diffusion coefficient normal to propagation direction
+    real                                  :: dss      ! waveage-dependent diffusion coefficient in propagation direction
+    real                                  :: dxx      ! xx-component of diffusion coefficient in Cartesian coordinates
+    real                                  :: dxy      ! xy-component of diffusion coefficient in Cartesian coordinates
+    real                                  :: dyy      ! yy-component of diffusion coefficient in Cartesian coordinates
+    double precision                      :: x0       ! x-coordinate of the centroid of present cell
+    double precision                      :: x1       ! x-coordinate of the centroid of next cell
+    double precision                      :: y0       ! y-coordinate of the centroid of present cell
+    double precision                      :: y1       ! y-coordinate of the centroid of next cell
+    !
+    type(celltype), dimension(:), pointer :: cell     ! datastructure for cells with their attributes
+    type(verttype), dimension(:), pointer :: vert     ! datastructure for vertices with their attributes
+!
+!   Structure
+!
+!   Description of the pseudo code
+!
+!   Source text
+!
+    if (ltrace) call strace (ient,'SwanGSECorr')
+    !
+    ! point to vertex and cell objects
+    !
+    vert => gridobject%vert_grid
+    cell => gridobject%cell_grid
+    !
+    ivert = vs(1)
+    !
+    if ( vert(ivert)%atti(VMARKER) == 1 ) return    ! no GSE correction in boundary vertex
+    !
+    cslat = cos(DEGRAD*(vert(ivert)%attr(VERTY) + YOFFS))
+    !
+    do is = isslow, isstop
+       !
+       ! calculate waveage-dependent diffusion coefficients in polar coordinates
+       !
+       if ( is == 1 ) then
+          dcg = abs(cgo(is+1,1)-cgo(is,1))
+       elseif ( is == isstop ) then
+          dcg = abs(cgo(is,1)-cgo(is-1,1))
+       else
+          dcg = 0.5 * abs(cgo(is+1,1)-cgo(is-1,1))
+       endif
+       !
+       dss = dcg**2*WAVAGE/12.
+       dnn = (cgo(is,1)*DDIR)**2 * WAVAGE/12.
+       !
+       do iddum = idcmin(is), idcmax(is)
+          id = mod ( iddum - 1 + MDC , MDC ) + 1
+          !
+          ! calculate diffusion coefficients in Cartesian coordinates
+          !
+          dxx = dss*spcdir(id,4) + dnn*spcdir(id,6)
+          dyy = dss*spcdir(id,6) + dnn*spcdir(id,4)
+          dxy = (dss-dnn)*spcdir(id,5)
+          !
+          ! compute contribution to the diffusion terms in present vertex
+          !
+          carea = 0d0
+          dgxdx = 0.
+          dgydy = 0.
+          !
+          ! loop over cells around considered vertex
+          !
+          do jc = 1, vert(ivert)%noc
+             !
+             ! get present cell and its vertices
+             !
+             icell = vert(ivert)%cell(jc)%atti(CELLID)
+             !
+             v(1) = cell(icell)%atti(CELLV1)
+             v(2) = cell(icell)%atti(CELLV2)
+             v(3) = cell(icell)%atti(CELLV3)
+             !
+             ! determine centroid and area of present cell
+             !
+             x0    = cell(icell)%attr(CELLCX)
+             y0    = cell(icell)%attr(CELLCY)
+             area0 = cell(icell)%attr(CELLAREA)
+             !
+             ! determine derivatives of action density inside present cell
+             !
+             dac2dx = 0.5*( ac2(id,is,v(1))*(ycugrd(v(2))-ycugrd(v(3))) + &
+                            ac2(id,is,v(2))*(ycugrd(v(3))-ycugrd(v(1))) + &
+                            ac2(id,is,v(3))*(ycugrd(v(1))-ycugrd(v(2))) )/area0
+             !
+             dac2dy = 0.5*( ac2(id,is,v(1))*(xcugrd(v(3))-xcugrd(v(2))) + &
+                            ac2(id,is,v(2))*(xcugrd(v(1))-xcugrd(v(3))) + &
+                            ac2(id,is,v(3))*(xcugrd(v(2))-xcugrd(v(1))) )/area0
+             !
+             ! in case of spherical coordinates, transform back to Cartesian coordinates
+             !
+             if ( KSPHER > 0 ) then
+                !
+                dac2dx = dac2dx/(cslat * LENDEG)
+                dac2dy = dac2dy/LENDEG
+                !
+             endif
+             !
+             ! determine diffusion gradients in centroid of present cell
+             !
+             dgx0 = dxx*dac2dx + dxy*dac2dy
+             dgy0 = dxy*dac2dx + dyy*dac2dy
+             !
+             ! get next cell in counterclockwise direction
+             !
+             jcell = vert(ivert)%cell(jc)%atti(NEXTCELL)
+             !
+             v(1) = cell(jcell)%atti(CELLV1)
+             v(2) = cell(jcell)%atti(CELLV2)
+             v(3) = cell(jcell)%atti(CELLV3)
+             !
+             ! determine centroid and area of next cell
+             !
+             x1    = cell(jcell)%attr(CELLCX)
+             y1    = cell(jcell)%attr(CELLCY)
+             area1 = cell(jcell)%attr(CELLAREA)
+             !
+             ! determine derivatives of action density inside next cell
+             !
+             dac2dx = 0.5*( ac2(id,is,v(1))*(ycugrd(v(2))-ycugrd(v(3))) + &
+                            ac2(id,is,v(2))*(ycugrd(v(3))-ycugrd(v(1))) + &
+                            ac2(id,is,v(3))*(ycugrd(v(1))-ycugrd(v(2))) )/area1
+             !
+             dac2dy = 0.5*( ac2(id,is,v(1))*(xcugrd(v(3))-xcugrd(v(2))) + &
+                            ac2(id,is,v(2))*(xcugrd(v(1))-xcugrd(v(3))) + &
+                            ac2(id,is,v(3))*(xcugrd(v(2))-xcugrd(v(1))) )/area1
+             !
+             ! in case of spherical coordinates, transform back to Cartesian coordinates
+             !
+             if ( KSPHER > 0 ) then
+                !
+                dac2dx = dac2dx/(cslat * LENDEG)
+                dac2dy = dac2dy/LENDEG
+                !
+             endif
+             !
+             ! determine diffusion gradients in centroid of next cell
+             !
+             dgx1 = dxx*dac2dx + dxy*dac2dy
+             dgy1 = dxy*dac2dx + dyy*dac2dy
+             !
+             ! compute contribution to area of centroid dual
+             !
+             carea = carea + x0*y1 - x1*y0
+             !
+             ! compute contribution to x-gradient of x-diffusion gradient in centroid dual
+             !
+             dgxdx = dgxdx + ( dgx0 + dgx1 ) * real( y1 - y0 )
+             !
+             ! compute contribution to y-gradient of y-diffusion gradient in centroid dual
+             !
+             dgydy = dgydy + ( dgy0 + dgy1 ) * real( x0 - x1 )
+             !
+          enddo
+          !
+          ! add diffusion terms to the right-hand side of the action balance equation
+          !
+          if ( carea > 0d0 ) then
+             !
+             dgxdx = dgxdx/real(carea)
+             dgydy = dgydy/real(carea)
+             !
+             ! in case of spherical coordinates, transform back to Cartesian coordinates
+             !
+             if ( KSPHER > 0 ) then
+                !
+                dgxdx = dgxdx/(cslat * LENDEG)
+                dgydy = dgydy/LENDEG
+                !
+             endif
+             !
+             rhs  (id,is  ) = rhs  (id,is  ) + dgxdx + dgydy
+             trac0(id,is,1) = trac0(id,is,1) - dgxdx - dgydy
+             !
+          endif
+          !
+       enddo
+       !
+    enddo
+    !
+end subroutine SwanGSECorr
--- swanhcat.ftn	2008-05-25 14:59:50.000000000 +0200
+++ swanhcat.ftn	2009-03-03 17:50:34.000000000 +0100
@@ -16,7 +16,7 @@
 !     hcat.nml file.
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -45,6 +45,7 @@
       logical ydivide, verbose, haloflag, exists, stomp, nonstat
       real, allocatable :: x(:,:),y(:,:)
       integer, allocatable :: locations(:),index(:), source(:)
+      integer, allocatable :: ownlocat(:)
       namelist /arg_nml/ basefile,halo,verbose,stomp
 
  33   FORMAT('-',I3.3,A)
@@ -159,7 +160,7 @@
          ENDIF
       ENDDO
 
-!.....Specifiy the number of location for each hotfile including
+!.....Specify the number of location for each hotfile including
 !.....the one being created at location(0)
       nonstat = .FALSE. !nonstationary unless time is found
       allocate(locations(0:numfiles))
@@ -297,11 +298,22 @@
       ENDDO
       locations(0)=index(0)
 
+!.....Specify the number of own location for each hotfile
+      allocate(ownlocat(0:numfiles))
+      ownlocat=0
+      DO j=1,numfiles
+         index(0)=1
+         DO WHILE (index(0)<=locations(0))
+            IF (source(index(0))==j) ownlocat(j)=ownlocat(j)+1
+            index(0)=index(0)+1
+         END DO
+      END DO
+
 !.....Write total number of locations to basefile
       WRITE (10,37) locations(0),TRIM(RLINE)
 
 !.....Write grid points to basefile
-      DO i=1,index(0)
+      DO i=1,locations(0)
          WRITE (10,36) x(0,i),y(0,i)
       ENDDO
 
@@ -435,7 +447,7 @@
                IF(source(index(0))==j) THEN
 !.................IF source equals 1 there are no halo points at the beginning
                   IF(j .NE. 1) THEN
-                     DO WHILE (index(j-1)<=locations(j-1))
+                     DO WHILE (index(j-1)<=ownlocat(j-1))
                         READ (10+j,35,IOSTAT=ios) RLINE
                         IF (RLINE(1:6).EQ.'FACTOR' .and. ios==0)THEN
                            DO l=1,freq+1
--- swanimp.tex	2008-05-25 14:59:38.000000000 +0200
+++ swanimp.tex	2009-03-03 17:50:34.000000000 +0100
@@ -15,7 +15,7 @@
 \end{center}
 \vfill
 \begin{center}
-{\Large\bf SWAN Cycle III version 40.72}
+{\Large\bf SWAN Cycle III version 40.72AB}
 \end{center}
 
 \cleardoublepage
@@ -37,14 +37,14 @@
              & & The Netherlands \\
              & & \\
 e-mail       &:& swan-info-citg@tudelft.nl \\
-home page    &:& \hl{http://www.fluidmechanics.tudelft.nl/swan/index.htm}
+home page    &:& \hl{http://www.swan.tudelft.nl}
 \end{tabular}
 \end{table}
 
 \vfill
 
 \noindent
-Copyright (c) 2008 Delft University of Technology.
+Copyright (c) 2009 Delft University of Technology.
 \\[2ex]
 \noindent
 Permission is granted to copy, distribute and/or modify this document
@@ -80,10 +80,10 @@
     \item some test cases.
 \end{itemize}
 All of the material can be found on the following SWAN web page\\
-\hl{http://www.fluidmechanics.tudelft.nl/swan/download/info.htm}.
+\hl{http://www.swan.tudelft.nl}.
 \\[2ex]
 \noindent
-On the SWAN home page \hl{http://www.fluidmechanics.tudelft.nl/swan/index.htm},
+On the SWAN home page \hl{http://www.swan.tudelft.nl},
 general information is given about the functionalities, physics and limitations
 of SWAN. Moreover, the modification history of SWAN is given. Finally,
 information on support, links to the related web pages and various
@@ -168,6 +168,8 @@
                           \> $\,\,$ SwanPropvelS.ftn90 \\
                           \> $\,\,$ SwanTranspAc.ftn90 \\
                           \> $\,\,$ SwanTranspX.ftn90 \\
+                          \> $\,\,$ SwanDiffPar.ftn90 \\
+                          \> $\,\,$ SwanGSECorr.ftn90 \\
                           \> $\,\,$ SwanInterpolatePoint.ftn90 \\
                           \> $\,\,$ SwanInterpolateAc.ftn90 \\
                           \> $\,\,$ SwanInterpolateOutput.ftn90 \\
@@ -436,7 +438,7 @@
   \item If desirable, you may clean-up the generated object files and modules by
         typing {\tt make~clean}. If you want to go back to the original state
         with respect to the source code, i.e. removing everything that has been
-        generated by the Makefile, just type {\tt make~allclean}.
+        generated by the Makefile, just type {\tt make~clobber}.
   \item If you are unable to install SWAN using the Makefile and Perl scripts for
         whatever reason, see Section \ref{sec:manual} for instructions on manual
         installation.
@@ -771,7 +773,7 @@
 a (parallel MPI) run can be either from a single (concatenated) hotfile or from multiple hotfiles.
 In the latter case, the number of processors must be equal to the number of generated hotfiles. If
 appropriate, the single hotfile can be created from a set of multiple hotfiles using the program
-{\tt hcat.exe} as available from SWAN version 40.72AB. This executable is generated from the Fortran
+{\tt hcat.exe} as available from SWAN version 40.51A. This executable is generated from the Fortran
 program {\tt swanhcat.ftn}. A self-contained input file {\tt hcat.nml} is provided. This file contains,
 e.g. the (basis) name of the hotfile. To concatenate multiple hotfiles into single hotfile just execute
 {\tt hcat.exe}.
--- SwanInitCompGrid.ftn90	2008-05-25 14:59:44.000000000 +0200
+++ SwanInitCompGrid.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanInterpolateAc.ftn90	2009-03-03 17:49:37.000000000 +0100
+++ SwanInterpolateAc.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanInterpolateOutput.ftn90	2009-03-03 17:49:37.000000000 +0100
+++ SwanInterpolateOutput.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanInterpolatePoint.ftn90	2008-05-25 14:59:45.000000000 +0200
+++ SwanInterpolatePoint.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanIntgratSpc.ftn90	2008-05-25 14:59:45.000000000 +0200
+++ SwanIntgratSpc.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- swanmain.ftn	2009-03-03 17:49:37.000000000 +0100
+++ swanmain.ftn	2009-03-03 17:50:34.000000000 +0100
@@ -44,7 +44,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -150,7 +150,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -418,8 +418,11 @@
 !       --- initialisation of necessary grids for depth,
 !           current, wind and friction
 
-        IF(.NOT.ALLOCATED(COMPDA))
-     &                        ALLOCATE(COMPDA(MCGRD,MCMVAR),STAT=ISTAT)   40.41 40.31
+        IF (ALOCMP.AND.ALLOCATED(COMPDA)) DEALLOCATE(COMPDA)              40.97
+        IF (.NOT.ALLOCATED(COMPDA)) THEN                                  40.97
+           ALLOCATE(COMPDA(MCGRD,MCMVAR),STAT=ISTAT)                      40.97 40.41 40.31
+           ALOCMP = .FALSE.                                               40.97
+        END IF                                                            40.97
         IF ( ISTAT.NE.0 ) THEN                                            40.41
            CHARS(1) = NUMSTR(ISTAT,RNAN,'(I6)')                           40.41
            CALL TXPBLA(CHARS(1),IF1,IL1)                                  40.41
@@ -593,7 +596,7 @@
 
 !TIMG            CALL SWTSTA(9)                                                40.30
             CALL SWOUTP ( AC2   , SPCSIG, SPCDIR, COMPDA, XYTST ,         40.31
-     &                    KGRPNT, XCGRID, YCGRID, KGRBND, OURQT )         40.51 40.31
+     &                    KGRPNT, XCGRID, YCGRID, OURQT )                 40.51 40.31
 !TIMG            CALL SWTSTO(9)                                                40.30
             IF (STPNOW()) RETURN                                          40.30
 !
@@ -634,7 +637,10 @@
       IF ( PARLL ) THEN                                                   40.30
          ALLOCATE (BLKND(MXC*MYC))                                        40.51 40.41
          BLKND = REAL(INODE)                                              40.41
-         IF ( INODE.EQ.MASTER ) ALLOCATE(BLKNDC(MXCGL*MYCGL))             40.51 40.41
+         IF ( INODE.EQ.MASTER ) THEN
+            ALLOCATE(BLKNDC(MXCGL*MYCGL))                                 40.51 40.41
+            BLKNDC = 0.                                                   40.96
+         END IF
          CALL SWCOLLECT ( BLKNDC, BLKND, .TRUE. )                         40.51 40.41
          IF (STPNOW()) RETURN                                             40.41
          IF ( INODE.EQ.MASTER ) THEN                                      40.30
@@ -695,7 +701,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -850,6 +856,7 @@
       VERNUM = 40.72
       WRITE (VERTXT, '(F5.2)') VERNUM                                     40.03
       CALL BUGFIX ('A')
+      CALL BUGFIX ('B')
 !
       CALL OCPINI ('swaninit', .TRUE.,INERR)                              34.01
       IF (INERR.GT.0) RETURN                                              34.01
@@ -2219,7 +2226,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2545,7 +2552,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2928,7 +2935,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -3322,7 +3329,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -3725,7 +3732,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -3854,7 +3861,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -3983,7 +3990,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -4150,7 +4157,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -4249,7 +4256,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -4380,6 +4387,11 @@
          ENDIF                                                            40.80
       ENDIF                                                               40.80
 !
+!     --- reset PWIND(9) and PWIND(17) if RHO has a user value            40.96
+!
+      PWIND( 9) = PWIND(16)/RHO                                           40.96
+      PWIND(17) = RHO                                                     40.96
+!
 !     -----------------------------------------------------------------
 !
 !     *** WARNINGS AND ERROR MESSAGES ***
@@ -4600,7 +4612,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -5128,7 +5140,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -5723,7 +5735,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -6146,7 +6158,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -6467,7 +6479,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -6763,7 +6775,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- swanout1.ftn	2009-03-03 17:49:37.000000000 +0100
+++ swanout1.ftn	2009-03-03 17:50:34.000000000 +0100
@@ -20,8 +20,7 @@
      &                   SPCSIG          ,SPCDIR  ,                       30.72
      &                   COMPDA          ,XYTST   ,
      &                   KGRPNT          ,XCGRID  ,                       30.72
-     &                   YCGRID          ,KGRBND  ,                       40.30
-     &                   OURQT                    )                       40.51 40.30
+     &                   YCGRID          ,OURQT   )                       40.51 40.30
 !                                                                      *
 !***********************************************************************
 !
@@ -48,7 +47,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -231,8 +230,7 @@
 ! 13. Source text
 !
       INTEGER   VOQR(NMOVAR)   ,BKC           ,                           40.31
-     &          XYTST(*)       ,KGRPNT(MXC,MYC),IERR          ,           30.81 30.21
-     &          KGRBND(*)                                                 40.00
+     &          XYTST(*)       ,KGRPNT(MXC,MYC),IERR                      30.81 30.21
 !
       REAL      AC2(MDC,MSC,MCGRD) ,
      &          COMPDA(MCGRD,MCMVAR)
@@ -339,7 +337,7 @@
      &               VOQ(1+MIP)          ,VOQ(1+2*MIP)        ,           40.31 30.90
      &               VOQ(1+3*MIP)        ,KGRPNT              ,           40.31 30.90
      &               XCGRID              ,YCGRID              ,           30.21
-     &               KGRBND              ,CROSS               )           40.86 40.00
+     &               CROSS                                    )           40.86 40.00
 !
 !       Compute wave-induced force on unstructured grid                   40.80
 !
@@ -483,7 +481,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -784,7 +782,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1022,7 +1020,7 @@
      &                    YP        ,XC         ,
      &                    YC        ,KGRPNT     ,
      &                    XCGRID    ,YCGRID     ,                         30.21
-     &                    KGRBND    ,CROSS      )                         40.86 40.00
+     &                    CROSS                 )                         40.86 40.00
 !                                                                      *
 !***********************************************************************
 !
@@ -1047,7 +1045,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1103,7 +1101,6 @@
 !
 !  4. Argument variables
 !
-      INTEGER, INTENT(INOUT) :: KGRBND(*)                                 40.02
       INTEGER OPI(2)                                                      40.31
 !
 !     XCGRID: input  Coordinates of computational grid in x-direction     30.72
@@ -1370,7 +1367,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2421,7 +2418,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2684,7 +2681,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -4185,7 +4182,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -4405,7 +4402,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- swanout2.ftn	2009-03-03 17:49:37.000000000 +0100
+++ swanout2.ftn	2009-03-03 17:50:34.000000000 +0100
@@ -36,7 +36,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -309,7 +309,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -535,7 +535,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -683,7 +683,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -837,8 +837,8 @@
 !               Else if Rtype = 'TABD'
 !                    Then Format = '(E12.4)'
 !                         make lfield = 12
-!                    Else Format = '(F11.X)'
-!                         make lfield = 11
+!                    Else Format = '(F13.X)'
+!                         make lfield = 13
 !                         determine number of decimals and write into
 !                         Format
 !                    ---------------------------------------------------
@@ -891,15 +891,15 @@
               OUTLIN(1:9) = OUT_COMMENT // '        '                     40.13
               LINKAR = 9
             ELSE
-              OUTLIN(1:1) = OUT_COMMENT                                   40.13
-              LINKAR = 2                                                  40.13
+              OUTLIN(1:3) = OUT_COMMENT                                   40.13
+              LINKAR = 4                                                  40.13
             ENDIF
             DO  JVAR = 1, NVAR
               IVTYPE = IVTYP(JVAR)                                        40.31
               IF (IVTYPE.EQ.40) THEN
                 LFIELD = 18
               ELSE
-                LFIELD = 11
+                LFIELD = 13
               ENDIF
               IF (OVSVTY(IVTYPE).LE.2) THEN
                 OUTLIN(LINKAR:LINKAR+LFIELD) =
@@ -920,15 +920,15 @@
               OUTLIN(1:9) = OUT_COMMENT // '        '                     40.13
               LINKAR = 9
             ELSE
-              OUTLIN(1:1) = OUT_COMMENT                                   40.13
-              LINKAR = 2                                                  40.13
+              OUTLIN(1:3) = OUT_COMMENT                                   40.13
+              LINKAR = 4                                                  40.13
             ENDIF
             DO  JVAR = 1, NVAR
               IVTYPE = IVTYP(JVAR)                                        40.31
               IF (IVTYPE.EQ.40) THEN
                 LFIELD = 18
               ELSE
-                LFIELD = 11
+                LFIELD = 13
               ENDIF
               DO ISTR = LEN(OVUNIT(IVTYPE)), 1, -1
                 IF (OVUNIT(IVTYPE)(ISTR:ISTR) .NE. ' ') THEN
@@ -1032,8 +1032,8 @@
               FSTR = FLT_TABLE                                            40.13
               LFIELD = FLD_TABLE                                          40.13
             ELSE
-              FSTR = '(F11.X)'
-              LFIELD = 11
+              FSTR = '(F13.X)'
+              LFIELD = 13
 !             NUMDEC is number of decimals in the table for each output quantity
               NUMDEC = MAX (0, 6-NINT(LOG10(ABS(OVHEXP(IVTYPE)))))        40.03
               IF (NUMDEC.GT.9) NUMDEC = 9                                 40.00
@@ -1076,7 +1076,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1197,7 +1197,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1541,7 +1541,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1838,7 +1838,7 @@
  800  IERR = 1
 !
  900  CONTINUE
-      IF (OPTG.EQ.5) DEALLOCATE(ACL)                                      40.80
+      DEALLOCATE(ACL)                                                     40.90 40.80
       RETURN
 ! * end of subroutine SWCMSP *
       END
@@ -1865,7 +1865,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- swanparll.ftn	2009-03-03 17:49:37.000000000 +0100
+++ swanparll.ftn	2009-03-03 17:50:34.000000000 +0100
@@ -11,6 +11,8 @@
 !     SWBROADC
 !     SWGATHER
 !     SWREDUCE
+!     SWREDUCI
+!     SWREDUCR
 !     SWSTRIP
 !     SWPARTIT
 !     SWBLADM
@@ -48,7 +50,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -222,7 +224,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -343,7 +345,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -470,7 +472,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -605,7 +607,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -746,7 +748,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -880,7 +882,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1054,7 +1056,6 @@
 !
 !****************************************************************
 !
-!MPI      USE MPI
       USE OCPCOMM4                                                        40.41
       USE M_PARALL                                                        40.31
 !
@@ -1072,7 +1073,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1094,11 +1095,13 @@
 !
 !     40.30: Marcel Zijlema
 !     40.41: Marcel Zijlema
+!     40.96: Marcel Zijlema
 !
 !  1. Updates
 !
 !     40.30, Feb. 03: New subroutine
 !     40.41, Oct. 04: common blocks replaced by modules, include files removed
+!     40.96, Dec. 08: call SWREDUCI/R instead of passing startaddress of the array
 !
 !  2. Purpose
 !
@@ -1106,10 +1109,6 @@
 !     array (IPTR) of type ITYPE to collect values from
 !     all processes
 !
-!  3. Method
-!
-!     Wrapper for MPI_ALLREDUCE
-!
 !  4. Argument variables
 !
 !     ILEN        length of array to be collect
@@ -1125,17 +1124,131 @@
 !
 !  6. Local variables
 !
+!     IENT  :     number of entries
+!
+      INTEGER IENT
+!
+!  8. Subroutines used
+!
+!     STRACE           Tracing routine for debugging
+!
+!  9. Subroutines calling
+!
+!     SWCOMP
+!
+! 10. Error messages
+!
+!     ---
+!
+! 11. Remarks
+!
+!     ---
+!
+! 12. Structure
+!
+!     Performs a global reduction of data across all nodes
+!
+! 13. Source text
+!
+      SAVE IENT
+      DATA IENT/0/
+      IF (LTRACE) CALL STRACE (IENT,'SWREDUCE')
+
+!     --- if not parallel, return
+      IF (.NOT.PARLL) RETURN
+!
+!     --- actual reduction of field array based on its type
+      IF ( ITYPE.EQ.SWINT ) THEN
+         CALL SWREDUCI ( IPTR, ILEN, ITYPRD )
+      ELSE IF ( ITYPE.EQ.SWREAL ) THEN
+         CALL SWREDUCR ( IPTR, ILEN, ITYPRD )
+      END IF
+
+      RETURN
+      END
+!****************************************************************
+!
+      SUBROUTINE SWREDUCI ( IARR, ILEN, ITYPRD )
+!
+!****************************************************************
+!
+!MPI      USE MPI
+      USE OCPCOMM4
+      USE M_PARALL
+!
+      IMPLICIT NONE
+!
+!
+!   --|-----------------------------------------------------------|--
+!     | Delft University of Technology                            |
+!     | Faculty of Civil Engineering and Geosciences              |
+!     | Environmental Fluid Mechanics Section                     |
+!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
+!     |                                                           |
+!     | Programmer: Marcel Zijlema                                |
+!   --|-----------------------------------------------------------|--
+!
+!
+!     SWAN (Simulating WAves Nearshore); a third generation wave model
+!     Copyright (C) 2009  Delft University of Technology
+!
+!     This program is free software; you can redistribute it and/or
+!     modify it under the terms of the GNU General Public License as
+!     published by the Free Software Foundation; either version 2 of
+!     the License, or (at your option) any later version.
+!
+!     This program is distributed in the hope that it will be useful,
+!     but WITHOUT ANY WARRANTY; without even the implied warranty of
+!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
+!     GNU General Public License for more details.
+!
+!     A copy of the GNU General Public License is available at
+!     http://www.gnu.org/copyleft/gpl.html#SEC3
+!     or by writing to the Free Software Foundation, Inc.,
+!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
+!
+!
+!  0. Authors
+!
+!     40.96: Marcel Zijlema
+!
+!  1. Updates
+!
+!     40.96, Dec. 08: New subroutine
+!
+!  2. Purpose
+!
+!     Performs a global reduction of type ITYPRD on integer
+!     array IARR to collect values from all processes
+!
+!  3. Method
+!
+!     Wrapper for MPI_ALLREDUCE
+!
+!  4. Argument variables
+!
+!     IARR        integer array
+!     ILEN        length of array to be collect
+!     ITYPRD      type of reduction
+!
+      INTEGER ILEN, ITYPRD
+      INTEGER IARR(ILEN)
+!
+!  6. Local variables
+!
 !     CHARS :     character for passing info to MSGERR
 !     IENT  :     number of entries
 !     IERR  :     error value of MPI call
 !     IF    :     first non-character in string
 !     IL    :     last non-character in string
-!     IOPTR :     pointer to first element of output array
+!     ITEMP :     temporary array to store collected data
 !     MSGSTR:     string to pass message to call MSGERR
 !
-      INTEGER      IENT, IERR, IF, IL, IOPTR
+      INTEGER      IENT, IERR, IF, IL
       CHARACTER*20 INTSTR, CHARS
       CHARACTER*80 MSGSTR
+
+      INTEGER, ALLOCATABLE :: ITEMP(:)
 !
 !  8. Subroutines used
 !
@@ -1144,15 +1257,13 @@
 !MPI!                      distribute the result back to all processes
 !     MSGERR           Writes error message
 !     STRACE           Tracing routine for debugging
-!     SWCOPI           Makes copy of an integer array
-!     SWCOPR           Makes copy of a real array
 !TIMG!     SWTSTA           Start timing for a section of code
 !TIMG!     SWTSTO           Stop timing for a section of code
 !     TXPBLA           Removes leading and trailing blanks in string
 !
 !  9. Subroutines calling
 !
-!     SWCOMP
+!     SWREDUCE
 !
 ! 10. Error messages
 !
@@ -1164,20 +1275,19 @@
 !
 ! 12. Structure
 !
-!MPI!     Performs a global reduction of data across all nodes
+!MPI!     Performs a global reduction of integers across all nodes
 !MPI!     with command MPI_ALLREDUCE
 !
 ! 13. Source text
 !
       SAVE IENT
       DATA IENT/0/
-      IF (LTRACE) CALL STRACE (IENT,'SWREDUCE')
+      IF (LTRACE) CALL STRACE (IENT,'SWREDUCI')
 
-!     --- if not parallel, return
-      IF (.NOT.PARLL) RETURN
+      ALLOCATE(ITEMP(ILEN))
 
 !TIMG      CALL SWTSTA(202)
-!MPI      CALL MPI_ALLREDUCE ( IPTR, IOPTR, ILEN, ITYPE,
+!MPI      CALL MPI_ALLREDUCE ( IARR, ITEMP, ILEN, SWINT,
 !MPI     &                     ITYPRD, MPI_COMM_WORLD, IERR )
 !MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
 !MPI         CHARS = INTSTR(IERR)
@@ -1187,13 +1297,149 @@
 !MPI         CALL MSGERR ( 4, MSGSTR )
 !MPI         RETURN
 !MPI      END IF
-      IF ( ITYPE.EQ.SWINT ) THEN
-         CALL SWCOPI ( IOPTR, IPTR, ILEN )
-      ELSE IF ( ITYPE.EQ.SWREAL ) THEN
-         CALL SWCOPR ( IOPTR, IPTR, ILEN )
-      END IF
+      IARR = ITEMP
 !TIMG      CALL SWTSTO(202)
 
+      DEALLOCATE(ITEMP)
+
+      RETURN
+      END
+!****************************************************************
+!
+      SUBROUTINE SWREDUCR ( ARR, ILEN, ITYPRD )
+!
+!****************************************************************
+!
+!MPI      USE MPI
+      USE OCPCOMM4
+      USE M_PARALL
+!
+      IMPLICIT NONE
+!
+!
+!   --|-----------------------------------------------------------|--
+!     | Delft University of Technology                            |
+!     | Faculty of Civil Engineering and Geosciences              |
+!     | Environmental Fluid Mechanics Section                     |
+!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
+!     |                                                           |
+!     | Programmer: Marcel Zijlema                                |
+!   --|-----------------------------------------------------------|--
+!
+!
+!     SWAN (Simulating WAves Nearshore); a third generation wave model
+!     Copyright (C) 2009  Delft University of Technology
+!
+!     This program is free software; you can redistribute it and/or
+!     modify it under the terms of the GNU General Public License as
+!     published by the Free Software Foundation; either version 2 of
+!     the License, or (at your option) any later version.
+!
+!     This program is distributed in the hope that it will be useful,
+!     but WITHOUT ANY WARRANTY; without even the implied warranty of
+!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
+!     GNU General Public License for more details.
+!
+!     A copy of the GNU General Public License is available at
+!     http://www.gnu.org/copyleft/gpl.html#SEC3
+!     or by writing to the Free Software Foundation, Inc.,
+!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
+!
+!
+!  0. Authors
+!
+!     40.96: Marcel Zijlema
+!
+!  1. Updates
+!
+!     40.96, Dec. 08: New subroutine
+!
+!  2. Purpose
+!
+!     Performs a global reduction of type ITYPRD on real
+!     array ARR to collect values from all processes
+!
+!  3. Method
+!
+!     Wrapper for MPI_ALLREDUCE
+!
+!  4. Argument variables
+!
+!     ARR         real array
+!     ILEN        length of array to be collect
+!     ITYPRD      type of reduction
+!
+      INTEGER ILEN, ITYPRD
+      REAL    ARR(ILEN)
+!
+!  6. Local variables
+!
+!     CHARS :     character for passing info to MSGERR
+!     IENT  :     number of entries
+!     IERR  :     error value of MPI call
+!     IF    :     first non-character in string
+!     IL    :     last non-character in string
+!     MSGSTR:     string to pass message to call MSGERR
+!     TEMP  :     temporary array to store collected data
+!
+      INTEGER      IENT, IERR, IF, IL
+      CHARACTER*20 INTSTR, CHARS
+      CHARACTER*80 MSGSTR
+
+      REAL, ALLOCATABLE :: TEMP(:)
+!
+!  8. Subroutines used
+!
+!     INTSTR           Converts integer to string
+!MPI!     MPI_ALLREDUCE    Combines values from all processes and
+!MPI!                      distribute the result back to all processes
+!     MSGERR           Writes error message
+!     STRACE           Tracing routine for debugging
+!TIMG!     SWTSTA           Start timing for a section of code
+!TIMG!     SWTSTO           Stop timing for a section of code
+!     TXPBLA           Removes leading and trailing blanks in string
+!
+!  9. Subroutines calling
+!
+!     SWREDUCE
+!
+! 10. Error messages
+!
+!     ---
+!
+! 11. Remarks
+!
+!     ---
+!
+! 12. Structure
+!
+!MPI!     Performs a global reduction of reals across all nodes
+!MPI!     with command MPI_ALLREDUCE
+!
+! 13. Source text
+!
+      SAVE IENT
+      DATA IENT/0/
+      IF (LTRACE) CALL STRACE (IENT,'SWREDUCR')
+
+      ALLOCATE(TEMP(ILEN))
+
+!TIMG      CALL SWTSTA(202)
+!MPI      CALL MPI_ALLREDUCE ( ARR, TEMP, ILEN, SWREAL,
+!MPI     &                     ITYPRD, MPI_COMM_WORLD, IERR )
+!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
+!MPI         CHARS = INTSTR(IERR)
+!MPI         CALL TXPBLA(CHARS,IF,IL)
+!MPI         MSGSTR = 'MPI produces some internal error - '//
+!MPI     &            'return code is '//CHARS(IF:IL)
+!MPI         CALL MSGERR ( 4, MSGSTR )
+!MPI         RETURN
+!MPI      END IF
+      ARR = TEMP
+!TIMG      CALL SWTSTO(202)
+
+      DEALLOCATE(TEMP)
+
       RETURN
       END
 !****************************************************************
@@ -1218,7 +1464,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1270,8 +1516,9 @@
 !     MYC         maximum counter of gridpoints in y-direction
 !     NPART       number of parts to be created
 !
-      INTEGER IDIR, MXC, MYC, NPART
-      INTEGER IPOWN(*), IWORK(2,*)
+      INTEGER   IDIR, MXC, MYC, NPART
+      INTEGER   IPOWN(*)
+      INTEGER*8 IWORK(2,*)
 !
 !  6. Local variables
 !
@@ -1420,7 +1667,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1482,8 +1729,9 @@
 !     NACTP :     total number of active gridpoints
 !     NPCUM :     cumulative number of gridpoints
 !
-      INTEGER   I, ICNT, IDIR, IENT, IX, IY, NACTP, NPCUM
-      INTEGER   IWORK(2,NPROC)
+      INTEGER   I, IDIR, IENT, IX, IY
+      INTEGER*8 ICNT, NACTP, NPCUM
+      INTEGER*8 IWORK(2,NPROC)
 !
 !  8. Subroutines used
 !
@@ -1586,7 +1834,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1928,7 +2176,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2095,7 +2343,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2285,7 +2533,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2451,7 +2699,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2618,7 +2866,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2794,25 +3042,48 @@
             IOFF2 = 0
 
             DO IP = 0, NPROC-1
-               IF ( IARRC(1,IP).EQ.1 ) THEN
-                  MXFGL = 1
-               ELSE
-                  MXFGL = IARRC(1,IP) + IHALOX
-               END IF
-               IF ( IARRC(3,IP).EQ.1 ) THEN
-                  MYFGL = 1
-               ELSE
-                  MYFGL = IARRC(3,IP) + IHALOY
-               END IF
-               IF ( IARRC(2,IP).EQ.MXCGL ) THEN
-                  MXLGL = MXCGL
-               ELSE
-                  MXLGL = IARRC(2,IP) - IHALOX
-               END IF
-               IF ( IARRC(4,IP).EQ.MYCGL ) THEN
-                  MYLGL = MYCGL
+               IF ( MXCGL.GT.MYCGL ) THEN
+                  IF ( IARRC(1,IP).EQ.1 .AND. IP.EQ.0 ) THEN
+                     MXFGL = 1
+                  ELSE
+                     MXFGL = IARRC(1,IP) + IHALOX
+                  END IF
+                  IF ( IARRC(3,IP).EQ.1 ) THEN
+                     MYFGL = 1
+                  ELSE
+                     MYFGL = IARRC(3,IP) + IHALOY
+                  END IF
+                  IF ( IARRC(2,IP).EQ.MXCGL .AND. IP.EQ.NPROC-1 ) THEN
+                     MXLGL = MXCGL
+                  ELSE
+                     MXLGL = IARRC(2,IP) - IHALOX
+                  END IF
+                  IF ( IARRC(4,IP).EQ.MYCGL ) THEN
+                     MYLGL = MYCGL
+                  ELSE
+                     MYLGL = IARRC(4,IP) - IHALOY
+                  END IF
                ELSE
-                  MYLGL = IARRC(4,IP) - IHALOY
+                  IF ( IARRC(1,IP).EQ.1 ) THEN
+                     MXFGL = 1
+                  ELSE
+                     MXFGL = IARRC(1,IP) + IHALOX
+                  END IF
+                  IF ( IARRC(3,IP).EQ.1 .AND. IP.EQ.0 ) THEN
+                     MYFGL = 1
+                  ELSE
+                     MYFGL = IARRC(3,IP) + IHALOY
+                  END IF
+                  IF ( IARRC(2,IP).EQ.MXCGL ) THEN
+                     MXLGL = MXCGL
+                  ELSE
+                     MXLGL = IARRC(2,IP) - IHALOX
+                  END IF
+                  IF ( IARRC(4,IP).EQ.MYCGL .AND. IP.EQ.NPROC-1 ) THEN
+                     MYLGL = MYCGL
+                  ELSE
+                     MYLGL = IARRC(4,IP) - IHALOY
+                  END IF
                END IF
 
                ILEN = IARRC(2,IP)-IARRC(1,IP)+1
@@ -2879,7 +3150,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -3130,7 +3401,7 @@
      &                      YP                  ,XC                  ,
      &                      YC                  ,KGRPGL              ,
      &                      XGRDGL              ,YGRDGL              ,
-     &                      KGRBGL              ,CROSS               )    40.86
+     &                      CROSS               )                         40.86
                MXC    = ITMP1
                MYC    = ITMP2
                MCGRD  = ITMP3
@@ -3237,7 +3508,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -3567,7 +3838,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -3894,7 +4165,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- swanpgr.tex	2008-05-25 14:59:38.000000000 +0200
+++ swanpgr.tex	2009-03-03 17:50:34.000000000 +0100
@@ -37,14 +37,14 @@
              & & The Netherlands \\
              & & \\
 e-mail       &:& swan-info-citg@tudelft.nl \\
-home page    &:& \hl{http://www.fluidmechanics.tudelft.nl/swan/index.htm}
+home page    &:& \hl{http://www.swan.tudelft.nl}
 \end{tabular}
 \end{table}
 
 \vfill
 
 \noindent
-Copyright (c) 2008 Delft University of Technology.
+Copyright (c) 2009 Delft University of Technology.
 \\[2ex]
 \noindent
 Permission is granted to copy, distribute and/or modify this document
@@ -687,7 +687,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2005  Delft University of Technology
+!     Copyright (C) 2002  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -790,7 +790,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2007  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -977,7 +977,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2007  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanPointinMesh.ftn90	2008-05-25 14:59:46.000000000 +0200
+++ SwanPointinMesh.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- swanpre1.ftn	2009-03-03 17:49:37.000000000 +0100
+++ swanpre1.ftn	2009-03-03 17:50:34.000000000 +0100
@@ -50,7 +50,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -730,7 +730,7 @@
           PROPSN = 1                                                      40.02
           PROPSS = 1                                                      40.02
         ELSE IF (KEYWIS ('GSE')) THEN                                     40.02
-          IF (PROPSN.NE.3) THEN                                           40.02
+          IF (OPTG.NE.5 .AND. PROPSN.NE.3) THEN                           41.00 40.02
             CALL MSGERR(2,                                                40.02
      &      'Anti-GSE only allowed for S&L scheme.')                      40.02
           ENDIF                                                           40.02
@@ -1753,6 +1753,7 @@
         JSETUP = MCMVAR + 1                                               32.02
         JDPSAV = MCMVAR + 2                                               32.02
         MCMVAR = MCMVAR + 2                                               40.41 32.02
+        ALOCMP = .TRUE.                                                   40.97
         GOTO 100                                                          32.02
       ENDIF                                                               32.02
 !
@@ -1767,11 +1768,6 @@
 ! =============================================================           40.21
 !                                                                         40.21
       IF (KEYWIS ('DIFFR')) THEN                                          40.21
-        IF (OPTG.EQ.5) THEN                                               40.80
-           CALL MSGERR(2, ' No diffraction for unstructured grid')        40.80
-           IDIFFR = 0                                                     40.80
-           GOTO 100                                                       40.80
-        ENDIF                                                             40.80
         CALL ININTG ('IDIFFR', IDIFFR, 'STA', 1)                          40.21
         CALL INREAL ('SMPAR', PDIFFR(1), 'STA', 0.0)                      40.21
         CALL INREAL ('SMNUM', PDIFFR(2), 'STA', 0.)                       40.21
@@ -1952,8 +1948,14 @@
 ! ============================================================
 !
 !           | -> CONstant [alpha] [gamma]                                      |
-! BREaking <                                                                    >
+!           |                                                                  |
 !           |    VARiable [alpha] [gammin] [gammax] [gamneg] [coeff1] [coeff2] |
+!           |                                                                  |
+! BREaking <     RUEssink [alpha] [a] [b]                                       >
+!           |                                                                  |
+!           |    TG       [alpha] [gamma] [pown]                               |
+!           |                                                                  |
+!           !    BIPhase  [alpha] [pown] [bref]                                |
 !
 ! ============================================================
 !
@@ -1963,7 +1965,7 @@
           ISURF = 1
           CALL INREAL ('ALPHA', PSURF(1), 'STA', 1.0)
           CALL INREAL ('GAMMA', PSURF(2), 'STA', 0.73)
-        ELSE IF (KEYWIS('VAR')) THEN                                      28/FEB
+        ELSE IF (KEYWIS('VAR') .OR. KEYWIS('NEL')) THEN                   41.03
           ISURF = 2
           CALL INREAL ('ALPHA',  PSURF(1), 'STA', 1.5)
           CALL INREAL ('GAMMIN', PSURF(4), 'STA', 0.55)
@@ -1971,6 +1973,21 @@
           CALL INREAL ('GAMNEG', PSURF(6), 'STA', 0.73)
           CALL INREAL ('COEFF1', PSURF(7), 'STA', 0.88)
           CALL INREAL ('COEFF2', PSURF(8), 'STA', 0.012)
+        ELSE IF (KEYWIS('RUE')) THEN                                      41.03
+          ISURF = 3
+          CALL INREAL ('ALPHA', PSURF(1), 'STA', 1.0)
+          CALL INREAL ('A'    , PSURF(4), 'STA', 0.76)
+          CALL INREAL ('B'    , PSURF(5), 'STA', 0.29)
+        ELSE IF (KEYWIS('TG')) THEN                                       41.03
+          ISURF = 4
+          CALL INREAL ('ALPHA',  PSURF(1), 'STA', 1.0)
+          CALL INREAL ('GAMMA',  PSURF(2), 'STA', 0.42)
+          CALL INREAL ('POWN' ,  PSURF(4), 'STA', 4.0)
+        ELSE IF (KEYWIS('BIP')) THEN                                      41.03
+          ISURF = 5
+          CALL INREAL ('ALPHA',  PSURF(1), 'STA', 1.0)
+          CALL INREAL ('POWN' ,  PSURF(4), 'STA', 2.4)
+          CALL INREAL ('BREF' ,  PSURF(5), 'STA', -1.5708)
         ENDIF
         GOTO 100
       ENDIF
@@ -2022,6 +2039,7 @@
           JCDRAG = MCMVAR+3                                               30.22
           JTAUW  = MCMVAR+4                                               30.22
           MCMVAR = MCMVAR+4                                               30.22
+          ALOCMP = .TRUE.                                                 40.97
           CALL INREAL ('PWTAIL', PWTAIL(1), 'STA', 5.)                    30.70
           IF (PWTAIL(1).LE.1.) CALL MSGERR (3, 'Incorrect PWTAIL')        30.70
           PWTAIL(3) = PWTAIL(1) + 1.                                      30.70
@@ -2062,17 +2080,19 @@
 !
 !     ------------------------------------------------------------------
 !
-!     FRICTION  setting parameters for bottom friction
+!     FRIC  setting parameters for bottom friction
 !
 ! ===============================================
 !
-!                   | -> JONSWAP       [cfjon]
+!                   |              | -> CONstant [cfjon]
+!                   | -> JONswap  <
+!                   |              |    VARiable [cfj1] [cfj2] [dsp1] [dsp2]   (NOT documented)
 !                   |
-!   FRICTION       <     COLLINS       [cfw]    &
+!                   |    COLLins   [cfw]    &
+!   FRICtion       <
+!                   |              [cfc]      (NOT documented)
 !                   |
-!                   |                  [cfc]    (NOT documented)
-!                   |
-!                   |    MADSEN        [kn]
+!                   |    MADsen    [kn]
 !
 ! ===============================================
 !
@@ -2080,8 +2100,18 @@
         IBOT = 1
         CALL INKEYW ('STA','JON')
         IF (KEYWIS('JON')) THEN
-          IBOT = 1
-          CALL INREAL('CFJON',PBOT(3),'UNC',0.)
+          CALL INKEYW ('STA', 'CON')
+          IF (KEYWIS('VAR')) THEN                                         41.04
+            IBOT = 4
+            CALL INREAL ('CFJ1', PBOT(6), 'STA', 0.038)
+            CALL INREAL ('CFJ2', PBOT(7), 'STA', 0.067)
+            CALL INREAL ('DSP1', PBOT(8), 'STA', 10.)
+            CALL INREAL ('DSP2', PBOT(9), 'STA', 30.)
+          ELSE
+            CALL IGNORE ('CON')
+            IBOT = 1
+            CALL INREAL('CFJON',PBOT(3),'UNC',0.)
+          ENDIF
         ELSE IF (KEYWIS('COLL')) THEN                                     20.68
           IBOT = 2
           CALL INREAL('CFW',PBOT(2),'UNC',0.)
@@ -2234,6 +2264,7 @@
           JCDRAG = MCMVAR+3                                               30.72
           JTAUW  = MCMVAR+4                                               30.72
           MCMVAR = MCMVAR+4                                               30.72
+          ALOCMP = .TRUE.                                                 40.97
 !         if [pwtail] is not changed, make it 5
           IF (EQREAL(PWTAIL(1),4.)) THEN                                  40.00
             PWTAIL(1) = 5.                                                40.00
@@ -2571,7 +2602,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2989,7 +3020,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -3139,6 +3170,7 @@
         MCMVAR = MCMVAR + 2                                               40.00
         JFRC2  = MCMVAR - 1                                               40.00
         JFRC3  = MCMVAR                                                   40.00
+        ALOCMP = .TRUE.                                                   40.97
       ELSE IF (KEYWIS ('WI')) THEN
         LWINDR = 2                                                        30.10
         IWIND  = LWINDM                                                   30.10
@@ -3166,6 +3198,7 @@
           MCMVAR = MCMVAR + 2                                             40.03
           JASTD2 = MCMVAR - 1                                             40.03
           JASTD3 = MCMVAR                                                 40.03
+          ALOCMP = .TRUE.                                                 40.97
         ENDIF                                                             40.03
       ELSE
         CALL  WRNKEY
@@ -3311,7 +3344,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -3478,7 +3511,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -3740,7 +3773,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -3960,7 +3993,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -4325,7 +4358,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -4389,8 +4422,8 @@
 ! 12. Structure
 !
 !     PVALID = .TRUE.
-!     If gridcounter in x-direction is less then 1 or greater than mxc or
-!     If gridcounter in y-direction is less then 1 or greater than myc or
+!     If gridcounter in x-direction is less then 1 or greater than MXC or
+!     If gridcounter in y-direction is less then 1 or greater than MYC or
 !     If gridpoint is exception point PVALID = .FALSE.
 !
 ! 13. Source text
@@ -4435,7 +4468,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -4607,7 +4640,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -4722,7 +4755,7 @@
 !
       DATA      IENT /0/
       CALL STRACE (IENT,'SEPARAREA')
-!      IF (MOD(IENT,2).EQ.0) THEN
+!
 !     In case there are 3 wet neighbouring points and validbp(ix,iy,.)
 !
       IF(.NOT.PVALID(IX-1,IY+1,KGRPNT) .AND.
@@ -4777,7 +4810,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -5245,7 +5278,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- swanpre2.ftn	2009-03-03 17:49:37.000000000 +0100
+++ swanpre2.ftn	2009-03-03 17:50:34.000000000 +0100
@@ -36,7 +36,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -202,7 +202,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1039,7 +1039,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1256,62 +1256,77 @@
                 IF (IVTYPE.EQ.50 .AND. JPBOT.LE.1) THEN                   40.65
                    MCMVAR = MCMVAR+1                                      40.65
                    JPBOT  = MCMVAR                                        40.65
+                   ALOCMP = .TRUE.                                        40.97
                 ENDIF                                                     40.65
                 IF (IVTYPE.EQ.52 .AND. JBOTLV.LE.1) THEN                  40.65
                    MCMVAR = MCMVAR+1                                      40.65
                    JBOTLV = MCMVAR                                        40.65
+                   ALOCMP = .TRUE.                                        40.97
                 ENDIF                                                     40.65
                 IF (IVTYPE.EQ.54 .AND. JDSXB.LE.1) THEN                   40.65
                    MCMVAR = MCMVAR+1                                      40.65
                    JDSXB  = MCMVAR                                        40.65
+                   ALOCMP = .TRUE.                                        40.97
                 ENDIF                                                     40.65
                 IF (IVTYPE.EQ.55 .AND. JDSXS.LE.1) THEN                   40.65
                    MCMVAR = MCMVAR+1                                      40.65
                    JDSXS  = MCMVAR                                        40.65
+                   ALOCMP = .TRUE.                                        40.97
                 ENDIF                                                     40.65
                 IF (IVTYPE.EQ.56 .AND. JDSXW.LE.1) THEN                   40.65
                    MCMVAR = MCMVAR+1                                      40.65
                    JDSXW  = MCMVAR                                        40.65
+                   ALOCMP = .TRUE.                                        40.97
                 ENDIF                                                     40.65
                 IF (IVTYPE.EQ.60 .AND. JGENR.LE.1) THEN                   40.85
                    MCMVAR = MCMVAR+1                                      40.85
                    JGENR  = MCMVAR                                        40.85
+                   ALOCMP = .TRUE.                                        40.97
                 ENDIF                                                     40.85
                 IF (IVTYPE.EQ.61 .AND. JGSXW.LE.1) THEN                   40.85
                    MCMVAR = MCMVAR+1                                      40.85
                    JGSXW  = MCMVAR                                        40.85
+                   ALOCMP = .TRUE.                                        40.97
                 ENDIF                                                     40.85
                 IF (IVTYPE.EQ.62 .AND. JREDS.LE.1) THEN                   40.85
                    MCMVAR = MCMVAR+1                                      40.85
                    JREDS  = MCMVAR                                        40.85
+                   ALOCMP = .TRUE.                                        40.97
                 ENDIF                                                     40.85
                 IF (IVTYPE.EQ.63 .AND. JRSXQ.LE.1) THEN                   40.85
                    MCMVAR = MCMVAR+1                                      40.85
                    JRSXQ  = MCMVAR                                        40.85
+                   ALOCMP = .TRUE.                                        40.97
                 ENDIF                                                     40.85
                 IF (IVTYPE.EQ.64 .AND. JRSXT.LE.1) THEN                   40.85
                    MCMVAR = MCMVAR+1                                      40.85
                    JRSXT  = MCMVAR                                        40.85
+                   ALOCMP = .TRUE.                                        40.97
                 ENDIF                                                     40.85
                 IF (IVTYPE.EQ.65 .AND. JTRAN.LE.1) THEN                   40.85
                    MCMVAR = MCMVAR+1                                      40.85
                    JTRAN  = MCMVAR                                        40.85
+                   ALOCMP = .TRUE.                                        40.97
                 ENDIF                                                     40.85
                 IF (IVTYPE.EQ.66 .AND. JTSXG.LE.1) THEN                   40.85
                    MCMVAR = MCMVAR+1                                      40.85
                    JTSXG  = MCMVAR                                        40.85
+                   ALOCMP = .TRUE.                                        40.97
                 ENDIF                                                     40.85
                 IF (IVTYPE.EQ.67 .AND. JTSXT.LE.1) THEN                   40.85
                    MCMVAR = MCMVAR+1                                      40.85
                    JTSXT  = MCMVAR                                        40.85
+                   ALOCMP = .TRUE.                                        40.97
                 ENDIF                                                     40.85
                 IF (IVTYPE.EQ.68 .AND. JTSXS.LE.1) THEN                   40.85
                    MCMVAR = MCMVAR+1                                      40.85
                    JTSXS  = MCMVAR                                        40.85
+                   ALOCMP = .TRUE.                                        40.97
                 ENDIF                                                     40.85
                 IF (IVTYPE.EQ.69 .AND. JRADS.LE.1) THEN                   40.85
                    MCMVAR = MCMVAR+1                                      40.85
                    JRADS  = MCMVAR                                        40.85
+                   ALOCMP = .TRUE.                                        40.97
                 ENDIF                                                     40.85
              ENDIF
              GOTO 70
@@ -1430,62 +1445,77 @@
             IF (IVTYPE.EQ.50 .AND. JPBOT.LE.1) THEN                       40.65
                MCMVAR = MCMVAR+1                                          40.65
                JPBOT  = MCMVAR                                            40.65
+               ALOCMP = .TRUE.                                            40.97
             ENDIF                                                         40.65
             IF (IVTYPE.EQ.52 .AND. JBOTLV.LE.1) THEN                      40.65
                MCMVAR = MCMVAR+1                                          40.65
                JBOTLV = MCMVAR                                            40.65
+               ALOCMP = .TRUE.                                            40.97
             ENDIF                                                         40.65
             IF (IVTYPE.EQ.54 .AND. JDSXB.LE.1) THEN                       40.65
                MCMVAR = MCMVAR+1                                          40.65
                JDSXB  = MCMVAR                                            40.65
+               ALOCMP = .TRUE.                                            40.97
             ENDIF                                                         40.65
             IF (IVTYPE.EQ.55 .AND. JDSXS.LE.1) THEN                       40.65
                MCMVAR = MCMVAR+1                                          40.65
                JDSXS  = MCMVAR                                            40.65
+               ALOCMP = .TRUE.                                            40.97
             ENDIF                                                         40.65
             IF (IVTYPE.EQ.56 .AND. JDSXW.LE.1) THEN                       40.65
                MCMVAR = MCMVAR+1                                          40.65
                JDSXW  = MCMVAR                                            40.65
+               ALOCMP = .TRUE.                                            40.97
             ENDIF                                                         40.65
             IF (IVTYPE.EQ.60 .AND. JGENR.LE.1) THEN                       40.85
                MCMVAR = MCMVAR+1                                          40.85
                JGENR  = MCMVAR                                            40.85
+               ALOCMP = .TRUE.                                            40.97
             ENDIF                                                         40.85
             IF (IVTYPE.EQ.61 .AND. JGSXW.LE.1) THEN                       40.85
                MCMVAR = MCMVAR+1                                          40.85
                JGSXW  = MCMVAR                                            40.85
+               ALOCMP = .TRUE.                                            40.97
             ENDIF                                                         40.85
             IF (IVTYPE.EQ.62 .AND. JREDS.LE.1) THEN                       40.85
                MCMVAR = MCMVAR+1                                          40.85
                JREDS  = MCMVAR                                            40.85
+               ALOCMP = .TRUE.                                            40.97
             ENDIF                                                         40.85
             IF (IVTYPE.EQ.63 .AND. JRSXQ.LE.1) THEN                       40.85
                MCMVAR = MCMVAR+1                                          40.85
                JRSXQ  = MCMVAR                                            40.85
+               ALOCMP = .TRUE.                                            40.97
             ENDIF                                                         40.85
             IF (IVTYPE.EQ.64 .AND. JRSXT.LE.1) THEN                       40.85
                MCMVAR = MCMVAR+1                                          40.85
                JRSXT  = MCMVAR                                            40.85
+               ALOCMP = .TRUE.                                            40.97
             ENDIF                                                         40.85
             IF (IVTYPE.EQ.65 .AND. JTRAN.LE.1) THEN                       40.85
                MCMVAR = MCMVAR+1                                          40.85
                JTRAN  = MCMVAR                                            40.85
+               ALOCMP = .TRUE.                                            40.97
             ENDIF                                                         40.85
             IF (IVTYPE.EQ.66 .AND. JTSXG.LE.1) THEN                       40.85
                MCMVAR = MCMVAR+1                                          40.85
                JTSXG  = MCMVAR                                            40.85
+               ALOCMP = .TRUE.                                            40.97
             ENDIF                                                         40.85
             IF (IVTYPE.EQ.67 .AND. JTSXT.LE.1) THEN                       40.85
                MCMVAR = MCMVAR+1                                          40.85
                JTSXT  = MCMVAR                                            40.85
+               ALOCMP = .TRUE.                                            40.97
             ENDIF                                                         40.85
             IF (IVTYPE.EQ.68 .AND. JTSXS.LE.1) THEN                       40.85
                MCMVAR = MCMVAR+1                                          40.85
                JTSXS  = MCMVAR                                            40.85
+               ALOCMP = .TRUE.                                            40.97
             ENDIF                                                         40.85
             IF (IVTYPE.EQ.69 .AND. JRADS.LE.1) THEN                       40.85
                MCMVAR = MCMVAR+1                                          40.85
                JRADS  = MCMVAR                                            40.85
+               ALOCMP = .TRUE.                                            40.97
             ENDIF                                                         40.85
             CALL INKEYW ('STA', ' ')                                      40.00
             IF (KEYWIS('UNIT')) THEN
@@ -1717,7 +1747,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1887,7 +1917,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2005,7 +2035,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2137,7 +2167,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -3052,29 +3082,35 @@
                      IARR1(K) = IP                                        40.92
                   ENDIF
                ENDDO
-               ALLOCATE(IARR2(K))                                         40.92
-               IARR2(1:K) = IARR1(1:K)                                    40.92
-               ISH = 0                                                    40.92
-               DO IPP = 2, K                                              40.92
-                  IF ( IARR2(IPP)/=IARR2(IPP-1)+IXI ) THEN                40.92
-                     ISH = IPP-1                                          40.92
-                     EXIT                                                 40.92
-                  ENDIF                                                   40.92
-               ENDDO                                                      40.92
-               IARR2 = CSHIFT(IARR2,ISH)                                  40.92
                !
-               DO IPP = 1, K                                              40.92
-                  IP = IARR2(IPP)                                         40.92
-                  IX = blist(IP,JBG)                                      40.92 40.80
-                  KOUNTR = KOUNTR + 1                                     40.80
-                  vert(IX)%atti(VBC) = 1                                  40.80
-                  ALLOCATE(TMP)                                           40.80
-                  TMP%JX = IX                                             40.80
-                  NULLIFY(TMP%NEXTXY)                                     40.80
-                  CURR%NEXTXY => TMP                                      40.80
-                  CURR => TMP                                             40.80
-               ENDDO                                                      40.80
-               DEALLOCATE(IARR1,IARR2)                                    40.92
+               IF ( K/=0 ) THEN                                           40.92
+                  !
+                  ALLOCATE(IARR2(K))                                      40.92
+                  IARR2(1:K) = IARR1(1:K)                                 40.92
+                  ISH = 0                                                 40.92
+                  DO IPP = 2, K                                           40.92
+                     IF ( IARR2(IPP)/=IARR2(IPP-1)+IXI ) THEN             40.92
+                        ISH = IPP-1                                       40.92
+                        EXIT                                              40.92
+                     ENDIF                                                40.92
+                  ENDDO                                                   40.92
+                  IARR2 = CSHIFT(IARR2,ISH)                               40.92
+                  !
+                  DO IPP = 1, K                                           40.92
+                     IP = IARR2(IPP)                                      40.92
+                     IX = blist(IP,JBG)                                   40.92 40.80
+                     KOUNTR = KOUNTR + 1                                  40.80
+                     vert(IX)%atti(VBC) = 1                               40.80
+                     ALLOCATE(TMP)                                        40.80
+                     TMP%JX = IX                                          40.80
+                     NULLIFY(TMP%NEXTXY)                                  40.80
+                     CURR%NEXTXY => TMP                                   40.80
+                     CURR => TMP                                          40.80
+                  ENDDO                                                   40.80
+                  DEALLOCATE(IARR2)                                       40.92
+                  !
+               ENDIF                                                      40.92
+               DEALLOCATE(IARR1)                                          40.92
                !
             ENDDO                                                         40.92
             !
@@ -3364,7 +3400,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -3855,7 +3891,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -4444,7 +4480,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -4905,7 +4941,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -5207,7 +5243,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -5369,7 +5405,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanPrepComp.ftn90	2008-05-25 14:59:46.000000000 +0200
+++ SwanPrepComp.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanPrintGridInfo.ftn90	2008-05-25 14:59:47.000000000 +0200
+++ SwanPrintGridInfo.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanPropvelS.ftn90	2008-05-25 14:59:47.000000000 +0200
+++ SwanPropvelS.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -16,7 +16,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -37,10 +37,12 @@
 !   Authors
 !
 !   40.80: Marcel Zijlema
+!   41.02: Marcel Zijlema
 !
 !   Updates
 !
-!   40.80, July 2007: New subroutine
+!   40.80,     July 2007: New subroutine
+!   41.02, February 2009: adaption of velocities in case of diffraction
 !
 !   Purpose
 !
@@ -52,6 +54,7 @@
     use swcomm2
     use swcomm3
     use swcomm4
+    use m_diffr
     use SwanGriddata
     use SwanGridobjects
     use SwanCompdata
@@ -380,6 +383,7 @@
           !
           cs(1) =  kwave(is,1) * spcsig(is) / sinh (2.* kd)
           cs(5) = -cgo(is,1) * kwave(is,1)
+          if ( IDIFFR /= 0 ) cs(5) = cs(5)*DIFPARAM(iv1)
           !
           cs(10) = cs(1) * cs(2)
           cs(11) = cs(1) * (cs(3)+cs(4))
@@ -450,7 +454,8 @@
              id = mod ( iddum - 1 + MDC , MDC ) + 1
              !
              cad(id,is) = esin(id)*cd(8) - ecos(id)*cd(9)
-             if ( ICUR /= 0 ) cad(id,is) = cad(id,is) + sincos(id)*cd(10) + sinsin(id)*cd(6) - coscos(id)*cd(7)
+             if ( IDIFFR /= 0 ) cad(id,is) = cad(id,is)*DIFPARAM(iv1) - DIFPARDX(iv1)*cgo(is,1)*esin(id) + DIFPARDY(iv1)*cgo(is,1)*ecos(id)
+             if ( ICUR   /= 0 ) cad(id,is) = cad(id,is) + sincos(id)*cd(10) + sinsin(id)*cd(6) - coscos(id)*cd(7)
              !
           enddo
           !
@@ -461,7 +466,8 @@
           id = mod ( idcmin(is) - 2 + MDC , MDC ) + 1  ! this direction belongs to sweep below considered sweep
           !
           cad(id,is) = esin(id)*cd(17) - ecos(id)*cd(18)
-          if ( ICUR /= 0 ) cad(id,is) = cad(id,is) + sincos(id)*cd(19) + sinsin(id)*cd(15) - coscos(id)*cd(16)
+          if ( IDIFFR /= 0 ) cad(id,is) = cad(id,is)*DIFPARAM(iv1) - DIFPARDX(iv1)*cgo(is,1)*esin(id) + DIFPARDY(iv1)*cgo(is,1)*ecos(id)
+          if ( ICUR   /= 0 ) cad(id,is) = cad(id,is) + sincos(id)*cd(19) + sinsin(id)*cd(15) - coscos(id)*cd(16)
           !
           cd(26) = cd( 1) * cd(20)
           cd(27) = cd( 1) * cd(21)
@@ -470,7 +476,8 @@
           id = mod ( idcmax(is) + MDC , MDC ) + 1      ! this direction belongs to sweep above considered sweep
           !
           cad(id,is) = esin(id)*cd(26) - ecos(id)*cd(27)
-          if ( ICUR /= 0 ) cad(id,is) = cad(id,is) + sincos(id)*cd(28) + sinsin(id)*cd(24) - coscos(id)*cd(25)
+          if ( IDIFFR /= 0 ) cad(id,is) = cad(id,is)*DIFPARAM(iv1) - DIFPARDX(iv1)*cgo(is,1)*esin(id) + DIFPARDY(iv1)*cgo(is,1)*ecos(id)
+          if ( ICUR   /= 0 ) cad(id,is) = cad(id,is) + sincos(id)*cd(28) + sinsin(id)*cd(24) - coscos(id)*cd(25)
           !
        enddo
        !
--- SwanPropvelX.ftn90	2008-05-25 14:59:47.000000000 +0200
+++ SwanPropvelX.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -32,10 +32,12 @@
 !   Authors
 !
 !   40.80: Marcel Zijlema
+!   41.02: Marcel Zijlema
 !
 !   Updates
 !
-!   40.80, July 2007: New subroutine
+!   40.80,     July 2007: New subroutine
+!   41.02, February 2009: adaption of velocities in case of diffraction
 !
 !   Purpose
 !
@@ -45,6 +47,7 @@
 !
     use ocpcomm4
     use swcomm3
+    use m_diffr
     use SwanGriddata
     use SwanCompdata
 !
@@ -88,6 +91,17 @@
           enddo
        enddo
        !
+       ! adapt the celerities in case of diffraction
+       !
+       if ( IDIFFR /= 0 .and. PDIFFR(3) /= 0. ) then
+          do is = 1, MSC
+             do id = 1 ,MDC
+                cax(id,is,ic) = cax(id,is,ic)*DIFPARAM(ivert)
+                cay(id,is,ic) = cay(id,is,ic)*DIFPARAM(ivert)
+             enddo
+          enddo
+       endif
+       !
        ! ambient currents added
        !
        if ( ICUR /= 0 )  then
--- SwanReadADCGrid.ftn90	2008-05-25 14:59:47.000000000 +0200
+++ SwanReadADCGrid.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanReadEasymeshGrid.ftn90	2008-05-25 14:59:47.000000000 +0200
+++ SwanReadEasymeshGrid.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanReadGrid.ftn90	2008-05-25 14:59:47.000000000 +0200
+++ SwanReadGrid.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanReadTriangleGrid.ftn90	2008-05-25 14:59:48.000000000 +0200
+++ SwanReadTriangleGrid.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- swanser.ftn	2008-05-25 14:59:50.000000000 +0200
+++ swanser.ftn	2009-03-03 17:50:34.000000000 +0100
@@ -68,7 +68,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -202,7 +202,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -329,7 +329,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -428,7 +428,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -582,7 +582,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -755,7 +755,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -826,7 +826,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -900,7 +900,6 @@
 !     FIRST = True
 !     For ix=1 to MXC-1 do
 !         For iy=1 to MYC-1 do
-!             Inmesh = True
 !             For iside=1 to 4 do
 !                 Case iside=
 !                 1: K1 = KGRPNT(ix,iy), K2 = KGRPNT(ix+1,iy),
@@ -1057,7 +1056,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1242,7 +1241,7 @@
 !         --- First attempt, to find XC,YC with Newton-Raphson method
           MXITNR = 5
           CALL NEWTON  (XP, YP, XCGRID, YCGRID, KGRPNT,                   40.00
-     &                  MXITNR ,ITER, XC ,YC ,FINDXY, KGRBND)             40.41 40.02
+     &                  MXITNR ,ITER, XC ,YC ,FINDXY)                     40.41 40.02
           IF ((ITEST .GE. 150 .OR. INTES .GE. 20) .AND. FINDXY) THEN      40.02
             WRITE(PRINTF,25) XP+XOFFS ,YP+YOFFS ,XC ,YC                   40.03
           ENDIF
@@ -1293,7 +1292,7 @@
               YC = YC0 - 0.2
             ENDIF
             CALL NEWTON  (XP, YP, XCGRID, YCGRID, KGRPNT,                 40.00
-     &                    MXITNR ,ITER, XC ,YC ,FINDXY, KGRBND)           40.41 40.02
+     &                    MXITNR ,ITER, XC ,YC ,FINDXY)                   40.41 40.02
             IF (FINDXY) THEN
               IF (ITEST .GE. 150 .OR. INTES .GE. 20) THEN
                 WRITE(PRINTF,25) XP+XOFFS ,YP+YOFFS ,XC ,YC               40.00
@@ -1404,7 +1403,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1572,7 +1571,7 @@
 !***********************************************************************
 !                                                                      *
       SUBROUTINE NEWTON (XP, YP, XCGRID, YCGRID, KGRPNT,                  40.00
-     &                   MXITNR, ITER, XC, YC, FIND, KGRBND)              40.41 40.02
+     &                   MXITNR, ITER, XC, YC, FIND)                      40.41 40.02
 !                                                                      *
 !***********************************************************************
 !
@@ -1592,7 +1591,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1649,11 +1648,10 @@
 !
 !  4. Argument variables
 !
-! i   KGRBND: Grid adresses of the boundary points                        40.02
 ! i   KGRPNT: Grid adresses                                               40.00
 ! i   MXITNR: Maximum number of iterations                                30.82
 !
-      INTEGER KGRBND(*), KGRPNT(MXC,MYC), MXITNR, ITER                    40.41 40.02
+      INTEGER KGRPNT(MXC,MYC), MXITNR, ITER                               40.41 40.02
 !
 !   o XC    : X-coordinate in computational coordinates                   30.82
 ! i   XCGRID: Coordinates of computational grid in x-direction            30.72
@@ -1812,7 +1810,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1965,7 +1963,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2087,7 +2085,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2303,7 +2301,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2502,7 +2500,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2719,7 +2717,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -3005,7 +3003,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -3487,7 +3485,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -3838,7 +3836,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -4270,7 +4268,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -4672,7 +4670,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -4797,7 +4795,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -4895,7 +4893,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -4995,7 +4993,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -5183,7 +5181,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -5333,7 +5331,9 @@
             GOTO 70
           ENDIF
           IF (X1A.GT.X2HI) GOTO 300
+          IF (X1B.LT.X2LO) GOTO 300
           IF (X1A.LT.X2LO .AND. X1A+PERIOD.LT.X2HI) TWICE = .TRUE.
+          IF (X1B.GT.X2HI .AND. X1B-PERIOD.GT.X2LO) TWICE = .TRUE.
         ENDIF
 !
 !       loop over cells in X2
@@ -5377,8 +5377,14 @@
 !
 !       Cell in X1 covers both ends of sector boundary
         IF (TWICE) THEN
-          X1A = X1A + PERIOD
-          X1B = X1B + PERIOD
+          IF (X1A.LT.X2LO) THEN
+             X1A = X1A + PERIOD
+             X1B = X1B + PERIOD
+          ENDIF
+          IF (X1B.GT.X2HI) THEN
+             X1A = X1A - PERIOD
+             X1B = X1B - PERIOD
+          ENDIF
           TWICE = .FALSE.
           GOTO 100
         ENDIF
@@ -5481,7 +5487,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -5617,7 +5623,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -5733,7 +5739,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -5902,7 +5908,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -6064,7 +6070,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -6397,7 +6403,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -6499,7 +6505,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -6591,7 +6597,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -6673,7 +6679,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -6788,7 +6794,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -6903,7 +6909,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -7044,7 +7050,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- SwanSweepSel.ftn90	2008-05-25 14:59:48.000000000 +0200
+++ SwanSweepSel.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -14,7 +14,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- swantech.tex	2009-03-03 17:49:37.000000000 +0100
+++ swantech.tex	2009-03-03 17:55:03.000000000 +0100
@@ -17,7 +17,7 @@
 \end{center}
 \vfill
 \begin{center}
-{\Large\bf SWAN Cycle III version 40.72A}
+{\Large\bf SWAN Cycle III version 40.72AB}
 \end{center}
 
 \cleardoublepage
@@ -39,14 +39,14 @@
              & & The Netherlands \\
              & & \\
 e-mail       &:& swan-info-citg@tudelft.nl \\
-home page    &:& \hl{http://www.fluidmechanics.tudelft.nl/swan/index.htm}
+home page    &:& \hl{http://www.swan.tudelft.nl}
 \end{tabular}
 \end{table}
 
 \vfill
 
 \noindent
-Copyright (c) 2008 Delft University of Technology.
+Copyright (c) 2009 Delft University of Technology.
 \\[2ex]
 \noindent
 Permission is granted to copy, distribute and/or modify this document
@@ -105,7 +105,7 @@
 \\[2ex]
 \noindent
 SWAN is developed at Delft University of Technology and is freely available from\\
-\hl{http://www.fluidmechanics.tudelft.nl/swan/index.htm}. It is used by many goverment
+\hl{http://www.swan.tudelft.nl}. It is used by many goverment
 authorities, research institutes and consultants worldwide. The feedback has widely
 indicated the reliability of SWAN in different experiment and field cases.
 \\[2ex]
@@ -481,7 +481,10 @@
 $ds=(c_x \cos \theta + c_y \sin \theta)dt$ yields
 \begin{equation}
   \frac{d\theta}{dt} = -\frac{c_x \cos \theta + c_y \sin \theta}{R} \cos \theta \tan \varphi
+  \label{eq:dthdt}
 \end{equation}
+This term (\ref{eq:dthdt}) accounts for the change of propagation direction relative to true North when travelling
+along a great circle.
 This holds for deep water and without currents. Hence,
 \begin{equation}
   \tilde{c}_\theta = c_\theta - \frac{c_x \cos \theta + c_y \sin \theta}{R} \cos \theta \tan \varphi
@@ -592,17 +595,20 @@
 \noindent
 A number of alternative whitecapping expressions have been proposed to improve the accuracy of SWAN. These range
 from alternative calibrations of the Komen~{\it et~al} (1984) expression, e.g. Rogers~{\it et~al} (2003), to
-alternative ways of calculating mean spectral steepness, e.g. Van Vledder and Hurdle (2002). In SWAN, two
-alternatives are presented.
+alternative ways of calculating mean spectral steepness, e.g. Van Vledder and Hurdle (2002).
+In SWAN, another alternative is presented.
+%In SWAN, two alternatives are presented.
 \nocite{Vle02H}
 \\[2ex]
 \noindent
-An alternative formulation for whitecapping is based on the Cumulative Steepness Method as described
-in Hurdle and Van Vledder (2004). With this method dissipation due to whitecapping depends on the steepness of
-the wave spectrum at and below a particular frequency.
-\\[2ex]
-\noindent
-A second alternative of the whitecapping expression is based on Alves and Banner (2003). This expression is based
+%An alternative formulation for whitecapping is based on the Cumulative Steepness Method as described
+%in Hurdle and Van Vledder (2004). With this method dissipation due to whitecapping depends on the steepness of
+%the wave spectrum at and below a particular frequency.
+%\\[2ex]
+%\noindent
+%A second
+An
+alternative of the whitecapping expression is based on Alves and Banner (2003). This expression is based
 on experimental findings that whitecapping dissipation appears to be related to the nonlinear hydrodynamics within
 wave groups. This yields a dissipation term that primarily depends on quantities that are local in the frequency
 spectrum, as opposed to ones that are distributed over the spectrum, as in the expression of Komen~{\it et~al} (1984).
@@ -669,11 +675,13 @@
 Boussinesq model that the high-frequency levels are insensitive to such frequency dependency because
 an increased dissipation at high frequencies is compensated approximately by increased nonlinear
 energy transfer (but they did find the frequency dependency to be relevant in time domain). The value of
-$D_{\rm tot}$ depends critically on the breaking parameter $\gamma = H_{\max}/d$
+$D_{\rm tot}$ depends critically on the breaker parameter $\gamma = H_{\max}/d$
 (in which $H_{\max}$ is the maximum possible individual wave height in the local water depth $d$).
-In SWAN, both a constant value and a variable value are available. The constant value is $\gamma=0.73$ found as
+In SWAN, both a constant value and a variable value are available. Examples of a variable breaker parameter
+can be found in Nelson (1987) and Ruessink et al. (2003). (Both are implemented in SWAN.)
+The constant value is $\gamma=0.73$ found as
 the mean value of the data set of Battjes and Stive (1985).
-\nocite{Bat78J,Tho83G,Bat92B,Vin94SD,Arc94RC,Eld96B,Eld95B,Mas92K,Elg97GRHG,Che97G,Bat85S}
+\nocite{Bat78J,Tho83G,Bat92B,Vin94SD,Arc94RC,Eld96B,Eld95B,Mas92K,Elg97GRHG,Che97G,Nel87,Rue03WS,Bat85S}
 \\[2ex]
 \noindent
 \underline{Nonlinear wave-wave interactions}\\[2ex]
@@ -729,7 +737,7 @@
 waves and is therefore not suitable in many practical applications of wind waves. The breakthrough in the
 development came with the work of Eldeberky and Battjes (1995) who transformed the amplitude part of
 the Boussinesq model of Madsen and S{\o}rensen (1993) into an energy density formulation and who
-parameterized the bi-phase of the waves on the basis of laboratory observations (Battjes and Beji, 1992;
+parameterized the biphase of the waves on the basis of laboratory observations (Battjes and Beji, 1992;
 Arcilla~{\it et~al}., 1994). A discrete triad approximation (DTA) for co-linear waves was subsequently
 obtained by considering only the dominant self-self interactions. Their model has been verified with flume
 observations of long-crested, random waves breaking over a submerged bar (Beji and Battjes, 1993) and
@@ -887,32 +895,32 @@
 Cycle~4; Komen~{\it et~al}., 1994).
 \\[2ex]
 \noindent
-\underline{Whitecapping: CSM formulation}\\[2ex]
-An alternative formulation for whitecapping is based on the Cumulative Steepness Method as described
-in Hurdle and Van Vledder (2004). With this method dissipation due to whitecapping depends on the steepness of
-the wave spectrum at and below a particular frequency. It is defined as (directionally dependent):
-\begin{equation}
-  S_{\rm st} (\sigma,\theta) = A_m \int_{0}^{\sigma} \int_{0}^{2\pi} k^2 | \cos (\theta - \theta')|^m E
-  (\sigma,\theta)d\sigma d\theta
-  \label{eq2-5}
-\end{equation}
-with $A_m$ the normalisation coefficient as determined by
-\begin{equation}
-  \int_{0}^{2\pi} A_m \cos^m (\theta) d\theta = 1
-\end{equation}
-In expression (\ref{eq2-5}) the coefficient $m$ controls the directional dependence. It is expected that this
-coefficient will be order 1 if the straining mechanism is dominant, $m$ is more than 10 if other mechanism play
-a role (e.g. instability that occurs when vertical acceleration in the waves becomes greater than gravity).
-Default in SWAN is $m=2$. The alternative whitecapping source term is given by
-\begin{equation}
-  S_{\rm wc}^{\rm st} = -C_{\rm wc}^{\rm st} \left( S_{\rm st} (\sigma,\theta) \right)^p E (\sigma,\theta)
-  \label{eq2-6}
-\end{equation}
-with $C_{wc}^{st}$ a tuneable coefficient and $p$ a parameter controlling the proportionality of the dissipation
-rate on the steepness. In SWAN, $p=1$ is assumed.
-\nocite{Hur04V}
-\\[2ex]
-\noindent
+%\underline{Whitecapping: CSM formulation}\\[2ex]
+%An alternative formulation for whitecapping is based on the Cumulative Steepness Method as described
+%in Hurdle and Van Vledder (2004). With this method dissipation due to whitecapping depends on the steepness of
+%the wave spectrum at and below a particular frequency. It is defined as (directionally dependent):
+%\begin{equation}
+%  S_{\rm st} (\sigma,\theta) = A_m \int_{0}^{\sigma} \int_{0}^{2\pi} k^2 | \cos (\theta - \theta')|^m E
+%  (\sigma,\theta)d\sigma d\theta
+%  \label{eq2-5}
+%\end{equation}
+%with $A_m$ the normalisation coefficient as determined by
+%\begin{equation}
+%  \int_{0}^{2\pi} A_m \cos^m (\theta) d\theta = 1
+%\end{equation}
+%In expression (\ref{eq2-5}) the coefficient $m$ controls the directional dependence. It is expected that this
+%coefficient will be order 1 if the straining mechanism is dominant, $m$ is more than 10 if other mechanism play
+%a role (e.g. instability that occurs when vertical acceleration in the waves becomes greater than gravity).
+%Default in SWAN is $m=2$. The alternative whitecapping source term is given by
+%\begin{equation}
+%  S_{\rm wc}^{\rm st} = -C_{\rm wc}^{\rm st} \left( S_{\rm st} (\sigma,\theta) \right)^p E (\sigma,\theta)
+%  \label{eq2-6}
+%\end{equation}
+%with $C_{wc}^{st}$ a tuneable coefficient and $p$ a parameter controlling the proportionality of the dissipation
+%rate on the steepness. In SWAN, $p=1$ is assumed.
+%\nocite{Hur04V}
+%\\[2ex]
+%\noindent
 \underline{Whitecapping: saturation-based model}\\[2ex]
 The whitecapping formulation used in SWAN is an adapted form of the expression of Alves and Banner (2003), which is based
 on the apparent relationship between wave groups and whitecapping dissipation. This adaption is due to the fact that it
@@ -984,9 +992,9 @@
   \label{eq3-17}
 \end{equation}
 in which $C_{\rm b}$ is a bottom friction coefficient that generally depends on the bottom orbital motion
-represented by $U_{rms}$:
+represented by $U_{\rm rms}$:
 \begin{equation}
-  U^2_{rms} = \int_{0}^{2\pi} \int_{0}^{\infty} \frac{\sigma^2}{g^2 \sinh^2 kd} E(\sigma,\theta) d\sigma d\theta
+  U^2_{\rm rms} = \int_{0}^{2\pi} \int_{0}^{\infty} \frac{\sigma^2}{g^2 \sinh^2 kd} E(\sigma,\theta) d\sigma d\theta
   \label{eq3-18}
 \end{equation}
 Hasselmann~{\it et~al}. (1973) found from the results of the JONSWAP experiment
@@ -1006,7 +1014,7 @@
 model the bottom friction factor is a function of the bottom roughness height and the actual wave
 conditions. Their bottom friction coefficient is given by:
 \begin{equation}
-  C_{\rm b} = f_w \frac{g}{\sqrt{2}} U_{rms}
+  C_{\rm b} = f_w \frac{g}{\sqrt{2}} U_{\rm rms}
   \label{eq3-19}
 \end{equation}
 in which $f_w$ is a non-dimensional friction factor estimated by using the formulation of Jonsson (1966) cf.
@@ -1084,6 +1092,42 @@
 of 0.73. From a compilation of a large number of experiments Kaminsky and Kraus (1993) have found
 breaker parameters in the range of 0.6 to 1.59 with an average of 0.79.
 \nocite{Arc90L,Kam93K,Nel87,Nel94}
+\\[2ex]
+\noindent
+An alternative to the bore-based model of Battjes and Janssen (1978) is proposed by Thornton and Guza (1983).
+This model can be regarded as an alteration of Battjes and Janssen with respect to the description of the wave
+height probability density function. The total dissipation due to depth-induced breaking is formulated as
+\begin{equation}
+  D_{\rm tot} = - \frac{B^3 \tilde{\sigma}}{8 \pi d} \int_{0}^{\infty} H^3 \, p_b(H)\,dH
+  \label{eq:tg83}
+\end{equation}
+in which $B$ is a proportionality coefficient and $p_b(H)$ is the probability density function of breaking
+waves times the fraction of breakers, $Q_b$. Based on field observations, the wave heights in the surf zone
+are assumed to remain Rayleigh distributed, even after breaking. This implies that all waves will break, not
+only the highest as assumed by Battjes and Janssen (1978). The function $p_b(H)$ is obtained by multiplying
+the Rayleigh wave height probability density function $p(H)$, given by
+\begin{equation}
+  p(H) = \frac{2H}{H^2_{\rm rms}}\exp \left( - \left( \frac{H}{H_{\rm rms}} \right)^2 \right)
+\end{equation}
+by a weighting function $W(H)$ defined so that $0 \leq W(H) \leq 1$, to yield
+\begin{equation}
+  p_b(H) = W(H)\,p(H)
+\end{equation}
+Thornton and Guza (1983) proposed the following weighting function in which the fraction of breaking waves
+is independent of the wave height:
+\begin{equation}
+  W(H) = Q_b = \left( \frac{H_{\rm rms}}{\gamma d} \right)^n
+\end{equation}
+with a calibration parameter $n$ and a breaker index $\gamma$ (not to be confused with the Battjes and Janssen breaker index!).
+The integral in expression (\ref{eq:tg83}) can then be simplified, as follows:
+\begin{equation}
+  \int_{0}^{\infty} H^3 \, p_b(H)\,dH = Q_b\,\int_{0}^{\infty} H^3 \, p(H)\,dH = \frac{3}{4}\sqrt{\pi}\,Q_b\,H^3_{\rm rms}
+\end{equation}
+Hence,
+\begin{equation}
+  D_{\rm tot} = - \frac{3 B^3 \tilde{\sigma}}{32 \sqrt{\pi} d} \, Q_b\, H^3_{\rm rms}
+  \label{eq:tg83_2}
+\end{equation}
 
 \subsection{Nonlinear wave-wave interactions ($S_{\rm nl}$)} \label{sec:waveint}
 
@@ -1352,7 +1396,7 @@
   S^-_{\rm nl3} (\sigma,\theta) = -2 S^+_{\rm nl3} (2\sigma,\theta)
   \label{eq3-35}
 \end{equation}
-in which $\alpha _{\rm EB}$ is a tunable proportionality coefficient. The bi-phase $\beta$ is approximated with
+in which $\alpha _{\rm EB}$ is a tunable proportionality coefficient. The biphase $\beta$ is approximated with
 \begin{equation}
   \beta = -\frac{\pi}{2} + \frac{\pi}{2} \tanh (\frac{0.2}{Ur})
   \label{eq3-36}
@@ -1834,7 +1878,14 @@
   D_{xx} = D_{ss} \cos^2\theta + D_{nn} \sin^2 \theta \,, \quad D_{yy} = D_{ss} \sin^2\theta + D_{nn} \cos^2 \theta \,, \quad
   D_{xy} = (D_{ss}-D_{nn}) \cos \theta \sin \theta
 \end{equation}
-The diffusion terms are computed at the time level $i_t-1$. The diffusion terms are computed as follows
+The diffusion terms are given by
+\begin{equation}
+  -D_{xx}\frac{\partial^2 N}{\partial x^2}
+  -2D_{xy}\frac{\partial^2 N}{\partial x \partial y}
+  -D_{yy}\frac{\partial^2 N}{\partial y^2}
+  \label{eq:gsediff}
+\end{equation}
+and are computed at the time level $i_t-1$ as follows
 \begin{equation}
   D_{xx} \left( \frac{(N)_{i_x+1} - 2(N)_{i_x} + (N)_{i_x-1}}{\Delta x^2} \right)^{i_t-1}_{i_y, i_{\sigma}, i_{\theta}}
 \end{equation}
@@ -2508,7 +2559,7 @@
 \noindent
 After the transmission coefficient has been calculated, it is used in
 the propagation terms of the action balance equation.
-In curvilinear coordinates, the propagation terms
+In curvi-linear coordinates, the propagation terms
 (including time-derivative, but ignoring dependence on
 $\sigma$ and $\theta$ temporarily) read:
 \begin{eqnarray}
@@ -2715,57 +2766,71 @@
 
 \section{Computation of breaking source term}
 
-The surf breaking dissipation reads
+The surf breaking dissipation of Battjes and Janssen (1978) reads
 \begin{equation}
-  D_{\rm tot} = -\beta Q_{b} \frac{H^2_m}{8}
+  D_{\rm tot} = -\alpha_{\rm BJ} Q_{b} \tilde{\sigma} \frac{H^2_{\rm max}}{8\pi}
 \end{equation}
 The surf breaking source term for each spectral bin $i$ is
 \begin{equation}
-  S_{\rm br} = \frac{D_{\rm tot}}{E_{\rm tot}} E_i = -\beta \frac{Q_b}{B} E_i
+  S_{i} = \frac{D_{\rm tot}}{E_{\rm tot}} E_i = {\tilde D}\, E_i
+  \label{eq:st1}
 \end{equation}
-with
+with the normalized total dissipation
 \begin{equation}
-  B = \frac{8E_{\rm tot}}{H^2_m}
+  {\tilde D} = -\frac{\alpha_{\rm BJ} \tilde{\sigma} Q_{b}}{\pi {\cal B}} < 0
 \end{equation}
-In the numerical process, we need old values (from previous iteration) to express $S$ in $E$. This is linearized
-in the numerical process by using the derivative of $Q_b$ (which is the factor varying the most) with respect
-to $B$. The two are related by
+and
 \begin{equation}
-  1- Q_b + B \ln Q_b = 0
-  \label{eq:breakfract}
+  {\cal B} = \frac{8E_{\rm tot}}{H^2_{\rm max}} = \left( \frac{H_{\rm rms}}{\gamma d} \right)^2
 \end{equation}
-The derivative is found by differentiating this with respect to $B$:
+Since, the source term is strongly nonlinear in $E$ (since ${\tilde D}$ depends on $E$ through ${\cal B}$), we apply the Newton linearisation
+to approximate the source term at iteration level $n+1$, as follows:
 \begin{equation}
-  -{Q'}_{b} + \ln Q_b + \frac{B}{Q_b} {Q'}_{b} = 0
+   S^{n+1}_{i} \approx {\tilde D} E_i^{n} + \left( \frac{\partial S}{\partial E} \right)_i^n (E_i^{n+1} - E_i^n)
+   \label{eq:newlin}
 \end{equation}
-Hence,
+In SWAN, this approximation has been slightly adapted for reasons of numerical stability; the first term in the right-hand
+side, ${\tilde D} E_i^{n}$, is replaced by ${\tilde D} E_i^{n+1}$. This preserves positivity of energy density $E$, if the following
+inequality holds
 \begin{equation}
-  {Q'}_{b} = \frac{\ln Q_b}{1 - B/Q_b} = \frac{Q_b}{B} \frac{Q_b - 1}{Q_b - B}
+  \frac{\partial S}{\partial E} < 0
 \end{equation}
-using Eq. (\ref{eq:breakfract}).
-The approximation of the source term is based on
+We derive an expression for this derivative as follows. From (\ref{eq:st1}), we have
 \begin{equation}
-  S_{\rm br,i}^{+} = S_{\rm br,i}^{-} + \frac{d S_{\rm br}}{d E_{\rm tot}} (E_{i}^{+} - E_{i}^{-})
+  \frac{\partial S}{\partial E} |_i = \frac{\partial {\tilde D}}{\partial E}|_i E_i + {\tilde D}
 \end{equation}
-Here,
+The normalized dissipation ${\tilde D}$ is a function of ${\cal B}$ which is proportional to $E$, so
 \begin{equation}
-  \frac{d S_{\rm br}}{d E_{\rm tot}} = -\beta \frac{H_m^2}{8} \frac{dQ_b}{dB} \frac{B}{E_{\rm tot}} =
-  -\beta \frac{dQ_b}{dB}
+  \frac{\partial S}{\partial E} |_i = \frac{\partial {\tilde D}}{\partial {\cal B}}|_i {\cal B}_i + {\tilde D}
+  \label{eq:st4}
 \end{equation}
-So,
+Since, $Q_b$ is a function of ${\cal B}$, we get (using the quotient rule)
 \begin{equation}
-  S_{\rm br,i}^{+} =-\beta \frac{Q_b}{B}E_{i}^{-} - \beta \frac{Q_b}{B} \frac{Q_b - 1}{Q_b - B} (E_{i}^{+} - E_{i}^{-})
+  \frac{\partial S}{\partial E} |_i = -\frac{\alpha_{\rm BJ} \tilde{\sigma}}{\pi}\, \frac{\partial Q_b}{\partial {\cal B}}
+  \label{eq:st2}
 \end{equation}
-or
+Since,
 \begin{equation}
-  S_{\rm br,i}^{+} =\beta \frac{Q_b}{B} \frac{1 - B}{B - Q_b} E_{i}^{-} + \beta \frac{Q_b}{B} \frac{1 - Q_b}{Q_b - B} E_{i}^{+}
+  1- Q_b + {\cal B} \ln Q_b = 0
+  \label{eq:breakfract}
+\end{equation}
+the derivative of $Q_b$ is found by differentiating this with respect to ${\cal B}$:
+\begin{equation}
+  -{Q'}_{b} + \ln Q_b + \frac{\cal B}{Q_b} {Q'}_{b} = 0
+\end{equation}
+Hence,
+\begin{equation}
+  {Q'}_{b} = \frac{\ln Q_b}{1 - {\cal B}/Q_b} = \frac{Q_b}{\cal B} \frac{Q_b - 1}{Q_b - {\cal B}}
+\end{equation}
+using Eq. (\ref{eq:breakfract}). Now, ${Q'}_{b} > 0$, because $0 < Q_b < 1$ and ${\cal B} > Q_b$.
+Substitution in (\ref{eq:st2}) gives
+\begin{equation}
+  \frac{\partial S}{\partial E} |_i = {\tilde D} \, \frac{Q_b - 1}{Q_b - {\cal B}} |_i < 0
+\end{equation}
+Finally, the approximation of the source term reads
+\begin{equation}
+  S^{n+1}_{i} = {\tilde D} \left( 1 + \frac{Q_b - 1}{Q_b - {\cal B}} \right)_i^n E_i^{n+1} - {\tilde D} \, \frac{Q_b - 1}{Q_b - {\cal B}} |^n_i E_i^n
 \end{equation}
-In order to prevent oscillations the factor $\beta \frac{Q_b}{B} \frac{1 - B}{B - Q_b}$ was enhanced to
-$\beta \frac{Q_b}{B} \frac{1 - Q_b}{B - Q_b}$. This does not influence the limit of the process.
-\\[2ex]
-In the SWAN program, $S$ is dissipation in terms of action, and $E$ is used for action density. This does not matter since, the ratio
-between dissipation in terms of energy and dissipation in terms of action is the same as the ratio between energy density and
-action density.
 
 \chap{Wave boundary and initial conditions} \label{ch:bc}
 
@@ -2902,11 +2967,10 @@
 Appendix A, the design of the modifications of the code is presented in Appendix
 B. After each iteration performed in SWAN new values of the setup are being
 calculated and added to the depth, so that the SWAN model incorporates the effect
-of setup on the wave field. An output quantity SETUP is added so that the user can
+of setup on the wave field. An output quantity {\tt SETUP} is added so that the user can
 be informed about the magnitude and distribution of the wave setup.
 \\[2ex]
 \noindent
-% stuk van Guus
 \section{Analysis and Results}
 \subsection{Discretization of the 2D setup equation}
 \label{sec1G}
@@ -3436,6 +3500,92 @@
 Hence, SWAN is stable for any time step.
 Because of this nice property, this solution technique is tailored to unstructured grids.
 
+\section{Description of an unstructured grid}
+\label{sec:defunstr}
+
+\subsection{Definitions}
+
+We distinguish between two types of grids, namely structured and unstructured grids. A two-dimensional structured grid may contain quadrilaterals. These
+can be recti-linear or curvi-linear. The number
+of cells that meet each other in an internal vertex is always 4. In unstructured meshes this restriction is abandoned. Moreover, 2D unstructured grids usually
+consist of triangles or a combination of triangles and quadrilaterals, a so-called hybrid grid. The unstructured meshes that we consider in SWAN consist solely
+of triangles, also called cells. The edges of the triangles are called faces.
+
+\subsection{Relations between number of cells, vertices and faces}
+
+For a two-dimensional triangular mesh, the number of cells $C$, the number of boundary faces $E_b$ and internal faces $E_i$ are related according to:
+\begin{equation}
+  E_b + 2 E_i = 3C
+\end{equation}
+The total number of faces $E = E_i + E_b$. With $V$ the number of vertices and $H$ the number of holes ('islands'), we have the following Euler's relation
+for a triangulation:
+\begin{equation}
+  C + V  - E = 1 - H
+\end{equation}
+Usually, $E_b << E_i$ and the number of holes $H$ is negligibly small, so
+\begin{equation}
+  C \approx 2V\, , \quad E \approx 3V
+\end{equation}
+There are approximately twice as many cells as vertices
+in a triangular mesh.
+Therefore, it is an optimal choice to locate the action density in vertices
+as the number of unknowns is minimal on a given grid.
+Concerning the
+time-consuming evaluation of the physical processes representing the wave energy generation,
+dissipation and redistribution, this allows SWAN to save a considerable amount of computing time.
+
+\subsection{Conditions imposed to the grid}
+
+In order to avoid badly shaped grids, the grids must satisfy the following
+properties:
+\begin{itemize}
+  \item The number of cells that meet at each vertex in the interior of the mesh must be at least 4 and at most 10.
+  \item The angles inside each triangle must be smaller than a certain value. Let $\vec{a}$ and $\vec{b}$ be the tangential vectors of two faces of a
+        triangle, then the angle $\phi$ between these two faces equals
+        \begin{equation}
+           \cos \phi = \frac{\vec{a} \cdot \vec{b}}{|\vec{a}||\vec{b}|}
+        \end{equation}
+        For safety, we do not allow for angles with $\cos \phi < -0.8$ or, equivalently, $\phi > 143^o$.
+\end{itemize}
+
+\section{Some notes on grid generation}
+
+We briefly outline some issues related to grid generation from a practical point of view.
+The process of grid generation can be difficult and time consuming. A common approach is proceeding from coarse to fine grid through refinement in various ways.
+Generally, one would like to have an optimal grid in which areas where the bathymetry or evolution of the waves change
+rapidly require a higher resolution than areas where the physics or depth changes less.
+This goes around by having an indication how to determine the refinement based on
+bathymetry or geometric variations through preliminary evaluations.
+To facilitate this procedure,
+a mesh generation package BatTri (Bilgili and Smith, 2003) is used.
+This grid generator is a public-domain, graphical Matlab interface to Triangle (Shewchuk, 1996).
+Triangle is a freely-distributed, two-dimensional Delaunay triangulator.
+Generation of the meshes presented in this paper was accomplished using BatTri.
+\nocite{Bil03S,She96}
+
+An important key ingredient for the preparation of the grid for the model domain is bathymetry data.
+Boundary nodes, segments and holes
+can be created from this data with the
+use of the mesh editing options of BatTri.
+After checking and improving grid quality, the final information on nodes and segments is forced into the
+triangulation of the domain. Triangle is called within BatTri to generate
+the actual grid. This triangulation includes only acute triangles.
+
+BatTri provides many pre-defined depth-dependent contraints for further mesh refinement.
+From a numerical point of view, mesh refinement is often directly related to properly resolve the shape of the wave,
+i.e. to keep the wavelength to grid size ratio relatively large.
+When wavelength decreases in shallower water, the grid size must decrease as well. Therefore, this criterion, called the $h$-refinement,
+has the effect of using smaller cells in shallow water and larger
+cells in deeper water. Here, $h$ refers to the water depth. Another useful criterion is known as the topographic length scale constraint,
+when one try to keep the ratio $\Delta h/h$ less than one. Here, $\Delta h$ equals the difference between the maximum depth of a triangle and the minimum depth and
+$h$ is the average depth. This
+criterion addresses the bathymetric slope and cells with a high value of $\Delta h/h$ indicate areas of steep bottom topography that will need to be
+more finely resolved.
+When refining the grid, one must balance the need to fully meet the refinement criteria with the desire to keep the triangle sizes from becoming too small.
+Thus, these criteria
+are generally imposed along with a minimum area constraint.
+The refinement process is repeated iteratively until a final grid with the appropriate resolution is obtained.
+
 \section{Numerical method}
 \label{sec:nummeth}
 
@@ -3504,7 +3654,7 @@
   \vec{e}^{(1)} = \frac{1}{D} ( e^2_{(2)},-e^1_{(2)} )^{\top}\, ,\, \, \vec{e}^{(2)} = \frac{1}{D} (-e^2_{(1)}, e^1_{(1)} )^{\top}\, , \, \,
   D = e^2_{(2)} e^1_{(1)} - e^2_{(1)} e^1_{(2)} \, .
 \end{equation}
-This space discretization is first order accurate and conserves action. The upwind difference scheme (\ref{eq:spacedisc}) is employed for two reasons.
+This space discretization is of lowest order accurate and conserves action. The upwind difference scheme (\ref{eq:spacedisc}) is employed for two reasons.
 First, it enforces the propagation of wave action to follow the characteristics.
 Second, it is monotone (i.e. guaranteeing $N > 0$ everywhere) and compact (i.e. operating on one triangle only), while sufficiently accurate for nearshore applications.
 Given the action densities $N^n_2$ and $N^n_3$ at vertices 2 and 3 of triangle $\triangle$123,
@@ -3568,43 +3718,178 @@
 This numerical process is iterated until an {\em a priori} convergence condition is satisfied. Here, the curvature-based stopping criteria
 as outlined in Section~\ref{sec:stop} will be applied.
 
-\section{Some notes on grid generation}
+\section{Interpolation at user-defined locations}
 
-We briefly outline some issues related to grid generation from a practical point of view.
-The process of grid generation can be difficult and time consuming. A common approach is proceeding from coarse to fine grid through refinement in various ways.
-Generally, one would like to have an optimal grid in which areas where the bathymetry or evolution of the waves change
-rapidly require a higher resolution than areas where the physics or depth changes less.
-This goes around by having an indication how to determine the refinement based on
-bathymetry or geometric variations through preliminary evaluations.
-To facilitate this procedure,
-a mesh generation package BatTri (Bilgili and Smith, 2003) is used.
-This grid generator is a public-domain, graphical Matlab interface to Triangle (Shewchuk, 1996).
-Triangle is a freely-distributed, two-dimensional Delaunay triangulator.
-Generation of the meshes presented in this paper was accomplished using BatTri.
-\nocite{Bil03S,She96}
+All the quantities deals with in SWAN are located at the vertices. Hence, due to the user-defined locations of the wave parameters, interpolations are required.
+Let parameter $\varphi_j = \varphi(\vec{x}_j)$ and Cartesian coordinates $\vec{x}_j = (x_j, y_j)$, with $j \in \{1,2,3\}$ indicating the vertices of cell $i$, be given.
+The vertices 1, 2 and 3 are ordered in a counterclockwise direction, see Figure~\ref{fig:gsunstruc}. The associated 3 edges are denoted as 12, 23 and 31.
 
-An important key ingredient for the preparation of the grid for the model domain is bathymetry data.
-Boundary nodes, segments and holes
-can be created from this data with the
-use of the mesh editing options of BatTri.
-After checking and improving grid quality, the final information on nodes and segments is forced into the
-triangulation of the domain. Triangle is called within BatTri to generate
-the actual grid. This triangulation includes only acute triangles.
+Linear interpolation, with $\vec{x}_0$ inside cell $i$ and $\varphi_0 = \varphi(\vec{x}_0)$, is given by
+\begin{equation}
+  \varphi(\vec{x}) = \varphi_0 + \nabla \varphi \cdot (\vec{x} - \vec{x}_0)
+\end{equation}
+where $\nabla \varphi$ is a constant vector inside cell $i$. We apply Green-Gauss reconstruction, i.e.,
+\begin{equation}
+  \nabla \varphi \approx \frac{1}{A_i} \int_{\triangle i} \nabla \varphi d\Omega = \frac{1}{A_i} \oint_{\partial \triangle i}
+  \varphi \vec{n} d\Gamma \approx \frac{1}{A_i} \sum_e \varphi_e \vec{n}_e
+\end{equation}
+where $A_i$ is the area of cell $i$ and the summation runs over the 3 edges $e \in \{12,23,31\}$ of cell $i$.
+The values $\varphi_e$ at edges are taken as averages:
+\begin{equation}
+   \varphi_{12} = \frac{1}{2} (\varphi_1+\varphi_2)\, , \quad
+   \varphi_{23} = \frac{1}{2} (\varphi_2+\varphi_3)\, , \quad
+   \varphi_{31} = \frac{1}{2} (\varphi_3+\varphi_1)\
+\end{equation}
+Furthermore, $\vec{n}_e$ is the outward pointing normal
+at edge $e$ and is obtained by rotating the edge over 90$^o$ in the clockwise
+direction. Hence,
+\begin{equation}
+  \vec{n}_{12} = R\vec{t}_{12}\, , \quad  R = \left( \begin{array}{ll} 0 & 1 \\ -1 & 0 \end{array} \right) \, , \quad \vec{t}_{12} = \vec{x}_2 - \vec{x}_1
+\end{equation}
+We also need the following identity
+\begin{equation}
+   \vec{n}_{12} + \vec{n}_{23} + \vec{n}_{31} = 0
+\end{equation}
+It is not difficult to show that
+\begin{eqnarray}
+  \nabla \varphi &=& \frac{1}{2A_i} \left [ \vec{n}_{12} (\varphi_1 - \varphi_3) + \vec{n}_{31} (\varphi_1 - \varphi_2) \right ] \nonumber \\
+                 &=& -\frac{1}{2A_i} \left [ \varphi_1 \vec{n}_{23} + \varphi_2 \vec{n}_{31} + \varphi_3 \vec{n}_{12} \right ]
+  \label{eq:gradcell}
+\end{eqnarray}
+or
+\begin{equation}
+  \frac{\partial \varphi}{\partial x} = \frac{1}{2A_i} \left [ \varphi_1 (y_2 - y_3) + \varphi_2 (y_3 - y_1) + \varphi_3 (y_1 - y_2) \right]
+  \label{eq:xgradcell}
+\end{equation}
+and
+\begin{equation}
+  \frac{\partial \varphi}{\partial y} = \frac{1}{2A_i} \left [ \varphi_1 (x_3 - x_2) + \varphi_2 (x_1 - x_3) + \varphi_3 (x_2 - x_1) \right]
+  \label{eq:ygradcell}
+\end{equation}
+The area $A_i$ of cell $i$ is given by $|\vec{t}_{12} \cdot \vec{n}_{13}|/2$. Hence, with
+\begin{equation}
+  \vec{n}_{13} = R\vec{t}_{13} = \left( \begin{array}{l} y_3 - y_1 \\ x_1 - x_3 \end{array} \right)
+\end{equation}
+we have
+\begin{equation}
+  A_i = \frac{1}{2} | (x_2 - x_1)(y_3 - y_1) - (x_3 - x_1)(y_2 - y_1) |
+\end{equation}
+\\[2ex]
+\noindent
+Alternatively, we may interpolate using the following relation
+\begin{equation}
+  \varphi(\vec{x}) = \sum_k \varphi_k \lambda_k(\vec{x}) = \varphi_1 \lambda_1 + \varphi_2 \lambda_2 +\varphi_3 \lambda_3
+\end{equation}
+where $\lambda_k$ is a linear shape function with the following properties:
+\begin{enumerate}
+  \item $\lambda_k$ is linear per cell and
+  \item $\lambda_k(\vec{x}_j) = \delta_{kj}$ with $\delta_{kj}$ is the Kronecker delta.
+\end{enumerate}
+We choose the following shape function
+\begin{equation}
+  \lambda_k(\vec{x}) = a_0^k + a_x^k x + a_y^k y
+\end{equation}
+and the coefficients $a$ follow from solving
+\begin{equation}
+  \left(
+     \begin{array}{lll}
+        1 & x_1 & y_1 \\
+        1 & x_2 & y_2 \\
+        1 & x_3 & y_3
+     \end{array}
+  \right)
+  \left(
+     \begin{array}{lll}
+        a_0^1 & a_0^2 & a_0^3 \\
+        a_x^1 & a_x^2 & a_x^3 \\
+        a_y^1 & a_y^2 & a_y^3
+     \end{array}
+  \right)
+  =
+  \left(
+     \begin{array}{lll}
+        1 & 0 & 0 \\
+        0 & 1 & 0 \\
+        0 & 0 & 1
+     \end{array}
+  \right)
+\end{equation}
 
-BatTri provides many pre-defined depth-dependent contraints for further mesh refinement.
-From a numerical point of view, mesh refinement is often directly related to properly resolve the shape of the wave,
-i.e. to keep the wavelength to grid size ratio relatively large.
-When wavelength decreases in shallower water, the grid size must decrease as well. Therefore, this criterion, called the $h$-refinement,
-has the effect of using smaller cells in shallow water and larger
-cells in deeper water. Here, $h$ refers to the water depth. Another useful criterion is known as the topographic length scale constraint,
-when one try to keep the ratio $\Delta h/h$ less than one. Here, $\Delta h$ equals the difference between the maximum depth of a triangle and the minimum depth and
-$h$ is the average depth. This
-criterion addresses the bathymetric slope and cells with a high value of $\Delta h/h$ indicate areas of steep bottom topography that will need to be
-more finely resolved.
-When refining the grid, one must balance the need to fully meet the refinement criteria with the desire to keep the triangle sizes from becoming too small.
-Thus, these criteria
-are generally imposed along with a minimum area constraint.
-The refinement process is repeated iteratively until a final grid with the appropriate resolution is obtained.
+\section{Computation of wave-induced force}
+
+{\tt FORCE} is the wave-driven stress, i.e. the force per unit surface driving the wave-driven current, expressed in N/m$^2$, is defined
+as the gradient of the radiation stresses:
+\begin{equation}
+  S_{xx} = \rho g \int \lfloor n \cos^2\theta + n - \frac{1}{2} \rfloor E d\sigma d\theta
+\end{equation}
+\begin{equation}
+  S_{xy} = S_{yx} = \rho g \int n \sin\theta \cos\theta E d\sigma d\theta
+\end{equation}
+\begin{equation}
+  S_{yy} = \rho g \int \lfloor n \sin^2\theta + n - \frac{1}{2} \rfloor E d\sigma d\theta
+\end{equation}
+with $n$ the ratio of group velocity and phase velocity.
+The force is then
+\begin{equation}
+  F_x = -\frac{\partial S_{xx}}{\partial x} - \frac{\partial S_{xy}}{\partial y}
+\end{equation}
+and
+\begin{equation}
+  F_y = -\frac{\partial S_{yx}}{\partial x} - \frac{\partial S_{yy}}{\partial y}
+\end{equation}
+\\[2ex]
+\noindent
+In order to compute the force in all internal vertices of the unstructured mesh, we consider a control volume (CV) as depicted in Figure~\ref{fig:cv}.
+\begin{figure}[htb]
+   \centerline{
+      \epsfig{file=CV.eps,height=7cm}
+              }
+      \caption{Control volume (centroid dual) of vertex is shaded. Some notation is introduced.}
+      \label{fig:cv}
+\end{figure}
+This CV is called centroid dual and is constructed by joining the centroids neighbouring the vertex under consideration. The set of CVs must fill
+the whole computational domain and must also be non-overlapping.
+In the following we use the numbering from Figure~\ref{fig:cv}.
+Let $\varphi$ be one of the
+radiation stresses $S_{xx}$, $S_{xy}$ and $S_{yy}$. The gradient of $\varphi$ is computed as follows:
+\begin{equation}
+  \nabla \varphi \approx \frac{1}{A_{\rm CV}} \sum_e \varphi_e \vec{n}_e
+\end{equation}
+where $A_{\rm CV}$ is the area of the CV and the summation runs over the associated edges $e$ of this CV.
+The values $\varphi_e$ at edges of the centroid dual are taken as averages, i.e. $(\varphi_0+\varphi_1)/2$, $(\varphi_1+\varphi_2)/2$, etc.
+Moreover, the value of the radiation stresses
+inside each triangle is simply the average of the radiation stresses in the associated vertices of the cell. Now, the derivatives of $\varphi$ inside
+CV are
+\begin{equation}
+  \frac{\partial \varphi}{\partial x} = \frac{1}{2A_{\rm CV}} \sum_{i=0}^{n-1} (\varphi_i+\varphi_{i+1})\,(y_{i+1}-y_i)
+  \label{eq:xgradcv}
+\end{equation}
+and
+\begin{equation}
+  \frac{\partial \varphi}{\partial y} = \frac{1}{2A_{\rm CV}} \sum_{i=0}^{n-1} (\varphi_i+\varphi_{i+1})\,(x_i-x_{i+1})
+  \label{eq:ygradcv}
+\end{equation}
+with $n$ the number of surrounding cells of the considered vertex and $\varphi_n = \varphi_0$, $x_n = x_0$ and $y_n = y_0$.
+The area of the CV is given by
+\begin{equation}
+  A_{\rm CV} = \frac{1}{2} \sum_{i=0}^{n-1} (x_i\,y_{i+1} - x_{i+1}\,y_i)
+\end{equation}
+
+\section{Calculation of diffusion-like terms}
+
+There are situations in which the following diffusion-like term need to be computed on unstructured meshes:
+\begin{equation}
+  \nabla \cdot \left ( \kappa \nabla \varphi \right )
+\end{equation}
+in vertices with $\kappa$ a space-varying diffusion coefficient (a tensor) and $\varphi$  a scalar defined in vertices. In SWAN, these are
+\begin{itemize}
+  \item the alleviation of the garden-sprinkler effect; see Eq. (\ref{eq:gsediff}), and
+  \item the computation of the diffraction parameter; see Eq. (\ref{eq3-44}).
+\end{itemize}
+We consider the centroid dual as shown in Figure~\ref{fig:cv}. The calculation consists of 3 steps. First, we compute the gradient of $\varphi$ inside
+each surrounding cell using expressions (\ref{eq:xgradcell}) and (\ref{eq:ygradcell}). Next, this gradient is multiplied with the appropriate diffusion
+coefficient $\kappa$ as given in the centroid. Finally, we compute the gradient of $\kappa \nabla \varphi$ inside the CV according to
+Eqs. (\ref{eq:xgradcv}) and (\ref{eq:ygradcv}).
 
 \chap{The overall solution algorithm} \label{ch:concl}
 This chapter is under preparation.
@@ -4137,7 +4422,7 @@
 Inst. Eng. Aust}., {\bf 29}, 157-161
 
 \bibitem{Nel94}
-Nelson, R.C., 1994: Depth limited wave heights in very flat regions, {\it Coastal Engineering}, {\bf 23}, 43-59
+Nelson, R.C., 1994: Depth limited design wave heights in very flat regions, {\it Coastal Engineering}, {\bf 23}, 43-59
 
 \bibitem{Nel97}
 Nelson, R.C., 1997: Height limits in top down and bottom up wave environments, {\it Coastal Engineering},
@@ -4244,6 +4529,10 @@
 Rogers, W.E., P.A. Hwang and D.W. Wang, 2003: Investigation of wave growth and decay in the SWAN model:
 three regional-scale applications, {\it J. Phys. Oceanogr.}, {\bf 33}, 366-389.
 
+\bibitem{Rue03WS}
+Ruessink, B.G., D.J.R. Walstra and H.N. Southgate, 2003: Calibration and verification of a parametric wave model
+on barred beaches, {\it Coastal Engineering}, {\bf 48}, 139-149
+
 \bibitem{Sak83KI}
 Sakai, T., M. Koseki and Y. Iwagaki, 1983: Irregular wave refraction due to current, {\it J. of Hydr. Eng.}, ASCE,
 {\bf 109}, 9, 1203-1215
@@ -4442,92 +4731,97 @@
 \begin{theindex}
 \addcontentsline{toc}{chapter}{Index}
   \item ambient, \hyperpage{2}, \hyperpage{7}, \hyperpage{11, 12}, 
-		\hyperpage{29}, \hyperpage{91}
+		\hyperpage{30}, \hyperpage{97}
 
   \indexspace
 
-  \item bathymetry, \hyperpage{24}, \hyperpage{77}, \hyperpage{81, 82}, 
-		\hyperpage{93}
-  \item bottom, \hyperpage{1--3}, \hyperpage{12, 13}, \hyperpage{15}, 
-		\hyperpage{22--24}, \hyperpage{30}, \hyperpage{80--82}, 
-		\hyperpage{88}, \hyperpage{90}, \hyperpage{93--98}
+  \item bathymetry, \hyperpage{24}, \hyperpage{77}, \hyperpage{79}, 
+		\hyperpage{99}
+  \item bottom, \hyperpage{1--3}, \hyperpage{12}, \hyperpage{14, 15}, 
+		\hyperpage{22--24}, \hyperpage{31}, \hyperpage{79}, 
+		\hyperpage{82, 83}, \hyperpage{94}, \hyperpage{96}, 
+		\hyperpage{99--104}
   \item boundary, \hyperpage{3, 4}, \hyperpage{13}, \hyperpage{15}, 
 		\hyperpage{43}, \hyperpage{59, 60}, \hyperpage{62, 63}, 
 		\hyperpage{65--69}, \hyperpage{78}, \hyperpage{80}, 
-		\hyperpage{92}, \hyperpage{94}
-  \item breaking, \hyperpage{1}, \hyperpage{3}, \hyperpage{12, 13}, 
-		\hyperpage{15--17}, \hyperpage{21}, \hyperpage{23}, 
-		\hyperpage{42}, \hyperpage{56}, \hyperpage{85, 86}, 
-		\hyperpage{88, 89}, \hyperpage{93}, \hyperpage{97}
+		\hyperpage{82}, \hyperpage{98}, \hyperpage{100}
+  \item breaking, \hyperpage{1}, \hyperpage{3}, \hyperpage{12}, 
+		\hyperpage{14--17}, \hyperpage{21}, \hyperpage{23--25}, 
+		\hyperpage{42}, \hyperpage{56, 57}, \hyperpage{91, 92}, 
+		\hyperpage{94, 95}, \hyperpage{99}, 
+		\hyperpage{103, 104}
 
   \indexspace
 
   \item Cartesian, \hyperpage{3}, \hyperpage{10--12}, \hyperpage{64}, 
-		\hyperpage{79}
+		\hyperpage{81}, \hyperpage{83}
   \item co-ordinate, \hyperpage{3}, \hyperpage{7}, \hyperpage{10--12}, 
-		\hyperpage{27}, \hyperpage{48--50}, \hyperpage{66}, 
+		\hyperpage{27, 28}, \hyperpage{48--50}, \hyperpage{66}, 
 		\hyperpage{68}
   \item coastal, \hyperpage{1, 2}, \hyperpage{15}, \hyperpage{45}, 
-		\hyperpage{73}, \hyperpage{77}, \hyperpage{86, 87}, 
-		\hyperpage{89}, \hyperpage{91--93}, \hyperpage{95, 96}
+		\hyperpage{73}, \hyperpage{77}, \hyperpage{92, 93}, 
+		\hyperpage{95}, \hyperpage{97--99}, 
+		\hyperpage{101, 102}
   \item convergence, \hyperpage{36, 37}, \hyperpage{41}, 
-		\hyperpage{43--47}, \hyperpage{75, 76}, \hyperpage{81}, 
-		\hyperpage{95}, \hyperpage{99}
+		\hyperpage{43--47}, \hyperpage{75, 76}, \hyperpage{83}, 
+		\hyperpage{101}, \hyperpage{105}
   \item Courant, \hyperpage{40}
   \item current, \hyperpage{2, 3}, \hyperpage{7}, \hyperpage{10--13}, 
-		\hyperpage{15}, \hyperpage{29}, \hyperpage{32, 33}, 
-		\hyperpage{37, 38}, \hyperpage{40}, \hyperpage{43}, 
-		\hyperpage{49}, \hyperpage{55}, \hyperpage{61}, 
-		\hyperpage{74}, \hyperpage{76}, \hyperpage{88}, 
-		\hyperpage{91}, \hyperpage{95--98}
+		\hyperpage{15}, \hyperpage{30}, \hyperpage{32}, 
+		\hyperpage{34}, \hyperpage{37, 38}, \hyperpage{40}, 
+		\hyperpage{43}, \hyperpage{49}, \hyperpage{55}, 
+		\hyperpage{61}, \hyperpage{74}, \hyperpage{76}, 
+		\hyperpage{85}, \hyperpage{94}, \hyperpage{97}, 
+		\hyperpage{101--104}
   \item curvi-linear, \hyperpage{2, 3}, \hyperpage{48, 49}, 
-		\hyperpage{62}, \hyperpage{68, 69}, \hyperpage{87}, 
-		\hyperpage{91}
+		\hyperpage{51}, \hyperpage{62}, \hyperpage{68, 69}, 
+		\hyperpage{77}, \hyperpage{93}, \hyperpage{97}
 
   \indexspace
 
-  \item dam, \hyperpage{16}, \hyperpage{29--31}, \hyperpage{91}, 
+  \item dam, \hyperpage{16}, \hyperpage{30--32}, \hyperpage{97}, 
+		\hyperpage{105}
+  \item diffraction, \hyperpage{3}, \hyperpage{10}, \hyperpage{30}, 
+		\hyperpage{32--34}, \hyperpage{87}, \hyperpage{92}, 
 		\hyperpage{98, 99}
-  \item diffraction, \hyperpage{3}, \hyperpage{10}, \hyperpage{30--33}, 
-		\hyperpage{86}, \hyperpage{92, 93}
-  \item diffusion, \hyperpage{38--40}
-  \item dissipation, \hyperpage{1}, \hyperpage{3}, \hyperpage{12--16}, 
-		\hyperpage{20, 21}, \hyperpage{23, 24}, 
-		\hyperpage{56, 57}, \hyperpage{85--87}, \hyperpage{90}, 
-		\hyperpage{95}, \hyperpage{97, 98}
+  \item diffusion, \hyperpage{38--40}, \hyperpage{86, 87}
+  \item dissipation, \hyperpage{1}, \hyperpage{3}, \hyperpage{12}, 
+		\hyperpage{14--16}, \hyperpage{20--24}, 
+		\hyperpage{56, 57}, \hyperpage{78}, \hyperpage{91--93}, 
+		\hyperpage{96}, \hyperpage{101}, \hyperpage{103, 104}
 
   \indexspace
 
-  \item filter, \hyperpage{13}, \hyperpage{18}, \hyperpage{32}
-  \item flow, \hyperpage{7}, \hyperpage{61}, \hyperpage{87}, 
-		\hyperpage{92}, \hyperpage{94}, \hyperpage{99}
-  \item force, \hyperpage{33}, \hyperpage{49, 50}, \hyperpage{61--63}, 
-		\hyperpage{68, 69}, \hyperpage{79}, \hyperpage{82}, 
-		\hyperpage{88}
+  \item filter, \hyperpage{13}, \hyperpage{18}, \hyperpage{33}
+  \item flow, \hyperpage{7}, \hyperpage{61}, \hyperpage{93}, 
+		\hyperpage{98}, \hyperpage{100}, \hyperpage{105}
+  \item force, \hyperpage{34}, \hyperpage{49, 50}, \hyperpage{61--63}, 
+		\hyperpage{68, 69}, \hyperpage{79}, \hyperpage{81}, 
+		\hyperpage{85}, \hyperpage{94}
   \item frequency, \hyperpage{7--9}, \hyperpage{11}, \hyperpage{13, 14}, 
-		\hyperpage{16}, \hyperpage{18--22}, \hyperpage{24, 25}, 
-		\hyperpage{28, 29}, \hyperpage{35--37}, \hyperpage{40}, 
+		\hyperpage{16--23}, \hyperpage{25, 26}, 
+		\hyperpage{28--30}, \hyperpage{35--37}, \hyperpage{40}, 
 		\hyperpage{43--46}, \hyperpage{53--55}, 
-		\hyperpage{59, 60}, \hyperpage{85}, \hyperpage{93}, 
-		\hyperpage{99}
-  \item friction, \hyperpage{1}, \hyperpage{3}, \hyperpage{12, 13}, 
-		\hyperpage{15}, \hyperpage{18}, \hyperpage{22, 23}, 
-		\hyperpage{90}, \hyperpage{92, 93}, \hyperpage{95}, 
-		\hyperpage{97, 98}
+		\hyperpage{59, 60}, \hyperpage{91}, \hyperpage{99}, 
+		\hyperpage{105}
+  \item friction, \hyperpage{1}, \hyperpage{3}, \hyperpage{12--15}, 
+		\hyperpage{18}, \hyperpage{22, 23}, \hyperpage{96}, 
+		\hyperpage{98, 99}, \hyperpage{101}, 
+		\hyperpage{103, 104}
 
   \indexspace
 
-  \item garden-sprinkler, \hyperpage{39, 40}
+  \item garden-sprinkler, \hyperpage{39, 40}, \hyperpage{87}
 
   \indexspace
 
-  \item harbour, \hyperpage{3}, \hyperpage{33}
+  \item harbour, \hyperpage{3}, \hyperpage{34}
 
   \indexspace
 
-  \item initial, \hyperpage{13}, \hyperpage{15}, \hyperpage{28}, 
+  \item initial, \hyperpage{13}, \hyperpage{15}, \hyperpage{29}, 
 		\hyperpage{59, 60}
-  \item island, \hyperpage{77}
+  \item island, \hyperpage{77, 78}
 
   \indexspace
 
@@ -4541,20 +4835,21 @@
 
   \indexspace
 
-  \item obstacle, \hyperpage{3}, \hyperpage{29--31}, \hyperpage{33}, 
-		\hyperpage{38, 39}, \hyperpage{50--52}
+  \item obstacle, \hyperpage{3}, \hyperpage{30--34}, \hyperpage{38, 39}, 
+		\hyperpage{50--52}
   \item ocean, \hyperpage{1}, \hyperpage{7, 8}, \hyperpage{11}, 
-		\hyperpage{29}, \hyperpage{33}, \hyperpage{38}, 
-		\hyperpage{40}, \hyperpage{86}, \hyperpage{90--94}, 
-		\hyperpage{96}, \hyperpage{98}
+		\hyperpage{30}, \hyperpage{34}, \hyperpage{38}, 
+		\hyperpage{40}, \hyperpage{92}, \hyperpage{96--100}, 
+		\hyperpage{102}, \hyperpage{104, 105}
 
   \indexspace
 
   \item propagation, \hyperpage{1--4}, \hyperpage{10--12}, 
-		\hyperpage{32}, \hyperpage{37}, \hyperpage{39}, 
+		\hyperpage{32, 33}, \hyperpage{37}, \hyperpage{39}, 
 		\hyperpage{41, 42}, \hyperpage{48, 49}, \hyperpage{51}, 
-		\hyperpage{73, 74}, \hyperpage{77}, \hyperpage{79, 80}, 
-		\hyperpage{86--89}, \hyperpage{95}, \hyperpage{97}
+		\hyperpage{73, 74}, \hyperpage{77}, \hyperpage{79}, 
+		\hyperpage{81, 82}, \hyperpage{92--95}, 
+		\hyperpage{101}, \hyperpage{103}
 
   \indexspace
 
@@ -4562,64 +4857,64 @@
 
   \indexspace
 
-  \item recti-linear, \hyperpage{69}
-  \item reflection, \hyperpage{3}, \hyperpage{29}, \hyperpage{31}, 
-		\hyperpage{33}, \hyperpage{35}, \hyperpage{52}
+  \item recti-linear, \hyperpage{69}, \hyperpage{77}
+  \item reflection, \hyperpage{3}, \hyperpage{30}, \hyperpage{32, 33}, 
+		\hyperpage{35}, \hyperpage{52}
   \item refraction, \hyperpage{1}, \hyperpage{3}, \hyperpage{11}, 
-		\hyperpage{31, 32}, \hyperpage{43}, \hyperpage{59}, 
-		\hyperpage{80}, \hyperpage{86, 87}, \hyperpage{92, 93}, 
-		\hyperpage{96}, \hyperpage{98}
+		\hyperpage{32}, \hyperpage{43}, \hyperpage{59}, 
+		\hyperpage{82}, \hyperpage{92, 93}, \hyperpage{98, 99}, 
+		\hyperpage{102}, \hyperpage{104}
   \item regular, \hyperpage{2}, \hyperpage{4}, \hyperpage{7}, 
-		\hyperpage{15}, \hyperpage{30}, \hyperpage{80}, 
-		\hyperpage{96}, \hyperpage{98}
+		\hyperpage{15}, \hyperpage{30}, \hyperpage{82}, 
+		\hyperpage{102}, \hyperpage{104}
 
   \indexspace
 
-  \item set-up, \hyperpage{3, 4}, \hyperpage{24}, \hyperpage{33}, 
-		\hyperpage{61}, \hyperpage{69}, \hyperpage{85}
+  \item set-up, \hyperpage{3, 4}, \hyperpage{24}, \hyperpage{34}, 
+		\hyperpage{61}, \hyperpage{69}, \hyperpage{91}
   \item shoaling, \hyperpage{1}, \hyperpage{3}, \hyperpage{15}, 
-		\hyperpage{89}
+		\hyperpage{25}, \hyperpage{95}
   \item SORDUP, \hyperpage{38}, \hyperpage{40}, \hyperpage{73}
   \item spherical, \hyperpage{3}, \hyperpage{11, 12}, 
-		\hyperpage{80, 81}
-  \item stability, \hyperpage{21}, \hyperpage{35, 36}, 
-		\hyperpage{39, 40}, \hyperpage{43}, \hyperpage{45, 46}, 
-		\hyperpage{78}, \hyperpage{80}
+		\hyperpage{82, 83}
+  \item stability, \hyperpage{35, 36}, \hyperpage{40}, \hyperpage{43}, 
+		\hyperpage{45, 46}, \hyperpage{56}, \hyperpage{58}, 
+		\hyperpage{80}, \hyperpage{82}
   \item stationary, \hyperpage{2, 3}, \hyperpage{7}, \hyperpage{36--38}, 
 		\hyperpage{40}, \hyperpage{44}, \hyperpage{49}, 
-		\hyperpage{60, 61}, \hyperpage{87, 88}, \hyperpage{91}, 
-		\hyperpage{99}
-  \item steepness, \hyperpage{14}, \hyperpage{19--21}, \hyperpage{24}, 
-		\hyperpage{29}
-  \item swell, \hyperpage{14, 15}, \hyperpage{21}, \hyperpage{23}, 
-		\hyperpage{87}, \hyperpage{90}, \hyperpage{98}
+		\hyperpage{60, 61}, \hyperpage{93, 94}, \hyperpage{97}, 
+		\hyperpage{105}
+  \item steepness, \hyperpage{14}, \hyperpage{20}, \hyperpage{24}, 
+		\hyperpage{30}
+  \item swell, \hyperpage{14, 15}, \hyperpage{20}, \hyperpage{22}, 
+		\hyperpage{93}, \hyperpage{96}, \hyperpage{104}
 
   \indexspace
 
   \item triads, \hyperpage{16}
-  \item triangular, \hyperpage{78}
+  \item triangular, \hyperpage{78}, \hyperpage{80}
 
   \indexspace
 
-  \item unstructured, \hyperpage{4}, \hyperpage{77}, \hyperpage{80}, 
-		\hyperpage{86}
+  \item unstructured, \hyperpage{4}, \hyperpage{77, 78}, \hyperpage{82}, 
+		\hyperpage{85, 86}, \hyperpage{92}
 
   \indexspace
 
   \item WAM, \hyperpage{1}, \hyperpage{13, 14}, \hyperpage{18--20}, 
-		\hyperpage{25}, \hyperpage{35, 36}, \hyperpage{44, 45}, 
-		\hyperpage{55}, \hyperpage{59}, \hyperpage{90, 91}, 
-		\hyperpage{96}, \hyperpage{98}
-  \item WAVEWATCH, \hyperpage{1}, \hyperpage{59}, \hyperpage{97}
-  \item whitecapping, \hyperpage{3}, \hyperpage{12--14}, \hyperpage{16}, 
-		\hyperpage{19--22}, \hyperpage{29}, \hyperpage{90}, 
-		\hyperpage{92}, \hyperpage{97, 98}
+		\hyperpage{26}, \hyperpage{35, 36}, \hyperpage{44, 45}, 
+		\hyperpage{55}, \hyperpage{59}, \hyperpage{96, 97}, 
+		\hyperpage{102}, \hyperpage{104}
+  \item WAVEWATCH, \hyperpage{1}, \hyperpage{59}, \hyperpage{103}
+  \item whitecapping, \hyperpage{3}, \hyperpage{12}, \hyperpage{14}, 
+		\hyperpage{16}, \hyperpage{19--22}, \hyperpage{30}, 
+		\hyperpage{96}, \hyperpage{98}, \hyperpage{103, 104}
   \item wind, \hyperpage{1--4}, \hyperpage{7}, \hyperpage{12--15}, 
-		\hyperpage{17--22}, \hyperpage{29}, \hyperpage{35}, 
+		\hyperpage{17--22}, \hyperpage{30}, \hyperpage{35}, 
 		\hyperpage{37--42}, \hyperpage{44, 45}, \hyperpage{50}, 
 		\hyperpage{60}, \hyperpage{68}, \hyperpage{77}, 
-		\hyperpage{79}, \hyperpage{85--87}, \hyperpage{90}, 
-		\hyperpage{92--99}
+		\hyperpage{81}, \hyperpage{91--93}, \hyperpage{96}, 
+		\hyperpage{98--105}
 
 \end{theindex}
 
--- SwanTranspAc.ftn90	2009-03-03 17:49:37.000000000 +0100
+++ SwanTranspAc.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -1,8 +1,9 @@
 subroutine SwanTranspAc ( amat  , rhs   , leakcf, ac2   , ac1   , &
-                          cax   , cay   , cad   , cas   , anybin, &
-                          rdx   , rdy   , obredf, spcsig, idcmin, &
-                          idcmax, iscmin, iscmax, iddlow, iddtop, &
-                          isslow, isstop, trac0 , trac1 )
+                          cgo   , cax   , cay   , cad   , cas   , &
+                          anybin, rdx   , rdy   , spcsig, spcdir, &
+                          obredf, idcmin, idcmax, iscmin, iscmax, &
+                          iddlow, iddtop, isslow, isstop, trac0 , &
+                          trac1 )
 !
 !   --|-----------------------------------------------------------|--
 !     | Delft University of Technology                            |
@@ -15,7 +16,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -39,8 +40,9 @@
 !
 !   Updates
 !
-!   40.80, August 2007: New subroutine
-!   40.85, August 2008: add tranport terms for output purposes
+!   40.80,   August 2007: New subroutine
+!   40.85,   August 2008: add tranport terms for output purposes
+!   41.00, February 2009: add GSE correction
 !
 !   Purpose
 !
@@ -80,11 +82,18 @@
     real, dimension(MDC,MSC), intent(in)        :: cas    ! wave transport velocity in sigma-direction
     real, dimension(MDC,MSC,ICMAX), intent(in)  :: cax    ! wave transport velocity in x-direction
     real, dimension(MDC,MSC,ICMAX), intent(in)  :: cay    ! wave transport velocity in y-direction
+    real, dimension(MSC,ICMAX), intent(in)      :: cgo    ! group velocity
     real, dimension(MDC,MSC), intent(out)       :: leakcf ! leak coefficient
     real, dimension(MDC,MSC,2), intent(in)      :: obredf ! action reduction coefficient based on transmission
     real, dimension(2), intent(in)              :: rdx    ! first component of contravariant base vector rdx(b) = a^(b)_1
     real, dimension(2), intent(in)              :: rdy    ! second component of contravariant base vector rdy(b) = a^(b)_2
     real, dimension(MDC,MSC), intent(out)       :: rhs    ! right-hand side of system of equations in (sigma,theta) space
+    real, dimension(MDC,6), intent(in)          :: spcdir ! (*,1): spectral direction bins (radians)
+                                                          ! (*,2): cosine of spectral directions
+                                                          ! (*,3): sine of spectral directions
+                                                          ! (*,4): cosine^2 of spectral directions
+                                                          ! (*,5): cosine*sine of spectral directions
+                                                          ! (*,6): sine^2 of spectral directions
     real, dimension(MSC), intent(in)            :: spcsig ! relative frequency bins
     real, dimension(MDC,MSC,MTRNP), intent(out) :: trac0  ! explicit part of propagation in present vertex for output purposes
     real, dimension(MDC,MSC,MTRNP), intent(out) :: trac1  ! implicit part of propagation in present vertex for output purposes
@@ -120,6 +129,10 @@
     call SwanTranspX ( amat   , rhs  , ac2   , ac1   , cax   , cay   , &
                        rdx    , rdy  , obredf, idcmin, idcmax, isslow, &
                        isstop , trac0, trac1 )
+    !
+    ! add GSE correction, if appropriate
+    !
+    if ( WAVAGE > 0. ) call SwanGSECorr ( rhs, ac2, cgo, spcdir, idcmin, idcmax, isslow, isstop, trac0 )
 !TIMG    call SWTSTO(140)
     !
     ! compute transport in theta space
--- SwanTranspX.ftn90	2009-03-03 17:49:37.000000000 +0100
+++ SwanTranspX.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -13,7 +13,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -42,7 +42,7 @@
 !
 !   Purpose
 !
-!   Computes transport in x-y space using the first order upwind scheme
+!   Computes transport in x-y space using the lowest order upwind scheme
 !
 !   Modules used
 !
@@ -125,7 +125,7 @@
        do iddum = idcmin(is), idcmax(is)
           id = mod ( iddum - 1 + MDC , MDC ) + 1
           !
-          ! compute the contributions based on the first order upwind scheme
+          ! compute the contributions based on the lowest order upwind scheme
           !
           asum = (rdx(1)+rdx(2))*cax(id,is,1) + (rdy(1)+rdy(2))*cay(id,is,1)
           !
--- swanuse.tex	2009-03-03 17:49:37.000000000 +0100
+++ swanuse.tex	2009-03-03 17:52:17.000000000 +0100
@@ -32,7 +32,7 @@
 \end{center}
 \vfill
 \begin{center}
-{\Large\bf SWAN Cycle III version 40.72A}
+{\Large\bf SWAN Cycle III version 40.72AB}
 \end{center}
 
 \cleardoublepage
@@ -54,14 +54,14 @@
              & & The Netherlands \\
              & & \\
 e-mail       &:& swan-info-citg@tudelft.nl \\
-home page    &:& \hl{http://www.fluidmechanics.tudelft.nl/swan/index.htm}
+home page    &:& \hl{http://www.swan.tudelft.nl}
 \end{tabular}
 \end{table}
 
 \vfill
 
 \noindent
-Copyright (c) 2008 Delft University of Technology.
+Copyright (c) 2009 Delft University of Technology.
 \\[2ex]
 \noindent
 Permission is granted to copy, distribute and/or modify this document
@@ -89,14 +89,13 @@
 \chap{Introduction} \label{ch:intro}
 \pagenumbering{arabic}
 
-The information about the SWAN package is distributed over five different documents. This User Manual describes
+The information about the SWAN package is distributed over four different documents. This User Manual describes
 the complete input and usage of the SWAN package. The Implementation Manual explains the installation procedure
-of SWAN on a single- or multi-processor machine with shared or distributed memory. The System documentation outlines
-the internals of the program and discusses program maintenance. The Programming rules is meant for programmers
-who want to develop SWAN. The Scientific/Technical documentation discusses the physical and mathematical details
-and the discretizations that are used in the SWAN program. The mapping of these numerical techniques in SWAN code
-is also discussed.
-\nocite{Impman,Sysdoc,Progrul,Techdoc}
+of SWAN on a single- or multi-processor machine with shared or distributed memory. The Programming rules is meant
+for programmers who want to develop SWAN. The Scientific/Technical documentation discusses the physical and
+mathematical details and the discretizations that are used in the SWAN program. The mapping of these numerical
+techniques in SWAN code is also discussed.
+\nocite{Impman,Progrul,Techdoc}
 \\[2ex]
 \noindent
 In Chapter~\ref{ch:defin} some general definitions and remarks concerning the usage of SWAN, the treatment of grids,
@@ -736,12 +735,16 @@
                         & Snyder {\it et al}. (1981)             &          &          & $\times$ \\
 \cline{2-5}
                         & Janssen (1989, 1991)                   &          &          & $\times$ \\
+\cline{2-5}
+                        & Yan (1987) (modified)                  &          &          & $\times$ \\
 \hline
 Whitecapping            & Holthuijsen and De Boer (1988)         & $\times$ & $\times$ &          \\
 \cline{2-5}
                         & Komen {\it et al}. (1984)              &          &          & $\times$ \\
 \cline{2-5}
                         & Janssen (1991)                         &          &          & $\times$ \\
+\cline{2-5}
+                        & Alves and Banner (2003)                &          &          & $\times$ \\
 \hline
 Quadruplets             & Hasselmann {\it et al}. (1985)         &          &          & $\times$ \\
 \hline
@@ -895,7 +898,7 @@
 {\tt GEN1}     \> SWAN runs in first generation mode\\
 {\tt GEN2}     \> SWAN runs in second generation mode\\
 {\tt GEN3}     \> SWAN runs in third generation mode\\
-{\tt WCAPPING} \> activates cumulative steepness method for whitecapping\\
+%{\tt WCAPPING} \> activates cumulative steepness method for whitecapping\\
 {\tt QUAD}     \> controls the computation of quadruplets\\
 {\tt BREAKING} \> activates dissipation by depth-induced wave breaking\\
 {\tt FRICTION} \> activates dissipation by bottom friction\\
@@ -1743,7 +1746,8 @@
                        =6: \> as {\tt [idla]}=5 but a new column in the map need not start on a new line\+\\
                               in the file.\-\\
                        \poptabs
-                       Default: {\tt [idla]}=1.\-\\
+                       Default: {\tt [idla]}=1.\\
+                       ONLY MEANT FOR STRUCTURED GRIDS.\-\\
 {\tt [nhedf]}       \> is the number of header lines at the start of the file. The text in the header\+\\
                        lines is reproduced in the print file created by SWAN (see Section~\ref{sec:prtfile}). The\\
                        file may start with more header lines than {\tt [nhedf]} because the start of the\\
@@ -1895,7 +1899,7 @@
           <                                                    >
            |                                                  |
            |           | PAR  < [len] [hs] [per] [dir] [dd] > |
-           |VARiable  <                                       |
+           | VARiable <                                       |
            |           | FILE < [len] 'fname' [seq] >         |
 \end{verbatim}
 \linecmd
@@ -1965,7 +1969,7 @@
                         ONLY MEANT FOR STRUCTURED GRIDS.\-\\
 {\tt [k]}            \> index of boundary vertex of the segment. This can be obtained in a grid\+\\
                         generator file ({\tt fort.14}, {\tt .node} and {\tt.n} files of ADCIRC, Triangle and\\
-                        Easymesh, respectively).\\
+                        Easymesh, respectively). The order must be counterclockwise!\\
                         ONLY MEANT FOR UNSTRUCTURED MESHES.\-\\
 {\tt CONSTANT}       \> with this option the wave spectra are constant along the side or segment.\\
 {\tt VARIABLE}       \> with this option the wave spectra can vary along the side or segment. The\+\\
@@ -2411,26 +2415,26 @@
                  Default: {\tt [a]} = 0.0015.\-\\
 \end{tabbing}
 
-\idxcmd{WCAPPING}
-\linecmd
-\begin{verbatim}
-WCAPping CSM [cst] [pow]
-\end{verbatim}
-\linecmd
-
-\noindent
-With this command the user wants to choose the Cumulative Steepness Method (CSM) for approximating whitecapping (see Scientific/Technical documentation)
-and {\bf not} the formulation of Komen {\it et al.} (1984) and {\bf not} Janssen (1991a).
-\begin{tabbing}
-xxxxxxxxxxxx\= \kill
-{\tt [cst]} \> the tuneable coefficient $C_{wc}^{st}$\+\\
-               Default: {\tt [cst]} = 4.0.\-\\
-{\tt [pow]} \> power $m$.\+\\
-               Default: {\tt [pow]} = 2.0.\-\\
-\end{tabbing}
-Note that the CSM method in SWAN is still in its experimental phase. Its results are promising, but the method still suffers some numerical
-problems.
-
+%\idxcmd{WCAPPING}
+%\linecmd
+%\begin{verbatim}
+%WCAPping CSM [cst] [pow]
+%\end{verbatim}
+%\linecmd
+%
+%\noindent
+%With this command the user wants to choose the Cumulative Steepness Method (CSM) for approximating whitecapping (see Scientific/Technical documentation)
+%and {\bf not} the formulation of Komen {\it et al.} (1984) and {\bf not} Janssen (1991a).
+%\begin{tabbing}
+%xxxxxxxxxxxx\= \kill
+%{\tt [cst]} \> the tuneable coefficient $C_{wc}^{st}$\+\\
+%               Default: {\tt [cst]} = 4.0.\-\\
+%{\tt [pow]} \> power $m$.\+\\
+%               Default: {\tt [pow]} = 2.0.\-\\
+%\end{tabbing}
+%Note that the CSM method in SWAN is still in its experimental phase. Its results are promising, but the method still suffers some numerical
+%problems.
+%
 \idxcmd{QUADRUPL}
 \linecmd
 \begin{verbatim}
@@ -2478,7 +2482,7 @@
 \idxcmd{BREAKING}
 \linecmd
 \begin{verbatim}
-BREaking CONSTANT [alpha] [gamma]
+BREaking  CONstant [alpha] [gamma]           
 \end{verbatim}
 \linecmd
 
@@ -2500,9 +2504,11 @@
 \idxcmd{FRICTION}
 \linecmd
 \begin{verbatim}
-           | -> JONswap [cfjon]
-           |
-FRICTION  <     COLLins [cfw]
+           |             | -> CONstant [cfjon]
+           | -> JONswap <
+           |             |    VARiable
+FRICtion  <
+           |    COLLins [cfw]
            |
            |    MADsen  [kn]
 \end{verbatim}
@@ -2513,32 +2519,44 @@
 \\[2ex]
 \noindent
 In SWAN three different formulations are available, i.e., that of Hasselmann et al. (1973, JONSWAP),
-Collins (1972), and Madsen et al. (1988). The default option is: {\tt JONSWAP}.
-\begin{tabbing}
-xxxxxxxxxxxx\= \kill
-{\tt JONSWAP} \> indicates that the semi-empirical expression derived from the JONSWAP results\+\\
-                 for bottom friction dissipation (Hasselmann et al., 1973, JONSWAP) should be\\
-                 activated. This option is default.\-\\
-{\tt [cfjon]} \> coefficient of the JONSWAP formulation. {\tt [cfjon]} is equal to 0.038$m^2 s^{-3}$ for\+\\
-                 swell conditions (Hasselmann et al., 1973) and equal to 0.067$m^2 s^{-3}$ for wind\\
-                 sea conditions.\\
-                 Default: {\tt [cfjon]} = 0.067.\-\\
-{\tt COLLINS} \> indicates that the expression of Collins (1972) should be activated.\\
-{\tt [cfw]}   \> Collins bottom friction coefficient.\+\\
-                 Default: {\tt [cfw]} = 0.015.\\
-                 Note that {\tt [cfw]} is allowed to vary over the computational region; in that\\
-                 case use the commands {\tt INPGRID FRICTION} and {\tt READINP FRICTION} to define\\
-                 and read the friction data. The command {\tt FRICTION} is still required to define\\
-                 the type of friction expression. The value of {\tt [cfw]} in this command is then\\
-                 not required (it will be ignored).\-\\
-{\tt MADSEN}  \> indicates that the expression of Madsen et al. (1988) should be activated.\\
-{\tt [kn]}    \> equivalent roughness length scale of the bottom (in m).\+\\
-                 Default: {\tt [kn]} = 0.05.\\
-                 Note that {\tt [kn]} is allowed to vary over the computational region; in that case\\
-                 use the commands {\tt INPGRID FRICTION} and {\tt READINP FRICTION} to define and read\\
-                 the friction data. This command {\tt FRICTION} is still required to define the type of\\
-                 friction expression. The value of {\tt [kn]} in this command is then not required\\
-                 (it will be ignored).\-\\
+Collins (1972), and Madsen et al. (1988). The JONSWAP formulation can be applied with a constant friction
+coefficient or with a varying friction coefficient that depends on the frequency-dependent directional
+spreading. In the latter formulation the friction coefficient is set to 0.038$m^2 s^{-3}$ for frequencies
+with directional spreading lower than or equal to 10$^o$ and it is set to 0.067$m^2 s^{-3}$ for frequencies
+with a directional spreading higher than or equal to 30$^o$. For intermediate values of the directional
+spreading the friction coefficient varies linearly between these two values.
+\\[2ex]
+\noindent
+The default option is: {\tt JONSWAP} with a constant friction coefficient.
+\begin{tabbing}
+xxxxxxxxxxxx\= \kill
+{\tt JONSWAP}  \> indicates that the semi-empirical expression derived from the JONSWAP results\+\\
+                  for bottom friction dissipation (Hasselmann et al., 1973, JONSWAP) should be\\
+                  activated. This option is default.\-\\
+{\tt CONSTANT} \> this default option indicates that the JONSWAP coefficient is constant.\\
+{\tt [cfjon]}  \> coefficient of the JONSWAP formulation. {\tt [cfjon]} is equal to 0.038$m^2 s^{-3}$ for\+\\
+                  swell conditions (Hasselmann et al., 1973) and equal to 0.067$m^2 s^{-3}$ for wind\\
+                  sea conditions.\\
+                  Default: {\tt [cfjon]} = 0.067.\-\\
+{\tt VARIABLE} \> indicates that the frequency-dependent JONSWAP formulation is used, with the\+\\
+                  friction coefficient depending on the frequency-dependent directional spreading\\
+                  as described above.\-\\
+{\tt COLLINS}  \> indicates that the expression of Collins (1972) should be activated.\\
+{\tt [cfw]}    \> Collins bottom friction coefficient.\+\\
+                  Default: {\tt [cfw]} = 0.015.\\
+                  Note that {\tt [cfw]} is allowed to vary over the computational region; in that\\
+                  case use the commands {\tt INPGRID FRICTION} and {\tt READINP FRICTION} to define\\
+                  and read the friction data. The command {\tt FRICTION} is still required to define\\
+                  the type of friction expression. The value of {\tt [cfw]} in this command is then\\
+                  not required (it will be ignored).\-\\
+{\tt MADSEN}   \> indicates that the expression of Madsen et al. (1988) should be activated.\\
+{\tt [kn]}     \> equivalent roughness length scale of the bottom (in m).\+\\
+                  Default: {\tt [kn]} = 0.05.\\
+                  Note that {\tt [kn]} is allowed to vary over the computational region; in that case\\
+                  use the commands {\tt INPGRID FRICTION} and {\tt READINP FRICTION} to define and read\\
+                  the friction data. This command {\tt FRICTION} is still required to define the type of\\
+                  friction expression. The value of {\tt [kn]} in this command is then not required\\
+                  (it will be ignored).\-\\
 \end{tabbing}
 
 \idxcmd{TRIAD}
@@ -2587,7 +2605,7 @@
 \idxcmd{OBSTACLE}
 \linecmd
 \begin{verbatim}
-          |  TRANSm [trcoef]                          |
+          | -> TRANSm [trcoef]                        |
 OBSTacle <                                            |
           |       | -> GODA [hgt] [alpha] [beta]       >                        &
           |  DAM <                                    |
@@ -2684,7 +2702,7 @@
 {\tt [reflc]}    \> constant reflection coefficient, formulated in terms of wave height, i.e. ratio\+\\
                     of reflected significant wave height over incoming significant wave height.\\
                     Restriction: $0 \leq$ {\tt [reflc]} $\leq 1$.\\
-                    Default: {\tt [reflc]}=1, if the keyword REFL is present.\\
+                    Default: {\tt [reflc]}=1, if the keyword REFL is present, otherwise {\tt [reflc]}=0.\\
                     NOTE: the program checks if the criterion $0 \leq$ {\tt [reflc]}$^2 + ${\tt [trcoef]}$^2 \leq 1$\\
                     is fulfilled.\-\\
 {\tt RSPEC}      \> indicates specular reflection which is the default. The angle of reflection\+\\
@@ -2741,9 +2759,6 @@
 \linecmd
 
 \noindent
-CANNOT BE USED IN CASE OF \underline{UNSTRUCTURED GRIDS}.
-\\[2ex]
-\noindent
 If this optional command is given, the diffraction is included in the wave computation.
 But the diffraction approximation in SWAN does not properly handle diffraction in harbours
 or in front of reflecting obstacles (see Scientific/Technical documentation). Behind breakwaters with a
@@ -2764,7 +2779,7 @@
         \end{displaymath}
         For $a = 0.2$ (recommended), the final width of the filter is $\varepsilon_x = \frac{1}{2} \sqrt{3n} \Delta x$
         (in $x-$direction and similarly in $y-$direction) and $n$ is the number of repetitions (see Scientific/Technical
-        documentation, Eq. (2.100)).
+        documentation, Eq. (2.100)). Note that this smoothing option can not be applied in case of unstructured meshes.
 \end{enumerate}
 \begin{tabbing}
 xxxxxxxxxxxx\= \kill
@@ -2837,15 +2852,16 @@
 \begin{itemize}
   \item the BSBT scheme (stationary and nonstationary) instead of the default S\&L scheme (in case
         of nonstationary cases) or the default SORDUP scheme (in case of stationary cases) or
-  \item the wave age (for the default nonstationary S\&L scheme).
+  \item the alleviation of the garden-sprinkler effect (GSE).
 \end{itemize}
 \begin{tabbing}
 xxxxxxxxxxxx\= \kill
 {\tt BSBT}      \> the BSBT scheme will be used in the computations.\\
-{\tt GSE}       \> garden-sprinkler-effect is to be counteracted in the S\&L propagation scheme\+\\
-                   (default for nonstationary computations) by adding a diffusion term to the basic\\
-                   equation. This may affect the numerical stability of SWAN (see Scientific/\\
-                   Technical documentation).\-\\
+{\tt GSE}       \> garden-sprinkler effect is to be counteracted in the S\&L propagation scheme\+\\
+                   (default for nonstationary regular grid computations) or in the propagation\\
+                   scheme for unstructured grids by adding a diffusion term to the basic equation.\\
+                   This may affect the numerical stability of SWAN (see Scientific/Technical\\
+                   documentation).\-\\
 {\tt [waveage]} \> the time interval used to determine the diffusion which counteracts the so-called\+\\
                    garden-sprinkler effect. The default value of {\tt [waveage]} is zero, i.e. no added\\
                    diffusion. The value of {\tt [waveage]} should correspond to the travel time of\\
@@ -2857,10 +2873,10 @@
         higher order schemes (S\&L and SORDUP) it is important to use a gradually varying grid, otherwise there
         may be a severe loss of accuracy. If sharp transitions in the grid cannot be avoided it is safer to use
         the BSBT scheme.
-  \item In the computation with unstructured meshes, a first order upwind scheme will be employed. This
+  \item In the computation with unstructured meshes, a lowest order upwind scheme will be employed. This
         scheme is very robust but rather diffusive. This may only be significant for the case when swell
         waves propagate over relative large distances (in the order of thousands of kilometers) within the
-        model domain.
+        model domain. However and most fortunately, in such a case, this will alleviate the garden-sprinkler effect.
 \end{itemize}
 
 \idxcmd{NUMERIC}
@@ -3658,10 +3674,21 @@
 {\tt VEL}    \> current velocity (vector; in m/s).\\
 {\tt FRCOEF} \> friction coefficient (equal to {\tt [cfw]} or {\tt [kn]} in command {\tt FRICTION}).\\
 {\tt WIND}   \> wind velocity (vector; in m/s).\\
-{\tt PROPAGAT} \> total energy propagation\+\\
-                (in W/m$^2$ or m$^2$/s, depending on command {\tt SET}).\-\\
-{\tt PROPXY} \> energy propagation in geographic space\+\\
+\end{tabbing}
+\begin{picture}(16,0.1)
+  \put(0,0){\line(1,0){15}}
+  \put(0,0.04){\line(1,0){15}}
+\end{picture}
+\\
+\noindent
+Integral over all frequencies and directions of the absolute value of each term in\\
+the energy balance equation:
+\begin{tabbing}
+ xxxxxxxxxxxx\= \kill
+{\tt PROPAGAT} \> sum of {\tt PROPXY}, {\tt PROPTHETA} and {\tt PROPSIGMA}\+\\
                 (in W/m$^2$ or m$^2$/s, depending on command {\tt SET}).\-\\
+{\tt PROPXY} \> energy propagation in geographic space; sum of $x-$ and $y-$direction\+\\
+                terms (in W/m$^2$ or m$^2$/s, depending on command {\tt SET}).\-\\
 {\tt PROPTHETA} \> energy propagation in theta space\+\\
                 (in W/m$^2$ or m$^2$/s, depending on command {\tt SET}).\-\\
 {\tt PROPSIGMA} \> energy propagation in sigma space\+\\
@@ -3684,8 +3711,15 @@
                 (in W/m$^2$ or m$^2$/s, depending on command {\tt SET}).\-\\
 {\tt DISWCAP} \> energy dissipation due to whitecapping\+\\
                 (in W/m$^2$ or m$^2$/s, depending on command {\tt SET}).\-\\
-{\tt RADSTR} \> radiation stress\+\\
+{\tt RADSTR} \> energy transfer between waves and currents due to radiation stress\+\\
                 (in W/m$^2$ or m$^2$/s, depending on command {\tt SET}).\-\\
+\end{tabbing}
+\begin{picture}(16,0.1)
+  \put(0,0){\line(1,0){15}}
+  \put(0,0.04){\line(1,0){15}}
+\end{picture}
+\begin{tabbing}
+ xxxxxxxxxxxx\= \kill
 {\tt QB}     \> fraction of breaking waves due to depth-induced breaking.\\
 {\tt TRANSP} \> transport of energy (vector; in W/m or m$^3$/s, depending on command {\tt SET}).\\
 {\tt FORCE}  \> wave-induced force per unit surface area (vector; in N/m$^2$).\\
@@ -4209,7 +4243,7 @@
                         (Kuik~{\it et~al}. (1988); this is the standard definition for WAVEC buoys\\
                         integrated over all frequencies):\\
                         \\
-$({\rm DSPR} \frac{\pi}{180})^2 = 2\left( 1 - \sqrt{\left[ \left( \frac{\int\sin\theta E(\sigma,\theta)d\sigma d\theta}{\int E(\sigma)d\sigma} \right)^2 + \left( \frac{\int \cos\theta E(\sigma,\theta)d\sigma d\theta}{\int E(\sigma)d\sigma} \right)^2 \right]}\right)$ \\
+$({\rm DSPR} \frac{\pi}{180})^2 = 2\left( 1 - \sqrt{\left[ \left( \frac{\int\sin\theta E(\sigma,\theta)d\sigma d\theta}{\int E(\sigma,\theta)d\sigma d\theta} \right)^2 + \left( \frac{\int \cos\theta E(\sigma,\theta)d\sigma d\theta}{\int E(\sigma,\theta)d\sigma d\theta} \right)^2 \right]}\right)$ \\
                         \-\\
 QP                   \> The peakedness of the wave spectrum, defined as\+\\
                         \\
@@ -4278,7 +4312,7 @@
                         whitecapping and depth-induced
                         surf breaking (in W/${\rm m}^2$ or ${\rm m}^2$/s,\\
                         depending on the command {\tt SET}).\-\\
-RADSTR               \> Radiation stress per unit time, defined as\+\\
+RADSTR               \> Work done by the radiation stress per unit time, defined as\+\\
                         \\
 $\int\limits_{0}\limits^{2\pi} \int\limits_{\sigma_{\mbox{\tiny low}}}\limits^{\sigma_{\mbox{\tiny high}}} |S_{\rm tot} - \frac{\partial E}{\partial t} - \nabla_{\vec{x}} \cdot [({\vec{c}}_g + \vec{U}) E] - \nabla_{(\sigma,\theta)} \cdot ({\vec{c}}_{(\sigma,\theta)} E)| d\sigma d\theta$ \\
                         \\
@@ -4646,7 +4680,7 @@
 !         | WESTHuysen                  |
 !
 !
-!         | ->KOMen    [cds2] [stpm] [powst] [delta] [powk]  |
+!         | -> KOMen   [cds2] [stpm] [powst] [delta] [powk]  |
 !         |                                                  |
 !         |   JANSsen  [cds1]  [delta] [pwtail]              |
 !         |                                                  |
@@ -4666,15 +4700,29 @@
 !   MDIA LAMbda < [lambda] > <                                 >
 !                             | CNL4_12 < [cnl4_1] [cnl4_2] > |
 !
-!   BREaking CONSTANT [alpha] [gamma]
+!        | -> CON [alpha] [gamma]                                      |
+!        |                                                             |
+!        |    VAR [alpha] [gammin] [gammax] [gamneg] [coeff1] [coeff2] |
+!        |                                                             |
+!   BRE <     RUE [alpha] [a] [b]                                       >
+!        |                                                             |
+!        |    TG  [alpha] [gamma] [pown]                               |
+!        |                                                             |
+!        !    BIP [alpha] [pown] [bref]                                |
+!
+!              |             | -> CONstant [cfjon]
+!              | -> JONswap <
+!              |             |    VARiable [cfj1] [cfj2] [dsp1] [dsp2]
+!   FRICtion  <
+!              |    COLLins [cfw] [cfc]
+!              |
+!              |    MADsen  [kn]
 !
-!   FRICTION  JONSWAP [cfjon] | COLLINS [cfw] | MADSEN [kn]
-!
-!   TRIAD   [trfac] [cutfr] [urcrit] [urslim]
+!   TRIad   [trfac] [cutfr] [urcrit] [urslim]
 !
 !   LIMiter [ursell] [qb]
 !
-!               |  TRANSm [trcoef]                          |
+!               | -> TRANSm [trcoef]                        |
 !   OBSTacle   <                                            |
 !               |       | -> GODA [hgt] [alpha] [beta]       >              &
 !               |  DAM <                                    |
@@ -4801,7 +4849,7 @@
 
 \begin{verbatim}
 SWAN   1                                Swan standard spectral file, version
-$ Data produced by SWAN version 40.72A
+$ Data produced by SWAN version 40.72AB
 $ Project:'projname'     ;   run number:'runnum'
 TIME                                    time-dependent data
      1                                  time coding option
@@ -4907,7 +4955,7 @@
 
 \begin{verbatim}
 SWAN   1                                Swan standard spectral file, version
-$ Data produced by SWAN version 40.72A
+$ Data produced by SWAN version 40.72AB
 $ Project:'projname'     ;   run number:'runnum'
 LOCATIONS                               locations in x-y-space
      1                                  number of locations
@@ -5087,18 +5135,15 @@
 
 \bibitem{Implman}
 {SWAN -- Implementation manual}. Delft University of Technology, Environmental Fluid Mechanics Section, available from
-\hl{http://www.fluidmechanics.tudelft.nl/swan/index.htm} (Version 40.72, May 2008).
+\hl{http://www.swan.tudelft.nl} (Version 40.72AB, March 2009).
 
 \bibitem{Progrul}
 {SWAN -- Programming rules}. Delft University of Technology, Environmental Fluid Mechanics Section, available from
-\hl{http://www.fluidmechanics.tudelft.nl/swan/index.htm} (Version 1.3, August 2006).
-
-\bibitem{Sysdoc}
-{SWAN -- System documentation}. Delft University of Technology, Environmental Fluid Mechanics Section, to be available.
+\hl{http://www.swan.tudelft.nl} (Version 1.3, August 2006).
 
 \bibitem{Techdoc}
 {SWAN -- Scientific and Technical documentation}. Delft University of Technology, Environmental Fluid Mechanics Section, available from
-\hl{http://www.fluidmechanics.tudelft.nl/swan/index.htm} (Version 40.72A, October 2008).
+\hl{http://www.swan.tudelft.nl} (Version 40.72AB, March 2009).
 
 \end{thebibliography}
 
@@ -5110,27 +5155,28 @@
 
   \item bathymetry, \hyperpage{5}, \hyperpage{10, 11}, \hyperpage{64}, 
 		\hyperpage{85}
-  \item BLOCK, \hyperpage{74}
+  \item BLOCK, \hyperpage{75}
   \item bottom, \hyperpage{3}, \hyperpage{8--13}, \hyperpage{16}, 
 		\hyperpage{18, 19}, \hyperpage{21, 22}, \hyperpage{24}, 
 		\hyperpage{34, 35}, \hyperpage{38--40}, \hyperpage{53}, 
-		\hyperpage{55, 56}, \hyperpage{58}, \hyperpage{66, 67}, 
-		\hyperpage{69}, \hyperpage{79, 80}, \hyperpage{85, 86}, 
-		\hyperpage{92}, \hyperpage{94}
+		\hyperpage{55, 56}, \hyperpage{58, 59}, 
+		\hyperpage{66, 67}, \hyperpage{69}, \hyperpage{79, 80}, 
+		\hyperpage{86}, \hyperpage{92}, \hyperpage{94}
   \item BOUND SHAPE, \hyperpage{41}
   \item boundary, \hyperpage{1}, \hyperpage{3}, \hyperpage{6--14}, 
 		\hyperpage{22}, \hyperpage{26}, \hyperpage{28}, 
 		\hyperpage{30}, \hyperpage{32}, \hyperpage{42--51}, 
-		\hyperpage{53}, \hyperpage{61}, \hyperpage{66}, 
-		\hyperpage{70--72}, \hyperpage{83, 84}, \hyperpage{87}
+		\hyperpage{53}, \hyperpage{62}, \hyperpage{66}, 
+		\hyperpage{70}, \hyperpage{72}, \hyperpage{83, 84}, 
+		\hyperpage{87}
   \item BOUNDNEST1, \hyperpage{46}
   \item BOUNDNEST2, \hyperpage{47}
   \item BOUNDNEST3, \hyperpage{49}
   \item BOUNDSPEC, \hyperpage{42}
-  \item BREAKING, \hyperpage{55}
+  \item BREAKING, \hyperpage{54}
   \item breaking, \hyperpage{11}, \hyperpage{16--18}, \hyperpage{22}, 
-		\hyperpage{53}, \hyperpage{55--57}, \hyperpage{61}, 
-		\hyperpage{79}, \hyperpage{85, 86}, \hyperpage{92}
+		\hyperpage{53--55}, \hyperpage{57}, \hyperpage{61}, 
+		\hyperpage{80}, \hyperpage{86}, \hyperpage{92}
 
   \indexspace
 
@@ -5138,13 +5184,13 @@
 		\hyperpage{26, 27}, \hyperpage{29}, \hyperpage{36}, 
 		\hyperpage{41}, \hyperpage{44}, \hyperpage{46}, 
 		\hyperpage{48--51}, \hyperpage{59}, \hyperpage{67--71}, 
-		\hyperpage{73, 74}, \hyperpage{78}, \hyperpage{85}, 
+		\hyperpage{74}, \hyperpage{78}, \hyperpage{86}, 
 		\hyperpage{90}, \hyperpage{94}, \hyperpage{106}, 
 		\hyperpage{108, 109}
   \item CGRID, \hyperpage{28}
   \item co-ordinate, \hyperpage{32}
   \item coastal, \hyperpage{3--5}, \hyperpage{15}
-  \item COMPUTE, \hyperpage{86}
+  \item COMPUTE, \hyperpage{87}
   \item convergence, \hyperpage{6, 7}, \hyperpage{63}
   \item COORDINATES, \hyperpage{27}
   \item Courant, \hyperpage{14}
@@ -5152,15 +5198,15 @@
 		\hyperpage{14}, \hyperpage{16}, \hyperpage{19--21}, 
 		\hyperpage{26}, \hyperpage{32}, \hyperpage{34, 35}, 
 		\hyperpage{38}, \hyperpage{45}, \hyperpage{52}, 
-		\hyperpage{64, 65}, \hyperpage{67}, \hyperpage{79}, 
+		\hyperpage{64, 65}, \hyperpage{67}, \hyperpage{79, 80}, 
 		\hyperpage{83}, \hyperpage{90}
   \item CURVE, \hyperpage{68}
   \item curvi-linear, \hyperpage{3}, \hyperpage{6}, \hyperpage{8, 9}, 
 		\hyperpage{12--14}, \hyperpage{16}, \hyperpage{21, 22}, 
 		\hyperpage{28, 29}, \hyperpage{32}, \hyperpage{34, 35}, 
 		\hyperpage{38, 39}, \hyperpage{43}, \hyperpage{47, 48}, 
-		\hyperpage{50}, \hyperpage{62}, \hyperpage{65}, 
-		\hyperpage{67}
+		\hyperpage{50}, \hyperpage{62}, \hyperpage{66}, 
+		\hyperpage{68}
 
   \indexspace
 
@@ -5170,7 +5216,7 @@
 		\hyperpage{64}
   \item diffusion, \hyperpage{62}, \hyperpage{64, 65}
   \item dissipation, \hyperpage{22}, \hyperpage{52, 53}, \hyperpage{55}, 
-		\hyperpage{61}, \hyperpage{79}, \hyperpage{92}
+		\hyperpage{61}, \hyperpage{80}, \hyperpage{92}
 
   \indexspace
 
@@ -5182,15 +5228,15 @@
   \item frequency, \hyperpage{4}, \hyperpage{6}, \hyperpage{8}, 
 		\hyperpage{14, 15}, \hyperpage{26}, \hyperpage{30, 31}, 
 		\hyperpage{42}, \hyperpage{45}, \hyperpage{51, 52}, 
-		\hyperpage{54}, \hyperpage{56}, \hyperpage{61}, 
-		\hyperpage{64--66}, \hyperpage{72--74}, \hyperpage{78}, 
-		\hyperpage{83}, \hyperpage{85, 86}, \hyperpage{89}, 
-		\hyperpage{91}, \hyperpage{108, 109}
+		\hyperpage{54--56}, \hyperpage{61}, \hyperpage{64--66}, 
+		\hyperpage{72--74}, \hyperpage{79}, \hyperpage{83}, 
+		\hyperpage{86}, \hyperpage{89}, \hyperpage{91}, 
+		\hyperpage{108, 109}
   \item FRICTION, \hyperpage{55}
   \item friction, \hyperpage{9--12}, \hyperpage{16--19}, 
 		\hyperpage{21, 22}, \hyperpage{24}, \hyperpage{34, 35}, 
 		\hyperpage{38}, \hyperpage{53}, \hyperpage{55, 56}, 
-		\hyperpage{79, 80}, \hyperpage{85, 86}, \hyperpage{92}
+		\hyperpage{79, 80}, \hyperpage{86}, \hyperpage{92}
   \item Froude, \hyperpage{26}, \hyperpage{38}
 
   \indexspace
@@ -5204,14 +5250,14 @@
   \indexspace
 
   \item harbour, \hyperpage{5}, \hyperpage{60}
-  \item HOTFILE, \hyperpage{87}
+  \item HOTFILE, \hyperpage{88}
 
   \indexspace
 
   \item INITIAL, \hyperpage{50}
   \item initial, \hyperpage{8}, \hyperpage{10, 11}, \hyperpage{14}, 
 		\hyperpage{22}, \hyperpage{41}, \hyperpage{51}, 
-		\hyperpage{87}
+		\hyperpage{87, 88}
   \item INPGRID, \hyperpage{33}
   \item island, \hyperpage{4}, \hyperpage{6}, \hyperpage{11}
   \item ISOLINE, \hyperpage{69}
@@ -5222,13 +5268,13 @@
 
   \indexspace
 
-  \item latitude, \hyperpage{8}, \hyperpage{27}, \hyperpage{85}, 
+  \item latitude, \hyperpage{8}, \hyperpage{27}, \hyperpage{86}, 
 		\hyperpage{108}
   \item LIMITER, \hyperpage{56}
-  \item limiter, \hyperpage{5}, \hyperpage{7}, \hyperpage{56}, 
-		\hyperpage{64}
+  \item limiter, \hyperpage{5}, \hyperpage{7}, \hyperpage{57}, 
+		\hyperpage{64, 65}
   \item longitude, \hyperpage{8}, \hyperpage{27}, \hyperpage{49, 50}, 
-		\hyperpage{85}, \hyperpage{108}
+		\hyperpage{86}, \hyperpage{108}
 
   \indexspace
 
@@ -5237,9 +5283,9 @@
   \indexspace
 
   \item nautical, \hyperpage{7, 8}, \hyperpage{108, 109}
-  \item NESTOUT, \hyperpage{83}
+  \item NESTOUT, \hyperpage{84}
   \item NGRID, \hyperpage{70}
-  \item NUMERIC, \hyperpage{62}
+  \item NUMERIC, \hyperpage{63}
 
   \indexspace
 
@@ -5252,35 +5298,35 @@
 
   \indexspace
 
-  \item POINTS, \hyperpage{69}
+  \item POINTS, \hyperpage{70}
   \item PROJECT, \hyperpage{24}
-  \item PROP, \hyperpage{61}
+  \item PROP, \hyperpage{62}
   \item propagation, \hyperpage{6}, \hyperpage{12}, \hyperpage{22}, 
 		\hyperpage{57}, \hyperpage{61, 62}, \hyperpage{79}, 
 		\hyperpage{92}
 
   \indexspace
 
-  \item QUADRUPL, \hyperpage{54}
+  \item QUADRUPL, \hyperpage{53}
   \item quadruplets, \hyperpage{15, 16}, \hyperpage{22}, \hyperpage{54}, 
-		\hyperpage{56}, \hyperpage{79}, \hyperpage{92}
+		\hyperpage{57}, \hyperpage{79}, \hyperpage{92}
   \item QUANTITY, \hyperpage{72}
 
   \indexspace
 
-  \item RAY, \hyperpage{68}
+  \item RAY, \hyperpage{69}
   \item READGRID COORDINATES, \hyperpage{32}
   \item READGRID UNSTRUCTURED, \hyperpage{32}
   \item READINP, \hyperpage{37}
   \item recti-linear, \hyperpage{3}, \hyperpage{16}, \hyperpage{28, 29}, 
-		\hyperpage{38}, \hyperpage{67}
+		\hyperpage{38}, \hyperpage{68}
   \item reflection, \hyperpage{57--59}
   \item refraction, \hyperpage{6}, \hyperpage{11}, \hyperpage{61}, 
-		\hyperpage{64}, \hyperpage{69}
+		\hyperpage{64, 65}, \hyperpage{69}
   \item regular, \hyperpage{4}, \hyperpage{8, 9}, \hyperpage{12}, 
 		\hyperpage{16}, \hyperpage{22}, \hyperpage{28, 29}, 
-		\hyperpage{31, 32}, \hyperpage{34, 35}, 
-		\hyperpage{65, 66}
+		\hyperpage{31, 32}, \hyperpage{34, 35}, \hyperpage{62}, 
+		\hyperpage{66, 67}
 
   \indexspace
 
@@ -5291,23 +5337,22 @@
   \item SETUP, \hyperpage{59}
   \item shoaling, \hyperpage{18}
   \item SORDUP, \hyperpage{62}
-  \item SPECOUT, \hyperpage{82}
+  \item SPECOUT, \hyperpage{83}
   \item specular, \hyperpage{59}
   \item spherical, \hyperpage{3--5}, \hyperpage{8, 9}, \hyperpage{21}, 
 		\hyperpage{25}, \hyperpage{27--29}, \hyperpage{31}, 
 		\hyperpage{36}, \hyperpage{45}, \hyperpage{48--50}, 
-		\hyperpage{59}, \hyperpage{67--71}, \hyperpage{80}, 
-		\hyperpage{85}, \hyperpage{103}, \hyperpage{108}
+		\hyperpage{59, 60}, \hyperpage{67--71}, \hyperpage{81}, 
+		\hyperpage{86}, \hyperpage{103}, \hyperpage{108}
   \item stability, \hyperpage{14}, \hyperpage{62}
   \item stationary, \hyperpage{4}, \hyperpage{6}, \hyperpage{8}, 
 		\hyperpage{14}, \hyperpage{16}, \hyperpage{18}, 
 		\hyperpage{20, 21}, \hyperpage{26, 27}, \hyperpage{34}, 
 		\hyperpage{37--39}, \hyperpage{45, 46}, \hyperpage{51}, 
 		\hyperpage{53}, \hyperpage{62--64}, \hyperpage{80}, 
-		\hyperpage{86, 87}, \hyperpage{103}, \hyperpage{106}, 
+		\hyperpage{87, 88}, \hyperpage{103}, \hyperpage{106}, 
 		\hyperpage{108, 109}
-  \item steepness, \hyperpage{22}, \hyperpage{53}, \hyperpage{80}, 
-		\hyperpage{93}
+  \item steepness, \hyperpage{53}, \hyperpage{80}, \hyperpage{93}
   \item STOP, \hyperpage{88}
   \item swell, \hyperpage{12}, \hyperpage{15}, \hyperpage{55}, 
 		\hyperpage{62}, \hyperpage{73, 74}, \hyperpage{78}, 
@@ -5316,38 +5361,38 @@
   \indexspace
 
   \item TABLE, \hyperpage{81}
-  \item TEST, \hyperpage{84}
+  \item TEST, \hyperpage{85}
   \item TRIAD, \hyperpage{56}
-  \item triads, \hyperpage{16}, \hyperpage{56}, \hyperpage{79}, 
+  \item triads, \hyperpage{16}, \hyperpage{56}, \hyperpage{80}, 
 		\hyperpage{92}
-  \item triangular, \hyperpage{3}, \hyperpage{16}, \hyperpage{71}
+  \item triangular, \hyperpage{3}, \hyperpage{16}, \hyperpage{72}
 
   \indexspace
 
   \item unstructured, \hyperpage{3, 4}, \hyperpage{8, 9}, 
 		\hyperpage{12, 13}, \hyperpage{21}, \hyperpage{28, 29}, 
 		\hyperpage{32}, \hyperpage{34, 35}, \hyperpage{44}, 
-		\hyperpage{58}, \hyperpage{62--64}, \hyperpage{66}, 
-		\hyperpage{71}, \hyperpage{78}
+		\hyperpage{58}, \hyperpage{60}, \hyperpage{62}, 
+		\hyperpage{64}, \hyperpage{66}, \hyperpage{71}, 
+		\hyperpage{78}
 
   \indexspace
 
   \item WAM, \hyperpage{4--7}, \hyperpage{9, 10}, \hyperpage{15}, 
 		\hyperpage{18}, \hyperpage{22}, \hyperpage{37}, 
 		\hyperpage{46}, \hyperpage{48--50}, \hyperpage{54}, 
-		\hyperpage{81--84}, \hyperpage{87}
+		\hyperpage{81, 82}, \hyperpage{84}, \hyperpage{87}
   \item WAVEWATCH, \hyperpage{4--7}, \hyperpage{9, 10}, \hyperpage{15}, 
 		\hyperpage{22}, \hyperpage{46}, \hyperpage{49, 50}
-  \item WCAPPING, \hyperpage{53}
-  \item whitecapping, \hyperpage{7}, \hyperpage{16}, \hyperpage{22}, 
-		\hyperpage{53, 54}, \hyperpage{61}, \hyperpage{79}, 
-		\hyperpage{85, 86}, \hyperpage{92}
+  \item whitecapping, \hyperpage{7}, \hyperpage{16}, \hyperpage{53}, 
+		\hyperpage{61}, \hyperpage{80}, \hyperpage{86}, 
+		\hyperpage{92}
   \item WIND, \hyperpage{41}
   \item wind, \hyperpage{3}, \hyperpage{5}, \hyperpage{7--12}, 
 		\hyperpage{14--19}, \hyperpage{21}, \hyperpage{26}, 
 		\hyperpage{34, 35}, \hyperpage{39}, \hyperpage{41}, 
 		\hyperpage{51--53}, \hyperpage{55}, \hyperpage{61, 62}, 
-		\hyperpage{64, 65}, \hyperpage{79}, \hyperpage{85, 86}, 
+		\hyperpage{64, 65}, \hyperpage{79}, \hyperpage{86, 87}, 
 		\hyperpage{92}, \hyperpage{94}
 
 \end{theindex}
--- SwanVertlist.ftn90	2008-05-25 14:59:48.000000000 +0200
+++ SwanVertlist.ftn90	2009-03-03 17:50:34.000000000 +0100
@@ -11,7 +11,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- swmod1.ftn	2009-03-03 17:49:37.000000000 +0100
+++ swmod1.ftn	2009-03-03 17:50:34.000000000 +0100
@@ -26,7 +26,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -184,7 +184,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -306,7 +306,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -409,7 +409,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -552,7 +552,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1359,7 +1359,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1557,7 +1557,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -1805,6 +1805,7 @@
 !
 ! ALCP   [CALCUL] =-ALPC;
 !                  direction of user coordinates w.r.t. computational coordinates
+! ALOCMP [ FALSE]  if True, array COMPDA must be re-allocated
 ! ALPC   [    0.]  direction of x-axis of computational grid w.r.t. user coordinates
 ! COSPC  [CALCUL] =COS(ALPC)
 ! DX     [CALCUL] =XCLEN/MXS; mesh size in x-direction of computational grid
@@ -1868,6 +1869,7 @@
       REAL                YPC
       REAL                XCGMIN,      XCGMAX,      YCGMIN,      YCGMAX
       LOGICAL             FULCIR
+      LOGICAL ::          ALOCMP = .FALSE.                                40.97
 !$OMP THREADPRIVATE(IXCGRD,IYCGRD,KCGRD)
 !
 !     *** physical parameters ***
@@ -1975,6 +1977,9 @@
 !                 =0; for command OFF BRE, surf breaking is inactive
 !                 =1; for command BRE CON ..., surf breaking with constant parameter
 !                 =2; for command BRE VAR ..., surf breaking according to Nelson
+!                 =3; for command BRE RUE ..., surf breaking according to Ruessink
+!                 =4; for command BRE TG ... , surf breaking according to Thornton and Guza
+!                 =5; for command BRE BIP ..., surf breaking according to biphase scaling
 ! ITERMX [      ] maximum number of iterations:
 !                 is set equal to MXITST in case of stationary computations
 !                 is set equal to MXITNS in case of nonstationary computations
@@ -2370,7 +2375,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -2477,7 +2482,7 @@
 !                 1=first order (BSBT), 3=3rd order (S&L)
 ! PROPSL [CALCUL] indicates which numerical scheme is used locally
 ! REARTH [   6E6] radius of the earth
-! WAVAGE [    0.] indicates "wave age" parameter (used in countering Garden Sprinkler Effect
+! WAVAGE [    0.] indicates "wave age" parameter (used in counteracting garden-sprinkler effect
 !                 in subroutine SANDL)
 !
       INTEGER PROPSC,    KSPHER,    KREPTX
@@ -2526,7 +2531,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
--- swmod2.ftn	2008-05-25 14:59:50.000000000 +0200
+++ swmod2.ftn	2009-03-03 17:50:34.000000000 +0100
@@ -25,7 +25,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -146,7 +146,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -287,7 +287,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -400,7 +400,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -545,7 +545,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -657,7 +657,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -788,7 +788,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
@@ -967,7 +967,7 @@
 !
 !
 !     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2008  Delft University of Technology
+!     Copyright (C) 2009  Delft University of Technology
 !
 !     This program is free software; you can redistribute it and/or
 !     modify it under the terms of the GNU General Public License as
