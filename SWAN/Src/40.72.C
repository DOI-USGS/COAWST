1) Unstructured mesh computation is made more efficient, i.e. only a few
   sweeps per iteration or time step is needed.
2) Original fort.14 file can also be dealt with in SWAN.
   Boundary markers will be derived from the ADCIRC boundary information.
3) Block outputting for unstructured mesh cases is considerably optimized.
4) The default advanced stopping criterion (NUM STOPC) is based on the
   curvature of Hm0 only for reasons of robustness. However, the curvature
   of Tm01 can be included as an option.
5) Technical documentation is extended with useful information.
6) Bug fixes:
   *) never-ending sweep in unstructured mesh cases is prevented,
   *) assign reference point to deepest point in case of no boundary
      condition in unstructured mesh cases,
   *) 2 small corrections in collecting data for MPI parallel runs.

--- celldisc2.eps	1970-01-01 01:00:00.000000000 +0100
+++ celldisc2.eps	2009-08-05 21:33:07.000000000 +0200
@@ -0,0 +1,597 @@
+%!PS-Adobe-2.0 EPSF-2.0
+%%Title: celldisc2.fig
+%%Creator: fig2dev Version 3.2 Patchlevel 4
+%%CreationDate: Thu Jun 25 17:06:44 2009
+%%For: mzijlema@TUD11321 (Marcel Zijlema)
+%%BoundingBox: 0 0 503 324
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
+newpath 0 324 moveto 0 0 lineto 503 0 lineto 503 324 lineto closepath clip newpath
+-139.8 590.6 translate
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
+n 3675 7492 81 81 0 360 DrawEllipse gs /PC [[0.00 0.00 0.00] [0.00 0.00 0.00]] def
+15.00 15.00 sc P1 [16 0 0 -8 239.60 494.07] PATmp PATsp ef gr PATusp gs col0 s gr
+
+% Ellipse
+n 4965 4635 81 81 0 360 DrawEllipse gs /PC [[0.00 0.00 0.00] [0.00 0.00 0.00]] def
+15.00 15.00 sc P1 [16 0 0 -8 325.60 303.60] PATmp PATsp ef gr PATusp gs col0 s gr
+
+% Ellipse
+n 7350 6157 81 81 0 360 DrawEllipse gs /PC [[0.00 0.00 0.00] [0.00 0.00 0.00]] def
+15.00 15.00 sc P1 [16 0 0 -8 484.60 405.07] PATmp PATsp ef gr PATusp gs col0 s gr
+
+% Polyline
+gs  clippath
+2430 8220 m 2370 8220 l 2370 8371 l 2400 8251 l 2430 8371 l cp
+eoclip
+n 2400 9270 m
+ 2400 8235 l gs col0 s gr gr
+
+% arrowhead
+n 2430 8371 m 2400 8251 l 2370 8371 l 2430 8371 l  cp gs 0.00 setgray ef gr  col0 s
+% Polyline
+gs  clippath
+3450 9300 m 3450 9240 l 3299 9240 l 3419 9270 l 3299 9300 l cp
+eoclip
+n 2400 9270 m
+ 3435 9270 l gs col0 s gr gr
+
+% arrowhead
+n 3299 9300 m 3419 9270 l 3299 9240 l 3299 9300 l  cp gs 0.00 setgray ef gr  col0 s
+% Polyline
+n 7365 6123 m 4980 4638 l 3675 7518 l
+ 7365 6168 l gs col0 s gr 
+% Polyline
+30.000 slw
+gs  clippath
+8611 7009 m 8706 6856 l 8330 6624 l 8590 6890 l 8236 6777 l cp
+eoclip
+n 7431 6175 m
+ 8646 6925 l gs col0 s gr gr
+
+% arrowhead
+n 8236 6777 m 8590 6890 l 8330 6624 l  col0 s
+% Polyline
+gs  clippath
+9009 5696 m 8950 5526 l 8533 5671 l 8903 5638 l 8592 5841 l cp
+eoclip
+n 7421 6156 m
+ 8966 5616 l gs col0 s gr gr
+
+% arrowhead
+n 8592 5841 m 8903 5638 l 8533 5671 l  col0 s
+% Polyline
+7.500 slw
+gs  clippath
+9386 7449 m 9418 7398 l 9290 7317 l 9376 7407 l 9258 7367 l cp
+eoclip
+n 8565 6891 m
+ 9390 7416 l gs col0 s gr gr
+
+% arrowhead
+n 9258 7367 m 9376 7407 l 9290 7317 l  col0 s
+% Polyline
+gs  clippath
+9980 5276 m 9959 5219 l 9817 5272 l 9940 5259 l 9837 5328 l cp
+eoclip
+n 8891 5643 m
+ 9956 5253 l gs col0 s gr gr
+
+% arrowhead
+n 9837 5328 m 9940 5259 l 9817 5272 l  col0 s
+% Polyline
+15.000 slw
+ [90] 0 sd
+gs  clippath
+7753 7447 m 7866 7408 l 7773 7136 l 7794 7383 l 7659 7175 l cp
+eoclip
+n 7385 6184 m
+ 7805 7414 l gs col0 s gr gr
+ [] 0 sd
+% arrowhead
+n 7659 7175 m 7794 7383 l 7773 7136 l  col0 s
+% Polyline
+ [90] 0 sd
+gs  clippath
+8197 5073 m 8100 5003 l 7932 5237 l 8121 5078 l 8029 5307 l cp
+eoclip
+n 7375 6116 m
+ 8140 5051 l gs col0 s gr gr
+ [] 0 sd
+% arrowhead
+n 8029 5307 m 8121 5078 l 7932 5237 l  col0 s
+% Polyline
+7.500 slw
+n 7204 6219 m 7266 6419 l
+ 7426 6366 l gs col0 s gr 
+% Polyline
+n 7226 6029 m 7328 5864 l
+ 7466 5958 l gs col0 s gr 
+/Times-Roman ff 375.00 scf sf
+2220 8265 m
+gs 1 -1 sc (y) col0 sh gr
+/Times-Roman ff 375.00 scf sf
+3480 9375 m
+gs 1 -1 sc (x) col0 sh gr
+/Times-Roman ff 375.00 scf sf
+3375 7920 m
+gs 1 -1 sc (3) col0 sh gr
+/Times-Roman ff 375.00 scf sf
+4695 4485 m
+gs 1 -1 sc (2) col0 sh gr
+/Symbol ff 375.00 scf sf
+9285 7785 m
+gs 1 -1 sc (x) col0 sh gr
+/Symbol ff 375.00 scf sf
+9960 5175 m
+gs 1 -1 sc (h) col0 sh gr
+/Times-Bold ff 450.00 scf sf
+7890 6975 m
+gs 1 -1 sc (e) col0 sh gr
+/Times-Bold ff 450.00 scf sf
+8850 5985 m
+gs 1 -1 sc (e) col0 sh gr
+/Times-Bold ff 225.00 scf sf
+9015 6075 m
+gs 1 -1 sc (\(2\)) col0 sh gr
+/Times-Bold ff 225.00 scf sf
+7755 7590 m
+gs 1 -1 sc (\(1\)) col0 sh gr
+/Times-Bold ff 450.00 scf sf
+7830 5010 m
+gs 1 -1 sc (e) col0 sh gr
+/Times-Bold ff 225.00 scf sf
+7995 4800 m
+gs 1 -1 sc (\(2\)) col0 sh gr
+/Times-Roman ff 375.00 scf sf
+6885 6240 m
+gs 1 -1 sc (1) col0 sh gr
+/Times-Bold ff 450.00 scf sf
+7590 7785 m
+gs 1 -1 sc (e) col0 sh gr
+/Times-Bold ff 225.00 scf sf
+8055 7110 m
+gs 1 -1 sc (\(1\)) col0 sh gr
+% here ends figure;
+$F2psEnd
+rs
+end
+showpage
--- fsweep.eps	1970-01-01 01:00:00.000000000 +0100
+++ fsweep.eps	2009-08-05 21:33:07.000000000 +0200
@@ -0,0 +1,214 @@
+%!PS-Adobe-2.0 EPSF-2.0
+%%Title: fsweep.fig
+%%Creator: fig2dev Version 3.2 Patchlevel 4
+%%CreationDate: Mon Jul 20 16:32:41 2009
+%%For: mzijlema@TUD11321 (Marcel Zijlema)
+%%BoundingBox: 0 0 372 508
+%%Magnification: 1.0000
+%%EndComments
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
+newpath 0 508 moveto 0 0 lineto 372 0 lineto 372 508 lineto closepath clip newpath
+-100.2 584.7 translate
+1 -1 scale
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
+n 2847 5083 38 38 0 360 DrawEllipse gs col7 0.00 shd ef gr gs col0 s gr
+
+% Ellipse
+n 3510 5086 38 38 0 360 DrawEllipse gs col7 0.00 shd ef gr gs col0 s gr
+
+% Ellipse
+n 3510 5722 38 38 0 360 DrawEllipse gs col7 0.00 shd ef gr gs col0 s gr
+
+% Ellipse
+n 4458 8304 592 592 0 360 DrawEllipse gs col0 s gr
+
+% Polyline
+gs  clippath
+1980 1215 m 1920 1215 l 1920 1367 l 1950 1247 l 1980 1367 l cp
+eoclip
+n 1950 6750 m 1950 1275 l
+ 1950 1230 l gs col0 s gr gr
+
+% arrowhead
+n 1980 1367 m 1950 1247 l 1920 1367 l  col0 s
+% Polyline
+gs  clippath
+7485 6795 m 7485 6735 l 7333 6735 l 7453 6765 l 7333 6795 l cp
+eoclip
+n 1950 6765 m 7425 6765 l
+ 7470 6765 l gs col0 s gr gr
+
+% arrowhead
+n 7333 6795 m 7453 6765 l 7333 6735 l  col0 s
+% Polyline
+n 1950 2430 m 6930 2430 l
+ 6930 6765 l gs col0 s gr 
+% Polyline
+n 5865 2895 m 5865 3510 l
+ 5865 3525 l gs col0 s gr 
+% Polyline
+n 6045 2895 m 6045 3510 l
+ 6045 3525 l gs col0 s gr 
+% Polyline
+n 6225 2895 m 6225 3510 l
+ 6225 3525 l gs col0 s gr 
+% Polyline
+n 6360 3075 m 5745 3075 l
+ 5730 3075 l gs col0 s gr 
+% Polyline
+n 6360 3255 m 5745 3255 l
+ 5730 3255 l gs col0 s gr 
+% Polyline
+n 6375 3435 m 5760 3435 l
+ 5745 3435 l gs col0 s gr 
+% Polyline
+15.000 slw
+n 2850 5085 m
+ 3870 5085 l gs col0 s gr 
+% Polyline
+n 3510 4710 m
+ 3510 5730 l gs col0 s gr 
+% Polyline
+7.500 slw
+n 3821 5074 m 3812 4983 l 3776 4897 l 3700 4829 l 3626 4778 l 3516 4752 l
+ 3516 5074 l 3815 5074 l
+ cp gs col7 0.50 shd ef gr gs col0 s gr 
+% Polyline
+30.000 slw
+gs  clippath
+3206 5553 m 3078 5425 l 2767 5736 l 3086 5546 l 2895 5864 l cp
+eoclip
+n 2162 6470 m
+ 3132 5500 l gs col7 0.00 shd ef gr gs col0 s gr gr
+
+% arrowhead
+n 2895 5864 m 3086 5546 l 2767 5736 l 2895 5864 l  cp gs 0.00 setgray ef gr  col0 s
+% Polyline
+7.500 slw
+n 4470 7320 m
+ 4470 9270 l gs col0 s gr 
+% Polyline
+n 5430 8310 m
+ 3480 8310 l gs col0 s gr 
+% Polyline
+n 5051 8307 m 5033 8153 l 4986 8032 l 4882 7890 l 4758 7789 l 4640 7739 l
+ 4533 7715 l 4471 7712 l 4471 8307 l 5051 8310 l
+ cp gs col7 0.50 shd ef gr gs col0 s gr 
+/Times-Roman ff 375.00 scf sf
+1590 1620 m
+gs 1 -1 sc (y) col0 sh gr
+/Times-Roman ff 375.00 scf sf
+7290 7200 m
+gs 1 -1 sc (x) col0 sh gr
+/Times-Roman ff 375.00 scf sf
+5505 8490 m
+gs 1 -1 sc (0) col0 sh gr
+/Times-Roman ff 375.00 scf sf
+4560 7365 m
+gs 1 -1 sc (90) col0 sh gr
+/Times-Roman ff 300.00 scf sf
+4935 7187 m
+gs 1 -1 sc (o) col0 sh gr
+/Times-Roman ff 300.00 scf sf
+5693 8311 m
+gs 1 -1 sc (o) col0 sh gr
+% here ends figure;
+$F2psEnd
+rs
+showpage
--- SwanBpntlist.ftn90	2009-08-05 23:35:15.000000000 +0200
+++ SwanBpntlist.ftn90	2009-08-05 21:33:10.000000000 +0200
@@ -185,7 +185,8 @@
                 endif
              enddo
              !
-             if ( vn == 0 .or. vert(vn)%atti(VMARKER) /= 1 ) goto 10
+             if ( vn == 0 ) goto 10
+             if ( vert(vn)%atti(VMARKER) /= 1 ) goto 10
              firstvert = .false.
              !
           endif
--- swancom1.ftn	2009-08-05 23:35:15.000000000 +0200
+++ swancom1.ftn	2009-08-05 21:33:22.000000000 +0200
@@ -426,7 +426,7 @@
 !
 !     PBOT(1)   = CFC      0.005    (Collins equation)
 !     PBOT(2)   = CFW      0.01     (Collins equation)
-!     PBOT(3)   = GAMJNS   0.0038   (Jonswap equation)
+!     PBOT(3)   = GAMJNS   0.067    (Jonswap equation)
 !     PBOT(4)   = MF      -0.08     (Madsen equation)
 !     PBOT(5)   = KN       0.05     (bottom roughness)
 !
@@ -1227,7 +1227,7 @@
 !       initialise Ursell number to 0 for each iteration                  40.03
 !       this is done in parallel within OpenMP environment                40.31
 !
-        IF (ITRIAD.GT.0) THEN
+        IF (ITRIAD.GT.0 .OR. ISURF.EQ.5) THEN                             41.03
            DO IP = I1GRD,I2GRD                                            40.31
             COMPDA(IP,JURSEL) = 0.
           ENDDO
@@ -4541,7 +4541,7 @@
           CALL STRSSB (IDDLOW  ,IDDTOP  ,
      &                 IDCMIN  ,IDCMAX  ,ISSTOP  ,CAX     ,CAY     ,
      &                 CAS     ,AC2     ,SPCSIG  ,IMATRA  ,
-     &                 ANYBLK  ,RDX     ,RDY                       )      40.41 30.21
+     &                 ANYBLK  ,RDX     ,RDY     ,TRAC0            )      41.07 40.41 30.21
 !
         END IF
       END IF
@@ -4970,7 +4970,7 @@
 !     *** calculate Ursell number ***
 !     *** update only for first encounter in a sweep                      40.16
 !
-      IF ( ITRIAD.GT.0 ) THEN                                             40.41
+      IF ( ITRIAD.GT.0 .OR. ISURF.EQ.5 ) THEN                             41.03 40.41
          IF (( SWPDIR .EQ. 1) .OR.                                        40.16
      &       ( SWPDIR .EQ. 2 .AND. IXCGRD(1) .EQ. 1) .OR.                 40.16
      &       ( SWPDIR .EQ. 3 .AND. IYCGRD(1) .EQ. 1) .OR.                 40.16
@@ -6035,7 +6035,7 @@
 !
 !     PBOT(1)   = CFC      0.005    (Putnam and Collins equation)
 !     PBOT(2)   = CFW      0.01     (Putnam and Collins equation)
-!     PBOT(3)   = GAMJNS   0.0038   (Jonswap equation)
+!     PBOT(3)   = GAMJNS   0.067    (Jonswap equation)
 !     PBOT(4)   = MF      -0.08     (Madsen et al. equation)
 !     PBOT(5)   = KN       0.05     (Madsen et al. bottom roughness)
 !
@@ -6243,8 +6243,8 @@
 !
 !         *** calculate surf breaking source term (5 formulations) ***    41.03
 !
-          CALL SSURF (ETOT    ,HM      ,
-     &                QBLOC   ,SMEBRK  ,AC2     ,IMATRA  ,                30.81
+          CALL SSURF (ETOT    ,HM      ,QBLOC   ,SMEBRK  ,
+     &                SPCSIG  ,AC2     ,IMATRA  ,                         30.81
      &                IMATDA  ,IDCMIN  ,IDCMAX  ,PLWBRK  ,
      &                ISSTOP  ,DISSC0  ,DISSC1  )                         40.67 40.61 30.21
 !
@@ -6325,14 +6325,7 @@
 !       *** geographical gridpoint a continuous spectrum    ***
 !       *** is present, i.e., after first iteration         ***
 !
-        IF ( ICUR .EQ. 0 .AND. ITER .GE. 1 ) THEN
-!
-            CALL SWLTA ( AC2   , DEP2  , CGO   , SPCSIG,                  40.55
-     &                   KWAVE , IMATRA, IMATDA, REDC0 , REDC1 ,          40.85
-     &                   IDDLOW, IDDTOP, ISSTOP, IDCMIN, IDCMAX,
-     &                   HS    , SMEBRK, PLTRI , URSELL )
-!
-        ELSE IF ( ICUR .EQ. 1 .AND. ITER .GT. 1 ) THEN
+        IF ( ICUR .EQ. 0 .OR. ITER .GT. 1 ) THEN
 !
             CALL SWLTA ( AC2   , DEP2  , CGO   , SPCSIG,                  40.55
      &                   KWAVE , IMATRA, IMATDA, REDC0 , REDC1 ,          40.85
@@ -8675,8 +8668,8 @@
 
                IF ( ( HSCURV.LE.PNUMS(15) .AND.
      &               (HSREL.LE.PNUMS(1) .OR. HSABS.LE.PNUMS(2)) ) .AND.   40.93
-     &              ( TMCURV.LE.PNUMS(15) .AND.                           40.93
-     &               (TMREL.LE.PNUMS(1) .OR. TMABS.LE.PNUMS(2)) ) ) THEN  40.93
+     &              ( TMCURV.LE.PNUMS(16) .AND.                           40.93
+     &               (TMREL.LE.PNUMS(1) .OR. TMABS.LE.PNUMS(3)) ) ) THEN  40.93
                   IACCURt = IACCURt + 1
                END IF
 
--- swancom2.ftn	2009-08-05 23:35:15.000000000 +0200
+++ swancom2.ftn	2009-08-05 21:33:22.000000000 +0200
@@ -611,8 +611,8 @@
 !
 !****************************************************************
 !
-      SUBROUTINE SSURF (ETOT    ,HM               ,                       30.81
-     &                  QB      ,SMEBRK  ,AC2     ,IMATRA  ,              30.81
+      SUBROUTINE SSURF (ETOT    ,HM      ,QB      ,SMEBRK  ,              30.81
+     &                  SPCSIG  ,AC2     ,IMATRA  ,                       30.81
      &                  IMATDA  ,IDCMIN  ,IDCMAX  ,PLWBRK  ,              30.81
      &                  ISSTOP  ,DISSC0  ,DISSC1  )                       40.67 40.61 30.81 30.21
 !
@@ -664,6 +664,7 @@
 !     40.61: Marcel Zijlema
 !     40.67: Nico Booij
 !     41.03: Andre van der Westhuysen
+!     41.06: Gerbrant van Vledder
 !
 !  1. Updates
 !
@@ -676,6 +677,7 @@
 !     40.61, Sep. 06: introduce DISSRF variable for output purposes
 !     40.67, Jun. 07: more accurate computation of dissipation terms
 !     41.03, Feb. 09: extension to alternative surf breaking formula's
+!     41.06, Mar. 09: extension to frequency dependent surf breaking
 !
 !  2. Purpose
 !
@@ -851,6 +853,7 @@
      &         IMATDA(MDC,MSC)      ,
      &         IMATRA(MDC,MSC)      ,
      &         PLWBRK(MDC,MSC,NPTST)                                      40.00
+      REAL     SPCSIG(MSC)
 !
       REAL     ETOT,  HM,  QB, SMEBRK                                     30.81
 !
@@ -872,10 +875,14 @@
 !     WS      Wavebreaking source term coefficient = DTOT/ETOT
 !     SbrD    Derivative of source term for surf breaking (Sbr) to action density
 !
-      INTEGER          ID,      IDDUM,   IENT,   IS
-      REAL             BIPH,    DEPLOC,  URSLOC, WH
-      DOUBLE PRECISION BB,      DIS0,    SbrD,
-     &                 SURFA0,  SURFA1,  WS
+      INTEGER          ID,       IDDUM,   IENT,   IS
+      REAL             BIPH,     DEPLOC,  URSLOC, WH,
+     &                 PP,       FAC,     EPTOT,  ETOT0,
+     &                 ECS(MDC), FMIN,    FMAX,   FRFAC(MSC)
+      DOUBLE PRECISION BB,       DIS0,    SbrD,
+     &                 SURFA0,   SURFA1,  WS  ,
+     &                 TEMP1 ,   TEMP2
+      REAL             SwanIntgratSpc
 !
 !
 !  7. Common blocks used
@@ -1007,7 +1014,28 @@
 !     *** store the results for surf wave breaking  ***
 !     *** in the matrices IMATDA and IMATRA         ***
 !
+      FRFAC = 1.                                                          41.06
+      IF (IFRSRF.EQ.1) THEN                                               41.06
+         PP    = PSURF(11)
+         FMIN  = PI2*PSURF(12)
+         FMAX  = PI2*PSURF(13)
+         ECS   = 1.
+         ETOT0 = SwanIntgratSpc(0., FMIN, FMAX, SPCSIG, ECS, SPCSIG,
+     &                          ECS, 0., 0., AC2(1,1,KCGRD(1)), 1 )
+         EPTOT = SwanIntgratSpc(PP, FMIN, FMAX, SPCSIG, ECS, SPCSIG,
+     &                          ECS, 0., 0., AC2(1,1,KCGRD(1)), 1 )
+         FAC   = ETOT0/EPTOT
+         IF ( ETOT0.GT.1.E-8 ) THEN
+            DO IS = 1, ISSTOP
+               FRFAC(IS) = FAC*SPCSIG(IS)**PP
+            END DO
+         END IF
+      END IF
+      TEMP1 = SURFA0
+      TEMP2 = SURFA1
       DO 101 IS = 1, ISSTOP
+        SURFA0 = TEMP1*FRFAC(IS)                                          41.06
+        SURFA1 = TEMP2*FRFAC(IS)                                          41.06
         DO 100 IDDUM = IDCMIN(IS), IDCMAX(IS)
           ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
           IMATDA(ID,IS) = IMATDA(ID,IS) + REAL(SURFA1)                    30.82
--- swancom5.ftn	2009-08-05 23:35:15.000000000 +0200
+++ swancom5.ftn	2009-08-05 21:33:22.000000000 +0200
@@ -1193,6 +1193,7 @@
 !     40.30: Marcel Zijlema
 !     40.41: Marcel Zijlema
 !     40.61: John Warner
+!     41.06: Gerbrant van Vledder
 !
 !  1. Updates
 !
@@ -1215,6 +1216,7 @@
 !     40.30, Mar. 03: correcting indices of test point with offsets MXF, MYF
 !     40.41, Oct. 04: common blocks replaced by modules, include files removed
 !     40.61, Dec. 06: correction DO loop 60 (IDCMIN, IDCMAX -> IDDLOW,IDDTOP)
+!     41.06, Mar. 09: add option of limitation of velocity in theta-direction
 !
 !  2. Purpose
 !
@@ -1362,6 +1364,9 @@
 !        CADT..        aux. quantities to compute Ctheta
 !        DPDX, DPDY    depth gradient
 !        DUXDX,DUXDY,DUYDX,DUYDY  current velocity gradients
+!        FAC           a factor
+!        FRLIM         frequency range in which limit on Ctheta is applied
+!        PP            power of the frequency dependent limiter on refraction
 !
       REAL     VLSINH ,KD1   ,COEF
       REAL     RDXL(2),RDYL(2),XC1   ,YC1    ,DET    ,
@@ -1372,6 +1377,7 @@
      &         CADT1    ,CADT2(3) ,CADT3(3) ,
      &         CADT4(3) ,CADT5(3) ,CADT6(3) ,CADT7(3)
       REAL  :: DLOC1, DLOC2, DLOC3
+      REAL  :: FAC, FRLIM, PP                                             41.06
 !     local depths corrected in view of stability                         40.02
 !
 !  5. Parameter variables
@@ -1854,6 +1860,19 @@
         ENDDO
       ENDIF
 !
+!     limit Ctheta in some frequency range if requested
+!
+      IF ( INT(PNUMS(29)) .EQ. 1 ) THEN
+         FRLIM = PI2*PNUMS(26)
+         PP    =     PNUMS(27)
+         DO IS = 1, MSC
+            FAC = MIN(1.,(SPCSIG(IS)/FRLIM)**PP)
+            DO ID = 1, MDC
+               CAD(ID,IS,1) = FAC*CAD(ID,IS,1)
+            ENDDO
+         ENDDO
+      ENDIF
+!
 !     *** test output ***
 !
       IF (TESTFL .AND. ITEST.GE.140) THEN                                 40.00
@@ -2989,7 +3008,7 @@
      &        '100km or 1 deg on a side, as a rule of thumb)')            40.08
 !        --- Note that code will stop even without this line              40.08
 !            (at beginning of next time step).                            40.08
-!         IF(MAXERR.LT.2) STOP                                             40.08
+!         IF(MAXERR.LT.2) STOP                                            40.08
       END IF                                                              40.08
 !
       DO 200 IS = 1, ISSTOP
@@ -3465,7 +3484,7 @@
       SUBROUTINE STRSSB (IDDLOW  ,IDDTOP  ,
      &                   IDCMIN  ,IDCMAX  ,ISSTOP  ,CAX     ,CAY     ,
      &                   CAS     ,AC2     ,SPCSIG  ,IMATRA  ,
-     &                   ANYBLK  ,RDX     ,RDY                       )    40.41 30.21
+     &                   ANYBLK  ,RDX     ,RDY     ,TRAC0            )    41.07 40.41 30.21
 !
 !****************************************************************
 !
@@ -3507,11 +3526,13 @@
 !
 !     30.72: IJsbrand Haagsma
 !     40.41: Marcel Zijlema
+!     41.07: Marcel Zijlema
 !
 !  1. Updates
 !
 !     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
 !     40.41, Oct. 04: common blocks replaced by modules, include files removed
+!     41.07, Jul. 09: also central scheme blended with upwind scheme
 !
 !  2. Purpose
 !
@@ -3562,19 +3583,30 @@
 !     ---------------------------------- = -----------------------
 !                   DS                                DS
 !
-!                  /  CAS(IS+0.5) * AC2(IS+1)    IF CAS(IS+0.5) < 0
+!                  /  CAS(IS+0.5) * ( (1-0.5mu)*AC2(IS+1) + 0.5mu*AC2(IS) )    IF CAS(IS+0.5) < 0
 !     F(IS+0.5) =  |
-!                  \  CAS(IS+0.5) * AC2(IS)      IF CAS(IS+0.5) > 0
+!                  \  CAS(IS+0.5) * ( (1-0.5mu)*AC2(IS) + 0.5mu*AC2(IS+1) )    IF CAS(IS+0.5) > 0
 !
-!                  /  CAS(IS-0.5) * AC2(IS-1)    IF CAS(IS-0.5) > 0
+!                  /  CAS(IS-0.5) * ( (1-0.5mu)*AC2(IS-1) + 0.5mu*AC2(IS) )    IF CAS(IS-0.5) > 0
 !     F(IS-0.5) =  |
-!                  \  CAS(IS-0.5) * AC2(IS)      IF CAS(IS-0.5) < 0
+!                  \  CAS(IS-0.5) * ( (1-0.5mu)*AC2(IS) + 0.5mu*AC2(IS-1) )    IF CAS(IS-0.5) < 0
+!
+!     with
+!
+!           0 <= mu <= 1 a blending factor
+!
+!           mu = 0 corresponds to 1st order upwind scheme
+!           mu = 1 corresponds to 2nd order central scheme
+!
+!           default value, mu = 0.5
+!
+!    and
 !
-!     WITH:
 !           CAS(IS+0.5) = ( CAS(IS+1) + CAS(IS) ) / 2.
 !
 !           CAS(IS-0.5) = ( CAS(IS) + CAS(IS-1) ) / 2.
 !
+!
 !     ------------------------------------------------------------
 !     Courant-Friedlich-Levich criterion :
 !
@@ -3677,7 +3709,7 @@
 !
       REAL     FSA     ,FLEFT   ,FRGHT   ,DS      ,CFLMAX  ,CFLCEN  ,
      &         CAXCEN  ,CAYCEN  ,CASCEN  ,TX      ,TY      ,TS      ,
-     &         CASL    ,CASR
+     &         CASL    ,CASR    ,PN1     ,PN2
 !
       REAL     CAS(MDC,MSC,ICMAX)       ,
      &         CAX(MDC,MSC,ICMAX)       ,
@@ -3686,6 +3718,7 @@
      &         IMATRA(MDC,MSC)          ,
      &         RDX(10)                  ,                                 40.08
      &         RDY(10)                                                    40.08
+      REAL  :: TRAC0(MDC,MSC,MTRNP)                                       41.07
 !
       INTEGER  IDCMIN(MSC)              ,
      &         IDCMAX(MSC)
@@ -3696,7 +3729,12 @@
       DATA IENT/0/
       IF (LTRACE) CALL STRACE (IENT,'STRSSB')
 !
-!     *** initialization of array and CFLMAX value ***
+!     --- determine blending factor
+!
+      PN1 = 0.5*(1.+PNUMS(7))                                             41.07
+      PN2 = 0.5*(1.-PNUMS(7))                                             41.07
+!
+!     *** initialization of ANYBLK and CFLMAX value ***
 !
       DO IS = 1, MSC
         DO ID = 1, MDC
@@ -3745,16 +3783,16 @@
               IF ( CASR .LT. 0. ) THEN
                 FRGHT = CASR * AC2(ID,IS+1,KCGRD(1))                      30.21
               ELSE
-                FRGHT = CASR * AC2(ID,IS,KCGRD(1))                        30.21
+                FRGHT = CASR * AC2(ID,IS  ,KCGRD(1))                      30.21
               END IF
               FLEFT = 0.
             ELSE IF ( IS .EQ. MSC ) THEN
 !             *** for the last discrete point in frequency space ***
 !             *** an upwind scheme is used                       ***
               CASL = CAS(ID,IS-1,1)
-              CASR = CAS(ID,IS,1  )
+              CASR = CAS(ID,IS  ,1)
               IF ( CASL .LT. 0. ) THEN
-                FLEFT = CASL * AC2(ID,IS,KCGRD(1))                        30.21
+                FLEFT = CASL * AC2(ID,IS  ,KCGRD(1))                      30.21
               ELSE
                 FLEFT = CASL * AC2(ID,IS-1,KCGRD(1))                      30.21
               END IF
@@ -3770,14 +3808,18 @@
               CASL  = 0.5 * ( CAS(ID,IS,1) + CAS(ID,IS-1,1) )
               CASR  = 0.5 * ( CAS(ID,IS,1) + CAS(ID,IS+1,1) )
               IF ( CASL .LT. 0. ) THEN
-                FLEFT = CASL * AC2(ID,IS,KCGRD(1))                        30.21
+                FLEFT = CASL * ( PN1*AC2(ID,IS  ,KCGRD(1)) +              41.07
+     &                           PN2*AC2(ID,IS-1,KCGRD(1)) )              41.07 30.21
               ELSE
-                FLEFT = CASL * AC2(ID,IS-1,KCGRD(1))                      30.21
+                FLEFT = CASL * ( PN1*AC2(ID,IS-1,KCGRD(1)) +              41.07
+     &                           PN2*AC2(ID,IS  ,KCGRD(1)) )              41.07 30.21
               END IF
               IF ( CASR .LT. 0. ) THEN
-                FRGHT = CASR * AC2(ID,IS+1,KCGRD(1))                      30.21
+                FRGHT = CASR * ( PN1*AC2(ID,IS+1,KCGRD(1)) +              41.07
+     &                           PN2*AC2(ID,IS  ,KCGRD(1)) )              41.07 30.21
               ELSE
-                FRGHT = CASR * AC2(ID,IS,KCGRD(1))                        30.21
+                FRGHT = CASR * ( PN1*AC2(ID,IS  ,KCGRD(1)) +              41.07
+     &                           PN2*AC2(ID,IS+1,KCGRD(1)) )              41.07 30.21
               END IF
             END IF
 !
@@ -3786,6 +3828,7 @@
 !           *** all the terms are known, store in IMATRA ***
 !
             IMATRA(ID,IS) = IMATRA(ID,IS) - FSA
+            TRAC0(ID,IS,3) = TRAC0(ID,IS,3) + FSA                         41.07
           ENDIF
 !
 !         *** test output ***
--- SwanCompUnstruc.ftn90	2009-08-05 23:35:15.000000000 +0200
+++ SwanCompUnstruc.ftn90	2009-08-05 21:33:22.000000000 +0200
@@ -34,12 +34,14 @@
 !   40.80: Marcel Zijlema
 !   40.85: Marcel Zijlema
 !   41.02: Marcel Zijlema
+!   41.07: Marcel Zijlema
 !
 !   Updates
 !
 !   40.80,     July 2007: New subroutine
 !   40.85,   August 2008: add propagation, generation and redistribution terms for output purposes
 !   41.02, February 2009: implementation of diffraction
+!   41.07,   August 2009: bug fix: never-ending sweep is prevented
 !
 !   Purpose
 !
@@ -223,6 +225,7 @@
     logical                               :: swpfull   ! indicate whether all necessary sweeps are done or not
     !
     logical, dimension(:,:), allocatable  :: anybin    ! true if bin is active in considered sweep
+    logical, dimension(:,:), allocatable  :: anyblk    ! true if bin is blocked by a counter current based on a CFL criterion
     logical, dimension(:), allocatable    :: anywnd    ! true if wind input is active in considered bin
     logical, dimension(:,:), allocatable  :: groww     ! check for each frequency whether the waves are growing or not
                                                        ! in a spectral direction
@@ -313,6 +316,7 @@
     allocate(reflso(MDC,MSC))
     allocate( alimw(MDC,MSC))
     allocate( groww(MDC,MSC))
+    allocate(anyblk(MDC,MSC))
     !
     allocate( disc0(MDC,MSC,MDISP))
     allocate( disc1(MDC,MSC,MDISP))
@@ -414,7 +418,7 @@
        acnrms = -9999.
        !
        if ( IQUAD  > 2 ) memnl4 = 0.
-       if ( ITRIAD > 0 ) compda(:,JURSEL) = 0.
+       if ( ITRIAD > 0 .or. ISURF == 5 ) compda(:,JURSEL) = 0.
        !
        compda(:,JDISS) = 0.
        compda(:,JLEAK) = 0.
@@ -487,7 +491,7 @@
        !
        if ( IDIFFR /= 0 ) call SwanDiffPar ( ac2, compda(1,JDP2), spcsig )
        !
-       ! all vertices are set untagged except non-active ones and those where boundary conditions are given
+       ! all vertices are set untagged except non-active ones
        !
        vert(:)%fullupdated = .false.
        !
@@ -514,17 +518,6 @@
              !
           endif
           !
-          if ( vert(kvert)%atti(VBC) /= 0 ) then      ! boundary condition given in vertex
-             !
-             vert(kvert)%fullupdated = .true.
-             !
-             do jc = 1, vert(kvert)%noc
-                icell = vert(kvert)%cell(jc)%atti(CELLID)
-                vert(kvert)%updated(jc) = icell
-             enddo
-             !
-          endif
-          !
        enddo
        !
        swpdir = 0
@@ -575,7 +568,7 @@
             !
             ivert = vlist(kvert)
             !
-            if ( vert(ivert)%atti(VBC) == 0 .and. vert(ivert)%active ) then   ! this active vertex needs to be updated
+            if ( vert(ivert)%active ) then   ! this active vertex needs to be updated
                !
                ! determine whether the present vertex is a test point
                !
@@ -617,7 +610,7 @@
                   l = 0
                   do k = 1, 2
                      do j = 1, vert(vu(k))%noc
-                        if ( vert(vu(k))%updated(j) /= 0 ) then     ! this upwave vertex is geographically updated
+                        if ( vert(vu(k))%updated(j) /= 0 .or. vert(vu(k))%atti(VMARKER) /= 0 ) then   ! this upwave vertex is geographically updated or known
                            l = l + 1
                            exit
                         endif
@@ -776,8 +769,8 @@
                                          cgo   , cax   , cay   , cad   , cas   , &
                                          anybin, rdx   , rdy   , spcsig, spcdir, &
                                          obredf, idcmin, idcmax, iscmin, iscmax, &
-                                         iddlow, iddtop, isslow, isstop, trac0 , &
-                                         trac1 )
+                                         iddlow, iddtop, isslow, isstop, anyblk, &
+                                         trac0 , trac1 )
 !TIMG                     call SWTSTO(118)
                      !
                      ! compute the source part of the action balance equation
@@ -870,6 +863,17 @@
                                            idcmax     )
 !TIMG                              call SWTSTO(120)
                               !
+                           elseif (int(PNUMS(8)) == 2 ) then
+                              !
+                              ! explicit scheme in sigma space
+                              ! solve tridiagonal system of equations using Thomas' algorithm
+                              !
+!TIMG                              call SWTSTA(120)
+                              call SOLMT1  ( idcmin     , idcmax     , ac2        , rhs    , &
+                                             amat(1,1,1), amat(1,1,5), amat(1,1,4),          &
+                                             isstop     , anyblk     , iddlow     , iddtop )
+!TIMG                              call SWTSTO(120)
+                              !
                            endif
                            !
                         endif
@@ -956,7 +960,7 @@
           !
           nwetp = 0
           do ivert = 1, nverts
-             if ( vert(ivert)%atti(VBC) == 0 .and. vert(ivert)%active ) nwetp = nwetp +1
+             if ( vert(ivert)%active ) nwetp = nwetp +1
           enddo
           !
           mxnfr = maxval(nrscal)
@@ -975,7 +979,7 @@
           !
           nwetp = 0
           do ivert = 1, nverts
-             if ( vert(ivert)%atti(VBC) == 0 .and. vert(ivert)%active ) nwetp = nwetp +1
+             if ( vert(ivert)%active ) nwetp = nwetp +1
           enddo
           !
           mxnfl = maxval(nflim)
@@ -1101,6 +1105,7 @@
     deallocate(reflso)
     deallocate( alimw)
     deallocate( groww)
+    deallocate(anyblk)
     !
     deallocate( disc0)
     deallocate( disc1)
--- SwanConvAccur.ftn90	2009-08-05 23:35:15.000000000 +0200
+++ SwanConvAccur.ftn90	2009-08-05 21:33:13.000000000 +0200
@@ -119,7 +119,7 @@
     !
     do ivert = 1, nverts
        !
-       if ( vert(ivert)%atti(VBC) == 0 .and. vert(ivert)%active ) then
+       if ( vert(ivert)%active ) then
           !
           nwetp = nwetp + 1
           !
@@ -137,7 +137,7 @@
     !
     do ivert = 1, nverts
        !
-       if ( vert(ivert)%atti(VBC) == 0 .and. vert(ivert)%active ) then
+       if ( vert(ivert)%active ) then
           !
           ! determine whether the present vertex is a test point
           !
--- SwanConvStopc.ftn90	2009-08-05 23:35:15.000000000 +0200
+++ SwanConvStopc.ftn90	2009-08-05 21:33:13.000000000 +0200
@@ -124,7 +124,7 @@
     !
     do ivert = 1, nverts
        !
-       if ( vert(ivert)%atti(VBC) == 0 .and. vert(ivert)%active ) then
+       if ( vert(ivert)%active ) then
           !
           ! determine whether the present vertex is a test point
           !
--- swan.edt	2009-08-05 23:35:15.000000000 +0200
+++ swan.edt	2009-08-05 23:25:05.000000000 +0200
@@ -114,12 +114,15 @@
 !        |                                                             |
 !        |    VAR [alpha] [gammin] [gammax] [gamneg] [coeff1] [coeff2] |
 !        |                                                             |
-!   BRE <     RUE [alpha] [a] [b]                                       >
+!   BRE <     RUE [alpha] [a] [b]                                       >   &
 !        |                                                             |
 !        |    TG  [alpha] [gamma] [pown]                               |
 !        |                                                             |
 !        !    BIP [alpha] [pown] [bref]                                |
 !
+!     ( FREQDep [power] [fmin] [fmax] )
+!
+!
 !              |             | -> CONstant [cfjon]
 !              | -> JONswap <
 !              |             |    VARiable [cfj1] [cfj2] [dsp1] [dsp2]
@@ -153,17 +156,21 @@
 !   PROP   /  BSBT
 !          \  GSE  [waveage] SEC|MIN|HR|DAY
 !
-!             | -> ACCUR [drel] [dhoval] [dtoval] [npnts] |
-!   NUMeric (<                                             >                &
-!             |    STOPC [dabs] [drel] [curvat] [npnts]   |
+!             | -> ACCUR [drel] [dhoval] [dtoval] [npnts]               |
+!   NUMeric (<                                                           >  &
+!             |    STOPC [dabs] [drel] [curvat] [npnts] [dtabs] [curvt] |
 !
 !                    | -> STAT  [mxitst] [alfa] |
 !                   <                            >  [limiter]   )           &
 !                    | NONSTat  [mxitns]        |
 !
-!           ( DIRimpl [cdd] [cdlim]                             )           &
+!           ( DIRimpl [cdd] [cdlim]  WNUMber                    )           &
+!
+!           ( REFRLim [frlim] [power]                           )           &
 !
-!           ( SIGIMpl [css] [eps2] [outp] [niter]               )           &
+!              | -> SIGIMpl [css] [eps2] [outp] [niter]
+!           ( <                                                 )           &
+!              |    SIGEXpl [css] [cfl]
 !
 !           ( SETUP [eps2] [outp] [niter]                       )
 !
--- swanimp.tex	2009-08-05 23:35:15.000000000 +0200
+++ swanimp.tex	2009-08-05 21:33:07.000000000 +0200
@@ -15,7 +15,7 @@
 \end{center}
 \vfill
 \begin{center}
-{\Large\bf SWAN Cycle III version 40.72AB}
+{\Large\bf SWAN Cycle III version 40.72ABC}
 \end{center}
 
 \cleardoublepage
@@ -204,6 +204,9 @@
 can be done automatically or manually; see Chapter \ref{ch:instal}.
 \\[2ex]
 \noindent
+For any use of the SWAN source code in your environment, proper reference must be made to the origin of the software!
+\\[2ex]
+\noindent
 You are allow to make changes in the source code of SWAN, but Delft University of Technology
 will not support modified versions of SWAN. If you ever want your modifications to be
 implemented in the authorized version of SWAN (the version on the SWAN web page), you need
--- SwanInterpolateOutput.ftn90	2009-08-05 23:35:15.000000000 +0200
+++ SwanInterpolateOutput.ftn90	2009-08-05 21:33:16.000000000 +0200
@@ -1,4 +1,4 @@
-subroutine SwanInterpolateOutput ( foutp, x, y, finp, mip, excval )
+subroutine SwanInterpolateOutput ( foutp, x, y, finp, mip, kvert, excval )
 !
 !   --|-----------------------------------------------------------|--
 !     | Delft University of Technology                            |
@@ -33,11 +33,13 @@
 !
 !   40.80: Marcel Zijlema
 !   40.90: Nico Booij
+!   41.07: Marcel Zijlema
 !
 !   Updates
 !
 !   40.80, August 2007: New subroutine
 !   40.90,   June 2008: improved interpolation near obstacles
+!   41.07,   July 2009: optimization
 !
 !   Purpose
 !
@@ -63,6 +65,7 @@
 !
 !   Argument variables
 !
+    integer, dimension(mip), intent(in) :: kvert  ! vertex indices of output points
     integer, intent(in)                 :: mip    ! number of given points
     !
     real, intent(in)                    :: excval ! exception value for output quantity
@@ -137,9 +140,9 @@
     !
     pointloop: do ip = 1, mip
        !
-       ! find closest vertex for given point
+       ! assign vertex index of given point
        !
-       call SwanFindPoint ( x(ip), y(ip), ivert )
+       ivert = kvert(ip)
        !
        ! if point not found, go to next point
        !
--- swanmain.ftn	2009-08-05 23:35:15.000000000 +0200
+++ swanmain.ftn	2009-08-05 21:33:22.000000000 +0200
@@ -857,6 +857,7 @@
       WRITE (VERTXT, '(F5.2)') VERNUM                                     40.03
       CALL BUGFIX ('A')
       CALL BUGFIX ('B')
+      CALL BUGFIX ('C')
 !
       CALL OCPINI ('swaninit', .TRUE.,INERR)                              34.01
       IF (INERR.GT.0) RETURN                                              34.01
@@ -1053,6 +1054,10 @@
 !
       CSETUP = .TRUE.                                                     30.82
 !
+!     flag for frequency dependent surf breaking                          41.06
+!
+      IFRSRF = 0                                                          41.06
+!
 !     PSETUP(1) is currently unused, but can be used as setup nesting flag
 !     PSETUP(2) is the user defined correction for the level of the setup
 !
@@ -1156,6 +1161,16 @@
 !
       PNUMS(30) = 0.00                                                    40.23
 !
+!     --- parameters for limiting Ctheta
+!
+      PNUMS(26) = 0.2                                                     41.06
+      PNUMS(27) = 2.0                                                     41.06
+      PNUMS(29) = 0.0                                                     41.06
+!
+!     --- computation of Ctheta based on wave number
+!
+      PNUMS(32) = 0.                                                      41.07
+!
 !     *** (1) and (2): Komen et al. (1984) formulation ***
 !
       PWCAP(1)  = 2.36E-5
--- swanout1.ftn	2009-08-05 23:35:15.000000000 +0200
+++ swanout1.ftn	2009-08-05 21:33:22.000000000 +0200
@@ -1438,7 +1438,7 @@
 !     MIP     Int    input    number of output points
 !     XP, YP  real   outp     user coordinates of output point
 !     XC, YC  real   outp     comp. grid coordinates
-!     WX2, WY2  real   input    wind components
+!     WX2, WY2  real input    wind components
 !
 !  9. Subroutines calling
 !
@@ -1458,6 +1458,7 @@
       REAL       FORCE(nverts,2)                                          40.80
       INTEGER    VOQR(*), KGRPNT(MXC,MYC)
       INTEGER    IONOD(*)                                                 40.31
+      INTEGER, ALLOCATABLE :: KVERT(:)                                    41.07
       LOGICAL    OQPROC(*), EQREAL                                        30.72
       LOGICAL    CROSS(4,MIP)                                             40.86
       LOGICAL, ALLOCATABLE :: LTMP(:)                                     40.91
@@ -1495,6 +1496,17 @@
         ENDDO
       ENDIF
 !
+      IF (OPTG.EQ.5) THEN
+!
+!        ---find closest vertex for given point in case of unstructured grid
+!
+         ALLOCATE(KVERT(MIP))
+         DO IP = 1, MIP
+            CALL SwanFindPoint ( VOQ(IP,1), VOQ(IP,2), KVERT(IP) )        41.07
+         ENDDO
+!
+      ENDIF
+!
 !     depth
 !
  120  IF (OQPROC(4)) THEN
@@ -1507,7 +1519,7 @@
         ELSE                                                              40.80
            CALL SwanInterpolateOutput ( VOQ(1,VOQR(4)), VOQ(1,1),         40.80
      &                                  VOQ(1,2), COMPDA(1,JDP2),         40.80
-     &                                  MIP, OVEXCV(4) )                  40.80
+     &                                  MIP, KVERT, OVEXCV(4) )           40.80
         ENDIF                                                             40.80
       ENDIF
 !
@@ -1527,10 +1539,10 @@
           ELSE                                                            40.80
              CALL SwanInterpolateOutput ( VOQ(1,JVQX), VOQ(1,1),          40.80
      &                                    VOQ(1,2), COMPDA(1,JVX2),       40.80
-     &                                    MIP, OVEXCV(5) )                40.80
+     &                                    MIP, KVERT, OVEXCV(5) )         40.80
              CALL SwanInterpolateOutput ( VOQ(1,JVQY), VOQ(1,1),          40.80
      &                                    VOQ(1,2), COMPDA(1,JVY2),       40.80
-     &                                    MIP, OVEXCV(5) )                40.80
+     &                                    MIP, KVERT, OVEXCV(5) )         40.80
           ENDIF                                                           40.80
           DO IP = 1, MIP
             UXLOC = VOQ(IP,JVQX)
@@ -1557,7 +1569,7 @@
         ELSE                                                              40.80
            CALL SwanInterpolateOutput ( VOQ(1,VOQR(6)), VOQ(1,1),         40.80
      &                                  VOQ(1,2), COMPDA(1,JUBOT),        40.80
-     &                                  MIP, OVEXCV(6) )                  40.80
+     &                                  MIP, KVERT, OVEXCV(6) )           40.80
         ENDIF                                                             40.80
         KK = VOQR(6)
         RR = SQRT(2.)
@@ -1578,7 +1590,7 @@
         ELSE                                                              40.80
            CALL SwanInterpolateOutput ( VOQ(1,VOQR(34)), VOQ(1,1),        40.80
      &                                  VOQ(1,2), COMPDA(1,JUBOT),        40.80
-     &                                  MIP, OVEXCV(34) )                 40.80
+     &                                  MIP, KVERT, OVEXCV(34) )          40.80
         ENDIF                                                             40.80
       ENDIF
 !
@@ -1594,7 +1606,7 @@
            ELSE                                                           40.80
               CALL SwanInterpolateOutput ( VOQ(1,VOQR(50)), VOQ(1,1),     40.80
      &                                     VOQ(1,2), COMPDA(1,JPBOT),     40.80
-     &                                     MIP, OVEXCV(50) )              40.80
+     &                                     MIP, KVERT, OVEXCV(50) )       40.80
            ENDIF                                                          40.80
         ELSE
            DO IP = 1, MIP                                                 40.65
@@ -1614,7 +1626,7 @@
         ELSE                                                              40.80
            CALL SwanInterpolateOutput ( VOQ(1,VOQR(7)), VOQ(1,1),         40.80
      &                                  VOQ(1,2), COMPDA(1,JDISS),        40.80
-     &                                  MIP, OVEXCV(7) )                  40.80
+     &                                  MIP, KVERT, OVEXCV(7) )           40.80
         ENDIF                                                             40.80
         IF (INRHOG.EQ.1) THEN
           DO 152 IP = 1, MIP
@@ -1636,7 +1648,7 @@
            ELSE                                                           40.80
               CALL SwanInterpolateOutput ( VOQ(1,VOQR(54)), VOQ(1,1),     40.80
      &                                     VOQ(1,2), COMPDA(1,JDSXB),     40.80
-     &                                     MIP, OVEXCV(54) )              40.80
+     &                                     MIP, KVERT, OVEXCV(54) )       40.80
            ENDIF                                                          40.80
         ELSE
            DO IP = 1, MIP                                                 40.65
@@ -1663,7 +1675,7 @@
            ELSE                                                           40.80
               CALL SwanInterpolateOutput ( VOQ(1,VOQR(55)), VOQ(1,1),     40.80
      &                                     VOQ(1,2), COMPDA(1,JDSXS),     40.80
-     &                                     MIP, OVEXCV(55) )              40.80
+     &                                     MIP, KVERT, OVEXCV(55) )       40.80
            ENDIF                                                          40.80
         ELSE
            DO IP = 1, MIP                                                 40.65
@@ -1690,7 +1702,7 @@
            ELSE                                                           40.80
               CALL SwanInterpolateOutput ( VOQ(1,VOQR(56)), VOQ(1,1),     40.80
      &                                     VOQ(1,2), COMPDA(1,JDSXW),     40.80
-     &                                     MIP, OVEXCV(56) )              40.80
+     &                                     MIP, KVERT, OVEXCV(56) )       40.80
            ENDIF                                                          40.80
         ELSE
            DO IP = 1, MIP                                                 40.65
@@ -1717,7 +1729,7 @@
            ELSE
               CALL SwanInterpolateOutput ( VOQ(1,VOQR(60)), VOQ(1,1),
      &                                     VOQ(1,2), COMPDA(1,JGENR),
-     &                                     MIP, OVEXCV(60) )
+     &                                     MIP, KVERT, OVEXCV(60) )
            ENDIF
         ELSE
            DO IP = 1, MIP
@@ -1744,7 +1756,7 @@
            ELSE
               CALL SwanInterpolateOutput ( VOQ(1,VOQR(61)), VOQ(1,1),
      &                                     VOQ(1,2), COMPDA(1,JGSXW),
-     &                                     MIP, OVEXCV(61) )
+     &                                     MIP, KVERT, OVEXCV(61) )
            ENDIF
         ELSE
            DO IP = 1, MIP
@@ -1771,7 +1783,7 @@
            ELSE
               CALL SwanInterpolateOutput ( VOQ(1,VOQR(62)), VOQ(1,1),
      &                                     VOQ(1,2), COMPDA(1,JREDS),
-     &                                     MIP, OVEXCV(62) )
+     &                                     MIP, KVERT, OVEXCV(62) )
            ENDIF
         ELSE
            DO IP = 1, MIP
@@ -1798,7 +1810,7 @@
            ELSE
               CALL SwanInterpolateOutput ( VOQ(1,VOQR(63)), VOQ(1,1),
      &                                     VOQ(1,2), COMPDA(1,JRSXQ),
-     &                                     MIP, OVEXCV(63) )
+     &                                     MIP, KVERT, OVEXCV(63) )
            ENDIF
         ELSE
            DO IP = 1, MIP
@@ -1825,7 +1837,7 @@
            ELSE
               CALL SwanInterpolateOutput ( VOQ(1,VOQR(64)), VOQ(1,1),
      &                                     VOQ(1,2), COMPDA(1,JRSXT),
-     &                                     MIP, OVEXCV(64) )
+     &                                     MIP, KVERT, OVEXCV(64) )
            ENDIF
         ELSE
            DO IP = 1, MIP
@@ -1852,7 +1864,7 @@
            ELSE
               CALL SwanInterpolateOutput ( VOQ(1,VOQR(65)), VOQ(1,1),
      &                                     VOQ(1,2), COMPDA(1,JTRAN),
-     &                                     MIP, OVEXCV(65) )
+     &                                     MIP, KVERT, OVEXCV(65) )
            ENDIF
         ELSE
            DO IP = 1, MIP
@@ -1879,7 +1891,7 @@
            ELSE
               CALL SwanInterpolateOutput ( VOQ(1,VOQR(66)), VOQ(1,1),
      &                                     VOQ(1,2), COMPDA(1,JTSXG),
-     &                                     MIP, OVEXCV(66) )
+     &                                     MIP, KVERT, OVEXCV(66) )
            ENDIF
         ELSE
            DO IP = 1, MIP
@@ -1906,7 +1918,7 @@
            ELSE
               CALL SwanInterpolateOutput ( VOQ(1,VOQR(67)), VOQ(1,1),
      &                                     VOQ(1,2), COMPDA(1,JTSXT),
-     &                                     MIP, OVEXCV(67) )
+     &                                     MIP, KVERT, OVEXCV(67) )
            ENDIF
         ELSE
            DO IP = 1, MIP
@@ -1933,7 +1945,7 @@
            ELSE
               CALL SwanInterpolateOutput ( VOQ(1,VOQR(68)), VOQ(1,1),
      &                                     VOQ(1,2), COMPDA(1,JTSXS),
-     &                                     MIP, OVEXCV(68) )
+     &                                     MIP, KVERT, OVEXCV(68) )
            ENDIF
         ELSE
            DO IP = 1, MIP
@@ -1960,7 +1972,7 @@
            ELSE
               CALL SwanInterpolateOutput ( VOQ(1,VOQR(69)), VOQ(1,1),
      &                                     VOQ(1,2), COMPDA(1,JRADS),
-     &                                     MIP, OVEXCV(69) )
+     &                                     MIP, KVERT, OVEXCV(69) )
            ENDIF
         ELSE
            DO IP = 1, MIP
@@ -1986,7 +1998,7 @@
         ELSE                                                              40.80
            CALL SwanInterpolateOutput ( VOQ(1,VOQR(8)), VOQ(1,1),         40.80
      &                                  VOQ(1,2), COMPDA(1,JQB),          40.80
-     &                                  MIP, OVEXCV(8) )                  40.80
+     &                                  MIP, KVERT, OVEXCV(8) )           40.80
         ENDIF                                                             40.80
       ENDIF
 !
@@ -2006,10 +2018,10 @@
           ELSE                                                            40.80
              CALL SwanInterpolateOutput ( VOQ(1,JVQX), VOQ(1,1),          40.80
      &                                    VOQ(1,2), COMPDA(1,JWX2),       40.80
-     &                                    MIP, OVEXCV(26) )               40.80
+     &                                    MIP, KVERT, OVEXCV(26) )        40.80
              CALL SwanInterpolateOutput ( VOQ(1,JVQY), VOQ(1,1),          40.80
      &                                    VOQ(1,2), COMPDA(1,JWY2),       40.80
-     &                                    MIP, OVEXCV(26) )               40.80
+     &                                    MIP, KVERT, OVEXCV(26) )        40.80
           ENDIF                                                           40.80
           DO IP = 1, MIP
             UXLOC = VOQ(IP,JVQX)
@@ -2036,7 +2048,7 @@
         ELSE                                                              40.80
            CALL SwanInterpolateOutput ( VOQ(1,VOQR(30)), VOQ(1,1),        40.80
      &                                  VOQ(1,2), COMPDA(1,JDHS),         40.80
-     &                                  MIP, OVEXCV(30) )                 40.80
+     &                                  MIP, KVERT, OVEXCV(30) )          40.80
         ENDIF                                                             40.80
       ENDIF
 !
@@ -2051,7 +2063,7 @@
         ELSE                                                              40.80
            CALL SwanInterpolateOutput ( VOQ(1,VOQR(31)), VOQ(1,1),        40.80
      &                                  VOQ(1,2), COMPDA(1,JDTM),         40.80
-     &                                  MIP, OVEXCV(31) )                 40.80
+     &                                  MIP, KVERT, OVEXCV(31) )          40.80
         ENDIF                                                             40.80
       ENDIF
 !
@@ -2066,7 +2078,7 @@
         ELSE                                                              40.80
            CALL SwanInterpolateOutput ( VOQ(1,VOQR(9)), VOQ(1,1),         40.80
      &                                  VOQ(1,2), COMPDA(1,JLEAK),        40.80
-     &                                  MIP, OVEXCV(9) )                  40.80
+     &                                  MIP, KVERT, OVEXCV(9) )           40.80
         ENDIF                                                             40.80
         IF (INRHOG.EQ.1) THEN
           DO 202 IP = 1, MIP
@@ -2088,7 +2100,7 @@
           ELSE                                                            40.80
              CALL SwanInterpolateOutput ( VOQ(1,VOQR(35)), VOQ(1,1),      40.80
      &                                    VOQ(1,2), COMPDA(1,JUSTAR),     40.80
-     &                                    MIP, OVEXCV(35) )               40.80
+     &                                    MIP, KVERT, OVEXCV(35) )        40.80
           ENDIF                                                           40.80
         ELSE
           DO IP = 1, MIP                                                  31.02
@@ -2109,7 +2121,7 @@
           ELSE                                                            40.80
              CALL SwanInterpolateOutput ( VOQ(1,VOQR(36)), VOQ(1,1),      40.80
      &                                    VOQ(1,2), COMPDA(1,JZEL),       40.80
-     &                                    MIP, OVEXCV(36) )               40.80
+     &                                    MIP, KVERT, OVEXCV(36) )        40.80
           ENDIF                                                           40.80
         ELSE
           DO IP = 1, MIP                                                  31.02
@@ -2130,7 +2142,7 @@
           ELSE                                                            40.80
              CALL SwanInterpolateOutput ( VOQ(1,VOQR(37)), VOQ(1,1),      40.80
      &                                    VOQ(1,2), COMPDA(1,JTAUW),      40.80
-     &                                    MIP, OVEXCV(37) )               40.80
+     &                                    MIP, KVERT, OVEXCV(37) )        40.80
           ENDIF                                                           40.80
         ELSE
           DO IP = 1, MIP                                                  31.02
@@ -2151,7 +2163,7 @@
           ELSE                                                            40.80
              CALL SwanInterpolateOutput ( VOQ(1,VOQR(38)), VOQ(1,1),      40.80
      &                                    VOQ(1,2), COMPDA(1,JCDRAG),     40.80
-     &                                    MIP, OVEXCV(38) )               40.80
+     &                                    MIP, KVERT, OVEXCV(38) )        40.80
           ENDIF                                                           40.80
         ELSE
           DO IP = 1, MIP                                                  31.02
@@ -2172,7 +2184,7 @@
           ELSE                                                            40.80
              CALL SwanInterpolateOutput ( VOQ(1,VOQR(39)), VOQ(1,1),      40.80
      &                                    VOQ(1,2), COMPDA(1,JSETUP),     40.80
-     &                                    MIP, OVEXCV(39) )               40.80
+     &                                    MIP, KVERT, OVEXCV(39) )        40.80
           ENDIF                                                           40.80
         ELSE                                                              32.02
           DO IP = 1, MIP                                                  32.02
@@ -2188,10 +2200,10 @@
      &  VOQR(20)                                                          40.80
         CALL SwanInterpolateOutput ( VOQ(1,VOQR(20)), VOQ(1,1),           40.80
      &                               VOQ(1,2), FORCE(1,1),                40.80
-     &                               MIP, OVEXCV(20) )                    40.80
+     &                               MIP, KVERT, OVEXCV(20) )             40.80
         CALL SwanInterpolateOutput ( VOQ(1,VOQR(20)+1), VOQ(1,1),         40.80
      &                               VOQ(1,2), FORCE(1,2),                40.80
-     &                               MIP, OVEXCV(20) )                    40.80
+     &                               MIP, KVERT, OVEXCV(20) )             40.80
       ENDIF                                                               40.80
 !
 !     Ursell
@@ -2205,7 +2217,7 @@
         ELSE                                                              40.80
            CALL SwanInterpolateOutput ( VOQ(1,VOQR(45)), VOQ(1,1),        40.80
      &                                  VOQ(1,2), COMPDA(1,JURSEL),       40.80
-     &                                  MIP, OVEXCV(45) )                 40.80
+     &                                  MIP, KVERT, OVEXCV(45) )          40.80
         ENDIF                                                             40.80
       ENDIF
 !
@@ -2221,7 +2233,7 @@
            ELSE                                                           40.80
               CALL SwanInterpolateOutput ( VOQ(1,VOQR(46)), VOQ(1,1),     40.80
      &                                     VOQ(1,2), COMPDA(1,JASTD2),    40.80
-     &                                     MIP, OVEXCV(46) )              40.80
+     &                                     MIP, KVERT, OVEXCV(46) )       40.80
            ENDIF                                                          40.80
         ELSE
            DO IP = 1, MIP
@@ -2242,7 +2254,7 @@
           ELSE                                                            40.80
              CALL SwanInterpolateOutput ( VOQ(1,VOQR(49)), VOQ(1,1),      40.80
      &                                    VOQ(1,2), DIFPARAM(:),          40.80
-     &                                    MIP, OVEXCV(49) )               40.80
+     &                                    MIP, KVERT, OVEXCV(49) )        40.80
           ENDIF                                                           40.80
         ELSE
           DO IP = 1, MIP                                                  40.21
@@ -2279,7 +2291,7 @@
                vert(:)%active = .TRUE.                                    40.91
                CALL SwanInterpolateOutput ( VOQ(1,VOQR(27)), VOQ(1,1),    40.80
      &                                      VOQ(1,2), COMPDA(1,JFRC2),    40.80
-     &                                      MIP, OVEXCV(27) )             40.80
+     &                                      MIP, KVERT, OVEXCV(27) )      40.80
                vert(:)%active = LTMP(:)                                   40.91
                DEALLOCATE(LTMP)                                           40.91
             ENDIF                                                         40.80
@@ -2314,7 +2326,7 @@
               vert(:)%active = .TRUE.                                     40.91
               CALL SwanInterpolateOutput ( VOQ(1,VOQR(51)), VOQ(1,1),     40.80
      &                                     VOQ(1,2), COMPDA(1,JWLV2),     40.80
-     &                                     MIP, OVEXCV(51) )              40.80
+     &                                     MIP, KVERT, OVEXCV(51) )       40.80
               vert(:)%active = LTMP(:)                                    40.91
               DEALLOCATE(LTMP)                                            40.91
            ENDIF                                                          40.80
@@ -2349,7 +2361,7 @@
               vert(:)%active = .TRUE.                                     40.91
               CALL SwanInterpolateOutput ( VOQ(1,VOQR(52)), VOQ(1,1),     40.80
      &                                     VOQ(1,2), COMPDA(1,JBOTLV),    40.80
-     &                                     MIP, OVEXCV(52) )              40.80
+     &                                     MIP, KVERT, OVEXCV(52) )       40.80
               vert(:)%active = LTMP(:)                                    40.91
               DEALLOCATE(LTMP)                                            40.91
            ENDIF                                                          40.80
@@ -2384,15 +2396,15 @@
          IYE = MYC-IHALOY                                                 40.31
          IF ( LMYL ) IYE = MYC                                            40.41 40.31
          DO IP = 1, MIP                                                   40.31
-            IX = NINT(XC(IP)+3.) - 2                                      40.31
-            IY = NINT(YC(IP)+3.) - 2                                      40.31
-            IF ( IX.EQ.0 ) IX = 1                                         40.31
-            IF ( IY.EQ.0 ) IY = 1                                         40.31
+            IX = NINT(XC(IP)+100.) - 99                                   41.07 40.31
+            IY = NINT(YC(IP)+100.) - 99                                   41.07 40.31
             IF ( IX.GE.IXB .AND. IX.LE.IXE .AND.                          40.31
      &           IY.GE.IYB .AND. IY.LE.IYE ) IONOD(IP) = INODE            40.31
          END DO                                                           40.31
       END IF                                                              40.31
-
+!
+      IF (ALLOCATED(KVERT)) DEALLOCATE(KVERT)                             41.07
+!
       RETURN
       END
 !***********************************************************************
--- swanout2.ftn	2009-08-05 23:35:15.000000000 +0200
+++ swanout2.ftn	2009-08-05 21:33:22.000000000 +0200
@@ -162,15 +162,15 @@
       ENDIF
       IF (ITEST.GE.90) WRITE (PRTEST, 21)  RTYPE,NREF, OQI(3)             40.31 40.03
   21  FORMAT (' Test SWBLOK: RTYPE NREF NVAR ',A4,2(1X,I6))
+      FILENM = OUTP_FILES(OQI(2))                                         40.31 40.13
+      MATLAB = INDEX( FILENM, '.MAT' ).NE.0 .OR.                          40.41 40.30
+     &         INDEX (FILENM, '.mat' ).NE.0                               40.41 40.30
       IF (NREF.EQ.0) THEN
-        FILENM = OUTP_FILES(OQI(2))                                       40.31 40.13
         IOSTAT = -1                                                       20.75
         CALL FOR (NREF, FILENM, 'UF', IOSTAT)
         IF (STPNOW()) RETURN                                              34.01
         OQI(1) = NREF                                                     40.31 30.00
         OUTP_FILES(OQI(2)) = FILENM                                       40.41
-        MATLAB = INDEX( FILENM, '.MAT' ).NE.0 .OR.                        40.41 40.30
-     &           INDEX (FILENM, '.mat' ).NE.0                             40.41 40.30
         IF (MATLAB) THEN                                                  40.30
            CLOSE(NREF)                                                    40.30
            OPEN(UNIT=NREF, FILE=FILENM, FORM='UNFORMATTED',               40.30
--- swanparll.ftn	2009-08-05 23:35:15.000000000 +0200
+++ swanparll.ftn	2009-08-05 21:33:22.000000000 +0200
@@ -3730,14 +3730,13 @@
             END DO
          END IF
          IPLOOP : DO IP = 1, MIP                                          40.51
-            IF ( XC(IP).LT.-0.01 .OR. YC(IP).LT.-0.01 .OR.                40.51
-     &           XC(IP).GT.REAL(MXCGL-1)+0.01 .OR.                        40.51
-     &           YC(IP).GT.REAL(MYCGL-1)+0.01 ) THEN                      40.51
+            IXK = NINT(XC(IP)+100.)-99                                    41.07 40.51
+            IYK = NINT(YC(IP)+100.)-99                                    41.07 40.51
+            IF ( IXK.LT.1 .OR. IYK.LT.1 .OR. IXK.GT.MXCGL .OR.            41.07 40.51
+     &           IYK.GT.MYCGL ) THEN                                      41.07 40.51
                CALL WREXCV                                                40.51
                CYCLE IPLOOP                                               40.51
             END IF                                                        40.51
-            IXK = NINT(XC(IP)+3.)-2                                       40.51
-            IYK = NINT(YC(IP)+3.)-2                                       40.51
             IPROC = 1                                                     40.51
             PROCLOOP : DO                                                 40.51
               IF ( NINT(BLKND(IXK,IYK)).EQ.IPROC ) THEN                   40.51
@@ -4025,8 +4024,8 @@
      &              YC(IP).GT.REAL(MYCGL-1)+0.01 ) THEN                   40.51
                   IBLKN = NPROC+1                                         40.51
                ELSE                                                       40.51
-                  IXK   = NINT(XC(IP)+3.)-2                               40.51
-                  IYK   = NINT(YC(IP)+3.)-2                               40.51
+                  IXK   = NINT(XC(IP)+100.)-99                            41.07 40.51
+                  IYK   = NINT(YC(IP)+100.)-99                            41.07 40.51
                   IBLKN = NINT(BLKND(IXK,IYK))                            40.51
                END IF                                                     40.51
                EMPTY = .TRUE.
@@ -4091,14 +4090,13 @@
       END IF
 
       IPLOOP : DO IP = 1, MIP                                             40.51
-         IF ( XC(IP).LT.-0.01 .OR. YC(IP).LT.-0.01 .OR.                   40.51
-     &        XC(IP).GT.REAL(MXCGL-1)+0.01 .OR.                           40.51
-     &        YC(IP).GT.REAL(MYCGL-1)+0.01 ) THEN                         40.51
+         IXK = NINT(XC(IP)+100.)-99                                       41.07 40.51
+         IYK = NINT(YC(IP)+100.)-99                                       41.07 40.51
+         IF ( IXK.LT.1 .OR. IYK.LT.1 .OR. IXK.GT.MXCGL .OR.               41.07 40.51
+     &        IYK.GT.MYCGL ) THEN                                         41.07 40.51
             WRITE (NREF, '(A)') 'NODATA'                                  40.51
             CYCLE IPLOOP                                                  40.51
          END IF                                                           40.51
-         IXK = NINT(XC(IP)+3.)-2                                          40.51
-         IYK = NINT(YC(IP)+3.)-2                                          40.51
          IPROC = 1                                                        40.51
          PROCLOOP : DO                                                    40.51
            IF ( NINT(BLKND(IXK,IYK)).EQ.IPROC ) THEN                      40.51
@@ -4394,8 +4392,8 @@
          VOQ = OVEXCV(IVTYPE)
 
          DO IP = 1, MXK*MYK                                               40.51
-            IXK = NINT(XC(IP)+3.)-2                                       40.51
-            IYK = NINT(YC(IP)+3.)-2                                       40.51
+            IXK = NINT(XC(IP)+100.)-99                                    41.07 40.51
+            IYK = NINT(YC(IP)+100.)-99                                    41.07 40.51
             IPROC = 1                                                     40.51
             PROCLOOP1 : DO                                                40.51
               IF ( NINT(BLKND(IXK,IYK)).EQ.IPROC ) THEN                   40.51
@@ -4416,8 +4414,8 @@
          IF ( OVSVTY(IVTYPE).GE.3 ) THEN
 
             DO IP = 1, MXK*MYK                                            40.51
-               IXK = NINT(XC(IP)+3.)-2                                    40.51
-               IYK = NINT(YC(IP)+3.)-2                                    40.51
+               IXK = NINT(XC(IP)+100.)-99                                 41.07 40.51
+               IYK = NINT(YC(IP)+100.)-99                                 41.07 40.51
                IPROC = 1                                                  40.51
                PROCLOOP2 : DO                                             40.51
                  IF ( NINT(BLKND(IXK,IYK)).EQ.IPROC ) THEN                40.51
--- swanpgr.tex	2009-08-05 23:35:15.000000000 +0200
+++ swanpgr.tex	2009-08-05 21:33:07.000000000 +0200
@@ -679,10 +679,7 @@
 !     | Environmental Fluid Mechanics Section                     |
 !     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
 !     |                                                           |
-!     | Programmers: R.C. Ris, N. Booij,                          |
-!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
-!     |              M. Zijlema, E.E. Kriezi,                     |
-!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
+!     | Programmers: The SWAN team                                |
 !   --|-----------------------------------------------------------|--
 !
 !
@@ -874,10 +871,7 @@
 !     | Environmental Fluid Mechanics Section                     |
 !     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
 !     |                                                           |
-!     | Programmers: R.C. Ris, N. Booij,                          |
-!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
-!     |              M. Zijlema, E.E. Kriezi,                     |
-!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
+!     | Programmers: The SWAN team                                |
 !   --|-----------------------------------------------------------|--
 !
 !
--- swanpre1.ftn	2009-08-05 23:35:15.000000000 +0200
+++ swanpre1.ftn	2009-08-05 21:33:22.000000000 +0200
@@ -1610,12 +1610,13 @@
 !                   <                            >  [limiter]   )     &   40.03
 !                    | NONSTat  [mxitns]        |
 !
-!           ( DIRimpl [cdd] [cdlim]                                )  &
+!           ( DIRimpl [cdd] [cdlim]  WNUMber                       )  &
 !
+!           ( REFRLim [frlim] [power]          (NOT documented)    )  &
 !
 !           | -> SIGIMpl [css] [eps2] [outp] [niter]               |
 !          (<                                                      >) &
-!           |    SIGEXpl [cfl]                 (NOT documented)    |
+!           |    SIGEXpl [css] [cfl]           (NOT documented)    |
 !           |                                                      |
 !           |    FIL     [diffc]               (NOT documented)    |
 !
@@ -1680,6 +1681,17 @@
             IREFR = 0                                                     30.80
             CALL MSGERR(0, 'Refraction deactivated')                      40.02
           ENDIF
+          CALL INKEYW ('STA','    ')                                      41.07
+          IF (KEYWIS ('WNUM')) PNUMS(32) = 1.                             41.07
+        ENDIF
+!
+!       limit Ctheta if user want so                                      41.06
+!
+        CALL INKEYW ('STA','  ')
+        IF (KEYWIS ('REFRL')) THEN
+          PNUMS(29) = 1.
+          CALL INREAL ('FRLIM', PNUMS(26) ,'UNC', 0.)
+          CALL INREAL ('POWER', PNUMS(27), 'UNC', 0.)
         ENDIF
 !
 !       *** numerical scheme in frequency space :                  ***
@@ -1698,15 +1710,16 @@
 !         *** implicit solver ***
 !         ***  This is the default option   PNUMS(8) = 1.  ***
           PNUMS(8) = 1.
-          CALL INREAL ('CSS'   , PNUMS(7), 'UNC', 0.)
-          CALL INREAL ('EPS1'  , PNUMS(11) , 'UNC', 0.)
-          CALL INREAL ('EPS2'  , PNUMS(12) , 'UNC', 0.)
-          CALL INREAL ('OUTP'  , PNUMS(13) , 'UNC', 0.)
-          CALL INREAL ('NITER' , PNUMS(14), 'UNC', 0.)
+          CALL INREAL ('CSS'  , PNUMS (7), 'UNC', 0.)
+          CALL INREAL ('EPS1' , PNUMS(11), 'UNC', 0.)
+          CALL INREAL ('EPS2' , PNUMS(12), 'UNC', 0.)
+          CALL INREAL ('OUTP' , PNUMS(13), 'UNC', 0.)
+          CALL INREAL ('NITER', PNUMS(14), 'UNC', 0.)
         ELSE IF (KEYWIS('SIGEX') .OR. KEYWIS('EXP')) THEN                 30.20
 !
 !         *** 2) explicit scheme ***
           PNUMS(8) = 2.
+          CALL INREAL ('CSS'  , PNUMS (7), 'UNC', 0.)                     41.07
           CALL INREAL ('CFL'  , PNUMS(19), 'UNC', 0.)
 !
         ELSE IF (KEYWIS ('FIL')) THEN
@@ -1951,12 +1964,15 @@
 !           |                                                                  |
 !           |    VARiable [alpha] [gammin] [gammax] [gamneg] [coeff1] [coeff2] |
 !           |                                                                  |
-! BREaking <     RUEssink [alpha] [a] [b]                                       >
+! BREaking <     RUEssink [alpha] [a] [b]                                       > &
 !           |                                                                  |
 !           |    TG       [alpha] [gamma] [pown]                               |
 !           |                                                                  |
 !           !    BIPhase  [alpha] [pown] [bref]                                |
 !
+!
+!       ( FREQDep [power] [fmin] [fmax] ) (fmin and fmax not documented)
+!
 ! ============================================================
 !
       IF (KEYWIS ('BRE')) THEN
@@ -1989,6 +2005,14 @@
           CALL INREAL ('POWN' ,  PSURF(4), 'STA', 2.4)
           CALL INREAL ('BREF' ,  PSURF(5), 'STA', -1.5708)
         ENDIF
+!
+        CALL INKEYW ('STA', '  ')                                         41.06
+        IF (KEYWIS ('FREQD')) THEN                                        41.06
+           IFRSRF = 1                                                     41.06
+           CALL INREAL ('POWER', PSURF(11), 'STA', 2.0)                   41.06
+           CALL INREAL ('FMIN' , PSURF(12), 'STA', 0.0)                   41.06
+           CALL INREAL ('FMAX' , PSURF(13), 'STA', 1000.)                 41.06
+        END IF                                                            41.06
         GOTO 100
       ENDIF
 !
--- SwanPropvelS.ftn90	2009-08-05 23:35:15.000000000 +0200
+++ SwanPropvelS.ftn90	2009-08-05 21:33:19.000000000 +0200
@@ -38,11 +38,18 @@
 !
 !   40.80: Marcel Zijlema
 !   41.02: Marcel Zijlema
+!   41.06: Gerbrant van Vledder
+!   41.07: Marcel Zijlema
 !
 !   Updates
 !
 !   40.80,     July 2007: New subroutine
 !   41.02, February 2009: adaption of velocities in case of diffraction
+!   41.06,    March 2009: add option of limitation of velocity in theta-direction
+!   41.07,   August 2009: add option of alternative formula for computation of
+!                         wave transport velocity in theta-direction based on
+!                         (x,y)-derivatives of the wave number
+!                         (see Holthuijsen (2007), page 210, footnote 4)
 !
 !   Purpose
 !
@@ -102,31 +109,43 @@
     integer                               :: k        ! loop counter
     integer, dimension(3)                 :: v        ! vertices in considered cell
     !
-    real, dimension(28)                   :: cd       ! coefficients for computing cad
-    real, dimension(14)                   :: cs       ! coefficients for computing cas
+    real, dimension(3)                    :: cd       ! coefficients for computing cad
+    real, dimension(10)                   :: cs       ! coefficients for computing cas
+    real                                  :: dhdx     ! derivative of dep2 to x in considered sweep
+    real                                  :: dhdxm    ! derivative of dep2 to x in sweep below considered sweep
+    real                                  :: dhdxp    ! derivative of dep2 to x in sweep above considered sweep
+    real                                  :: dhdy     ! derivative of dep2 to y in considered sweep
+    real                                  :: dhdym    ! derivative of dep2 to y in sweep below considered sweep
+    real                                  :: dhdyp    ! derivative of dep2 to y in sweep above considered sweep
     real, dimension(3)                    :: dloc     ! local depth at vertices
-    real                                  :: dpdx     ! derivative of dep2 to x in considered sweep
-    real                                  :: dpdy     ! derivative of dep2 to y in considered sweep
-    real                                  :: dpdxm    ! derivative of dep2 to x in sweep below considered sweep
-    real                                  :: dpdym    ! derivative of dep2 to y in sweep below considered sweep
-    real                                  :: dpdxp    ! derivative of dep2 to x in sweep above considered sweep
-    real                                  :: dpdyp    ! derivative of dep2 to y in sweep above considered sweep
     real                                  :: duxdx    ! derivative of ux2 to x in considered sweep
-    real                                  :: duxdy    ! derivative of ux2 to y in considered sweep
-    real                                  :: duydx    ! derivative of uy2 to x in considered sweep
-    real                                  :: duydy    ! derivative of uy2 to y in considered sweep
     real                                  :: duxdxm   ! derivative of ux2 to x in sweep below considered sweep
-    real                                  :: duxdym   ! derivative of ux2 to y in sweep below considered sweep
-    real                                  :: duydxm   ! derivative of uy2 to x in sweep below considered sweep
-    real                                  :: duydym   ! derivative of uy2 to y in sweep below considered sweep
     real                                  :: duxdxp   ! derivative of ux2 to x in sweep above considered sweep
+    real                                  :: duxdy    ! derivative of ux2 to y in considered sweep
+    real                                  :: duxdym   ! derivative of ux2 to y in sweep below considered sweep
     real                                  :: duxdyp   ! derivative of ux2 to y in sweep above considered sweep
+    real                                  :: duydx    ! derivative of uy2 to x in considered sweep
+    real                                  :: duydxm   ! derivative of uy2 to x in sweep below considered sweep
     real                                  :: duydxp   ! derivative of uy2 to x in sweep above considered sweep
+    real                                  :: duydy    ! derivative of uy2 to y in considered sweep
+    real                                  :: duydym   ! derivative of uy2 to y in sweep below considered sweep
     real                                  :: duydyp   ! derivative of uy2 to y in sweep above considered sweep
+    real                                  :: fac      ! a factor
+    real                                  :: frlim    ! frequency range in which limit on velocity in theta-direction is applied
     real                                  :: kd       ! help variable, wave number times water depth
+    real                                  :: pp       ! power of the frequency dependent limiter on refraction
     real, dimension(2)                    :: rdxl     ! first component of local contravariant base vector rdx(b) = a^(b)_1
     real, dimension(2)                    :: rdyl     ! second component of local contravariant base vector rdy(b) = a^(b)_2
     !
+    real, dimension(MSC)                  :: arr      ! auxiliary array
+    real, dimension(MSC)                  :: dkdx     ! derivative of wave number to x in considered sweep
+    real, dimension(MSC)                  :: dkdxm    ! derivative of wave number to x in sweep below considered sweep
+    real, dimension(MSC)                  :: dkdxp    ! derivative of wave number to x in sweep above considered sweep
+    real, dimension(MSC)                  :: dkdy     ! derivative of wave number to y in considered sweep
+    real, dimension(MSC)                  :: dkdym    ! derivative of wave number to y in sweep below considered sweep
+    real, dimension(MSC)                  :: dkdyp    ! derivative of wave number to y in sweep above considered sweep
+    real, dimension(MSC,3)                :: kloc     ! local wave number
+    !
     type(celltype), dimension(:), pointer :: cell     ! datastructure for cells with their attributes
     type(verttype), dimension(:), pointer :: vert     ! datastructure for vertices with their attributes
 !
@@ -182,8 +201,8 @@
           !
        endif
        !
-       dpdx = rdx(1) * (dloc(1)-dloc(2)) + rdx(2) * (dloc(1)-dloc(3))
-       dpdy = rdy(1) * (dloc(1)-dloc(2)) + rdy(2) * (dloc(1)-dloc(3))
+       dhdx = rdx(1) * (dloc(1)-dloc(2)) + rdx(2) * (dloc(1)-dloc(3))
+       dhdy = rdy(1) * (dloc(1)-dloc(2)) + rdy(2) * (dloc(1)-dloc(3))
        !
     endif
     !
@@ -198,14 +217,40 @@
        !
     endif
     !
-    ! compute extra derivatives of depth and current belonging to neighbouring sweeps meant for refraction
+    ! compute the derivatives of the wave number for the considered sweep
+    !
+    if ( int(PNUMS(32)) == 1 ) then
+       !
+       ! compute wave numbers for all frequencies
+       !
+       call KSCIP1 (MSC, spcsig, dloc(1), kloc(1,1), arr, arr, arr)
+       call KSCIP1 (MSC, spcsig, dloc(2), kloc(1,2), arr, arr, arr)
+       call KSCIP1 (MSC, spcsig, dloc(3), kloc(1,3), arr, arr, arr)
+       !
+       dkdx(:) = rdx(1) * (kloc(:,1)-kloc(:,2)) + rdx(2) * (kloc(:,1)-kloc(:,3))
+       dkdy(:) = rdy(1) * (kloc(:,1)-kloc(:,2)) + rdy(2) * (kloc(:,1)-kloc(:,3))
+       !
+    endif
+    !
+    ! compute extra derivatives of depth, wave number and current belonging to neighbouring sweeps meant for refraction
     !
     if ( IREFR /= 0 ) then
        !
-       dpdxm = dpdx
-       dpdym = dpdy
-       dpdxp = dpdx
-       dpdyp = dpdy
+       if ( int(PNUMS(32)) == 0 ) then
+          !
+          dhdxm = dhdx
+          dhdym = dhdy
+          dhdxp = dhdx
+          dhdyp = dhdy
+          !
+       else
+          !
+          dkdxm = dkdx
+          dkdym = dkdy
+          dkdxp = dkdx
+          dkdyp = dkdy
+          !
+       endif
        !
        if ( ICUR /= 0 ) then
           !
@@ -259,8 +304,6 @@
        !
        if ( dloc(2) <= DEPMIN .or. dloc(3) <= DEPMIN ) goto 10
        !
-       ! compute the derivatives of the depth
-       !
        if ( IREFR == -1 ) then
           !
           ! limitation of depths in upwave vertices
@@ -270,10 +313,28 @@
           !
        endif
        !
-       dpdxm = rdxl(1) * (dloc(1)-dloc(2)) + rdxl(2) * (dloc(1)-dloc(3))
-       dpdym = rdyl(1) * (dloc(1)-dloc(2)) + rdyl(2) * (dloc(1)-dloc(3))
+       if ( int(PNUMS(32)) == 0 ) then
+          !
+          ! compute the derivatives of the depth in cell below considered sweep
+          !
+          dhdxm = rdxl(1) * (dloc(1)-dloc(2)) + rdxl(2) * (dloc(1)-dloc(3))
+          dhdym = rdyl(1) * (dloc(1)-dloc(2)) + rdyl(2) * (dloc(1)-dloc(3))
+          !
+       else
+          !
+          ! compute wave numbers for all frequencies
+          !
+          call KSCIP1 (MSC, spcsig, dloc(2), kloc(1,2), arr, arr, arr)
+          call KSCIP1 (MSC, spcsig, dloc(3), kloc(1,3), arr, arr, arr)
+          !
+          ! compute the derivatives of the wave number in cell below considered sweep
+          !
+          dkdxm(:) = rdxl(1) * (kloc(:,1)-kloc(:,2)) + rdxl(2) * (kloc(:,1)-kloc(:,3))
+          dkdym(:) = rdyl(1) * (kloc(:,1)-kloc(:,2)) + rdyl(2) * (kloc(:,1)-kloc(:,3))
+          !
+       endif
        !
-       ! compute the derivatives of the ambient current
+       ! compute the derivatives of the ambient current in cell below considered sweep
        !
        if ( ICUR /= 0 ) then
           !
@@ -325,8 +386,6 @@
        !
        if ( dloc(2) <= DEPMIN .or. dloc(3) <= DEPMIN ) goto 20
        !
-       ! compute the derivatives of the depth
-       !
        if ( IREFR == -1 ) then
           !
           ! limitation of depths in upwave vertices
@@ -336,10 +395,28 @@
           !
        endif
        !
-       dpdxp = rdxl(1) * (dloc(1)-dloc(2)) + rdxl(2) * (dloc(1)-dloc(3))
-       dpdyp = rdyl(1) * (dloc(1)-dloc(2)) + rdyl(2) * (dloc(1)-dloc(3))
+       if ( int(PNUMS(32)) == 0 ) then
+          !
+          ! compute the derivatives of the depth in cell above considered sweep
+          !
+          dhdxp = rdxl(1) * (dloc(1)-dloc(2)) + rdxl(2) * (dloc(1)-dloc(3))
+          dhdyp = rdyl(1) * (dloc(1)-dloc(2)) + rdyl(2) * (dloc(1)-dloc(3))
+          !
+       else
+          !
+          ! compute wave numbers for all frequencies
+          !
+          call KSCIP1 (MSC, spcsig, dloc(2), kloc(1,2), arr, arr, arr)
+          call KSCIP1 (MSC, spcsig, dloc(3), kloc(1,3), arr, arr, arr)
+          !
+          ! compute the derivatives of the wave number in cell above considered sweep
+          !
+          dkdxp(:) = rdxl(1) * (kloc(:,1)-kloc(:,2)) + rdxl(2) * (kloc(:,1)-kloc(:,3))
+          dkdyp(:) = rdyl(1) * (kloc(:,1)-kloc(:,2)) + rdyl(2) * (kloc(:,1)-kloc(:,3))
+          !
+       endif
        !
-       ! compute the derivatives of the ambient current
+       ! compute the derivatives of the ambient current in cell above considered sweep
        !
        if ( ICUR /= 0 ) then
           !
@@ -366,12 +443,8 @@
        !
        if ( ICUR /= 0 ) then
           !
-          cs(3) = ux2(iv1) * dpdx
-          cs(4) = uy2(iv1) * dpdy
-          cs(6) = duxdx
-          cs(7) = duxdy
-          cs(8) = duydx
-          cs(9) = duydy
+          cs(3) = ux2(iv1) * dhdx
+          cs(4) = uy2(iv1) * dhdy
           !
        endif
        !
@@ -385,17 +458,17 @@
           cs(5) = -cgo(is,1) * kwave(is,1)
           if ( IDIFFR /= 0 ) cs(5) = cs(5)*DIFPARAM(iv1)
           !
-          cs(10) = cs(1) * cs(2)
-          cs(11) = cs(1) * (cs(3)+cs(4))
-          cs(12) = cs(5) * cs(6)
-          cs(13) = cs(5) * (cs(7)+cs(8))
-          cs(14) = cs(5) * cs(9)
+          cs( 6) = cs(1) * cs(2)
+          cs( 7) = cs(1) * (cs(3)+cs(4))
+          cs( 8) = cs(5) * duxdx
+          cs( 9) = cs(5) * (duxdy+duydx)
+          cs(10) = cs(5) * duydy
           !
           do iddum = idcmin(is), idcmax(is)
              id = mod ( iddum - 1 + MDC , MDC ) + 1
              !
-             cas(id,is) = cs(10)
-             if ( ICUR /= 0 ) cas(id,is) = cas(id,is) + cs(11) + coscos(id)*cs(12) + sincos(id)*cs(13) + sinsin(id)*cs(14)
+             cas(id,is) = cs(6)
+             if ( ICUR /= 0 ) cas(id,is) = cas(id,is) + cs(7) + coscos(id)*cs(8) + sincos(id)*cs(9) + sinsin(id)*cs(10)
 !            !
           enddo
        enddo
@@ -406,78 +479,72 @@
     !
     if ( IREFR /= 0 ) then
        !
-       ! compute coefficients depending on depth
-       !
-       cd( 2) = dpdx
-       cd( 3) = dpdy
-       !
-       cd(11) = dpdxm
-       cd(12) = dpdym
-       !
-       cd(20) = dpdxp
-       cd(21) = dpdyp
-       !
-       ! compute coefficients depending on ambient currents
-       !
-       if ( ICUR /= 0 ) then
-          !
-          cd( 4) = duxdx
-          cd( 5) = duydy
-          cd( 6) = duydx
-          cd( 7) = duxdy
-          !
-          cd(13) = duxdxm
-          cd(14) = duydym
-          cd(15) = duydxm
-          cd(16) = duxdym
-          !
-          cd(22) = duxdxp
-          cd(23) = duydyp
-          cd(24) = duydxp
-          cd(25) = duxdyp
-          !
-       endif
-       !
        do is = 1, MSC
           !
           ! compute frequency-dependent coefficients
           !
-          kd = min(30.,kwave(is,1) * dep2(iv1))
-          !
-          cd(1) = spcsig(is) / sinh (2.* kd)
-          !
-          cd( 8) = cd(1) * cd(2)
-          cd( 9) = cd(1) * cd(3)
-          cd(10) = cd(4) - cd(5)
+          if ( int(PNUMS(32)) == 0 ) then
+             !
+             kd = min(30.,kwave(is,1) * dep2(iv1))
+             !
+             cd(1) = spcsig(is) / sinh (2.* kd)
+             !
+             cd(2) = cd(1) * dhdx
+             cd(3) = cd(1) * dhdy
+             !
+          else
+             !
+             cd(1) = -cgo(is,1) / kwave(is,1)
+             !
+             cd(2) = cd(1) * dkdx(is)
+             cd(3) = cd(1) * dkdy(is)
+             !
+          endif
           !
           do iddum = idcmin(is), idcmax(is)
              id = mod ( iddum - 1 + MDC , MDC ) + 1
              !
-             cad(id,is) = esin(id)*cd(8) - ecos(id)*cd(9)
+             cad(id,is) = esin(id)*cd(2) - ecos(id)*cd(3)
              if ( IDIFFR /= 0 ) cad(id,is) = cad(id,is)*DIFPARAM(iv1) - DIFPARDX(iv1)*cgo(is,1)*esin(id) + DIFPARDY(iv1)*cgo(is,1)*ecos(id)
-             if ( ICUR   /= 0 ) cad(id,is) = cad(id,is) + sincos(id)*cd(10) + sinsin(id)*cd(6) - coscos(id)*cd(7)
+             if ( ICUR   /= 0 ) cad(id,is) = cad(id,is) + sincos(id)*(duxdx-duydy) + sinsin(id)*duydx - coscos(id)*duxdy
              !
           enddo
           !
-          cd(17) = cd( 1) * cd(11)
-          cd(18) = cd( 1) * cd(12)
-          cd(19) = cd(13) - cd(14)
+          if ( int(PNUMS(32)) == 0 ) then
+             !
+             cd(2) = cd(1) * dhdxm
+             cd(3) = cd(1) * dhdym
+             !
+          else
+             !
+             cd(2) = cd(1) * dkdxm(is)
+             cd(3) = cd(1) * dkdym(is)
+             !
+          endif
           !
           id = mod ( idcmin(is) - 2 + MDC , MDC ) + 1  ! this direction belongs to sweep below considered sweep
           !
-          cad(id,is) = esin(id)*cd(17) - ecos(id)*cd(18)
+          cad(id,is) = esin(id)*cd(2) - ecos(id)*cd(3)
           if ( IDIFFR /= 0 ) cad(id,is) = cad(id,is)*DIFPARAM(iv1) - DIFPARDX(iv1)*cgo(is,1)*esin(id) + DIFPARDY(iv1)*cgo(is,1)*ecos(id)
-          if ( ICUR   /= 0 ) cad(id,is) = cad(id,is) + sincos(id)*cd(19) + sinsin(id)*cd(15) - coscos(id)*cd(16)
+          if ( ICUR   /= 0 ) cad(id,is) = cad(id,is) + sincos(id)*(duxdxm-duydym) + sinsin(id)*duydxm - coscos(id)*duxdym
           !
-          cd(26) = cd( 1) * cd(20)
-          cd(27) = cd( 1) * cd(21)
-          cd(28) = cd(22) - cd(23)
+          if ( int(PNUMS(32)) == 0 ) then
+             !
+             cd(2) = cd(1) * dhdxp
+             cd(3) = cd(1) * dhdyp
+             !
+          else
+             !
+             cd(2) = cd(1) * dkdxp(is)
+             cd(3) = cd(1) * dkdyp(is)
+             !
+          endif
           !
           id = mod ( idcmax(is) + MDC , MDC ) + 1      ! this direction belongs to sweep above considered sweep
           !
-          cad(id,is) = esin(id)*cd(26) - ecos(id)*cd(27)
+          cad(id,is) = esin(id)*cd(2) - ecos(id)*cd(3)
           if ( IDIFFR /= 0 ) cad(id,is) = cad(id,is)*DIFPARAM(iv1) - DIFPARDX(iv1)*cgo(is,1)*esin(id) + DIFPARDY(iv1)*cgo(is,1)*ecos(id)
-          if ( ICUR   /= 0 ) cad(id,is) = cad(id,is) + sincos(id)*cd(28) + sinsin(id)*cd(24) - coscos(id)*cd(25)
+          if ( ICUR   /= 0 ) cad(id,is) = cad(id,is) + sincos(id)*(duxdxp-duydyp) + sinsin(id)*duydxp - coscos(id)*duxdyp
           !
        enddo
        !
@@ -496,6 +563,24 @@
           !
        endif
        !
+       ! limit velocity in some frequency range if requested
+       !
+       if ( int(PNUMS(29)) == 1 ) then
+          !
+          frlim = PI2*PNUMS(26)
+          pp    =     PNUMS(27)
+          !
+          do is = 1, MSC
+             !
+             fac = min(1.,(spcsig(is)/frlim)**pp)
+             !
+             do id = 1, MDC
+                cad(id,is) = fac*cad(id,is)
+             enddo
+          enddo
+          !
+       endif
+       !
     endif
     !
 end subroutine SwanPropvelS
--- SwanReadADCGrid.ftn90	2009-08-05 23:35:15.000000000 +0200
+++ SwanReadADCGrid.ftn90	2009-08-05 21:33:19.000000000 +0200
@@ -32,10 +32,12 @@
 !   Authors
 !
 !   40.80: Marcel Zijlema
+!   41.07: Casey Dietrich
 !
 !   Updates
 !
 !   40.80, December 2007: New subroutine
+!   41.07,   August 2009: use ADCIRC boundary info to mark all boundary vertices
 !
 !   Purpose
 !
@@ -63,12 +65,15 @@
     integer                 :: ii       ! auxiliary integer
     integer                 :: iostat   ! I/O status in call FOR
     integer                 :: istat    ! indicate status of allocation
+    integer                 :: itype    ! ADCIRC boundary type
     integer                 :: ivert    ! vertex index
+    integer                 :: ivert1   ! another vertex index
     integer                 :: j        ! loop counter
     integer                 :: k        ! loop counter
     integer                 :: n1       ! auxiliary integer
     integer                 :: n2       ! another auxiliary integer
     integer                 :: ndsd     ! unit reference number of file
+    integer                 :: nopbc    ! number of open boundaries in ADCIRC
     integer                 :: vm       ! boundary marker
     character(80)           :: line     ! auxiliary textline
     logical                 :: stpnow   ! indicate whether program must be terminated or not
@@ -127,36 +132,51 @@
        if ( ii/=j ) call msgerr ( 1, 'numbering of triangles is not sequential in grid file fort.14 ' )
     enddo
     !
-    ! skip part containing ADCIRC boundary information (not relevant to SWAN)
+    if(.not.allocated(vmark)) allocate (vmark(nverts), stat = istat)
+    if ( istat /= 0 ) then
+       call msgerr ( 4, 'Allocation problem in SwanReadADCGrid: array vmark ' )
+       goto 900
+    endif
+    vmark = 0
     !
-    read(ndsd, *, end=950, err=910) n1
+    ! read ADCIRC boundary information and store boundary markers
+    !
+    read(ndsd, *, end=950, err=910) nopbc
     read(ndsd, *, end=950, err=910) idum
-    do j = 1, n1
+    do j = 1, nopbc
+       vm = j
        read(ndsd, *, end=950, err=910) n2
        do k = 1, n2
-          read(ndsd, *, end=950, err=910) idum
+           read(ndsd, *, end=950, err=910) ivert
+           vmark(ivert) = vm
        enddo
     enddo
     !
     read(ndsd, *, end=950, err=910) n1
     read(ndsd, *, end=950, err=910) idum
     do j = 1, n1
-       read(ndsd, *, end=950, err=910) n2, idum
-       do k = 1, n2
-          read(ndsd, *, end=950, err=910) idum
-       enddo
+       vm = nopbc + j
+       read(ndsd, *, end=950, err=910) n2, itype
+       if ( itype /= 4 .and. itype /= 24 ) then
+          do k = 1, n2
+             read(ndsd, *, end=950, err=910) ivert
+             vmark(ivert) = vm
+          enddo
+       else
+          do k = 1, n2
+             read(ndsd, *, end=950, err=910) ivert, ivert1
+             vmark(ivert ) = vm
+             vmark(ivert1) = vm
+          enddo
+       endif
     enddo
     !
-    if(.not.allocated(vmark)) allocate (vmark(nverts), stat = istat)
-    if ( istat /= 0 ) then
-       call msgerr ( 4, 'Allocation problem in SwanReadADCGrid: array vmark ' )
-       goto 900
-    endif
-    vmark = 0
+    ! read and store SWAN boundary markers, if appropriate
     !
-    ! read and store boundary markers
+    read(ndsd, *, end=920, err=910) n1
+    !
+    vmark = 0
     !
-    read(ndsd, *, end=950, err=910) n1
     read(ndsd, *, end=950, err=910) idum
     do j = 1, n1
        read(ndsd, *, end=950, err=910) n2, vm
@@ -168,7 +188,7 @@
     !
     ! close file fort.14
     !
-    close(ndsd)
+ 920 close(ndsd)
     !
  900 return
     !
--- SwanSweepSel.ftn90	2009-08-05 23:35:15.000000000 +0200
+++ SwanSweepSel.ftn90	2009-08-05 21:33:20.000000000 +0200
@@ -309,10 +309,10 @@
           return
        endif
        !
-       if ( isslow /= 1 ) then
-          call msgerr ( 4, 'inconsistency found in SwanSweepSel: isslow <> 1 ' )
-          return
-       endif
+       !if ( isslow /= 1 ) then
+       !   call msgerr ( 4, 'inconsistency found in SwanSweepSel: isslow <> 1 ' )
+       !   return
+       !endif
        !
        isslow = 1
        !
--- swantech.tex	2009-08-05 23:35:15.000000000 +0200
+++ swantech.tex	2009-08-05 21:55:52.000000000 +0200
@@ -17,7 +17,7 @@
 \end{center}
 \vfill
 \begin{center}
-{\Large\bf SWAN Cycle III version 40.72AB}
+{\Large\bf SWAN Cycle III version 40.72ABC}
 \end{center}
 
 \cleardoublepage
@@ -216,7 +216,7 @@
 \noindent
 We further want to acknowledge all contributors who helped us to improve SWAN, reported bugs, and tested SWAN
 thoroughly: Tim Campbell, John Cazes, IJsbrand Haagsma, Agnieszka Herman, Jim Kaihatu, Kees Kassels, Annette Kieftenburg,
-Ekaterini Kriezi, Roberto Padilla-Hernandez, Erick Rogers, Kees Vuik, Andre van der Westhuijsen and Marcel Zijlema.
+Ekaterini Kriezi, Roberto Padilla-Hernandez, Erick Rogers, Kees Vuik, Andre van der Westhuysen and Marcel Zijlema.
 \\[2ex]
 \noindent
 Many thanks are due to Gerbrant van Vledder and Noriaki Hashimoto who provided the source code for exact
@@ -228,8 +228,8 @@
 Rijkswaterstaat (as part of the Ministry of Transport, Public Works and Water Management, The Netherlands).
 \\[2ex]
 \noindent
-We are finally grateful to all those other people working on the Public Domain Software without which a project
-like SWAN would be unthinkable: Linux, Intel, GNU F95, \LaTeX, MPICH and many others.
+We are finally grateful to all those other people working on the Public Domain Software without which the development
+of SWAN would be unthinkable: Linux, Intel, GNU F95, \LaTeX, MPICH2, Perl and many others.
 
 \chap{Governing equations} \label{ch:goveq}
 
@@ -669,19 +669,14 @@
   \label{eq2-8}
 \end{equation}
 in which $E_{\rm tot}$ is the total wave energy and $D_{\rm tot}~<~0$ is the rate of dissipation of
-the total energy due to wave breaking according to Battjes and Janssen (1978). Adding a quadratic dependency
-on frequency as suggested by Mase and Kirby (1992) and supported by Elgar~{\it et~al}. (1997) seems to have no
-noticeable effect on the SWAN results. Chen and Guza (1997) inferred from observations and simulations with a
-Boussinesq model that the high-frequency levels are insensitive to such frequency dependency because
-an increased dissipation at high frequencies is compensated approximately by increased nonlinear
-energy transfer (but they did find the frequency dependency to be relevant in time domain). The value of
-$D_{\rm tot}$ depends critically on the breaker parameter $\gamma = H_{\max}/d$
+the total energy due to wave breaking according to Battjes and Janssen (1978).
+The value of $D_{\rm tot}$ depends critically on the breaker parameter $\gamma = H_{\max}/d$
 (in which $H_{\max}$ is the maximum possible individual wave height in the local water depth $d$).
 In SWAN, both a constant value and a variable value are available. Examples of a variable breaker parameter
 can be found in Nelson (1987) and Ruessink et al. (2003). (Both are implemented in SWAN.)
 The constant value is $\gamma=0.73$ found as
 the mean value of the data set of Battjes and Stive (1985).
-\nocite{Bat78J,Tho83G,Bat92B,Vin94SD,Arc94RC,Eld96B,Eld95B,Mas92K,Elg97GRHG,Che97G,Nel87,Rue03WS,Bat85S}
+\nocite{Bat78J,Tho83G,Bat92B,Vin94SD,Arc94RC,Eld96B,Eld95B,Mas92K,Elg97GRHG,Che97GE,Nel87,Rue03WS,Bat85S}
 \\[2ex]
 \noindent
 \underline{Nonlinear wave-wave interactions}\\[2ex]
@@ -997,10 +992,11 @@
   U^2_{\rm rms} = \int_{0}^{2\pi} \int_{0}^{\infty} \frac{\sigma^2}{g^2 \sinh^2 kd} E(\sigma,\theta) d\sigma d\theta
   \label{eq3-18}
 \end{equation}
-Hasselmann~{\it et~al}. (1973) found from the results of the JONSWAP experiment
-$C_{\rm b}~=~C_{\rm JON}~=~0.038$m$^{2}$s$^{-3}$ for swell conditions. Bouws and Komen (1983) selected a bottom
-friction coefficient of $C_{\rm JON}~=~0.067$m$^{2}$s$^{-3}$ for
-fully developed wave conditions in shallow water. Both values are available in SWAN.
+Hasselmann~{\it et~al}. (1973) found $C_{\rm b}=C_{\rm JON}=0.038$m$^{2}$s$^{-3}$ which is in
+agreement with the JONSWAP result for swell dissipation. However, Bouws and Komen (1983) suggest a
+value of $C_{\rm JON}~=~0.067$m$^{2}$s$^{-3}$ for depth-limited wind-sea conditions in the North Sea. This value is derived
+from revisiting the energy balance equation employing an alternative deep water dissipation.
+Both values are available in SWAN.
 \\[2ex]
 \noindent
 The expression of Collins (1972) is based on a conventional formulation for periodic waves with the
@@ -1729,7 +1725,7 @@
 of a much smaller (pseudo) time step (Ferziger and Peri\'{c}, 1999).
 Consequently, a limiter may no longer be needed. Although this
 approach may be suitable to SWAN, it slows down convergence significantly.
-In this paper, we propose a new method that finds a
+Here, we propose a new method that finds a
 compromise between fast convergence on the one hand and minimizing the role of the limiter in the
 energetic part of the
 spectrum on the other. The key idea to achieve this is to link the extent of updating to the wave
@@ -1739,7 +1735,7 @@
 \nocite{Fer99P}
 \\[2ex]
 \noindent
-The second objective of this paper concerns the formulation and the use of termination criteria required by
+The second objective concerns the formulation and the use of termination criteria required by
 the iteration procedure in SWAN.
 In principle, the iterative process should be stopped if the convergence error defined as the
 difference between the current iterate and the stationary solution is smaller
@@ -1786,16 +1782,22 @@
 \end{eqnarray}
 where $n$ is a time-level with $\Delta t$ a time step. Note that locations in between consecutive
 counters are reflected with the half-indices.
-\\[2ex]
-\noindent
+
+\subsection{Discretization in geographical space}
+
 Since, the unknown $N$ and the propagation velocities are only given in points $(i,j,l,m)$,
-further approximation is needed. In the present paper, we employ a first order upwind scheme
-in geographical space, since it is sufficient accurate for nearshore applications and fully
-monotone, i.e. it can not to give rise to spurious oscillations. It should be noted, however,
-that in applications at oceanic scales, a higher order upwind scheme should be employed. In
+further approximation is needed. A first order upwind scheme
+in geographical space may be employed, since it fully
+monotone, i.e. it can not to give rise to spurious oscillations.
+A disadvantage of this scheme is that it is numerically
+diffusive, which naturally degrades the accuracy of the model.
+This numerical diffusion is caused by gradients of wave action across geographic space, e.g. due to refraction by bathymetry or currents.
+However, in
 the current SWAN version, two alternatives to this scheme are implemented, namely the second
-order SORDUP and the third order Stelling/Leendertse schemes.
-See also Rogers~{\it et~al}. (2002) and Stelling and Leendertse (1992).
+order SORDUP and the third order Stelling/Leendertse schemes. These schemes produce far less numerical diffusion.
+\\[2ex]
+\noindent
+\underline {First order upwind scheme; BSBT}
 \\[2ex]
 \noindent
 The fluxes $c_x N$ at $(i+1/2,j,l,m)$ and $c_y N$ at $(i,j+1/2,l,m)$ are approximated in the
@@ -1837,6 +1839,7 @@
 \begin{equation}
   \left( \frac{1.5 (c_y N)_{i_y} - 2 (c_y N)_{i_y-1} + 0.5 (c_y N)_{i_y-2}}{\Delta y} \right)^{i_t,n}_{i_x, i_{\sigma}, i_{\theta}}
 \end{equation}
+See also Rogers~{\it et~al}. (2002).
 In the neighboorhood of open boundaries, land boundaries and obstacles (i.e., the last two grids adjoining such grid points for
 the SORDUP scheme), SWAN will revert to the first order upwind BSBT scheme. This scheme has a larger numerical diffusion but that
 is usually acceptable over the small distances involved.
@@ -1856,6 +1859,7 @@
   \left( \frac{\frac{5}{6} (c_y N)_{i_y} - \frac{5}{4} (c_y N)_{i_y-1} + \frac{1}{2} (c_y N)_{i_y-2} \frac{1}{12} (c_y N)_{i_y-3}}{\Delta y} \right)^{i_t,n}_{i_x, i_{\sigma}, i_{\theta}}
   + \left( \frac{(c_y N)_{i_y+1} - (c_y N)_{i_y-1}}{4 \Delta y} \right)^{i_t-1}_{i_x, i_{\sigma}, i_{\theta}}
 \end{equation}
+See also Stelling and Leendertse (1992).
 In the neighboorhood of open boundaries, land boundaries and obstacles (i.e., the last three grids adjoining such grid points for
 the Stelling and Leendertse scheme), SWAN will revert to the first order upwind BSBT scheme.
 \\[2ex]
@@ -1923,8 +1927,30 @@
 effect tends to be small on these scales and the diffusion can and should not be used to avoid the stability problem. For small-scale (local)
 applications, typically $\mu = {\cal O} (10-100)$. But such cases are usually treated as stationary and the SORDUP scheme should be used (no GSE
 correction is included in this scheme).
-\\[2ex]
-\noindent
+
+\subsection{Note on the choice of geographic propagation schemes}
+
+The main interest of the SWAN users is in simulating wind-generated waves and combined swell-sea cases in {\em coastal ocean waters},
+and it is particularly with the view to
+such computations that a simple and compact, but first order, scheme was implemented in SWAN.
+A substantial body of experience gathered over the past 15 years on the performance of both lower and higher upwind schemes in SWAN suggests that in many circumstances,
+the discretization of the propagation terms in geographical space is not a crucial issue. Many nearshore simulations have shown the solution for action density to
+be on the whole rather insensitive to the accuracy with which geographic propagation terms are approximated.
+This reflects the tendency for the level of wave action to be dictated by source terms, while the local changes of
+the energy field across geographical space is relatively weak.
+This is consonant with the established view that a certain amount of numerical diffusion can be safely tolerated in the numerical scheme for geographic propagation,
+as its impact on wave parameters is negligible
+(Rogers~{\it et~al}., 2002; WISE Group, 2007).
+Also, broad wave spectra will tend make numerical diffusion far less noticable in a wave field.
+\\[2ex]
+\noindent
+This would appear to suggest, however, that the use of higher order upwind schemes serves no useful purpose. This is probably not so since there might be some cases
+that are prone to diffusion, where the benefit of such schemes is obvious. One can think of a case of swell propagation over very long distances. While low-diffusive,
+higher order schemes did permit long-distance swell cases to be validated, the reduced diffusion was found to pose a serious difficulty as the garden sprinkler
+effect becomes more visible, see e.g. WISE Group (2007).
+
+\subsection{Discretization in spectral space}
+
 The fluxes in the spectral space $(\sigma,\theta)$, as given in (\ref{eq:actdisc}), should not be approximated
 with the first order upwind scheme since, it turns out to be very diffusive for frequencies near the blocking
 frequency\footnote{Waves can be blocked by the current at a relative high frequency.}. Central differences
@@ -1959,6 +1985,18 @@
 where the parameters $\mu$ and $\nu$ are still to be chosen. For all values $\mu \in [0,1]$ and $\nu \in [0,1]$,
 a blended form arises between first order upwind differencing ($\mu = \nu = 0$) and central differencing
 ($\mu = \nu = 1$).
+\\[2ex]
+\noindent
+If large gradients in action density in the directional space are present, numerical oscillations can
+arise. Negative action density (which is unphysical) is removed from the two-dimensional spectrum. To
+retain a conservative numerical scheme, first the total action density for a certain frequency in a
+quadrant is computed, as follows:
+\begin{equation}
+  N(\sigma) = \int_{\mbox{quadrant}} N(\sigma,\theta)\, d\theta
+\end{equation}
+then all negative action density is removed and the positive action density is multiplied by a factor
+determined by the total action density divided by the absolute value of the total negative action
+density.
 
 \section{Solution algorithm} \label{sec:sol}
 
@@ -2020,7 +2058,8 @@
 One recognizes that the subblocks on the main diagonal express coupling among the unknowns in
 the $(\sigma,\theta)-$space for each geographic grid point, whereas the off-diagonal subblocks
 represent coupling across geographical grid points. This system can be solved with a
-Gauss-Seidel technique in one step (Wesseling, 1992). Generally,
+Gauss-Seidel technique in one step (Wesseling, 1992). For an illustrative explanation of
+this technique, see Section~\ref{sec:simple}. Generally,
 the velocities $c_x$ and $c_y$ may have different
 signs in the geographical domain and hence, more steps are needed. However, it is well known
 that adapting the ordering of updates of the unknowns $N$ in geographical space to the propagation
@@ -2174,6 +2213,115 @@
 terms of Holthuijsen and De Boer (1988).
 \nocite{Hol88B}
 
+\section{An illustrative explanation of the sweeping approach}
+\label{sec:simple}
+
+In the absence of a current, the direction of propagation of the wave crest is equal to that of the wave
+energy. For this case, the propagation velocity of energy ($c_x$, $c_y$) is equal
+to the group velocity ($c_{g,x}$, $c_{g,y}$). In presence of a current this is not the case, since the propagation
+velocities $c_x$ and $c_y$ of energy are changed by the current. Considering the applied numerical procedure
+in SWAN, it is initially more convenient to explain the basic principles of the numerical procedure in
+the absence of a current than in the situation where a current is present. So, first, we shall focus on
+the sweeping technique in the absence of a current. After this, we shall discuss the numerical
+procedure in case a current is present.
+\\[2ex]
+\noindent
+The computational region is a rectangle covered with a rectangular grid. One of the axes (say
+the $x-$axis) is chosen arbitrary, for instance perpendicular to the coast.
+The state in a grid point ($x_i$,$y_j$) in an upwind stencil is determined by its up-wave grid points
+($x_{i-1}$,$y_j$) and ($x_i$,$y_{j-1}$). This stencil covers
+the propagation of action density within a sector of 0$^o -$90$^o$, in the entire geographic space; see Figure~\ref{fig:fsweep}.
+\begin{figure}[htb]
+   \centerline{
+      \epsfig{file=fsweep.eps,height=10cm}
+              }
+      \caption{Numerical scheme for wave propagation in geographic space with below the first quadrant for which the waves are propagated.}
+      \label{fig:fsweep}
+\end{figure}
+Hence, this procedure is called
+sweep 1 and encloses all wave energy propagation over the first quadrant. By rotating the stencil over
+90$^o$, the next quadrant 90$^o -$180$^o$ is propagated. Rotating the stencil twice more ensures propagation over
+all four quadrants (see Figure~\ref{fig:4sweep}). This allows waves to propagate from all directions. Hence, the method
+is characterized as a four-sweep technique.
+\\[2ex]
+\noindent
+The gain of such a stencil is that the propagation is unconditionally stable because the wave
+characteristics lie within the concerning quadrant. Propagation is not subjected to a CFL criterion.
+In cases with current refraction or bottom refraction, action density can shift from one quadrant to
+another. This is taken into account in the model by repeating the computations with converging results
+(iterative four-sweep technique). Typically, we choose a change of less then 1\% or so in significant wave
+height and mean wave period in all geographic grid points to terminate the iteration (see Section~\ref{sec:sol}).
+\\[2ex]
+\noindent
+The numerical procedure as described above remains in principle the same when a current ($U_x$,$U_y$) is present. The
+main difference is that the propagation velocities of energy are no longer equal to the group velocity
+of the waves but become equal to $c_x = c_{g,x}+U_x$ and $c_y = c_{g,y}+U_y$. To ensure an unconditionally
+stable propagation of action in geographical space in the presence of any current, it is first determined
+which spectral wave components of the spectrum can be propagated in one sweep. This implies that all wave
+components with $c_x>0$ and $c_y>0$ are propagated in the first sweep, components with $c_x<0$ and $c_y>0$
+in the second sweep, components with $c_x<0$ and $c_y<0$ in the third sweep, and finally, components $c_x>0$
+and $c_y<0$ in the fourth sweep.
+Since the group velocity of the waves decreases with increasing frequency, the higher frequencies
+are more influenced by the current. As a result, the sector boundaries in directional space for these
+higher frequencies change more compared to the sector boundaries for the lower frequencies. In general,
+four possible configurations do occur (see Figure~\ref{fig:sweepcur}). Consider, for instance, one fixed frequency
+propagating on a uniform current. The current propagates at an angle of 45$^o$ with the $x-$axis. The sign
+of the current vector and strength of the current are arbitrary. The shaded sectors in Figure~\ref{fig:sweepcur}
+indicate that all the wave components that are propagating in the direction within the shaded sector, are
+propagated in the first sweep ($c_x>0$, $c_y>0$).
+\begin{figure}[htb]
+   \centerline{
+      \epsfig{file=sweepcur.eps,height=10cm}
+              }
+      \caption{Four possible configurations and propagation velocities $c_x$, $c_y$ for a fixed frequency in the presence of
+               a current propagating at an angle of 45$^o$ with the $x-$axis.}
+      \label{fig:sweepcur}
+\end{figure}
+\\[2ex]
+\noindent
+The top-left panel (A) represents a situation in which both $c_x$ and $c_y$ are negative due to a strong
+opposing current, i.e. wave blocking occurs. None of the wave components is propagated within the first
+sweep. The top-right panel (B) represents a situation in which the current velocity is rather small.
+The sector boundaries in directional space are hardly changed by the current such that the sector
+boundaries are approximately the same as in the absence of a current. The bottom-left panel (C)
+reflects a following current that causes the propagation velocities of the wave components in two
+sectors to be larger than zero. In this specific case, all the waves of the shaded sectors are propagated
+within the first sweep. The bottom-right panel (D) represents a case with a strong following current
+for which all the action is take along with the current. For this case the fully 360$^o$ sector is
+propagated in the first sweep.
+\\[2ex]
+\noindent
+After it has been determined which wave components are propagated in one
+sweep, i.e., the sector boundaries in directional space have been determined for each frequency, the
+integration in frequency and directional space can be carried out for those wave components.
+
+\section{Implementation of DIA within the four-sweep technique}
+
+In SWAN, the quadruplets are integrated using the DIA, see Section~\ref{sec:quad}. As a
+consequence of the four-sweep technique, two different types of methods can be used to calculate the
+four wave-wave interactions.
+\begin{enumerate}
+\item
+The first method implies that the interactions are calculated in every
+iteration prior to the first sweep. This method ensures conservation of energy density, but has the
+disadvantage that the spectral source term $S_{\rm nl4}(\sigma,\theta)$ for every grid point in
+geographical space has to be stored in internal memory. Such an integration method increases the
+amount of required memory with a factor about 2.
+The source term is stored in memory and is then explicitly integrated for a particular sweep.
+\item
+The second method is slightly different, in which the interactions are
+calculated and integrated for every sweep separately. The calculation of the energy transfer for a specific
+quadrant requires also the calculation of the transfer within a sector of about 33$^o$ in the adjacent two
+quadrants. Calculating the interactions for one sweep increases the
+computation time for the quadruplets with a factor of about 1.66 (2 adjacent sectors of 33$^o$ times
+4 sweeps divided by 360$^o$). Contrary to the first method, the total rate of this energy shift is
+not stored in memory such that energy density is not conserved per sweep (and per iteration). However,
+this does not influence the converging results.
+Within this second method, two different numerical schemes are available, namely, a semi-implicit scheme
+and a fully explicit scheme. The use of a fully explicit scheme is recommended because of the computational
+efficiency. A semi-implicit scheme increases the computation time of the quadruplets with a factor of about 2.
+\end{enumerate}
+
 \section{Convergence-enhancing measures} \label{sec:stab}
 
 As explained in Section~\ref{sec:intnum},
@@ -2313,7 +2461,7 @@
 This
 behaviour is problematic when any form of stricter stopping criterion is developed
 based on $T_{m01}$. Therefore, in the improved termination criterion
-proposed in this paper, $T_{m01}$ has been abandoned as a convergence
+proposed, $T_{m01}$ has been abandoned as a convergence
 measure and only $H_{m0}$, which displays more monotonic behaviour near
 convergence, is retained.
 \\[2ex]
@@ -3321,7 +3469,52 @@
 \nocite{Bot85E}
 
 \chap{Iterative solvers} \label{ch:solver}
-This chapter is under preparation.
+
+\section{Strongly Implicit Procedure (SIP)}
+
+We want to solve the following linear system of equations
+\begin{equation}
+  A \, \vec{N} = \vec{b}
+\end{equation}
+where $A$ is some non-symmetric penta-diagonal matrix, $\vec{N}$ is the wave action vector to be solved and $\vec{b}$
+contains source terms and boundary values.
+
+The basis for the SIP method (Stone, 1968; Ferziger and Peri\'{c}, 1999) lies in the observation that an LU decomposition is an excellent
+general purpose solver, which unfortunately cannot take advantage of the sparseness of a
+matrix. Secondly, in an iterative method, if the matrix $M = LU$ is a good approximation to the
+matrix $A$, rapid convergence results. These observations lead to the idea
+of using an approximate LU factorization of $A$ as the iteration matrix $M$, i.e.:
+\begin{equation}
+   M = L\,U = A + K
+  \label{eq:lusip}
+\end{equation}
+where $L$ and $U$ are both sparse and $K$ is small. For non-symmetric matrices the incomplete LU
+(ILU) factorisation gives such an decomposition but unfortunately converges rather slowly. In
+the ILU method one proceeds as in a standard LU decomposition. However, for every element
+of the original matrix $A$ that is zero the corresponding elements in $L$ or $U$ is set to zero. This
+means that the product of $LU$ will contain more nonzero diagonals than the original matrix $A$.
+Therefore the matrix $K$ must contain these extra diagonals as well if Eq. (\ref{eq:lusip}) is to hold.
+
+Stone reasoned that if the equations approximate an elliptic partial differential equation the
+solution can be expected to be smooth. This means that the unknowns corresponding to
+the extra diagonals can be approximated by interpolation of the surrounding points. By
+allowing $K$ to have more non zero entries on all seven diagonals and using the interpolation
+mentioned above the SIP method constructs an LU factorization with the property that for a
+given approximate solution $\phi$ the product $K\phi \approx 0$ and thus the iteration matrix $M$ is close to
+$A$ by relation (\ref{eq:lusip}).
+
+To solve the system of equations the following iterations is performed,
+starting with an initial guess for the wave action vector ${\vec{N}}^0$ an iteration is performed solving:
+\begin{equation}
+   U\,{\vec{N}}^{s+1}  = L^{-1}\,K\,{\vec{N}}^s + L^{-1}\,\vec{b}
+\end{equation}
+Since the matrix $U$ is upper triangular this equation is efficiently solved by back substitution.
+An essential property which makes the method feasible is that the matrix $L$ is easily
+invertible. This iterative process is repeated $s=0,1,2,...$ until convergence is reached.
+
+\section{Successive Over Relaxation (SOR) technique}
+This section is under preparation.
+See also Botta and Ellenbroek (1985).
 
 \chap{Parallel implementation aspects} \label{ch:parall}
 
@@ -3560,7 +3753,6 @@
 a mesh generation package BatTri (Bilgili and Smith, 2003) is used.
 This grid generator is a public-domain, graphical Matlab interface to Triangle (Shewchuk, 1996).
 Triangle is a freely-distributed, two-dimensional Delaunay triangulator.
-Generation of the meshes presented in this paper was accomplished using BatTri.
 \nocite{Bil03S,She96}
 
 An important key ingredient for the preparation of the grid for the model domain is bathymetry data.
@@ -3638,27 +3830,76 @@
                in spectral space for which the waves are propagated.}
       \label{fig:gsunstruc}
 \end{figure}
-Considering a triangle $\triangle$123 where the faces towards vertex 1 are given by
+We want to find an approximation for the propagation term of Eq. (\ref{eq:waveeq2}). To this end, we employ some vector calculus.
+We consider a triangular cell as depicted in Figure~\ref{fig:celldisc}.
+\begin{figure}[htb]
+   \centerline{
+      \epsfig{file=celldisc2.eps,height=8cm}
+              }
+      \caption{A triangular cell with geometrical quantities used for discretization in geographical space. Definitions of these quantities are provided in the text.}
+      \label{fig:celldisc}
+\end{figure}
+In vertex 1, we apply a mapping from a local coordinate system $\vec{\xi}$ = ($\xi$,$\eta$) to the
+Cartesian one $\vec{x}$ = ($x$,$y$).
+Based on this transformation $\vec{x}(\vec{\xi})$, we have the following base vectors that are tangential to the coordinate lines $\xi$ and $\eta$, respectively,
 \begin{equation}
-  {\vec{e}}_{(1)}  = {\vec{x}}_1 - {\vec{x}}_2 \, , \quad {\vec{e}}_{(2)}  = {\vec{x}}_1 - {\vec{x}}_3
+  {\vec{e}}_{(1)}  = \frac{\partial \vec{x}}{\partial \xi} \, , \quad {\vec{e}}_{(2)}  = \frac{\partial \vec{x}}{\partial \eta} \, .
 \end{equation}
-with ${\vec{x}}_i = (x_i, y_i)$ the position vector of vertex $i$ in a Cartesian coordinate system.
-Next, the action densities at vertices 1, 2 and 3 are denoted by $N_1$, $N_2$ and $N_3$, respectively.
-With the help of vector analysis and after some algebra, the propagation term of Eq. (\ref{eq:waveeq}) may be approximated as follows
+The vectors
 \begin{equation}
-  \nabla_{\vec{x}} \cdot [\vec{c}_{\vec{x}} N] \approx c_x N|_2^1 e^{(1)}_1 + c_x N|_3^1 e^{(2)}_1 + c_y N|_2^1 e^{(1)}_2 + c_y N|_3^1 e^{(2)}_2
-  \label{eq:spacedisc}
+  {\vec{e}}^{(1)}  = \mbox{grad} \, \xi \, , \quad {\vec{e}}^{(2)}  = \mbox{grad} \, \eta
 \end{equation}
-where $c_x$ and $c_y$ are the $x-$ and $y-$components of the wave propagation vector $\vec{c}_{\vec{x}}$, respectively, and
+are normal to the coordinate surface of constant $\xi$ and $\eta$, respectively (see Figure~\ref{fig:celldisc}). Moreover, they are reciprocal to the base vectors, i.e.
+\begin{equation}
+  {\vec{e}}_{(\alpha)} \cdot {\vec{e}}^{(\beta)} = \delta_{\alpha}^{\beta} \, , \quad \alpha, \beta = \{1,2\} \, ,
+\end{equation}
+where $\delta_{\alpha}^{\beta}$ is Kronecker delta (which is unity if $\alpha$ = $\beta$, and zero otherwise). Using Cramer's rule, one can find
 \begin{equation}
   \vec{e}^{(1)} = \frac{1}{D} ( e^2_{(2)},-e^1_{(2)} )^{\top}\, ,\, \, \vec{e}^{(2)} = \frac{1}{D} (-e^2_{(1)}, e^1_{(1)} )^{\top}\, , \, \,
   D = e^2_{(2)} e^1_{(1)} - e^2_{(1)} e^1_{(2)} \, .
+  \label{eq:contravar}
 \end{equation}
+
+Next, we expand the propagation term of Eq. (\ref{eq:waveeq2}):
+\begin{equation}
+  \nabla_{\vec{x}} \cdot [\vec{c}_{\vec{x}} N] = \frac{\partial c_x N}{\partial x} + \frac{\partial c_y N}{\partial y} \, ,
+\end{equation}
+where $c_x$ and $c_y$ are the $x-$ and $y-$components of the wave propagation vector $\vec{c}_{\vec{x}}$, respectively.
+Using the chain rule, we obtain
+\begin{equation}
+  \nabla_{\vec{x}} \cdot [\vec{c}_{\vec{x}} N] = e^{(1)}_1 \frac{\partial c_x N}{\partial \xi} + e^{(2)}_1 \frac{\partial c_x N}{\partial \eta}
+  + e^{(1)}_2 \frac{\partial c_y N}{\partial \xi} + e^{(2)}_2 \frac{\partial c_y N}{\partial \eta} \, .
+  \label{eq:propterm}
+\end{equation}
+Further, we approximate the derivatives in Eq. (\ref{eq:propterm}). The most simplest one is a one-sided first order difference scheme, as follows
+\begin{eqnarray}
+  && \frac{\partial c_x N}{\partial \xi} \approx \frac{c_x N_1 - c_x N_2 }{\Delta \xi} \, , \quad
+  \frac{\partial c_x N}{\partial \eta} \approx \frac{c_x N_1 - c_x N_3 }{\Delta \eta} \, , \nonumber \\
+  && \nonumber \\
+  &&\frac{\partial c_y N}{\partial \xi} \approx \frac{c_y N_1 - c_y N_2 }{\Delta \xi} \, , \quad
+  \frac{\partial c_y N}{\partial \eta} \approx \frac{c_y N_1 - c_y N_3 }{\Delta \eta} \, ,
+  \label{eq:onedif}
+\end{eqnarray}
+where the action densities at vertices 1, 2 and 3 are denoted by $N_1$, $N_2$ and $N_3$, respectively.
+Here, we choose the mapping $\vec{x}(\vec{\xi})$ such that $\Delta \xi$ = $\Delta \eta$ = 1.
+The approximation is completed by substituting (\ref{eq:onedif}) in (\ref{eq:propterm}):
+\begin{equation}
+  \nabla_{\vec{x}} \cdot [\vec{c}_{\vec{x}} N] \approx c_x N|_2^1 e^{(1)}_1 + c_x N|_3^1 e^{(2)}_1 + c_y N|_2^1 e^{(1)}_2 + c_y N|_3^1 e^{(2)}_2 \, .
+  \label{eq:spacedisc}
+\end{equation}
+Note that the components of the vectors ${\vec{e}}^{(1)}$ and ${\vec{e}}^{(2)}$ in Eq. (\ref{eq:spacedisc}) are given by Eqs. (\ref{eq:contravar}), while
+the base vectors are calculated according to
+\begin{equation}
+  {\vec{e}}_{(1)}  = {\vec{x}}_1 - {\vec{x}}_2 \, , \quad {\vec{e}}_{(2)}  = {\vec{x}}_1 - {\vec{x}}_3
+  \label{eq:basevec}
+\end{equation}
+with ${\vec{x}}_i = (x_i, y_i)$ the position vector of vertex $i$ in a Cartesian coordinate system.
+
 This space discretization is of lowest order accurate and conserves action. The upwind difference scheme (\ref{eq:spacedisc}) is employed for two reasons.
 First, it enforces the propagation of wave action to follow the characteristics.
 Second, it is monotone (i.e. guaranteeing $N > 0$ everywhere) and compact (i.e. operating on one triangle only), while sufficiently accurate for nearshore applications.
 Given the action densities $N^n_2$ and $N^n_3$ at vertices 2 and 3 of triangle $\triangle$123,
-the wave action in vertex 1 is determined according to
+the wave action in vertex 1 is readily determined according to
 \begin{eqnarray}
   &&\left[ \frac{1}{\Delta t} + c_{x,1} \left( e^{(1)}_1 + e^{(2)}_1 \right) + c_{y,1} \left( e^{(1)}_2 + e^{(2)}_2 \right) \right] N_1^n = \nonumber \\
   &&\frac{N_1^{n-1}}{\Delta t}+\left( c_{x,2} e^{(1)}_1 + c_{y,2} e^{(1)}_2 \right) N^n_2 + \left( c_{x,3} e^{(2)}_1 + c_{y,3} e^{(2)}_2 \right) N^n_3 + F^n \, .
@@ -3907,10 +4148,6 @@
 Alves, J.H.G.M. and M.L. Banner, 2003: Performance of a saturation-based dissipation-rate source term
 in modelling the fetch-limited evolution of wind waves, {\it J. Phys. Oceanogr}., {\bf 33}, 1274-1298
 
-\bibitem{And98HJK}
-Andorka Gal, J.H., L.H. Holthuijsen, J.C.M. de Jong and A.T.M.M. Kieftenburg, 1998: Wave transformation
-near a quasi-1D coast, {\it 26th Int. Conf. Coastal Engng}., Copenhagen, 150-160
-
 \bibitem{Arc90L}
 Arcilla, A.S. and C.M. Lemos, 1990: {\it Surf-Zone Hydrodynamics}, Centro Internacional de M\'{e}todos
 Num\'{e}ricos en Ingenieria, Barcelona, 310 p.
@@ -4000,10 +4237,6 @@
 Booij, N., L.H. Holthuijsen and R.C. Ris, 1996: The "SWAN" wave model for shallow water, {\it Proc. 25th Int.
 Conf. Coastal Engng.}, Orlando, 668-676
 
-\bibitem{Boo97HDK}
-Booij, N., L.H. Holthuijsen, N. Doorn and A.T.M.M. Kieftenburg, 1997: Diffraction in a spectral wave model,
-{\it Proc. 3rd Int. Symp. Ocean Wave Measurement and Analysis, WAVES'97}, ASCE, 234-255
-
 \bibitem{Boo97HPa}
 Booij, N., L.H. Holthuijsen and R. Padilla-Hernandez, 1997: A nonstationary, parametric coastal wave
 model, {\it Conf. Coastal Dynamics '97}, Plymouth, 99-107
@@ -4071,8 +4304,8 @@
 Chen, Y. and H. Wang, 1983: Numerical model for nonstationary shallow water wave spectral
 transformations, {\it J. Geophys. Res}., {\bf 88}, 9851-9863
 
-\bibitem{Che97G}
-Chen, Y. and R.T. Guza, 1997: Modelling of breaking surface waves in shallow water, {\it J. Geophys. Res}.,
+\bibitem{Che97GE}
+Chen, Y., R.T. Guza, and S. Elgar, 1997: Modelling of breaking surface waves in shallow water, {\it J. Geophys. Res}.,
 102, C11, 25035-25046
 
 \bibitem{Chr94HR}
@@ -4662,6 +4895,9 @@
 Wilson, B.W., 1965: Numerical prediction of ocean waves in the North Atlantic for December 1959,
 {\it Deutsch. Hydrogr. Z.}, {\bf 18}, 3, p. 114-130
 
+\bibitem{WIS07}
+WISE~Group, 2007. Wave modelling - The state of the art. Progr. Oceanogr., {\bf 75}, 603-674
+
 \bibitem{Whi74}
 Whitham, G.B., 1974: {\it Linear and nonlinear waves}, Wiley, New York, 636 p.
 
@@ -4714,7 +4950,7 @@
 
 \bibitem{Zij98W}
 Zijlema, M. and Wesseling, P., 1998. Higher-order flux-limiting schemes for the finite volume
-  computation of incompressible flow. Int. J. Comput. Fluid Dyn., 9:89--109.
+  computation of incompressible flow. Int. J. Comput. Fluid Dyn., {\bf 9}, 89-109
 
 \bibitem{Zij05}
 Zijlema, M., 2005: Parallelization of a nearshore wind wave model for distributed memory architectures,
@@ -4723,198 +4959,208 @@
 Elsevier Science B.V., Amsterdam, The Netherlands, 207-214
 
 \bibitem{Zij05W}
-Zijlema, M. and A.J. van der Westhuysen: On convergence behaviour and numerical accuracy in stationary
-SWAN simulation of nearshore wind wave spectra. {\it Coastal Engineering}, {\bf 52}, 237-256.
+Zijlema, M. and A.J. van der Westhuysen, 2005: On convergence behaviour and numerical accuracy in stationary
+SWAN simulation of nearshore wind wave spectra. {\it Coastal Engineering}, {\bf 52}, 237-256
+
+\bibitem{Zij09}
+Zijlema, M., 2009: Parallel, unstructured mesh implementation for SWAN,
+{\it Proc. 31$^{{\rm th}}$ Int. Conf. Coastal Engineering}, ASCE, 470-482
 
 \end{thebibliography}
 
 \begin{theindex}
 \addcontentsline{toc}{chapter}{Index}
   \item ambient, \hyperpage{2}, \hyperpage{7}, \hyperpage{11, 12}, 
-		\hyperpage{30}, \hyperpage{97}
+		\hyperpage{29, 30}, \hyperpage{101}
 
   \indexspace
 
-  \item bathymetry, \hyperpage{24}, \hyperpage{77}, \hyperpage{79}, 
-		\hyperpage{99}
+  \item bathymetry, \hyperpage{23}, \hyperpage{38}, \hyperpage{81}, 
+		\hyperpage{83}, \hyperpage{103}
   \item bottom, \hyperpage{1--3}, \hyperpage{12}, \hyperpage{14, 15}, 
-		\hyperpage{22--24}, \hyperpage{31}, \hyperpage{79}, 
-		\hyperpage{82, 83}, \hyperpage{94}, \hyperpage{96}, 
-		\hyperpage{99--104}
+		\hyperpage{22, 23}, \hyperpage{30, 31}, 
+		\hyperpage{46, 47}, \hyperpage{83}, \hyperpage{87, 88}, 
+		\hyperpage{98}, \hyperpage{100}, \hyperpage{103--108}
   \item boundary, \hyperpage{3, 4}, \hyperpage{13}, \hyperpage{15}, 
-		\hyperpage{43}, \hyperpage{59, 60}, \hyperpage{62, 63}, 
-		\hyperpage{65--69}, \hyperpage{78}, \hyperpage{80}, 
-		\hyperpage{82}, \hyperpage{98}, \hyperpage{100}
+		\hyperpage{44}, \hyperpage{63, 64}, \hyperpage{66, 67}, 
+		\hyperpage{69--73}, \hyperpage{75}, \hyperpage{82}, 
+		\hyperpage{84}, \hyperpage{88}, \hyperpage{102, 103}
   \item breaking, \hyperpage{1}, \hyperpage{3}, \hyperpage{12}, 
-		\hyperpage{14--17}, \hyperpage{21}, \hyperpage{23--25}, 
-		\hyperpage{42}, \hyperpage{56, 57}, \hyperpage{91, 92}, 
-		\hyperpage{94, 95}, \hyperpage{99}, 
-		\hyperpage{103, 104}
+		\hyperpage{14--17}, \hyperpage{20, 21}, 
+		\hyperpage{23--25}, \hyperpage{43}, \hyperpage{60}, 
+		\hyperpage{62}, \hyperpage{95}, \hyperpage{98, 99}, 
+		\hyperpage{102}, \hyperpage{107, 108}
 
   \indexspace
 
-  \item Cartesian, \hyperpage{3}, \hyperpage{10--12}, \hyperpage{64}, 
-		\hyperpage{81}, \hyperpage{83}
+  \item Cartesian, \hyperpage{3}, \hyperpage{10--12}, \hyperpage{68}, 
+		\hyperpage{85}, \hyperpage{87, 88}
   \item co-ordinate, \hyperpage{3}, \hyperpage{7}, \hyperpage{10--12}, 
-		\hyperpage{27, 28}, \hyperpage{48--50}, \hyperpage{66}, 
-		\hyperpage{68}
-  \item coastal, \hyperpage{1, 2}, \hyperpage{15}, \hyperpage{45}, 
-		\hyperpage{73}, \hyperpage{77}, \hyperpage{92, 93}, 
-		\hyperpage{95}, \hyperpage{97--99}, 
-		\hyperpage{101, 102}
-  \item convergence, \hyperpage{36, 37}, \hyperpage{41}, 
-		\hyperpage{43--47}, \hyperpage{75, 76}, \hyperpage{83}, 
-		\hyperpage{101}, \hyperpage{105}
+		\hyperpage{27, 28}, \hyperpage{52--54}, \hyperpage{70}, 
+		\hyperpage{72}
+  \item coastal, \hyperpage{1, 2}, \hyperpage{15}, \hyperpage{40}, 
+		\hyperpage{49}, \hyperpage{77}, \hyperpage{81}, 
+		\hyperpage{96, 97}, \hyperpage{99}, 
+		\hyperpage{101--103}, \hyperpage{105, 106}
+  \item convergence, \hyperpage{36, 37}, \hyperpage{42}, \hyperpage{44}, 
+		\hyperpage{49--52}, \hyperpage{75, 76}, 
+		\hyperpage{79, 80}, \hyperpage{88}, \hyperpage{105}, 
+		\hyperpage{109}
   \item Courant, \hyperpage{40}
   \item current, \hyperpage{2, 3}, \hyperpage{7}, \hyperpage{10--13}, 
-		\hyperpage{15}, \hyperpage{30}, \hyperpage{32}, 
-		\hyperpage{34}, \hyperpage{37, 38}, \hyperpage{40}, 
-		\hyperpage{43}, \hyperpage{49}, \hyperpage{55}, 
-		\hyperpage{61}, \hyperpage{74}, \hyperpage{76}, 
-		\hyperpage{85}, \hyperpage{94}, \hyperpage{97}, 
-		\hyperpage{101--104}
-  \item curvi-linear, \hyperpage{2, 3}, \hyperpage{48, 49}, 
-		\hyperpage{51}, \hyperpage{62}, \hyperpage{68, 69}, 
-		\hyperpage{77}, \hyperpage{93}, \hyperpage{97}
+		\hyperpage{15}, \hyperpage{29, 30}, \hyperpage{32--34}, 
+		\hyperpage{37, 38}, \hyperpage{41}, \hyperpage{44--47}, 
+		\hyperpage{54}, \hyperpage{59}, \hyperpage{65}, 
+		\hyperpage{78}, \hyperpage{80}, \hyperpage{90}, 
+		\hyperpage{98}, \hyperpage{100, 101}, 
+		\hyperpage{105--108}
+  \item curvi-linear, \hyperpage{2, 3}, \hyperpage{52--54}, 
+		\hyperpage{56}, \hyperpage{66}, \hyperpage{72, 73}, 
+		\hyperpage{81}, \hyperpage{96}, \hyperpage{101}
 
   \indexspace
 
-  \item dam, \hyperpage{16}, \hyperpage{30--32}, \hyperpage{97}, 
-		\hyperpage{105}
+  \item dam, \hyperpage{16}, \hyperpage{30, 31}, \hyperpage{101}, 
+		\hyperpage{109}
   \item diffraction, \hyperpage{3}, \hyperpage{10}, \hyperpage{30}, 
-		\hyperpage{32--34}, \hyperpage{87}, \hyperpage{92}, 
-		\hyperpage{98, 99}
-  \item diffusion, \hyperpage{38--40}, \hyperpage{86, 87}
+		\hyperpage{32, 33}, \hyperpage{92}, \hyperpage{96}, 
+		\hyperpage{102, 103}
+  \item diffusion, \hyperpage{38--41}, \hyperpage{92}
   \item dissipation, \hyperpage{1}, \hyperpage{3}, \hyperpage{12}, 
 		\hyperpage{14--16}, \hyperpage{20--24}, 
-		\hyperpage{56, 57}, \hyperpage{78}, \hyperpage{91--93}, 
-		\hyperpage{96}, \hyperpage{101}, \hyperpage{103, 104}
+		\hyperpage{60, 61}, \hyperpage{82}, \hyperpage{95}, 
+		\hyperpage{97}, \hyperpage{100}, \hyperpage{105}, 
+		\hyperpage{107, 108}
 
   \indexspace
 
   \item filter, \hyperpage{13}, \hyperpage{18}, \hyperpage{33}
-  \item flow, \hyperpage{7}, \hyperpage{61}, \hyperpage{93}, 
-		\hyperpage{98}, \hyperpage{100}, \hyperpage{105}
-  \item force, \hyperpage{34}, \hyperpage{49, 50}, \hyperpage{61--63}, 
-		\hyperpage{68, 69}, \hyperpage{79}, \hyperpage{81}, 
-		\hyperpage{85}, \hyperpage{94}
+  \item flow, \hyperpage{7}, \hyperpage{65}, \hyperpage{97}, 
+		\hyperpage{102}, \hyperpage{104}, \hyperpage{109}
+  \item force, \hyperpage{33, 34}, \hyperpage{54}, \hyperpage{65--67}, 
+		\hyperpage{72, 73}, \hyperpage{83}, \hyperpage{87}, 
+		\hyperpage{90}, \hyperpage{98}
   \item frequency, \hyperpage{7--9}, \hyperpage{11}, \hyperpage{13, 14}, 
-		\hyperpage{16--23}, \hyperpage{25, 26}, 
-		\hyperpage{28--30}, \hyperpage{35--37}, \hyperpage{40}, 
-		\hyperpage{43--46}, \hyperpage{53--55}, 
-		\hyperpage{59, 60}, \hyperpage{91}, \hyperpage{99}, 
-		\hyperpage{105}
+		\hyperpage{16}, \hyperpage{18--21}, \hyperpage{23}, 
+		\hyperpage{25}, \hyperpage{28}, \hyperpage{30}, 
+		\hyperpage{35--37}, \hyperpage{41, 42}, \hyperpage{44}, 
+		\hyperpage{46--50}, \hyperpage{57--60}, 
+		\hyperpage{63, 64}, \hyperpage{95}, \hyperpage{103}, 
+		\hyperpage{109}
   \item friction, \hyperpage{1}, \hyperpage{3}, \hyperpage{12--15}, 
-		\hyperpage{18}, \hyperpage{22, 23}, \hyperpage{96}, 
-		\hyperpage{98, 99}, \hyperpage{101}, 
-		\hyperpage{103, 104}
+		\hyperpage{18}, \hyperpage{21, 22}, \hyperpage{100}, 
+		\hyperpage{102, 103}, \hyperpage{105}, 
+		\hyperpage{107, 108}
 
   \indexspace
 
-  \item garden-sprinkler, \hyperpage{39, 40}, \hyperpage{87}
+  \item garden-sprinkler, \hyperpage{39, 40}, \hyperpage{92}
 
   \indexspace
 
-  \item harbour, \hyperpage{3}, \hyperpage{34}
+  \item harbour, \hyperpage{3}, \hyperpage{33}
 
   \indexspace
 
-  \item initial, \hyperpage{13}, \hyperpage{15}, \hyperpage{29}, 
-		\hyperpage{59, 60}
-  \item island, \hyperpage{77, 78}
+  \item initial, \hyperpage{13}, \hyperpage{15}, \hyperpage{28}, 
+		\hyperpage{45}, \hyperpage{63, 64}, \hyperpage{76}
+  \item island, \hyperpage{81, 82}
 
   \indexspace
 
-  \item Jonswap, \hyperpage{59, 60}
+  \item Jonswap, \hyperpage{63, 64}
 
   \indexspace
 
   \item latitude, \hyperpage{11, 12}
-  \item limiter, \hyperpage{36}, \hyperpage{44--47}
+  \item limiter, \hyperpage{36}, \hyperpage{49, 50}, \hyperpage{52}
   \item longitude, \hyperpage{11, 12}
 
   \indexspace
 
-  \item obstacle, \hyperpage{3}, \hyperpage{30--34}, \hyperpage{38, 39}, 
-		\hyperpage{50--52}
+  \item obstacle, \hyperpage{3}, \hyperpage{30--33}, \hyperpage{38, 39}, 
+		\hyperpage{55--57}
   \item ocean, \hyperpage{1}, \hyperpage{7, 8}, \hyperpage{11}, 
-		\hyperpage{30}, \hyperpage{34}, \hyperpage{38}, 
-		\hyperpage{40}, \hyperpage{92}, \hyperpage{96--100}, 
-		\hyperpage{102}, \hyperpage{104, 105}
+		\hyperpage{29}, \hyperpage{33}, \hyperpage{40}, 
+		\hyperpage{96}, \hyperpage{99--105}, 
+		\hyperpage{108, 109}
 
   \indexspace
 
   \item propagation, \hyperpage{1--4}, \hyperpage{10--12}, 
-		\hyperpage{32, 33}, \hyperpage{37}, \hyperpage{39}, 
-		\hyperpage{41, 42}, \hyperpage{48, 49}, \hyperpage{51}, 
-		\hyperpage{73, 74}, \hyperpage{77}, \hyperpage{79}, 
-		\hyperpage{81, 82}, \hyperpage{92--95}, 
-		\hyperpage{101}, \hyperpage{103}
+		\hyperpage{32}, \hyperpage{37}, \hyperpage{39--43}, 
+		\hyperpage{45--47}, \hyperpage{53}, \hyperpage{56}, 
+		\hyperpage{77, 78}, \hyperpage{81}, \hyperpage{83}, 
+		\hyperpage{85--88}, \hyperpage{96}, \hyperpage{98, 99}, 
+		\hyperpage{105}, \hyperpage{107}
 
   \indexspace
 
-  \item quadruplets, \hyperpage{16, 17}, \hyperpage{25}
+  \item quadruplets, \hyperpage{16}, \hyperpage{25}, \hyperpage{48}
 
   \indexspace
 
-  \item recti-linear, \hyperpage{69}, \hyperpage{77}
+  \item recti-linear, \hyperpage{73}, \hyperpage{81}
   \item reflection, \hyperpage{3}, \hyperpage{30}, \hyperpage{32, 33}, 
-		\hyperpage{35}, \hyperpage{52}
+		\hyperpage{35}, \hyperpage{56}
   \item refraction, \hyperpage{1}, \hyperpage{3}, \hyperpage{11}, 
-		\hyperpage{32}, \hyperpage{43}, \hyperpage{59}, 
-		\hyperpage{82}, \hyperpage{92, 93}, \hyperpage{98, 99}, 
-		\hyperpage{102}, \hyperpage{104}
+		\hyperpage{32}, \hyperpage{38}, \hyperpage{44}, 
+		\hyperpage{46}, \hyperpage{63}, \hyperpage{87}, 
+		\hyperpage{96, 97}, \hyperpage{102, 103}, 
+		\hyperpage{106}, \hyperpage{108}
   \item regular, \hyperpage{2}, \hyperpage{4}, \hyperpage{7}, 
-		\hyperpage{15}, \hyperpage{30}, \hyperpage{82}, 
-		\hyperpage{102}, \hyperpage{104}
+		\hyperpage{15}, \hyperpage{30}, \hyperpage{87}, 
+		\hyperpage{106}, \hyperpage{108}
 
   \indexspace
 
-  \item set-up, \hyperpage{3, 4}, \hyperpage{24}, \hyperpage{34}, 
-		\hyperpage{61}, \hyperpage{69}, \hyperpage{91}
+  \item set-up, \hyperpage{3, 4}, \hyperpage{23}, \hyperpage{33, 34}, 
+		\hyperpage{65}, \hyperpage{73}, \hyperpage{95}
   \item shoaling, \hyperpage{1}, \hyperpage{3}, \hyperpage{15}, 
-		\hyperpage{25}, \hyperpage{95}
-  \item SORDUP, \hyperpage{38}, \hyperpage{40}, \hyperpage{73}
-  \item spherical, \hyperpage{3}, \hyperpage{11, 12}, 
-		\hyperpage{82, 83}
-  \item stability, \hyperpage{35, 36}, \hyperpage{40}, \hyperpage{43}, 
-		\hyperpage{45, 46}, \hyperpage{56}, \hyperpage{58}, 
-		\hyperpage{80}, \hyperpage{82}
+		\hyperpage{25}, \hyperpage{98, 99}
+  \item SORDUP, \hyperpage{38}, \hyperpage{40}, \hyperpage{77}
+  \item spherical, \hyperpage{3}, \hyperpage{11, 12}, \hyperpage{88}
+  \item stability, \hyperpage{35, 36}, \hyperpage{40}, \hyperpage{44}, 
+		\hyperpage{50}, \hyperpage{61, 62}, \hyperpage{84}, 
+		\hyperpage{87}
   \item stationary, \hyperpage{2, 3}, \hyperpage{7}, \hyperpage{36--38}, 
-		\hyperpage{40}, \hyperpage{44}, \hyperpage{49}, 
-		\hyperpage{60, 61}, \hyperpage{93, 94}, \hyperpage{97}, 
-		\hyperpage{105}
-  \item steepness, \hyperpage{14}, \hyperpage{20}, \hyperpage{24}, 
+		\hyperpage{40}, \hyperpage{48}, \hyperpage{53}, 
+		\hyperpage{64, 65}, \hyperpage{96}, \hyperpage{98}, 
+		\hyperpage{101}, \hyperpage{109}
+  \item steepness, \hyperpage{14}, \hyperpage{19, 20}, \hyperpage{23}, 
 		\hyperpage{30}
   \item swell, \hyperpage{14, 15}, \hyperpage{20}, \hyperpage{22}, 
-		\hyperpage{93}, \hyperpage{96}, \hyperpage{104}
+		\hyperpage{40, 41}, \hyperpage{97}, \hyperpage{100}, 
+		\hyperpage{108}
 
   \indexspace
 
   \item triads, \hyperpage{16}
-  \item triangular, \hyperpage{78}, \hyperpage{80}
+  \item triangular, \hyperpage{76}, \hyperpage{82}, \hyperpage{84--86}
 
   \indexspace
 
-  \item unstructured, \hyperpage{4}, \hyperpage{77, 78}, \hyperpage{82}, 
-		\hyperpage{85, 86}, \hyperpage{92}
+  \item unstructured, \hyperpage{4}, \hyperpage{81, 82}, \hyperpage{87}, 
+		\hyperpage{90}, \hyperpage{92}, \hyperpage{96}, 
+		\hyperpage{109}
 
   \indexspace
 
   \item WAM, \hyperpage{1}, \hyperpage{13, 14}, \hyperpage{18--20}, 
-		\hyperpage{26}, \hyperpage{35, 36}, \hyperpage{44, 45}, 
-		\hyperpage{55}, \hyperpage{59}, \hyperpage{96, 97}, 
-		\hyperpage{102}, \hyperpage{104}
-  \item WAVEWATCH, \hyperpage{1}, \hyperpage{59}, \hyperpage{103}
+		\hyperpage{26}, \hyperpage{35, 36}, \hyperpage{49}, 
+		\hyperpage{59}, \hyperpage{63}, \hyperpage{99}, 
+		\hyperpage{101}, \hyperpage{106}, \hyperpage{108}
+  \item WAVEWATCH, \hyperpage{1}, \hyperpage{63}, \hyperpage{107}
   \item whitecapping, \hyperpage{3}, \hyperpage{12}, \hyperpage{14}, 
-		\hyperpage{16}, \hyperpage{19--22}, \hyperpage{30}, 
-		\hyperpage{96}, \hyperpage{98}, \hyperpage{103, 104}
+		\hyperpage{16}, \hyperpage{19--21}, \hyperpage{30}, 
+		\hyperpage{100}, \hyperpage{102}, \hyperpage{107, 108}
   \item wind, \hyperpage{1--4}, \hyperpage{7}, \hyperpage{12--15}, 
-		\hyperpage{17--22}, \hyperpage{30}, \hyperpage{35}, 
-		\hyperpage{37--42}, \hyperpage{44, 45}, \hyperpage{50}, 
-		\hyperpage{60}, \hyperpage{68}, \hyperpage{77}, 
-		\hyperpage{81}, \hyperpage{91--93}, \hyperpage{96}, 
-		\hyperpage{98--105}
+		\hyperpage{17--22}, \hyperpage{29}, \hyperpage{35}, 
+		\hyperpage{37--41}, \hyperpage{43}, \hyperpage{45}, 
+		\hyperpage{48, 49}, \hyperpage{55}, \hyperpage{64}, 
+		\hyperpage{72}, \hyperpage{81}, \hyperpage{87}, 
+		\hyperpage{95--97}, \hyperpage{100}, 
+		\hyperpage{102--109}
 
 \end{theindex}
 
--- SwanTranspAc.ftn90	2009-08-05 23:35:15.000000000 +0200
+++ SwanTranspAc.ftn90	2009-08-05 21:33:21.000000000 +0200
@@ -2,8 +2,8 @@
                           cgo   , cax   , cay   , cad   , cas   , &
                           anybin, rdx   , rdy   , spcsig, spcdir, &
                           obredf, idcmin, idcmax, iscmin, iscmax, &
-                          iddlow, iddtop, isslow, isstop, trac0 , &
-                          trac1 )
+                          iddlow, iddtop, isslow, isstop, anyblk, &
+                          trac0 , trac1 )
 !
 !   --|-----------------------------------------------------------|--
 !     | Delft University of Technology                            |
@@ -43,6 +43,7 @@
 !   40.80,   August 2007: New subroutine
 !   40.85,   August 2008: add tranport terms for output purposes
 !   41.00, February 2009: add GSE correction
+!   41.07,     July 2009: add explicit scheme for sigma space
 !
 !   Purpose
 !
@@ -99,6 +100,7 @@
     real, dimension(MDC,MSC,MTRNP), intent(out) :: trac1  ! implicit part of propagation in present vertex for output purposes
     !
     logical, dimension(MDC,MSC), intent(in)     :: anybin ! true if bin is active in considered sweep
+    logical, dimension(MDC,MSC), intent(out)    :: anyblk ! true if bin is blocked by a counter current based on a CFL criterion
 !
 !   Local variables
 !
@@ -153,19 +155,22 @@
 !TIMG    call SWTSTA(141)
     if ( (DYNDEP .OR. ICUR /= 0) .and. ITFRE /= 0 ) then
        !
-       ! only implicit scheme
-       !
        if ( int(PNUMS(8)) == 1 ) then
           !
+          ! implicit scheme
+          !
           call STRSSI ( spcsig     , cas   , amat(1,1,2), amat(1,1,1), &
                         amat(1,1,3), anybin, rhs        , ac2        , &
                         iscmin     , iscmax, iddlow     , iddtop     , &
                         trac0      , trac1 )
           !
-       else
+       elseif ( int(PNUMS(8)) == 2 ) then
+          !
+          ! explicit scheme
           !
-          call msgerr (4, 'No explicit scheme in sigma space available ')
-          return
+          call STRSSB ( iddlow, iddtop, idcmin, idcmax, isstop, &
+                        cax   , cay   , cas   , ac2   , spcsig, &
+                        rhs   , anyblk, rdx   , rdy   , trac0 )
           !
        endif
        !
--- swanuse.tex	2009-08-05 23:35:15.000000000 +0200
+++ swanuse.tex	2009-08-05 21:41:58.000000000 +0200
@@ -32,7 +32,7 @@
 \end{center}
 \vfill
 \begin{center}
-{\Large\bf SWAN Cycle III version 40.72AB}
+{\Large\bf SWAN Cycle III version 40.72ABC}
 \end{center}
 
 \cleardoublepage
@@ -248,9 +248,9 @@
         encountered during SWAN computations, the current is locally reduced to sub-critical flow.
   \item If the \underline{water depth} is less than some user-provided limit, the depth is set at that limit
         (default is 0.05 m, see command {\tt SET}).
-  \item The \underline{user-imposed wave boundary conditions} may not be reproduced by SWAN employing
-        \underline {structured grids}, as SWAN replaces the {\it imposed} waves at the boundaries that
-        propagate into the computational area with the {\it computed} waves that move out of  the computational
+  \item The \underline{user-imposed wave boundary conditions} may not be reproduced by SWAN,
+        as SWAN replaces the {\it imposed} waves at the boundaries that
+        propagate into the computational area with the {\it computed} waves that move out of the computational
         area at the boundaries.
   \item SWAN may have \underline{convergence} problems. There are three iteration processes in SWAN:
         \begin{enumerate}
@@ -537,8 +537,7 @@
 that this output differs from the boundary conditions that are imposed by the user. The reason is that
 SWAN accepts only the user-imposed \underline{incoming} wave components and that it replaces the user-imposed
 \underline{outgoing} wave components with computed outgoing components (propagating to the boundary from
-the interior region). {\bf This is only the case for structured grids (both regular and curvi-linear ones)}.
-The user is informed by means of a warning in the output when the
+the interior region). The user is informed by means of a warning in the output when the
 computed significant wave height differs more than 10\%, say (10\% is default), from the user-imposed
 significant wave height (command {\tt BOUND}...). The actual value of this difference can be set by the
 user (see the {\tt SET} command). Note that this warning will not apply in the case of unstructured grids.
@@ -1091,15 +1090,16 @@
                    \poptabs
                    Default: {\tt [inrhog]} = 0.\-\\
 {\tt [hsrerr]}  \> the relative difference between the user imposed significant wave height and the\+\\
-                   significant wave height computed by SWAN (anywhere along the boundary of\\
-                   structured grid) above which a warning will be given. This relative difference\\
-                   is the difference normalized with the user provided significant wave height.\\
-                   This warning will be given for each boundary grid point where the problem occurs\\
+                   significant wave height computed by SWAN (anywhere along the computational\\
+                   grid boundary) above which a warning will be given. This relative difference\\
+                   is the difference normalized with the user provided significant wave height. This\\
+                   warning will be given for each boundary grid point where the problem occurs\\
                    (with its $x-$ and $y-$index number of the computational grid). The cause of the\\
                    difference is explained in Section~\ref{sec:boundc}. To supress these warnings (in particular\\
                    for nonstationary computations), set {\tt [hsrerr]} at a very high value or use\\
                    command {\tt OFF BNDCHK}.\\
-                   Default: {\tt [hsrerr]} = 0.10.\-\\
+                   Default: {\tt [hsrerr]} = 0.10.\\
+                   ONLY MEANT FOR STRUCTURED GRIDS.\-\\
 {\tt NAUTICAL}  \> indicates that the Nautical convention for wind and wave direction (SWAN input\+\\
                    and output) will be used instead of the default Cartesian convention.\\
                    For definition, see Section~\ref{sec:units} or Appendix~\ref{app:defvar}.\-\\
@@ -2482,7 +2482,7 @@
 \idxcmd{BREAKING}
 \linecmd
 \begin{verbatim}
-BREaking  CONstant [alpha] [gamma]           
+BREaking  CONstant [alpha] [gamma]
 \end{verbatim}
 \linecmd
 
@@ -2926,13 +2926,12 @@
 {\tt STOPC}    \> With this alternative option the user can influence the criterion for terminating\+\\
                   the iterative procedure in the SWAN computations (both stationary and\\
                   nonstationary). The criterion make use of the second derivative, or curvature,\\
-                  of the iteration curve of both the significant wave height and the mean period.\\
-                  As the solution of a simulation approaches full convergence, the curvature of the\\
-                  iteration curve will tend to zero. SWAN stops the process {\bf if} the absolute change\\
-                  in both $H_s$ and $T_{m01}$ from one iteration to the next is less than {\tt [dabs]} {\bf or} the\\
-                  relative change in $H_s$ and $T_{m01}$ from one iteration to the next is less than\\
-                  {\tt [drel]} {\bf and} the curvature of the iteration curve of $H_s$ normalized with $H_s$\\
-                  and that of $T_{m01}$ normalized with $T_{m01}$ is less than {\tt [curvat]}.\\
+                  of the iteration curve of the significant wave height. As the solution of a\\
+                  simulation approaches full convergence, the curvature of the iteration curve will\\
+                  tend to zero. SWAN stops the process {\bf if} the absolute change in $H_s$ from one\\
+                  iteration to the next is less than {\tt [dabs]} {\bf or} the relative change in $H_s$ from one\\
+                  iteration to the next is less than {\tt [drel]} {\bf and} the curvature of the iteration\\
+                  curve of $H_s$ normalized with $H_s$ is less than {\tt [curvat]}.\\
                   DEFAULT IN CASE OF UNSTRUCTURED GRIDS.\-\\
 {\tt [dabs]}   \> Default: {\tt [dabs]} = 0.00 [$-$] in case of structured grids; {\tt [dabs]} = 0.005 [$-$]\+\\
                   in case of unstructured grids.\-\\
@@ -2971,18 +2970,19 @@
                   upwind scheme and it is more diffusive and therefore preferable if (strong)\\
                   gradients in depth or current are present.\\
                   Default: {\tt [cdd]} = 0.5.\-\\
-{\tt [cdlim]}  \> If the spatial discretization of the bathymetry or the flow currents is too coarse,\+\\
-                  the waves may turn too far (more than 90 degrees, say) over one spatial grid step.\\
-                  The computational results will then be very inaccurate. In such a case SWAN can\\
-                  limit the maximum turning of the waves over one spatial grid to 90 degrees to\\
-                  obtain robust (but not necessarily correct results).\\
+{\tt [cdlim]}  \> If the spatial discretization of the bathymetry or flow currents is too coarse, the\+\\
+                  waves may turn too far (more than 90 degrees, say) over one spatial grid step.\\
+                  The computational results will then be very inaccurate. In such a case SWAN\\
+                  can limit the maximum turning of the waves over one spatial grid to 90 degrees\\
+                  or so to obtain robust, but not necessarily correct, results.\\
                   \pushtabs
                   xxxxxxxxxxxxxxx\=xxxxxxxxxxx \kill
-                  {\tt [cdlim]} $<$ 0 \> then no limiter is used (this is default)\\
+                  {\tt [cdlim]} $<$ 0 \> no limiter is used\\
                   {\tt [cdlim]} = 0   \> refraction is off (same effect as command {\tt OFF REFRAC})\\
-                  {\tt [cdlim]} = 4   \> waves turning limited to about 90$^o$ over one spatial grid\+\\
-                                         step.\-\-\\
+                  {\tt [cdlim]} $>$ 0 \> wave turning is limited. Based on experiences, a suggestion\+\\
+                  for this parameter is {\tt [cdlim]} = 4.\-\\
                   \poptabs
+                  Default: {\tt [cdlim]} = $-$1 (i.e. no limiter!).\-\\
 {\tt SIGIMPL}  \> controls the accuracy of computing the frequency shifting and the stopping criterion\+\\
                   and amount of output for the SIP solver (used in the computations in the presence\\
                   of currents or time varying depth).\-\\
@@ -4704,12 +4704,15 @@
 !        |                                                             |
 !        |    VAR [alpha] [gammin] [gammax] [gamneg] [coeff1] [coeff2] |
 !        |                                                             |
-!   BRE <     RUE [alpha] [a] [b]                                       >
+!   BRE <     RUE [alpha] [a] [b]                                       >   &
 !        |                                                             |
 !        |    TG  [alpha] [gamma] [pown]                               |
 !        |                                                             |
 !        !    BIP [alpha] [pown] [bref]                                |
 !
+!     ( FREQDep [power] [fmin] [fmax] )
+!
+!
 !              |             | -> CONstant [cfjon]
 !              | -> JONswap <
 !              |             |    VARiable [cfj1] [cfj2] [dsp1] [dsp2]
@@ -4743,17 +4746,21 @@
 !   PROP   /  BSBT
 !          \  GSE  [waveage] SEC|MIN|HR|DAY
 !
-!             | -> ACCUR [drel] [dhoval] [dtoval] [npnts] |
-!   NUMeric (<                                             >                &
-!             |    STOPC [dabs] [drel] [curvat] [npnts]   |
+!             | -> ACCUR [drel] [dhoval] [dtoval] [npnts]               |
+!   NUMeric (<                                                           >  &
+!             |    STOPC [dabs] [drel] [curvat] [npnts] [dtabs] [curvt] |
 !
 !                    | -> STAT  [mxitst] [alfa] |
 !                   <                            >  [limiter]   )           &
 !                    | NONSTat  [mxitns]        |
 !
-!           ( DIRimpl [cdd] [cdlim]                             )           &
+!           ( DIRimpl [cdd] [cdlim]  WNUMber                    )           &
 !
-!           ( SIGIMpl [css] [eps2] [outp] [niter]               )           &
+!           ( REFRLim [frlim] [power]                           )           &
+!
+!              | -> SIGIMpl [css] [eps2] [outp] [niter]
+!           ( <                                                 )           &
+!              |    SIGEXpl [css] [cfl]
 !
 !           ( SETUP [eps2] [outp] [niter]                       )
 !
@@ -4849,7 +4856,7 @@
 
 \begin{verbatim}
 SWAN   1                                Swan standard spectral file, version
-$ Data produced by SWAN version 40.72AB
+$ Data produced by SWAN version 40.72ABC
 $ Project:'projname'     ;   run number:'runnum'
 TIME                                    time-dependent data
      1                                  time coding option
@@ -4955,7 +4962,7 @@
 
 \begin{verbatim}
 SWAN   1                                Swan standard spectral file, version
-$ Data produced by SWAN version 40.72AB
+$ Data produced by SWAN version 40.72ABC
 $ Project:'projname'     ;   run number:'runnum'
 LOCATIONS                               locations in x-y-space
      1                                  number of locations
@@ -5135,7 +5142,7 @@
 
 \bibitem{Implman}
 {SWAN -- Implementation manual}. Delft University of Technology, Environmental Fluid Mechanics Section, available from
-\hl{http://www.swan.tudelft.nl} (Version 40.72AB, March 2009).
+\hl{http://www.swan.tudelft.nl} (Version 40.72ABC, August 2009).
 
 \bibitem{Progrul}
 {SWAN -- Programming rules}. Delft University of Technology, Environmental Fluid Mechanics Section, available from
@@ -5143,7 +5150,7 @@
 
 \bibitem{Techdoc}
 {SWAN -- Scientific and Technical documentation}. Delft University of Technology, Environmental Fluid Mechanics Section, available from
-\hl{http://www.swan.tudelft.nl} (Version 40.72AB, March 2009).
+\hl{http://www.swan.tudelft.nl} (Version 40.72ABC, August 2009).
 
 \end{thebibliography}
 
@@ -5153,15 +5160,15 @@
 
   \indexspace
 
-  \item bathymetry, \hyperpage{5}, \hyperpage{10, 11}, \hyperpage{64}, 
-		\hyperpage{85}
+  \item bathymetry, \hyperpage{5}, \hyperpage{10, 11}, \hyperpage{65}, 
+		\hyperpage{86}
   \item BLOCK, \hyperpage{75}
   \item bottom, \hyperpage{3}, \hyperpage{8--13}, \hyperpage{16}, 
 		\hyperpage{18, 19}, \hyperpage{21, 22}, \hyperpage{24}, 
 		\hyperpage{34, 35}, \hyperpage{38--40}, \hyperpage{53}, 
-		\hyperpage{55, 56}, \hyperpage{58, 59}, 
-		\hyperpage{66, 67}, \hyperpage{69}, \hyperpage{79, 80}, 
-		\hyperpage{86}, \hyperpage{92}, \hyperpage{94}
+		\hyperpage{55, 56}, \hyperpage{59}, \hyperpage{66, 67}, 
+		\hyperpage{69}, \hyperpage{79, 80}, \hyperpage{86}, 
+		\hyperpage{92}, \hyperpage{94}
   \item BOUND SHAPE, \hyperpage{41}
   \item boundary, \hyperpage{1}, \hyperpage{3}, \hyperpage{6--14}, 
 		\hyperpage{22}, \hyperpage{26}, \hyperpage{28}, 
@@ -5175,7 +5182,7 @@
   \item BOUNDSPEC, \hyperpage{42}
   \item BREAKING, \hyperpage{54}
   \item breaking, \hyperpage{11}, \hyperpage{16--18}, \hyperpage{22}, 
-		\hyperpage{53--55}, \hyperpage{57}, \hyperpage{61}, 
+		\hyperpage{53--55}, \hyperpage{57}, \hyperpage{62}, 
 		\hyperpage{80}, \hyperpage{86}, \hyperpage{92}
 
   \indexspace
@@ -5191,19 +5198,19 @@
   \item co-ordinate, \hyperpage{32}
   \item coastal, \hyperpage{3--5}, \hyperpage{15}
   \item COMPUTE, \hyperpage{87}
-  \item convergence, \hyperpage{6, 7}, \hyperpage{63}
+  \item convergence, \hyperpage{6, 7}, \hyperpage{64}
   \item COORDINATES, \hyperpage{27}
   \item Courant, \hyperpage{14}
   \item current, \hyperpage{3}, \hyperpage{5, 6}, \hyperpage{9--12}, 
 		\hyperpage{14}, \hyperpage{16}, \hyperpage{19--21}, 
-		\hyperpage{26}, \hyperpage{32}, \hyperpage{34, 35}, 
+		\hyperpage{26, 27}, \hyperpage{32}, \hyperpage{34, 35}, 
 		\hyperpage{38}, \hyperpage{45}, \hyperpage{52}, 
-		\hyperpage{64, 65}, \hyperpage{67}, \hyperpage{79, 80}, 
+		\hyperpage{65}, \hyperpage{67}, \hyperpage{79, 80}, 
 		\hyperpage{83}, \hyperpage{90}
   \item CURVE, \hyperpage{68}
   \item curvi-linear, \hyperpage{3}, \hyperpage{6}, \hyperpage{8, 9}, 
-		\hyperpage{12--14}, \hyperpage{16}, \hyperpage{21, 22}, 
-		\hyperpage{28, 29}, \hyperpage{32}, \hyperpage{34, 35}, 
+		\hyperpage{13, 14}, \hyperpage{16}, \hyperpage{21, 22}, 
+		\hyperpage{29}, \hyperpage{32}, \hyperpage{34, 35}, 
 		\hyperpage{38, 39}, \hyperpage{43}, \hyperpage{47, 48}, 
 		\hyperpage{50}, \hyperpage{62}, \hyperpage{66}, 
 		\hyperpage{68}
@@ -5215,20 +5222,20 @@
   \item diffraction, \hyperpage{5}, \hyperpage{22}, \hyperpage{60, 61}, 
 		\hyperpage{64}
   \item diffusion, \hyperpage{62}, \hyperpage{64, 65}
-  \item dissipation, \hyperpage{22}, \hyperpage{52, 53}, \hyperpage{55}, 
-		\hyperpage{61}, \hyperpage{80}, \hyperpage{92}
+  \item dissipation, \hyperpage{22}, \hyperpage{52--56}, \hyperpage{62}, 
+		\hyperpage{80}, \hyperpage{92}
 
   \indexspace
 
-  \item filter, \hyperpage{13, 14}, \hyperpage{60}
+  \item filter, \hyperpage{13, 14}, \hyperpage{60, 61}
   \item flow, \hyperpage{6}, \hyperpage{30, 31}, \hyperpage{34}, 
-		\hyperpage{36, 37}, \hyperpage{64}
+		\hyperpage{36, 37}, \hyperpage{65}
   \item force, \hyperpage{80}, \hyperpage{93}
-  \item FRAME, \hyperpage{66}
+  \item FRAME, \hyperpage{67}
   \item frequency, \hyperpage{4}, \hyperpage{6}, \hyperpage{8}, 
 		\hyperpage{14, 15}, \hyperpage{26}, \hyperpage{30, 31}, 
 		\hyperpage{42}, \hyperpage{45}, \hyperpage{51, 52}, 
-		\hyperpage{54--56}, \hyperpage{61}, \hyperpage{64--66}, 
+		\hyperpage{54--56}, \hyperpage{62}, \hyperpage{64--66}, 
 		\hyperpage{72--74}, \hyperpage{79}, \hyperpage{83}, 
 		\hyperpage{86}, \hyperpage{89}, \hyperpage{91}, 
 		\hyperpage{108, 109}
@@ -5237,11 +5244,11 @@
 		\hyperpage{21, 22}, \hyperpage{24}, \hyperpage{34, 35}, 
 		\hyperpage{38}, \hyperpage{53}, \hyperpage{55, 56}, 
 		\hyperpage{79, 80}, \hyperpage{86}, \hyperpage{92}
-  \item Froude, \hyperpage{26}, \hyperpage{38}
+  \item Froude, \hyperpage{26, 27}, \hyperpage{38}
 
   \indexspace
 
-  \item garden-sprinkler, \hyperpage{62}
+  \item garden-sprinkler, \hyperpage{62, 63}
   \item GEN1, \hyperpage{52}
   \item GEN2, \hyperpage{52}
   \item GEN3, \hyperpage{53}
@@ -5270,7 +5277,7 @@
 
   \item latitude, \hyperpage{8}, \hyperpage{27}, \hyperpage{86}, 
 		\hyperpage{108}
-  \item LIMITER, \hyperpage{56}
+  \item LIMITER, \hyperpage{57}
   \item limiter, \hyperpage{5}, \hyperpage{7}, \hyperpage{57}, 
 		\hyperpage{64, 65}
   \item longitude, \hyperpage{8}, \hyperpage{27}, \hyperpage{49, 50}, 
@@ -5309,7 +5316,7 @@
 
   \item QUADRUPL, \hyperpage{53}
   \item quadruplets, \hyperpage{15, 16}, \hyperpage{22}, \hyperpage{54}, 
-		\hyperpage{57}, \hyperpage{79}, \hyperpage{92}
+		\hyperpage{57}, \hyperpage{80}, \hyperpage{92}
   \item QUANTITY, \hyperpage{72}
 
   \indexspace
@@ -5318,21 +5325,20 @@
   \item READGRID COORDINATES, \hyperpage{32}
   \item READGRID UNSTRUCTURED, \hyperpage{32}
   \item READINP, \hyperpage{37}
-  \item recti-linear, \hyperpage{3}, \hyperpage{16}, \hyperpage{28, 29}, 
+  \item recti-linear, \hyperpage{3}, \hyperpage{16}, \hyperpage{29}, 
 		\hyperpage{38}, \hyperpage{68}
-  \item reflection, \hyperpage{57--59}
-  \item refraction, \hyperpage{6}, \hyperpage{11}, \hyperpage{61}, 
+  \item reflection, \hyperpage{58, 59}
+  \item refraction, \hyperpage{6}, \hyperpage{11}, \hyperpage{62}, 
 		\hyperpage{64, 65}, \hyperpage{69}
-  \item regular, \hyperpage{4}, \hyperpage{8, 9}, \hyperpage{12}, 
-		\hyperpage{16}, \hyperpage{22}, \hyperpage{28, 29}, 
-		\hyperpage{31, 32}, \hyperpage{34, 35}, \hyperpage{62}, 
-		\hyperpage{66, 67}
+  \item regular, \hyperpage{4}, \hyperpage{8, 9}, \hyperpage{16}, 
+		\hyperpage{22}, \hyperpage{28, 29}, \hyperpage{31, 32}, 
+		\hyperpage{34, 35}, \hyperpage{62}, \hyperpage{66, 67}
 
   \indexspace
 
   \item SET, \hyperpage{25}
   \item set-up, \hyperpage{4}, \hyperpage{6}, \hyperpage{16, 17}, 
-		\hyperpage{22}, \hyperpage{28}, \hyperpage{59, 60}, 
+		\hyperpage{22}, \hyperpage{28}, \hyperpage{60}, 
 		\hyperpage{65}
   \item SETUP, \hyperpage{59}
   \item shoaling, \hyperpage{18}
@@ -5349,13 +5355,13 @@
 		\hyperpage{14}, \hyperpage{16}, \hyperpage{18}, 
 		\hyperpage{20, 21}, \hyperpage{26, 27}, \hyperpage{34}, 
 		\hyperpage{37--39}, \hyperpage{45, 46}, \hyperpage{51}, 
-		\hyperpage{53}, \hyperpage{62--64}, \hyperpage{80}, 
+		\hyperpage{53}, \hyperpage{62--64}, \hyperpage{81}, 
 		\hyperpage{87, 88}, \hyperpage{103}, \hyperpage{106}, 
 		\hyperpage{108, 109}
   \item steepness, \hyperpage{53}, \hyperpage{80}, \hyperpage{93}
   \item STOP, \hyperpage{88}
-  \item swell, \hyperpage{12}, \hyperpage{15}, \hyperpage{55}, 
-		\hyperpage{62}, \hyperpage{73, 74}, \hyperpage{78}, 
+  \item swell, \hyperpage{12}, \hyperpage{15}, \hyperpage{56}, 
+		\hyperpage{63}, \hyperpage{73, 74}, \hyperpage{78}, 
 		\hyperpage{89}
 
   \indexspace
@@ -5370,18 +5376,17 @@
   \indexspace
 
   \item unstructured, \hyperpage{3, 4}, \hyperpage{8, 9}, 
-		\hyperpage{12, 13}, \hyperpage{21}, \hyperpage{28, 29}, 
+		\hyperpage{12}, \hyperpage{21}, \hyperpage{29}, 
 		\hyperpage{32}, \hyperpage{34, 35}, \hyperpage{44}, 
-		\hyperpage{58}, \hyperpage{60}, \hyperpage{62}, 
-		\hyperpage{64}, \hyperpage{66}, \hyperpage{71}, 
-		\hyperpage{78}
+		\hyperpage{58}, \hyperpage{61--64}, \hyperpage{66}, 
+		\hyperpage{71}, \hyperpage{78}
 
   \indexspace
 
   \item WAM, \hyperpage{4--7}, \hyperpage{9, 10}, \hyperpage{15}, 
 		\hyperpage{18}, \hyperpage{22}, \hyperpage{37}, 
 		\hyperpage{46}, \hyperpage{48--50}, \hyperpage{54}, 
-		\hyperpage{81, 82}, \hyperpage{84}, \hyperpage{87}
+		\hyperpage{81, 82}, \hyperpage{84, 85}, \hyperpage{87}
   \item WAVEWATCH, \hyperpage{4--7}, \hyperpage{9, 10}, \hyperpage{15}, 
 		\hyperpage{22}, \hyperpage{46}, \hyperpage{49, 50}
   \item whitecapping, \hyperpage{7}, \hyperpage{16}, \hyperpage{53}, 
@@ -5391,9 +5396,9 @@
   \item wind, \hyperpage{3}, \hyperpage{5}, \hyperpage{7--12}, 
 		\hyperpage{14--19}, \hyperpage{21}, \hyperpage{26}, 
 		\hyperpage{34, 35}, \hyperpage{39}, \hyperpage{41}, 
-		\hyperpage{51--53}, \hyperpage{55}, \hyperpage{61, 62}, 
-		\hyperpage{64, 65}, \hyperpage{79}, \hyperpage{86, 87}, 
-		\hyperpage{92}, \hyperpage{94}
+		\hyperpage{51--53}, \hyperpage{56}, \hyperpage{61}, 
+		\hyperpage{63}, \hyperpage{65}, \hyperpage{79}, 
+		\hyperpage{86, 87}, \hyperpage{92}, \hyperpage{94}
 
 \end{theindex}
 
--- SwanVertlist.ftn90	2009-08-05 23:35:15.000000000 +0200
+++ SwanVertlist.ftn90	2009-08-05 21:33:21.000000000 +0200
@@ -32,10 +32,12 @@
 !   Authors
 !
 !   40.80: Marcel Zijlema
+!   41.07: Casey Dietrich
 !
 !   Updates
 !
 !   40.80, July 2007: New subroutine
+!   41.07, July 2009: small fix (assign ref.point to deepest point in case of no b.c.)
 !
 !   Purpose
 !
@@ -43,11 +45,13 @@
 !
 !   Method
 !
-!   Sorting based on distance and with respect to vertices where boundary condition is given
+!   Sorting based on distance with respect to a reference point
+!   This reference point can be either a vertex with boundary condition or a deepest point
 !
 !   Modules used
 !
     use ocpcomm4
+    use m_genarr, only: DEPTH
     use SwanGriddata
     use SwanGridobjects
     use SwanCompdata
@@ -67,6 +71,7 @@
     !
     real                            :: d1         ! distance of a point to origin
     real                            :: d2         ! distance of another point to origin
+    real                            :: depmax     ! maximum depth found
     real                            :: rtmp       ! temporary stored real for swapping
     real                            :: x0         ! x-coordinate of reference point
     real                            :: y0         ! y-coordinate of reference point
@@ -93,11 +98,34 @@
        return
     endif
     !
-    ! determine reference point nearest to the origin
-    ! this point is one of the vertices where boundary condition is given
+    ! determine reference point
     !
-    kx = minloc(vert(:)%attr(VERTX), vert(:)%atti(VBC)/=0)
-    ky = minloc(vert(:)%attr(VERTY), vert(:)%atti(VBC)/=0)
+    if ( all(mask=vert(:)%atti(VBC)==0) ) then
+       !
+       ! if no boundary condition is given then find the vertex with the maximum depth
+       !
+       depmax = -999.
+       !
+       do j = 1, nverts
+          !
+          if ( DEPTH(j) > depmax ) then
+             !
+             depmax = DEPTH(j)
+             kx(1)  = j
+             ky(1)  = j
+             !
+          endif
+          !
+       enddo
+       !
+    else
+       !
+       ! reference point is one of the vertices nearest to the origin where boundary condition is given
+       !
+       kx = minloc(vert(:)%attr(VERTX), vert(:)%atti(VBC)/=0)
+       ky = minloc(vert(:)%attr(VERTY), vert(:)%atti(VBC)/=0)
+       !
+    endif
     !
     if ( kx(1) == ky(1) ) then
        x0 = vert(kx(1))%attr(VERTX)
--- sweepcur.eps	1970-01-01 01:00:00.000000000 +0100
+++ sweepcur.eps	2009-08-05 21:33:07.000000000 +0200
@@ -0,0 +1,225 @@
+%!PS-Adobe-2.0 EPSF-2.0
+%%Title: sweepcur.fig
+%%Creator: fig2dev Version 3.2 Patchlevel 4
+%%CreationDate: Mon Jul 20 16:50:14 2009
+%%For: mzijlema@TUD11321 (Marcel Zijlema)
+%%BoundingBox: 0 0 623 468
+%%Magnification: 1.0000
+%%EndComments
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
+newpath 0 468 moveto 0 0 lineto 623 0 lineto 623 468 lineto closepath clip newpath
+-57.6 551.6 translate
+1 -1 scale
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
+n 8984 3038 652 652 0 360 DrawEllipse gs col0 s gr
+
+% Ellipse
+n 2437 3780 652 652 0 360 DrawEllipse gs col0 s gr
+
+% Ellipse
+n 9780 6532 660 660 0 360 DrawEllipse gs col7 0.50 shd ef gr gs col0 s gr
+
+% Ellipse
+n 3750 6727 652 652 0 360 DrawEllipse gs col0 s gr
+
+% Polyline
+n 3255 4470 m
+ 3255 1560 l gs col0 s gr 
+% Polyline
+n 4725 3000 m
+ 1815 3000 l gs col0 s gr 
+% Polyline
+n 8985 4530 m
+ 8985 1620 l gs col0 s gr 
+% Polyline
+n 10440 3045 m
+ 7530 3045 l gs col0 s gr 
+% Polyline
+n 9000 8745 m
+ 9000 5835 l gs col0 s gr 
+% Polyline
+n 10500 7290 m
+ 7590 7290 l gs col0 s gr 
+% Polyline
+n 3225 8670 m
+ 3225 5760 l gs col0 s gr 
+% Polyline
+n 4770 7260 m
+ 1860 7260 l gs col0 s gr 
+% Polyline
+n 9640 3045 m 9629 2923 l 9593 2784 l 9534 2678 l 9469 2595 l 9404 2533 l
+ 9309 2471 l 9200 2423 l 9111 2394 l 8984 2385 l 8984 3045 l
+ 9637 3045 l
+ cp gs col7 0.50 shd ef gr gs col0 s gr 
+% Polyline
+n 4141 7254 m 4235 7171 l 4289 7103 l 4339 7020 l 4377 6922 l 4398 6836 l
+ 4407 6757 l 4398 6632 l 4389 6570 l 4363 6484 l 4312 6390 l
+ 4241 6292 l 4126 6194 l 4046 6144 l 3978 6118 l 3898 6094 l
+ 3827 6079 l 3747 6073 l 3635 6082 l 3507 6120 l 3410 6171 l
+ 3345 6215 l 3283 6271 l 3229 6334 l 3735 6718 l 4138 7251 l
+
+ cp gs col7 0.50 shd ef gr gs col0 s gr 
+% Polyline
+n 3735 6715 m 3223 7112 l 3247 7144 l 3277 7174 l 3309 7206 l 3339 7233 l
+ 3368 7257 l
+ cp gs col7 0.50 shd ef gr gs col0 s gr 
+/Times-Roman ff 375.00 scf sf
+2700 1575 m
+gs 1 -1 sc (C) col0 sh gr
+/Times-Roman ff 375.00 scf sf
+8460 1680 m
+gs 1 -1 sc (C) col0 sh gr
+/Times-Roman ff 375.00 scf sf
+8490 5745 m
+gs 1 -1 sc (C) col0 sh gr
+/Times-Roman ff 375.00 scf sf
+2670 5895 m
+gs 1 -1 sc (C) col0 sh gr
+/Times-Roman ff 375.00 scf sf
+4530 3480 m
+gs 1 -1 sc (C) col0 sh gr
+/Times-Roman ff 375.00 scf sf
+10335 3450 m
+gs 1 -1 sc (C) col0 sh gr
+/Times-Roman ff 375.00 scf sf
+10380 7710 m
+gs 1 -1 sc (C) col0 sh gr
+/Times-Roman ff 375.00 scf sf
+4650 7635 m
+gs 1 -1 sc (C) col0 sh gr
+/Times-Roman ff 300.00 scf sf
+4905 7725 m
+gs 1 -1 sc (x) col0 sh gr
+/Times-Roman ff 300.00 scf sf
+10650 7785 m
+gs 1 -1 sc (x) col0 sh gr
+/Times-Roman ff 300.00 scf sf
+10575 3510 m
+gs 1 -1 sc (x) col0 sh gr
+/Times-Roman ff 300.00 scf sf
+4785 3555 m
+gs 1 -1 sc (x) col0 sh gr
+/Times-Roman ff 300.00 scf sf
+2910 5940 m
+gs 1 -1 sc (y) col0 sh gr
+/Times-Roman ff 300.00 scf sf
+2940 1605 m
+gs 1 -1 sc (y) col0 sh gr
+/Times-Roman ff 300.00 scf sf
+8715 1710 m
+gs 1 -1 sc (y) col0 sh gr
+/Times-Roman ff 300.00 scf sf
+8745 5790 m
+gs 1 -1 sc (y) col0 sh gr
+/Times-Roman ff 450.00 scf sf
+975 1860 m
+gs 1 -1 sc (\(A\)) col0 sh gr
+/Times-Roman ff 450.00 scf sf
+6855 6180 m
+gs 1 -1 sc (\(D\)) col0 sh gr
+/Times-Roman ff 450.00 scf sf
+6855 1905 m
+gs 1 -1 sc (\(B\)) col0 sh gr
+/Times-Roman ff 450.00 scf sf
+915 6150 m
+gs 1 -1 sc (\(C\)) col0 sh gr
+% here ends figure;
+$F2psEnd
+rs
+showpage
--- swmod1.ftn	2009-08-05 23:35:15.000000000 +0200
+++ swmod1.ftn	2009-08-05 21:33:22.000000000 +0200
@@ -1608,13 +1608,13 @@
 !     MDISP  [   3] dimension for dissipation arrays                      40.67
 !     MGENR  [   1] dimension for generation arrays                       40.85
 !     MICMAX [  10] max. number of points in comput. stencil
-!     MNUMS  [  30] dimension of array PNUMS
+!     MNUMS  [  35] dimension of array PNUMS
 !     MQUAD  [  10] dimension of array PQUAD
 !     MREDS  [   2] dimension for redistribution arrays                   40.85
 !     MSETUP [   2] dimension of array PSETUP
 !     MSHAPE [   5] dimension of array PSHAPE
 !     MSPPAR [   5] dimension of array SPPARM
-!     MSURF  [  10] dimension of array PSURF
+!     MSURF  [  15] dimension of array PSURF                              41.06
 !     MTRNP  [   3] dimension for propagation arrays                      40.85
 !     MTRIAD [  10] dimension of array PTRIAD
 !     MWCAP  [  15] dimension of array PWCAP
@@ -1635,12 +1635,12 @@
       PARAMETER           (MREDS  =  2)                                   40.85
       PARAMETER           (MTRNP  =  3)                                   40.85
       PARAMETER           (MICMAX = 10)
-      PARAMETER           (MNUMS  = 30)
+      PARAMETER           (MNUMS  = 35)
       PARAMETER           (MQUAD  = 10)
       PARAMETER           (MSETUP =  2)
       PARAMETER           (MSHAPE =  5)
       PARAMETER           (MSPPAR =  5)
-      PARAMETER           (MSURF  = 10)
+      PARAMETER           (MSURF  = 15)                                   41.06
       PARAMETER           (MTRIAD = 10)
       PARAMETER           (MWCAP  = 15)
       PARAMETER           (MWIND  = 40)
@@ -1947,6 +1947,9 @@
 ! IDIFFR [     0] diffraction method                                      40.21
 !                 0= no diffraction                                       40.21
 !                 1= diffraction                                          40.21
+! IFRSRF [     0] indicates frequency dependent surf breaking             41.06
+!                 0=no frequency dependency                               41.06
+!                 1=frequency dependency                                  41.06
 ! IGEN   [     3] indicates the generation mode
 !                 =1; for command GEN1 ...,
 !                 =2; for command GEN2 ...,
@@ -2321,6 +2324,7 @@
       INTEGER             MXITST,      MXITNS,                   NCOR
       INTEGER             NSTATC,      NSTATM,      NUMOBS,      NCOMPT   40.41
       INTEGER             IDIFFR                                          40.21
+      INTEGER             IFRSRF                                          41.06
       REAL                DEPMIN,      PBOT(MBOT),  PNUMS(MNUMS)
       REAL                PSETUP(MSETUP),           PSHAPE(MSHAPE)
       REAL                PSURF(MSURF),             PTRIAD(MTRIAD)
