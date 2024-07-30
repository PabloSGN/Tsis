% mythesis v1.2.0 
% By G. Blazquez-Calero <gblazquez@iaa.es>

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{mythesis}[30/09/2023, My thesis class v1.2.0]
\RequirePackage[english]{babel}
% ifthen – Conditional commands in LaTeX documents
\RequirePackage{ifthen}
% calc – Simple arithmetic in LaTeX commands
\RequirePackage{calc}
% microtype – Subliminal refinements towards typographical perfection
% \AtEndOfClass{\RequirePackage{microtype}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{book}}
\ProcessOptions*
\LoadClass[a4paper]{book}
% seting level of numbering (default for "report" is 3). With ''-1'' you have non number also for chapters
\setcounter{secnumdepth}{4} 


%-------
% Fonts
%-------
% amsmath – AMS mathematical facilities for LaTeX
% amsfonts – TeX fonts from the American Mathematical Society
\RequirePackage{amsmath,amsfonts,amssymb}
\RequirePackage{libertine, newpxmath}%newtxmath}
% Other fonts:
% \RequirePackage{mathptmx}
% \RequirePackage{helvet}
% \RequirePackage{kpfonts}
% \RequirePackage{Xcharter}
% \RequirePackage{charter}
% \RequirePackage{fourier}
% \RequirePackage{fouriernc}
% \RequirePackage{tgschola}
% \RequirePackage{pxfonts, newpxmath}
% \RequirePackage[varg]{txfonts}
% \RequirePackage{courier}
% \RequirePackage{utopia}
% \RequirePackage{biolinum}
% \RequirePackage{iwona}
% \RequirePackage{paratype}
% \RequirePackage{montserrat}
% \RequirePackage{lmodern}
% \RequirePackage{tgpagella}
% \RequirePackage{fbb}
% \RequirePackage{beraserif}
% \RequirePackage{chancery}
% \RequirePackage{mathpazo}
% \RequirePackage[defaultfam]{courier}

%----------------
% Useful packages
%----------------
% ifpdf – Provides the \ifpdf conditional
\RequirePackage{ifpdf}
% latex-graphics – The LaTeX standard graphics bundle
\RequirePackage{graphicx,xcolor}
% Some color definition
\definecolor{myorange}{rgb}{0.7,0.4,0}
\definecolor{mygreen}{rgb}{0.1,0.4,0.1}
\definecolor{myblue}{rgb}{0.00,0.1,0.6}
\definecolor{myred}{rgb}{0.40,0.0,0.0}

% booktabs – Publication quality tables in LaTeX
% \RequirePackage{booktabs}
% geometry – Flexible and complete interface to document dimensions
\RequirePackage[%
	bottom=40mm,
	top=35mm,
 	inner=25mm,
 	outer=30mm,
	headheight=14.5pt,
%	textheight=550pt,
%	footskip=0pt,
]{geometry}
\RequirePackage[
	% labelfont={bf,sf},%
        % labelsep=period,%
        % justification=raggedright
]{caption}

% \usepackage{deluxetable}
%hyperref – Extensive support for hypertext in LaTeX
% \RequirePackage[colorlinks=true, allcolors=blue]{hyperref}
\RequirePackage[colorlinks=true, allcolors=black]{hyperref}
% fontenc – Standard package for selecting font encodings
\RequirePackage[T1]{fontenc}
% inputenc – Accept different input encodings
\RequirePackage[utf8]{inputenc}
% lineno – Line numbers on paragraphs
\RequirePackage{lineno}
% enumitem – Control layout of itemize, enumerate, description
\RequirePackage{enumitem}
% Removes spacing between items in lists (enumerate, itemize)
%\setlist{nolistsep} 
% layout – View the layout of a document
\RequirePackage{layout}

%Footnotes
\usepackage[symbol]{footmisc}

%-------------------------
% Headers and footers
%-------------------------
\RequirePackage{fancyhdr}  % custom headers/footers
\pagestyle{fancy}          % Enables the custom headers/footers
\fancyfoot{}
\fancyhead[LE,RO]{\textbf{\thepage}}
\fancyhead[ER]{\nouppercase{\rightmark}}
\fancyhead[OL]{\nouppercase{\leftmark}}


%------------------------------
% Section and subsections setup 
%------------------------------
\RequirePackage{titlesec}
% \titleformat{\chapter}[display]{
%  \normalfont\scshape\huge
% }{% Prefix
% 	\centering\Large\chaptertitlename \ \thechapter
% }{% vertical space in between
%  -15.0pt
% }{% The name of the chapter
% \centering\rule{0.5\linewidth}{1.0pt}\\\vspace{10pt}\centering
% }

% No rule and only chapter number
% \titleformat{\chapter}[display]{
%  \normalfont\scshape\huge
% }{% Prefix
%  \Huge\thechapter\centering
% }{% vertical space in between
% 0pt
% }{% The name of the chapter
% \vspace{5pt}\centering
% }[\vspace{25pt}]

% Chapter and the chapter title with rule
\titleformat{\chapter}[display]{
 \normalfont\scshape\huge
}{% Prefix
 \centering\Large\bfseries\chaptertitlename~\thechapter\\[-5pt]\rule{0.5\linewidth}{1.5pt}
% \Huge\thechapter\centering
}{% vertical space in between
0pt
}{% The name of the chapter
\vspace{15pt}\centering
}[\vspace{25pt}]


\titlespacing{\chapter}{0pt}{0pt}{50pt}
% \titlespacing\numberless{0pt}{14pt}{4pt}
% \titlespacing\section{0pt}{14pt}{4pt}
% \titlespacing\subsection{0pt}{14pt}{4pt}
% \titlespacing\subsubsection{0pt}{14pt}{4pt}
% \titleformat*{\section}{\bfseries\scshape}
% \titleformat*{\subsection}{\bfseries\scshape}
% \titleformat*{\subsubsection}{\bfseries\scshape}
% \titleformat{\chapter}[display]{%
%     \normalsize%
% }{% Prefix.
%     \centering\hspace{1pt}\small\scshape\color{darkgray}\chaptertitlename \ \thechapter%
% }{% Vertical space between.
%     -9pt%
% }{% Chapter name.
%     \centering\LARGE\bfseries\scshape%
% }

% \definecolor{ChapterColor}{RGB}{10,50,90}
% \RequirePackage{titlesec}
% \titleformat{\chapter}[display]
%   {\normalfont\fontsize{50}{0}\raggedleft\bfseries\color{ChapterColor}}{\thechapter}{20pt}{\Huge\raggedright\color{ChapterColor}}
% \titlespacing*{\chapter}{0pt}{10pt}{30pt} % the asterisc means no indentation
\titlespacing*{\section}{0pt}{30pt}{10pt}
\titlespacing*{\subsection}{0pt}{25pt}{10pt}


%------------------------
% Table of contents setup 
%------------------------
\usepackage{titletoc}
\newlength{\tocsep} 
% Sets the indentation of the sections in the table of contents
\setlength\tocsep{1.5pc} 
% Show only three levels in the table of contents section: sections, subsections and subsubsections
\setcounter{tocdepth}{3} 

%-----------
% Glossaries
%-----------
% \usepackage[nonumberlist,toc]{glossaries}
% Be careful because this glossary make use of supertable, which is incompatible with the deluxetable


%--------------------------------
% Set up of the bibliography style
%--------------------------------
\RequirePackage{natbib}
% My custom bibliography style, based on aajournal.bst
\bibliographystyle{mycustom}


%----------------------------------------------
% Journal Abrevations
%----------------------------------------------
\let\jnl@style=\rmfamily
\def\ref@jnl#1{{\jnl@style#1}}%
% style\newcommand\aj{\ref@jnl{AJ}}%        % Astronomical Journal
\newcommand\aj{\ref@jnl{Astron. J.}}%        % Astronomical Journal
% \newcommand\psj{\ref@jnl{PSJ}}%       % Planetary Science Journal
\newcommand\psj{\ref@jnl{Planetary Science Journal}}
% \newcommand\araa{\ref@jnl{ARA\&A}}%  % Annual Review of Astron and Astrophys
\newcommand\araa{\ref@jnl{Ann. Rev. Astron. Astrophys.}}%  % Annual Review of Astron and Astrophys
% \newcommand\apj{\ref@jnl{ApJ}}%    % Astrophysical Journal
\newcommand\apj{\ref@jnl{Astrophys. J.}}%    % Astrophysical Journal
% \newcommand\apjl{\ref@jnl{ApJL}}     % Astrophysical Journal, Letters
\newcommand\apjl{\ref@jnl{Astrophys. J. Lett.}}     % Astrophysical Journal, Letters
% \newcommand\apjs{\ref@jnl{ApJS}}%    % Astrophysical Journal, Supplement
\newcommand\apjs{\ref@jnl{Astrophys. J. Suppl.}}
\newcommand\ao{\ref@jnl{ApOpt}}%   % Applied Optics
\newcommand\apss{\ref@jnl{Ap\&SS}}%  % Astrophysics and Space Science
% \newcommand\aap{\ref@jnl{A\&A}}%     % Astronomy and Astrophysics
\newcommand\aap{\ref@jnl{Astron. Astrophys.}}%     % Astronomy and Astrophysics
% \newcommand\aapr{\ref@jnl{A\&A~Rv}}%  % Astronomy and Astrophysics Reviews
\newcommand\aapr{\ref@jnl{Astron. Astrophys. Rev.}}%  % Astronomy and Astrophysics Reviews
% \newcommand\aaps{\ref@jnl{A\&AS}}%    % Astronomy and Astrophysics, Supplement
\newcommand\aaps{\ref@jnl{Astron. Astrophys. Suppl.}}%    % Astronomy and Astrophysics, Supplement
\newcommand\azh{\ref@jnl{AZh}}%       % Astronomicheskii Zhurnal
\newcommand\baas{\ref@jnl{BAAS}}%     % Bulletin of the AAS
\newcommand\icarus{\ref@jnl{Icarus}}% % Icarus
\newcommand\jaavso{\ref@jnl{JAAVSO}}  % The Journal of the American Association of Variable Star Observers
\newcommand\jrasc{\ref@jnl{JRASC}}%   % Journal of the RAS of Canada
\newcommand\memras{\ref@jnl{MmRAS}}%  % Memoirs of the RAS
\newcommand\mnras{\ref@jnl{MNRAS}}%   % Monthly Notices of the RAS
\newcommand\pra{\ref@jnl{PhRvA}}% % Physical Review A: General Physics
\newcommand\prb{\ref@jnl{PhRvB}}% % Physical Review B: Solid State
\newcommand\prc{\ref@jnl{PhRvC}}% % Physical Review C
\newcommand\prd{\ref@jnl{PhRvD}}% % Physical Review D
\newcommand\pre{\ref@jnl{PhRvE}}% % Physical Review E
\newcommand\prl{\ref@jnl{PhRvL}}% % Physical Review Letters
\newcommand\pasp{\ref@jnl{PASP}}%     % Publications of the ASP
\newcommand\pasj{\ref@jnl{PASJ}}%     % Publications of the ASJ
\newcommand\qjras{\ref@jnl{QJRAS}}%   % Quarterly Journal of the RAS
\newcommand\skytel{\ref@jnl{S\&T}}%   % Sky and Telescope
\newcommand\solphys{\ref@jnl{SoPh}}% % Solar Physics
\newcommand\sovast{\ref@jnl{Soviet~Ast.}}% % Soviet Astronomy
\newcommand\ssr{\ref@jnl{SSRv}}% % Space Science Reviews
\newcommand\zap{\ref@jnl{ZA}}%       % Zeitschrift fuer Astrophysik
\newcommand\nat{\ref@jnl{Nature}}%  % Nature
\newcommand\iaucirc{\ref@jnl{IAUC}}% % IAU Cirulars
\newcommand\aplett{\ref@jnl{Astrophys.~Lett.}}%  % Astrophysics Letters
\newcommand\apspr{\ref@jnl{Astrophys.~Space~Phys.~Res.}}% % Astrophysics Space Physics Research
\newcommand\bain{\ref@jnl{BAN}}% % Bulletin Astronomical Institute of the Netherlands
\newcommand\fcp{\ref@jnl{FCPh}}%   % Fundamental Cosmic Physics
\newcommand\gca{\ref@jnl{GeoCoA}}% % Geochimica Cosmochimica Acta
\newcommand\grl{\ref@jnl{Geophys.~Res.~Lett.}}%  % Geophysics Research Letters
\newcommand\jcp{\ref@jnl{JChPh}}%     % Journal of Chemical Physics
\newcommand\jgr{\ref@jnl{J.~Geophys.~Res.}}%     % Journal of Geophysics Research
\newcommand\jqsrt{\ref@jnl{JQSRT}}%   % Journal of Quantitiative Spectroscopy and Radiative Trasfer
\newcommand\memsai{\ref@jnl{MmSAI}}% % Mem. Societa Astronomica Italiana
\newcommand\nphysa{\ref@jnl{NuPhA}}%     % Nuclear Physics A
\newcommand\physrep{\ref@jnl{PhR}}%       % Physics Reports
\newcommand\physscr{\ref@jnl{PhyS}}%        % Physica Scripta
\newcommand\planss{\ref@jnl{Planet.~Space~Sci.}}%  % Planetary Space Science
\newcommand\procspie{\ref@jnl{Proc.~SPIE}}%      % Proceedings of the SPIE
\newcommand\actaa{\ref@jnl{AcA}}%  % Acta Astronomica
\newcommand\caa{\ref@jnl{ChA\&A}}%  % Chinese Astronomy and Astrophysics
\newcommand\cjaa{\ref@jnl{ChJA\&A}}%  % Chinese Journal of Astronomy and Astrophysics
\newcommand\jcap{\ref@jnl{JCAP}}%  % Journal of Cosmology and Astroparticle Physics
\newcommand\na{\ref@jnl{New. Astron.}}%  % New Astronomy
\newcommand\nar{\ref@jnl{New. Astron. Rev.}}%  % New Astronomy Review
\newcommand\pasa{\ref@jnl{PASA}}%  % Publications of the Astron. Soc. of Australia
\newcommand\rmxaa{\ref@jnl{Revista Mexicana de Astronomia y Astrofisica}}%  % Revista Mexicana de Astronomia y Astrofisica
\newcommand\maps{\ref@jnl{M\&PS}}% Meteoritics and Planetary Science
\newcommand\aas{\ref@jnl{AAS Meeting Abstracts}}% American Astronomical Society Meeting Abstracts
\newcommand\dps{\ref@jnl{AAS/DPS Meeting Abstracts}}% American Astronomical Society/Division for Planetary Sciences Meeting Abstracts
\let\astap=\aap
\let\apjlett=\apjl
\let\apjsupp=\apjs
\let\applopt=\ao
% \newcommand\ulap[1]{\vbox\@to\z@{{\vss#1}}}%
% \newcommand\dlap[1]{\vbox\@to\z@{{#1\vss}}}%


%----------------------
% Useful commands
%----------------------
\newcommand{\thco}{\mbox{$^{13}$CO}}
\newcommand{\twco}{\mbox{$^{12}$CO}}
\newcommand{\coonezero}{\mbox{CO(1-0)}}
\newcommand{\cotwoone}{\mbox{CO(2-1)}}
\newcommand{\cothreetwo}{\mbox{CO(3-2)}}
\newcommand{\kmspc}{\mbox{km s$^{-1}$ pc$^{-1}$}}
\newcommand{\kms}{\mbox{km s$^{-1}$}}
\newcommand{\gcmcu}{\mbox{g cm$^{-3}$}}
\newcommand{\gcmsq}{\mbox{g cm$^{-2}$}}
\newcommand{\cmcu}{\mbox{cm$^{-3}$}}
\newcommand{\cmsq}{\mbox{cm$^{-2}$}}
\newcommand{\htp}{\mbox{H$_3^+$}}
\newcommand{\lsun}{\mbox{$L_{\sun}$}}
\newcommand{\msun}{\mbox{$M_{\sun}$}}
\newcommand{\msunyr}{\mbox{$M_{\sun}$ yr$^{-1}$}}
\newcommand{\mearth}{\mbox{$M_{\earth}$}}
\newcommand{\tbol}{\mbox{T$_{bol}$}}
\newcommand{\lbol}{\mbox{L$_{bol}$}}
\newcommand{\masyr}{\mbox{mas yr$^{-1}$}}
\newcommand{\jybeamkms}{\mbox{Jy beam$^{-1}$ \kms}}
\newcommand\farcs{\mbox{$.\!\!^{\prime\prime}$}}%
\newcommand\arcsec{\mbox{$^{\prime\prime}$}}% 
\newcommand\arcmin{\mbox{$^\prime$}}% 
\newcommand\fs{\mbox{$.\!\!^{\mathrm s}$}}% 
\let\la=\lesssim
\let\ga=\gtrsim 

\newcommand*{\Myitem}{\item[{\includegraphics[width=.3cm]{images/sun_item_1.png}}]}

%-----------
% Title page 
%-----------
\def\advisor#1{\gdef\@advisor{#1}}
\def\@advisor{\@latex@warning@no@line{No \noexpand\advisor given}}
\newcommand{\institution}[1]{\gdef\@institution{#1}}
\newcommand{\@institution}{\@latex@warning@no@line{No \noexpand\institution given}}
\newcommand{\logo}[1]{\gdef\@logo{#1}}
\newcommand{\@logo}{\@latex@warning@no@line{No \noexpand\logo given}}

\def\maketitle{%
 \thispagestyle{empty}
 \newgeometry{top=4cm,bottom=-0cm}
 \begin{center}
  \vspace{1.5cm}
  {\bfseries\LARGE \@title}\\
  \vspace{2cm}
  PhD dissertation by\\
  \vspace{0.25cm}
  {\large\bf \@author}\\
  \vspace{1cm}
  \@institution
  \vspace{1.5cm}
  \noindent
  A thesis submitted in fulfillment\\
  of the requirements of the degree of\\
  Doctor of Philosophy\\
  \vspace{1cm}
  \@date\\
  \vspace{1cm}
  PhD thesis supervised by\\
  \vspace{0.5cm}
  \@advisor
  \vspace{4.5cm}
  \@logo
 \end{center}
 \restoregeometry
}


%----
% Deluxetable. Downloaded from https://fits.gsfc.nasa.gov/standard30/deluxetable.sty
%----

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% deluxetable.sty:  Standalone definitions of aastex's deluxetable. 
% 
% Credit: This file kindly provided by Alexey Vikhlinin. 
% 
% USAGE: 
% 
% 1. Put \usepackage{deluxetable} in the preamble of your document. 
% 2. Copy this text to a local file called deluxetable.sty. 
% 
\NeedsTeXFormat{LaTeX2e}[1995/12/01]% 
\ProvidesPackage{deluxetable}% 

\newcommand\nodata{ ~$\cdots$~ }% 
\newcommand\ulap[1]{\vbox\@to\z@{{\vss#1}}}% 
\newcommand\dlap[1]{\vbox\@to\z@{{#1\vss}}}% 

\newcommand\let@environment[2]{% 
 \expandafter\let 
 \csname#1\expandafter\endcsname\csname#2\endcsname 
 \expandafter\let 
 \csname end#1\expandafter\endcsname\csname end#2\endcsname 
 }% 

\def\substitute@command#1#2{% 
 \ClassWarning{aastex}{% 
  Command \protect#1\space is deprecated in aastex. 
  Using \protect#2\space instead (please fix your document). 
 }% 
 #2% 
}% 
\def\remove@command#1{% 
 \ClassWarning{aastex}{% 
  Command \protect#1\space is deprecated in aastex. 
  Please remove it from your document% 
 }% 
}% 
\def\@provide#1{%Implicit #2 
 \@ifxundefined#1{\def#1}{\def\@tempa}% 
}% 
\def\prepdef#1#2{% 
 \@ifxundefined#1{\toks@{}}{\toks@\expandafter{#1}}% 
 \toks@ii{#2}% 
 \edef#1{\the\toks@ii\the\toks@}% 
}% 
\def\appdef#1#2{% 
 \@ifxundefined#1{\toks@{}}{\toks@\expandafter{#1}}% 
 \toks@ii{#2}% 
 \edef#1{\the\toks@\the\toks@ii}% 
}% 
\def\appgdef#1#2{% 
 \@ifxundefined#1{\toks@{}}{\toks@\expandafter{#1}}% 
 \toks@ii{#2}% 
 \xdef#1{\the\toks@\the\toks@ii}% 
}% 
\toksdef\toks@ii=\tw@ 
\long\def\@ifxundefined#1{\@ifx{\undefined#1}}% 
\def\@boolean#1#2{% 
  \long\def#1{% 
    #2% \if<something> 
      \expandafter\boolean@true 
    \else 
      \expandafter\boolean@false 
    \fi 
  }% 
}% 
\def\@boole@def#1#{\@boolean{#1}}% Implicit #2 
\@boole@def\@ifx#1{\ifx#1}% 
\@boole@def\@ifx@empty#1{\ifx\@empty#1}% 
\@boole@def\@ifnum#1{\ifnum#1}% 
\@boole@def\@ifdim#1{\ifdim#1}% 
\@boole@def\@ifvoid#1{\ifvoid#1}% 
\long\def\boolean@true#1#2{#1}% 
\long\def\boolean@false#1#2{#2}% 

\newcommand\tablenotemark[1]{% 
 \rlap{$^{\mathrm #1}$}% 
}% 
\def\@tablenotetext#1#2{% 
 \vspace{.5ex}% 
 {\noindent\llap{$^{#1}$}#2\par}% 
}% 
\def\@tablenotes#1{% 
 \par 
 \vspace{4.5ex}\footnoterule\vspace{.5ex}% 
 {\footnotesize#1}% 
}% 
\def\@tablenotes@pptt#1{% 
 \par 
 \vspace{3.2ex}\footnoterule\vspace{.4ex}% 
 {\footnotesize#1}% 
}% 
\AtBeginDocument{% 
 \let\tblnote@list\@empty 
}% 
\newcommand\tablenotetext[2]{% 
 \appgdef\tblnote@list{\@tablenotetext{#1}{#2}}% 
}% 
\def\spew@tblnotes{% 
 \@ifx@empty\tblnote@list{}{% 
  \@tablenotes{\tblnote@list}% 
  \global\let\tblnote@list\@empty 
 }% 
}% 
\prepdef\endtable{\spew@tblnotes}% 
\expandafter\prepdef\csname endtable*\endcsname{\spew@tblnotes}% 
\let\tableline=\hline 
\long\def\@makecaption#1#2{\vskip 2ex\noindent#1\ #2\par}% 
\newcommand\tablenum[1]{% 
 \def\thetable{#1}% 
 \let\@currentlabel\thetable 
 \addtocounter{table}{\m@ne}% 
}% 
\newcommand\dummytable{\refstepcounter{table}}% 
\newbox\pt@box 
\newdimen\pt@width 
\newcount\pt@line 
\newcount\pt@column 
\newcount\pt@nlines 
\newcount\pt@ncol 
\newcount\pt@page 
\newcommand\colhead[1]{\multicolumn{1}{c}{#1}\pt@addcol}% 
\newcommand\tabletypesize[1]{\def\pt@typesize{#1}}%AO-access to different type sizes in deluxetable 
\def\pt@typesize{}% 
\newcommand\tablecolumns[1]{% 
 \pt@column=#1\relax 
 \pt@ncol=#1\relax 
 \global\let\pt@addcol\@empty 
}% 
\newcommand\tablecaption[1]{\gdef\pt@caption{#1}}% 
\newcommand\tablehead[1]{% 
 \gdef\pt@head{% 
  \hline\hline 
  \relax\\[-1.7ex]% 
  #1\hskip\tabcolsep 
  \\[.7ex]% 
  \hline 
  \relax\\[-1.5ex]% 
 }% 
}% 
\def\tablehead@aj#1{% 
 \gdef\pt@head{% 
  #1\hskip\tabcolsep 
  \\[.7ex]% 
  \hline 
  \relax\\[-1.5ex]% 
 }% 
}% 
\newif\if@pt@rot 
\newcommand\rotate{\@pt@rottrue}% 
\newcommand\tabletail[1]{\gdef\pt@tail{#1}}% 
\newcommand\tablewidth[1]{\pt@width=#1\relax}% 
\newcommand\tableheadfrac[1]{\gdef\pt@headfrac{#1}}% 
\AtBeginDocument{% 
 \let\pt@caption\@empty 
 \let\pt@head\@empty 
 \let\pt@tail\@empty 
 \pt@width\textwidth 
 \def\pt@headfrac{.1}% 
}% 
\def\pt@calcnlines{% 
 \begingroup 
  \if@pt@rot\textheight\textwidth\fi 
  \pt@typesize%AO-type sizes in deluxetable 
  \@tempdima\pt@headfrac\textheight 
  \@tempdimb\textheight\advance\@tempdimb by-\@tempdima 
  \@tempdima\arraystretch\baselineskip 
  \divide\@tempdimb by\@tempdima 
  \global\pt@nlines\@tempdimb 
 \endgroup 
}% 
\def\pt@tabacol{% 
 \edef\@preamble{\@preamble\hskip\tabcolsep\tabskip\fill}% 
}% 
\newdimen\pt@tmpcapwidth 
\long\def\@makecaption@plano#1#2{% 
 \@ifdim{\pt@width>\z@}{% 
  \pt@tmpcapwidth\pt@width 
 }{% 
  \pt@tmpcapwidth\textwidth 
 }% 
 \center 
  \parbox{\pt@tmpcapwidth}{% 
   \@make@caption@text{#1}{#2}% 
  }% 
 \endcenter 
}% 
\def\fnum@ptable{Table \thetable}% 
\def\fnum@ptablecont{Table \thetable \,\,  {\rmfamily (cont'd)}}% 
\long\def\@make@caption@text#1#2{% 
 \center\rmfamily#1.\quad#2\endcenter 
}% 
\long\def\@makecaption@plano@cont#1#2{% 
 \center\rmfamily#1\endcenter 
 \vskip 2.5ex 
}% 
\def\fnum@ptablecont@pptt{Table \thetable \,\,  {\itshape (cont'd)}}% 
\long\def\@make@caption@text@pptt#1#2{% 
 \center\sc#1\\[.5ex]#2\endcenter 
}% 
\long\def\@makecaption@plano@cont@pptt#1#2{% 
 \center\sc#1\endcenter 
 \vskip 1.5ex 
}% 
\def\fnum@ptable@apj{TABLE \thetable}% 
\def\fnum@ptablecont@apj{TABLE \thetable \,\,  {\itshape (cont'd)}}% 
\long\def\@make@caption@text@apj#1#2{% 
 \center\Large\sc#1\\[.5ex]#2\endcenter 
}% 
\long\def\@makecaption@plano@cont@apj#1#2{% 
 \center\Large\sc#1\endcenter 
 \vskip 1.5ex 
}% 
\def\fnum@ptable@aj{Table \thetable.}% 
\def\fnum@ptablecont@aj{Table \thetable.\space{\rmfamily (continued)}}% 
\long\def\@make@caption@text@aj#1#2{% 
 \center\large{\sc#1}\space#2\endcenter 
}% 
\long\def\@makecaption@plano@cont@aj#1#2{% 
 \center\large{\sc#1}\endcenter 
 \vskip 1.5ex 
}% 
\newenvironment{deluxetable}[1]{% 
 \def\pt@format{\string#1}% 
 \let\@acol\pt@tabacol 
 \let\@tabularcr\@ptabularcr 
 \let\@tablenotetext\@ptablenotetext 
 \let\@tablenotes\@ptablenotes 
 \global\pt@ncol\z@ 
 \global\pt@column\z@ 
 \global\pt@page\@ne 
 \gdef\pt@addcol{\global\advance\pt@ncol\@ne}% 
}{% 
 \spew@tblnotes 
 \typeout@deluxetable 
 \endcenter 
 \end@plano@float 
}% 
\let@environment{planotable}{deluxetable}% 
\def\@plano@float{% Invoked by \startdata 
 \begingroup 
  \if@pt@rot\columnwidth\textheight\fi %AO-deluxetable \rotate 
  \@plano@float@{table}% 
}% 
\def\end@plano@float{% S/B invoked by \enddata; instead by \enddeluxetable 
  \end@plano@float@ 
 \endgroup 
}% 
\def\@plano@float@{\@float}%AO-extra layer of abstraction for float processing 
\def\end@plano@float@{\end@float}% 
\appdef\@endfloatbox{% 
 \if@pt@rot\global\setbox\@currbox\vbox{\rotatebox{90}{\box\@currbox}}\fi 
}% 
\def\typeout@deluxetable@mss{% 
 \typeout{% 
  Table \thetable\space has been set to width \the\pt@width 
 }% 
}% 
\def\typeout@deluxetable@ppt{% 
 \typeout{% 
  Page  \the\pt@page \space of table \thetable\space has been set to 
  width \the\pt@width\space with \the\pt@nlines\space lines per page 
 }% 
}% 
\let\typeout@deluxetable\typeout@deluxetable@mss 
\newcommand\startdata{% 
 \pt@calcnlines 
 \@ifdim{\pt@width>\z@}{% 
  \def\@halignto{\@to\pt@width}% 
 }{% 
  \def\@halignto{}% 
 }% 
 \let\fnum@table=\fnum@ptable 
 \let\@makecaption\@makecaption@plano 
 \global\pt@line\z@ 
 \start@pt@tabular 
 \after@startdata 
}% 
\def\start@pt@tabular{% 
 \@plano@float 
 \center 
 \expandafter\caption\expandafter{\pt@caption}% 
 \pt@typesize%AO-type sizes in deluxetable 
 \expandafter\@tabular\expandafter{\pt@format}% 
}% 
\def\set@pt@box#1{\setbox\pt@box}% 
\def\after@startdata{\pt@head}% 
\def\after@startdata@aj{\hline\hline\relax\\[-1.7ex]\pt@head}% 
\def\enddata{% 
 \crcr 
 \noalign{\vskip .7ex}% 
 \before@enddata 
 \endtabular 
 \setbox\pt@box\lastbox 
 \pt@width\wd\pt@box\box\pt@box 
}% 
\def\before@enddata{\hline}% 
\def\before@enddata@aj{\hline\hline}% 
\def\nl{\substitute@command\nl\\}% 
\def\nextline{\substitute@command\nextline\\}% 
\def\@ptabularcr{% 
 {\ifnum0=`}\fi % A klootch just in case the next token is & or \cr 
 \@ifstar{% 
  \@testopt{\@argptabularcr\@empty}\z@skip 
 }{% 
  \@testopt{\@argptabularcr\ptable@@split}\z@skip 
 }% 
}% 
\def\@argptabularcr#1[#2]{% 
 \ifnum0=`{\fi}% To undo the effect of the klootch. 
 \@ifdim{#2>\z@}{% 
  \unskip\@xargarraycr{#2}% 
 }{% 
  \@yargarraycr{#2}% 
 }\ptable@split#1% 
}% 
\def\ptable@split#1{% 
 \noalign{% 
  \global\advance\pt@line\@ne 
  \@ifnum{\pt@line<\pt@nlines}{}{% 
   \aftergroup#1% 
  }% 
 }% 
}% 
\def\ptable@@split{% 
  \before@suspendpt 
  \endtabular 
  \setbox\pt@box\lastbox 
  \pt@width\wd\pt@box\box\pt@box 
  \typeout@pt@nl 
  \global\advance\pt@page\@ne 
  \endcenter 
  \end@plano@float 
  \clearpage 
  \addtocounter{table}{\m@ne}% 
  \let\fnum@table=\fnum@ptablecont 
  \let\@makecaption\@makecaption@plano@cont 
  \global\pt@ncol=\pt@column%  Either 0 or value of \tablecolumns 
  \global\pt@line\z@ 
  \start@pt@tabular 
  \before@resumept 
  \pt@head 
}% 
\def\before@suspendpt{}% 
\def\before@resumept{}% 
\def\before@suspendpt@aj{\@tabularcr\noalign{\vskip .7ex}\hline}% 
\def\before@resumept@aj{\hline\relax\@tabularcr[-1.7ex]}% 
\def\typeout@pt@nl@ppt{% 
 \typeout{% 
  Page \the\pt@page\space of table \thetable\space has been set to 
  width \the\pt@width 
 }% 
}% 
\def\typeout@pt@nl@mss{% 
 \typeout{% 
  Page \the\pt@page\space of table \thetable\space has been set to 
  width \the\pt@width\space with \the\pt@nlines\space lines per page 
 }% 
}% 
\def\typeout@pt@nl@aj{% 
 \typeout{% 
  Table \thetable\space has been set to 
  width \the\pt@width\space with \the\pt@nlines\space lines per page 
 }% 
}% 
\let\typeout@pt@nl\typeout@pt@nl@mss 
\newcommand\tablevspace[1]{\substitute@command\tablevspace\\[#1]}% 
\newcommand\tablebreak{\substitute@command\tablebreak\\}%{\cr\ptable@@split}% 
\newcommand\cutinhead[1]{% 
 \noalign{\vskip 1.5ex}% 
 \hline 
 \@ptabularcr 
 \noalign{\vskip -4ex}% 
 \multicolumn{\pt@ncol}{c}{#1}% 
 \@ptabularcr 
 \noalign{\vskip .8ex}% 
 \hline 
 \@ptabularcr 
 \noalign{\vskip -2ex}% 
}% 
\def\cutinhead@ppt#1{% 
 \noalign{\vskip 1.5ex}% 
 \hline 
 \@ptabularcr 
 \noalign{\vskip -2ex}% Style Note: in apj, it is -1.5ex 
 \multicolumn{\pt@ncol}{c}{#1}% 
 \@ptabularcr 
 \noalign{\vskip .8ex}% 
 \hline 
 \@ptabularcr 
 \noalign{\vskip -2ex}% 
}% 
\newcommand\sidehead[1]{% 
 \noalign{\vskip 1.5ex}% 
 \multicolumn{\pt@ncol}{@{\hskip\z@}l}{#1}% 
 \@ptabularcr 
 \noalign{\vskip .5ex}% 
}% 
\def\@ptablenotetext#1#2{% 
 \vspace{0ex}% Style Note: in ppt, it is gone 
 {\parbox{\pt@width}{\hskip1em$^{\mathrm#1}$#2}\par}% 
}% 
\def\@ptablenotes#1{% 
 \par 
 \vspace{2ex}% 
 {\setlength\parskip{1.5ex}#1}% 
}% 
\def\@ptablenotes@apj#1{% 
 \par 
 \vspace{2ex}% 
 {#1}% 
}% 
\newcommand\tablerefs[1]{\appgdef\tblnote@list{\@tableref{#1}}}% 
\def\@tableref#1{% 
 \par 
 \vspace*{3ex}% 
 {\parbox{\pt@width}{\hskip1em\rmfamily References. --- #1}\par}% 
}% 
\newcommand\tablecomments[1]{\appgdef\tblnote@list{\@tablecom{#1}}}% 
\def\@tablecom#1{% 
 \par 
 \vspace*{3ex}% 
 {\parbox{\pt@width}{\hskip1em\rmfamily Note. --- #1}\par}% 
}% 
\def\@tableref@pptt#1{% 
 \par 
 \vspace*{3ex}{% 
 \parbox{\pt@width}{\hskip1em{\sc References.---}#1}\par}% 
}% 
\def\@tablecom@pptt#1{% 
 \vspace*{3ex}{% 
 \parbox{\pt@width}{\hskip1em{\sc Note.---}#1}\par}% 
}% 
\endinput 
%% 
%% ................................... End of file `deluxetable.sty'.
