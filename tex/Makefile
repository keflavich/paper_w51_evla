# latex Makefile
texpath=/usr/texbin/
PDFLATEX=${texpath}pdflatex -halt-on-error -synctex=1 --interaction=nonstopmode
SKIPERR=${texpath}pdflatex --interaction=nonstopmode
LATEX=${PDFLATEX}
BIBTEX=bibtex
DVIPS=dvips
PS2PDF=ps2pdf

all: w51_evla_h2co

w51_evla_h2co: 
	@rm -f w51_evla_h2co*.aux w51_evla_h2co*.bbl w51_evla_h2co*.blg w51_evla_h2co*.dvi w51_evla_h2co*.log w51_evla_h2co*.lot w51_evla_h2co*.lof
	${LATEX} w51_evla_h2co.tex
	${BIBTEX} w51_evla_h2co
	${LATEX} w51_evla_h2co.tex
	${BIBTEX} w51_evla_h2co
	${LATEX} w51_evla_h2co.tex
	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=w51_evla_h2co_compressed.pdf w51_evla_h2co.pdf

diff:
	python parse_macros.py w51_evla_h2co.tex w51_evla_h2co_fordiff.tex
	python parse_macros.py w51_evla_h2co_jun7.tex w51_evla_h2co_jun7_fordiff.tex
	latexdiff w51_evla_h2co_jun7_fordiff.tex w51_evla_h2co_fordiff.tex > diff.tex
	${SKIPERR} diff.tex
	${BIBTEX} diff
	${SKIPERR} diff.tex
	${BIBTEX} diff
	${SKIPERR} diff.tex
	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=diff_compressed.pdf diff.pdf
