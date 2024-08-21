pdflatex project.tex
bibtex project
bibtex suppone
bibtex supptwo
pdflatex project.tex
pdflatex project.tex
evince project.pdf &

## Cleanup
rm *.aux
rm *.log
rm *.bbl
rm *.blg
