doc: polyafit-manual.pdf

polyafit-manual.pdf:
	R CMD Rd2pdf --output=polyafit-manual.pdf --no-preview .

clean:
	rm -f polyafit-manual.pdf

clean-rdoc:
	rm -f man/*.Rd
