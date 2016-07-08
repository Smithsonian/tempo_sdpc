
PS_FILES = area_check.ps
PNG_FILES = $(patsubst %.ps,%.png,$(PS_FILES))
PDF_FILES = $(patsubst %.ps,%.pdf,$(PS_FILES))

ps-files:	$(PS_FILES)
pngs:	$(PNG_FILES)
pdfs:	$(PDF_FILES)

area_check.ps:
	./area_check.sl ../../tests

%.pdf:	%.ps
	ps2pdf $<
	pdfcrop $@ _tmp.pdf
	gs -dCompatibilityLevel=1.4 \
	    -dPDFSETTINGS=/screen \
	    -dCompressFonts=true \
	    -dSubsetFonts=true \
	    -dNOPAUSE \
	    -dBATCH \
	    -sDEVICE=pdfwrite \
	    -sOutputFile=$@ \
	    -c ".setpdfwrite <</NeverEmbed [ ]>> setdistillerparams" \
	    -f _tmp.pdf

%.png:	%.pdf
	convert -verbose -background white -flatten -density 300 -trim $< -quality 100 $@

clean:
	-/bin/rm $(PS_FILES) $(PNG_FILES) $(PDF_FILES) _tmp.pdf
