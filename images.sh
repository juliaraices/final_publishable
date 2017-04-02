for f in *.pdf; do
	convert -quality 150 ./"$f" -density 100 ./"${f%.pdf}.jpg"
done


for f in *.pdf; do
	convert -quality 150 ./"$f" -density 100 ./"${f%.pdf}.png"
done
