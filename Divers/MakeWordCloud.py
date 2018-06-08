# usage: python MakeWordCloud.py file_in file_out
# file_in: Word file
# file_out: Figure

from docx import Document
from wordcloud import WordCloud
import matplotlib.pyplot as plt
import sys

# read input file
file=sys.argv[1]
file_out=sys.argv[2]

# read document
document=Document(file)
text=''
for par in document.paragraphs:
    text=text+par.text+' '
    
    
# make word cloud
wordcloud=WordCloud(background_color='white').generate(text)
 
# figure + printing
plt.imshow(wordcloud)
plt.axis('off')
plt.savefig(file_out)