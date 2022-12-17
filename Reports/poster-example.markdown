---
title: Thorough research of this special topic regarding the influence of various factors
author: [Poster Author 1, Poster Author 2]
email: author@university.edu
institute: Institute, Dept., University
longinstitute: Institute, Department, University, Country
web: 'university.edu'
biblio-files: parsed-references.bib
posteroptions: width=90,height=116.3,scale=1.2 #,grid
headerheight: 13cm
# large, Large, LARGE, huge, Huge, veryHuge, VeryHuge, VERYHuge
titlefont: size=\veryHuge,series=\bfseries
authorfont: size=\huge
institutefont: size=\Large
captionfont: size=\small
knit: (function(input, encoding, make = TRUE) { source('templates/makefile-renderer.R', local = TRUE) })
---

%% smart
%% to=latex
%% template=templates/poster.tex
%% filter=templates/poster-filters.py
%% biblatex

[columns=2]

[column]

# Introduction

### Graphs

![two figures side-by-side]({width=0.5\linewidth,trim=0 0 0 1cm}presentation-examplefig,{width=0.5\linewidth,trim=-1cm 0 0 0}presentation-examplefig-magenta)

<!-- Comments -->
### Default lists

- Citations [@Macherey2006] and @Macherey2006
- references have a clickable link to Pubmed or Amazon
- Standard abreviations \\eg and \\ie for \eg and \ie
- Units like \pps{900}
- **Highlights** and *highlights*

### Numbered lists

1.  First paragraph
2.  Second paragraph
3.  Third paragraph

    Continued paragraphs

# Bla

\lipsum[1]

# Blub

[columns=2]

[column=0.67]

### Columns within columns within blocks within columns

[columns=2]

[column=0.59]

\lipsum[3]

[column]

\lipsum[4]

[/columns]

[column]

\lipsum[5]

[/columns]

\vskip0.7cm

[column]

# Big figure

![electrodes!]({width=.8\linewidth}presentation-examplefig-electrodes)

# Baz

\lipsum[6-7]

### Table

<!-- this is still latex :-) -->
\begin{table}
    \rowcolors{2}{kuldark!10}{kuldark!20}
    \begin{tabular}{lrrrllll}
            \rowcolor{kuldark!20}
                &     &                     &         &      &          &
                \multicolumn{2}{c}{\cellcolor{kuldark!20}Blub} \\
        Bla & Blub & Bla & Blub & Bla & Blub &
        Bla & Blub \\
        Bla & Blub & Bla & Blub & Bla & Blub & Blablublbaba & Blubblabalbal \\
        Blub & Bla & Blub & Bla & Blub & Bla & Blub & Bla \\
        Bla & Blub & Bla & Blub & Bla & Blub & Bla & Blub \\
        Blub & Bla & Blub & Bla & Blub & Bla & Blub & Bla \\
        Bla & Blub & Bla & Blub & Bla & Blub & Bla & Blub \\
        Blub & Bla & Blub & Bla & Blub & Bla & Blub & Bla \\
    \end{tabular}
    \caption{\lipsum[9]}
    \label{tab:blub}
\end{table}

[/columns]

[columns=2]

[column=0.4]

# Conclusions

\lipsum[12]

\vskip1.4cm

[column]

# Bibliography

\printbibliography

\vskip7.3cm

[/columns]
