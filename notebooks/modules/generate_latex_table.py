from pylatex import Document, Table, Tabular, NoEscape

# Create a new LaTeX document
doc = Document()

# Add necessary packages for text wrapping and centering
doc.preamble.append(NoEscape(r'\usepackage{array}'))       # For column width and alignment
doc.preamble.append(NoEscape(r'\usepackage{makecell}'))    # To enable line breaks within cells
doc.preamble.append(NoEscape(r'\usepackage{graphicx}'))    # For resizing the table to fit within the page
doc.preamble.append(NoEscape(r'\usepackage[a4paper,margin=1in]{geometry}'))  # Set page margins
doc.preamble.append(NoEscape(r'\usepackage{caption}'))     # For customizing captions

# Configure the caption package to remove table numbering
doc.preamble.append(NoEscape(r'\captionsetup[table]{labelformat=empty}'))

# Create a table
with doc.create(Table(position='!h')) as table:
    # Use the \caption* command to create an unnumbered caption
    table.append(NoEscape(r'\caption*{294 $\beta \beta$ events $\rightarrow$ \textasciitilde 40910 s (\textasciitilde 11h 20mins)}'))

    # Center the table horizontally
    table.append(NoEscape(r'\centering'))

    # Start the \resizebox command to automatically resize the table
    table.append(NoEscape(r'\resizebox{1.1\textwidth}{!}{%'))

    # Specify column widths and alignments with m{} for vertical centering
    column_format = '|m{4cm}|m{2cm}|m{2cm}|m{5cm}|m{3.5cm}|'
    with doc.create(Tabular(column_format)) as tabular:
        tabular.add_hline()
        
        # Centering the column names
        tabular.add_row([
            NoEscape(r'\makecell[c]{Function}'), 
            NoEscape(r'\makecell[c]{Timecost [s]}'),
            NoEscape(r'\makecell[c]{$N_{calls}$}'), 
            NoEscape(r'\makecell[c]{Function description}'),
            NoEscape(r'\makecell[c]{Mother function}')
        ])
        tabular.add_hline()

        tabular.add_row([
            NoEscape(r'\makecell[c]{\\ CreateSignalHDF5}'), 
            NoEscape(r'\makecell[c]{\\ \textasciitilde 40908}'),
            NoEscape(r'\makecell[c]{\\ 1}'),
            NoEscape(r'\makecell[c]{Takes the hits from\\ nexus as input and\\ returns a .h5 file with the\\ s2 waveforms generated}'),
            NoEscape(r'\makecell[c]{\\ None}')
        ])
        tabular.add_hline()

        tabular.add_row([
            NoEscape(r'\makecell[c]{\\ AddShapinAndSamplin}'), 
            NoEscape(r'\makecell[c]{\\ \textasciitilde 26418}'),
            NoEscape(r'\makecell[c]{294\\ (1 per event\\ non-empty)}'),
            NoEscape(r'\makecell[c]{Takes raw waveforms and\\ adds the sensor\textquotesingle s response\\ effects returning s2 waveforms\\ shaped and sampled}'),
            NoEscape(r'\makecell[c]{\\ CreateSignalHDF5}')
        ])
        tabular.add_hline()

        tabular.add_row([
            NoEscape(r'\makecell[c]{convolve}'), 
            NoEscape(r'\makecell[c]{\\ \textasciitilde 26075}'),
            NoEscape(r'\makecell[c]{31752\\ (108 per\\ event, 1 per\\ sensor)}'),
            NoEscape(r'\makecell[c]{Creates a convolution between\\ a delta-like signal and a\\ generic SiPM\textquotesingle s response\\ waveform for the shaping}'),
            NoEscape(r'\makecell[c]{\\ AddShapinAndSamplin}')
        ])
        tabular.add_hline()

        tabular.add_row([
            NoEscape(r'\makecell[c]{\\ FindS2}'), 
            NoEscape(r'\makecell[c]{\\ \textasciitilde 11339}'),
            NoEscape(r'\makecell[c]{\\ 31752}'),
            NoEscape(r'\makecell[c]{Finds the corresponding\\ s2 signal from the light\\ maps according to the\\ particle\textquotesingle s position}'),
            NoEscape(r'\makecell[c]{\\ CreateSignalHDF5}')
        ])
        tabular.add_hline()

        tabular.add_row([
            NoEscape(r'\makecell[c]{FindRotation}'), 
            NoEscape(r'\makecell[c]{\textasciitilde 5394}'),
            NoEscape(r'\makecell[c]{3105745740\\ (1 per ie$^-$,\\ 1 per sensor)}'),
            NoEscape(r'\makecell[c]{Finds the rotation based\\ on the particle\textquotesingle s position}'),
            NoEscape(r'\makecell[c]{FindS2}')
        ])
        tabular.add_hline()

        tabular.add_row([
            NoEscape(r'\makecell[c]{safe\_write\_to\_hdf}'), 
            NoEscape(r'\makecell[c]{\textasciitilde 2887}'),
            NoEscape(r'\makecell[c]{589\\ (2 per event\\ $+$1 to save\\ configuration)}'),
            NoEscape(r'\makecell[c]{Saves signal, event, and\\ configuration info in a .h5}'),
            NoEscape(r'\makecell[c]{CreateSignalHDF5}')
        ])
        tabular.add_hline()

        tabular.add_row([
            NoEscape(r'\makecell[c]{FindSensor}'), 
            NoEscape(r'\makecell[c]{\textasciitilde 1737}'),
            NoEscape(r'\makecell[c]{ 3105745740}'),
            NoEscape(r'\makecell[c]{Finds the sensor ID\\ given it\textquotesingle s rotation}'),
            NoEscape(r'\makecell[c]{FindS2}')
        ])
        tabular.add_hline()

        tabular.add_row([
            NoEscape(r'\makecell[c]{\\ AddDriftAndDiffusion}'), 
            NoEscape(r'\makecell[c]{\\ \textasciitilde 170}'),
            NoEscape(r'\makecell[c]{\\ 294}'),
            NoEscape(r'\makecell[c]{Given the energy hits from\\ nexus, it simulates the\\ drift and diffusion of\\ the ie$^-$}'),
            NoEscape(r'\makecell[c]{\\ CreateSignalHDF5}')
        ])
        tabular.add_hline()

    # End the \resizebox command
    table.append(NoEscape(r'}'))

# Save the document to a .tex file
doc.generate_tex('latex_table')

import subprocess

# Run the LaTeX compiler (make sure pdflatex is installed and available in your PATH)
subprocess.run(["pdflatex", "latex_table.tex"])




