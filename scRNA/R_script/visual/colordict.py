import matplotlib.pyplot as plt

# create a color dictionary for celltype_major
def celltype_major_colors():

    celltype_major_colors_dict = {'Naive': '#8dd3c7',
             'Effector memory': '#ffffb3',
             'Memory': '#80b1d3',
             'MHC II': '#8da0cb',
             'Resident memory': '#b3de69',
             'Exhausted': '#fb8072',
             'Interferon': '#bebada',
             'Cycling': '#d9d9d9',
             'Nk-like': '#fc8d62',
             'MAIT': '#ccebc5',
             'γδT': '#fccde5'}
    
    return celltype_major_colors_dict

# create a color dictionary for celltype_minor

def celltype_minor_colors():

    celltype_minor_colors_dict = {'Tn.c01.CCR7': '#8dd3c7',
             'Tem.c02.GZMK': '#ffffb3',
             'Tem.c03.HMGB2': '#fdb462',
             'Tem.c04.GZMH': '#ffed6f',
             'Tm.c05.NR4A1': '#80b1d3',
             'T.c06.MHCII': '#8da0cb',
             'Trm.co6.ZNF683': '#b3de69',
             'Tex.c08.CXCL13': '#fb8072',
             'T.c09.IFN': '#bebada',
             'T.c10.STMN1': '#d9d9d9',
             'T.c11.ASPN': '#66c2a5',
             'T.c12.NK-like': '#fc8d62',
             'MAIT.c13': '#ccebc5',
             'γδT.c14.GNLY': '#fccde5',
             'γδT.c15.TRDV1': '#bc80bd',
             'γδT.c16.TRDV2': '#e78ac3'}
    
    return celltype_minor_colors_dict

# create a color dictionary for celltype_Tex
def celltype_Tex_colors():

    celltype_Tex_colors_dict = {'Tex.c01.CCL4': '#4dbbd5',
             'Tex.c02.GZMH': '#3c5488',
             'Tex.c03.IL7R': '#f39b7f',
             'Tex.c04.CREM': '#00a087'}
    
    return celltype_Tex_colors

# create a color dictionary for donors
def donor_colors():
    donors = ['TSP1','TSP2','TSP3','TSP4','TSP5','TSP6','TSP7','TSP8','TSP9','TSP10','TSP11','TSP12','TSP13','TSP14','TSP15']
    
    import matplotlib.colors as pltcolors
    
    cmap = plt.cm.get_cmap("YlGnBu")
        
    donor_color_dict = {}
    j=1/len(donors)
    for d in donors:
        donor_color_dict[d] = pltcolors.to_hex(cmap(j))
        j+=1/len(donors)
        
    return donor_color_dict


# create a color dictionary for donors
def compartment_colors():
    
    compartments = ['endothelial', 'epithelial', 'immune', 'stromal', "germ line",'PNS']
    
    import matplotlib.colors as pltcolors
    
    cmap = plt.cm.get_cmap("YlOrRd")
        
    compartment_color_dict = {}
    j=1/len(compartments)
    for c in compartments:
        compartment_color_dict[c] = pltcolors.to_hex(cmap(j))
        j+=1/len(compartments)
        
    return compartment_color_dict


# create a color dictionary for methods
def method_colors():
    methods = ['10X','smartseq2']
    
    import matplotlib.colors as pltcolors
    
    method_color_dict = {}
    method_color_dict['10X'] = '#90ee90'
    method_color_dict['smartseq2'] = '#006400'
    
    return method_color_dict

# create a color dictionary for sex
def sex_colors():
    sexes = ['male','female','undisclosed']
    
    import matplotlib.colors as pltcolors
    
    sex_color_dict = {}
    sex_color_dict['female'] = '#f4cae4'
    sex_color_dict['male'] = '#cbd5e8'
#     sex_color_dict['undisclosed'] = '#e6f5c9'
    
    return sex_color_dict

# visualize the color dictionaries
def plot_colortable(colors, title, sort_colors=True, emptycols=0):

    cell_width = 212
    cell_height = 22
    swatch_width = 48
    margin = 12
    topmargin = 40

    # Sort colors by hue, saturation, value and name.
    by_hsv = [(v, k) for k, v in colors.items()]
    
    if sort_colors is True:
        by_hsv = sorted(by_hsv)
    names = [name for hsv, name in by_hsv]

    n = len(names)
    ncols = 4 - emptycols
    nrows = n // ncols + int(n % ncols > 0)

    width = cell_width * 4 + 2 * margin
    height = cell_height * nrows + margin + topmargin
    dpi = 72

    fig, ax = plt.subplots(figsize=(width / dpi, height / dpi), dpi=dpi)
    fig.subplots_adjust(margin/width, margin/height,
                        (width-margin)/width, (height-topmargin)/height)
    ax.set_xlim(0, cell_width * 4)
    ax.set_ylim(cell_height * (nrows-0.5), -cell_height/2.)
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    ax.set_axis_off()
    ax.set_title(title, fontsize=24, loc="left", pad=10)

    for i, name in enumerate(names):
        row = i % nrows
        col = i // nrows
        y = row * cell_height

        swatch_start_x = cell_width * col
        swatch_end_x = cell_width * col + swatch_width
        text_pos_x = cell_width * col + swatch_width + 7

        ax.text(text_pos_x, y, name, fontsize=14,
                horizontalalignment='left',
                verticalalignment='center')

        ax.hlines(y, swatch_start_x, swatch_end_x,
                  color=colors[name], linewidth=18)

    return fig