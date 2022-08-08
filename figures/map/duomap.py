import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
import numpy as np
import argparse
import os

matplotlib.use('Qt5Agg')
plt.rcParams['font.family'] = 'Arial'

pd.set_option('display.max_columns', 500)
font = {'family': 'Arial', 'weight': 'regular', 'size': 8}
matplotlib.rc('font', **font)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Generate maps",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--config", required=True, help="Configuration file in TSV format")
    args = parser.parse_args()

    config = args.config

    path = os.path.abspath(os.getcwd()) + '/'
    config = args.config


    # path = '/Users/Anderson/Desktop/map/'
    # config = path + 'config_states_sgtf.tsv'

    def load_table(file):
        df = ''
        if str(file).split('.')[-1] == 'tsv':
            separator = '\t'
            df = pd.read_csv(file, encoding='utf-8', sep=separator, dtype='str')
        elif str(file).split('.')[-1] == 'csv':
            separator = ','
            df = pd.read_csv(file, encoding='utf-8', sep=separator, dtype='str')
        elif str(file).split('.')[-1] in ['xls', 'xlsx']:
            df = pd.read_excel(file, index_col=None, header=0, sheet_name=0, dtype='str')
            df.fillna('', inplace=True)
        else:
            print('Wrong file format. Compatible file formats: TSV, CSV, XLS, XLSX')
            exit()
        return df


    params = load_table(config)
    params.fillna('', inplace=True)
    params = params.set_index('param')
    # print(params)

    backend = params.loc['backend', 'value']
    matplotlib.use(backend)

    display_header = params.loc['display_header', 'value']
    if display_header.lower() == 'yes':
        shape1, shape2 = pd.DataFrame(), pd.DataFrame()
        shape1_path = params.loc['shape1', 'value']
        if shape1_path not in ['', None]:
            shape1 = gpd.read_file(shape1_path)

        shape2_path = params.loc['shape2', 'value']
        if shape2_path not in ['', None]:
            shape2 = gpd.read_file(shape2_path)

        for dataframe in [shape1, shape2]:
            if not dataframe.empty:
                # dict_df = dataframe.head(1).to_dict('list')
                dict_df = dataframe.to_dict('list')
                print('\nExample of column available in shapefile:')
                for col, val in dict_df.items():
                    if col != 'geometry':
                        print('\t- ' + col + ' = ' + str(val[0]))
                    else:
                        print('\t- ' + col + ' = POLYGON (...)')
        exit()

    # Load sample metadata
    input1 = path + params.loc['input1', 'value']
    df1 = load_table(input1)
    df1.fillna('X', inplace=True)

    # load bubble pinpoints, if any is available
    dfP = pd.DataFrame()
    input2 = params.loc['input2', 'value']
    if input2 not in [None, '']:
        input2 = path + input2
        dfP = load_table(input2)
        dfP['lat'] = dfP['lat'].astype(float)
        dfP['long'] = dfP['long'].astype(float)

    id_geocol = params.loc['unique_id1', 'value']
    id_datacol = params.loc['unique_id2', 'value']

    map_type = params.loc['map_type', 'value']

    xvar = params.loc['xvar', 'value']
    xvar_bins = [float(bound) for bound in params.loc['xbins', 'value'].split(',')]

    if map_type == 'bivariate':
        yvar = params.loc['yvar', 'value']
        yvar_bins = [float(bound) for bound in params.loc['ybins', 'value'].split(',')]

    show_legend = params.loc['legend', 'value']
    size_factor = params.loc['size_factor', 'value']
    if size_factor not in ['', None]:
        size_factor = float(params.loc['size_factor', 'value'])

    if params.loc['alpha', 'value'] not in ['', None]:
        alpha_value = float(params.loc['alpha', 'value'])
    else:
        alpha_value = 0.5

    # plot dimensions
    plot_width = int(params.loc['figsize', 'value'].split(',')[0].strip())
    plot_heigth = int(params.loc['figsize', 'value'].split(',')[1].strip())


    # filter
    def filter_df(df, criteria):
        new_df = pd.DataFrame()
        for filter_value in sorted(criteria.split(',')):
            filter_value = filter_value.strip()
            if not filter_value.startswith('~'):
                df_filtered = df[df[filter_value.split(':')[0]].isin([filter_value.split(':')[1]])]
                new_df = new_df.append(df_filtered)

        for filter_value in sorted(criteria.split(',')):
            filter_value = filter_value.strip()
            if filter_value.startswith('~'):
                filter_value = filter_value[1:]
                if new_df.empty:
                    df = df[~df[filter_value.split(':')[0]].isin([filter_value.split(':')[1]])]
                    new_df = new_df.append(df)
                else:
                    new_df = new_df[~new_df[filter_value.split(':')[0]].isin([filter_value.split(':')[1]])]
        return new_df


    filters = params.loc['filter', 'value']
    df2 = pd.DataFrame()
    if filters not in ['', None]:
        df2 = filter_df(df1, filters)
    else:
        df2 = df1

    if not dfP.empty:
        if filters not in ['', None]:
            dfP = filter_df(dfP, filters)

    labels_x = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']
    labels_y = ['A', 'B', 'C', 'D', 'E', 'F']

    hexcolours = [hex.strip() for hex in params.loc['map_colours', 'value'].replace('"', '').split(',')]
    notfoundcolour = params.loc['missing_colour', 'value']

    # print(hexcolours)

    colour_scheme = {}
    c = 0
    if map_type == 'bivariate':
        for y in labels_y[:len(yvar_bins)]:
            for x in labels_x[:len(xvar_bins)]:
                label = y + x
                # print(label)
                colour_scheme[label] = hexcolours[c]
                c += 1
    elif map_type == 'choropleth':
        for x in labels_x[:len(xvar_bins)]:
            label = x
            try:
                colour_scheme[label] = hexcolours[c]
            except:
                print('\nA total of ' + str(len(xvar_bins)) + ' HEX colours are required.')
                exit()
            c += 1

    ticks = []


    def classify(df, column, bins, axis, lowest):
        origin = 0
        if lowest not in [''] and str(lowest)[-1].isdigit():
            origin = float(lowest)
        bins = [origin] + bins
        # print(bins)
        for idx, row in df.iterrows():
            value = float(df.loc[idx, column])
            for num, varbin in enumerate(bins):
                # print(num, len(bins) - 1)
                if num < len(bins) - 1:
                    start, end = bins[num], bins[num + 1]
                    # print(bins[num], bins[num + 1])
                    tick_label = str(start) + '-' + str(end)
                    if tick_label not in ticks:
                        ticks.append(tick_label)
                    if start < value <= end:
                        # print(start, '<', value, '<=', end, ' = ', labels_x[num])
                        if axis == 'x':
                            df.loc[idx, 'x_category'] = labels_x[num]
                        else:
                            df.loc[idx, 'y_category'] = labels_y[num]


    # add x and y categories
    lowestX = float(params.loc['xlowest', 'value']) - 0.01
    classify(df2, xvar, xvar_bins, 'x', lowestX)
    if map_type == 'bivariate':
        lowestY = float(params.loc['ylowest', 'value']) - 0.01
        classify(df2, yvar, yvar_bins, 'y', lowestY)
        df2['category'] = df2['y_category'].astype(str) + df2['x_category'].astype(str)
    elif map_type == 'choropleth':
        df2['category'] = df2['x_category'].astype(str)

    # # output converted dataframes
    # output = params.loc['output_file', 'value']
    # if output not in ['', None]:
    #     df2.to_csv(output, sep='\t', index=False)
    df2.set_index(id_datacol, inplace=True)

    # shapefiles
    shape1_path = params.loc['shape1', 'value']
    shape1 = gpd.read_file(shape1_path)

    # keep only specific polygons in geodf 1
    keep_only1 = params.loc['keep_only1', 'value']
    geodf = pd.DataFrame()
    if keep_only1 not in ['', None]:
        geodf = filter_df(shape1, keep_only1)
    else:
        geodf = shape1

    if 'macrorregioes' in shape1_path:
        geodf['cd_hlt_'] = geodf['cd_hlt_'].astype(str).apply(lambda x: x[:-3])

    shape2_path = params.loc['shape2', 'value']
    shape2 = pd.DataFrame()
    geodf2 = pd.DataFrame()
    if shape2_path not in ['', None]:
        shape2_path = params.loc['shape2', 'value']
        shape2 = gpd.read_file(shape2_path)

        # keep only specific polygons in geodf 2
        keep_only2 = params.loc['keep_only2', 'value']
        geodf2 = pd.DataFrame()
        if keep_only2 not in ['', None]:
            geodf2 = filter_df(shape2, keep_only2)
        else:
            geodf2 = shape2

    # plot specifications
    nrows = int(params.loc['nrows', 'value'])
    ncols = int(params.loc['ncols', 'value'])

    figsize = (plot_width, plot_heigth)
    fig, axes = plt.subplots(nrows, ncols, sharex=True, sharey=True, figsize=figsize)

    group = params.loc['groupby', 'value']

    # table legend
    cells = []
    legend_hex = []
    for b in reversed(range(len(xvar_bins))):
        cells.append(list(colour_scheme.keys())[b * len(xvar_bins):b * len(xvar_bins) + len(xvar_bins)])
        legend_hex.append(list(colour_scheme.values())[b * len(xvar_bins):b * len(xvar_bins) + len(xvar_bins)])

    cells = [l for l in cells if len(l) > 0]
    legend_hex = [l for l in legend_hex if len(l) > 0]

    # print(df2.head())
    # print(cells)
    # print(legend_hex)
    # print(df2.columns.tolist())

    for i, (name, df) in enumerate(df2.groupby(group)):
        upperleft_coord = [float(value.strip()) for value in params.loc['upperleft_coord', 'value'].split(',')]
        lowerright_coord = [float(value.strip()) for value in params.loc['lowerright_coord', 'value'].split(',')]

        xlim = ([upperleft_coord[1], lowerright_coord[1]])
        ylim = ([lowerright_coord[0], upperleft_coord[0]])

        if nrows >= 2 and ncols >= 2:  # nrows + ncols > 3:
            # print(i // ncols, i % ncols)
            ax = axes[i // ncols][i % ncols]
            # ax.set_xlim(xlim)
            # ax.set_ylim(ylim)
        elif nrows < 2 or ncols < 2:
            # print(i // ncols, i % ncols)
            ax = axes[i]
            # ax.set_xlim(xlim)
            # ax.set_ylim(ylim)

        # print('\n' + name)
        for shape, data in geodf.groupby(id_geocol):
            shape = str(shape)  # [:-3]
            # print(shape)

            # define the color for each group using the dictionary
            if shape in df.index:
                # print(shape)
                # print(df.index.tolist())
                value = df.loc[shape, 'category']
            else:
                value = 'NA'

            if value in colour_scheme:
                colour = colour_scheme[value]
            else:
                # print(shape, 'No colour')
                colour = notfoundcolour

            # ax.axis('off')

            # Plot each group using the color defined above
            stroke1 = params.loc['stroke1', 'value']
            stroke2 = params.loc['stroke2', 'value']

            color1, line1 = '#000000', 0.25
            color2, line2 = '#000000', 0.5
            for stroke in [stroke1, stroke2]:
                if stroke not in ['', None]:
                    color, line = [element.strip() for element in stroke.split(',')]
                    if stroke == stroke1:
                        if color != '':
                            color1 = color
                        if line != '':
                            line1 = float(line)
                    if stroke == stroke2:
                        if color != '':
                            color2 = color
                        if line != '':
                            line2 = float(line)

            data.plot(color=colour,
                      ax=ax,
                      edgecolor=color1,
                      linewidth=line1,
                      label=value)
            ax.set_axis_off()

        if shape2_path not in ['', None]:
            geodf2.plot(ax=ax, color='none', edgecolor=color2, linewidth=line2)

        # plot pinpoints
        if input2 not in [None, '']:
            pinpoints = dfP.groupby(group).get_group(name)
            # print('>>>' + name)
            # print(pinpoints)
            if 'size' not in pinpoints.columns:
                pinpoints['size'] = 1
                min_size = 1
            else:
                min_size = min(pinpoints['size'].astype(float))
                max_size = max(pinpoints['size'].astype(float))
            if 'hex_colour' not in pinpoints.columns:
                pinpoints['hex_colour'] = '#000000'

            if 'label' not in pinpoints.columns:
                pinpoints['label'] = ''

            size_norm = params.loc['size_norm', 'value']  # show bubbles in log scale

            for idx, row in pinpoints.iterrows():
                label = pinpoints.loc[idx, 'label']
                point = (pinpoints.loc[idx, 'long'], pinpoints.loc[idx, 'lat'])
                hex = pinpoints.loc[idx, 'hex_colour']
                size = float(pinpoints.loc[idx, 'size'])
                if size_norm == 'yes':
                    # bubble_size = np.log10(size - min_size / max_size - min_size) * -1
                    bubble_size = (size - min_size / max_size - min_size) * size_factor
                    print(bubble_size)
                else:
                    bubble_size = (size / min_size) * size_factor

                draw_circle = plt.Circle(point, radius=bubble_size, linewidth=0, facecolor=hex, edgecolor=hex,
                                         fill=True, alpha=alpha_value)
                ax.add_artist(draw_circle)
                ax.text(point[0] + 0.05, point[1] + 0.05, s=label, horizontalalignment='left', fontsize=6,
                        color='#000000', zorder=30)

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        # top label
        ax.text(upperleft_coord[1] + 0.05, upperleft_coord[0] + 1,
                s=name, horizontalalignment='left', fontsize=10, color='#000000', zorder=30)

        # table legend
        axins = ax.inset_axes([-0.2, -0.1, 0.6, 0.6])
        axins.axis('tight')
        axins.axis('off')

        # inner legends
        if map_type == 'bivariate':
            counts = df['category'].value_counts().to_dict()
            newdata = []
            for row in cells:
                new_row = []
                for cat in row:
                    if cat in counts:
                        new_row.append(counts[cat])
                    else:
                        new_row.append('-')
                newdata.append(new_row)

            table = axins.table(cellText=newdata, cellLoc='center', loc="center",
                                colWidths=[0.1 for x in range(len(xvar_bins))], fontsize=6)

            table.scale(1, plot_heigth / (nrows * 4))

            for key, cell in table.get_celld().items():
                cell.set_linewidth(0.2)

    # legend_elements = [Line2D([0], [0], marker='o', color='w', label='Group A', markerfacecolor='#4393E5', markersize=10),
    #                    Line2D([0], [0], marker='o', color='w', label='Group B', markerfacecolor='#43BAE5', markersize=10),
    #                    Line2D([0], [0], marker='o', color='w', label='Group C', markerfacecolor='#7AE6EA', markersize=10)]

    add_legend = False
    if str(params.loc['legend', 'value']).lower() == 'true':
        add_legend = True

    if add_legend == True:
        if nrows >= 2 and ncols >= 2:
            ax = axes[nrows - 1][ncols - 1]
        else:
            ax = axes[-1]

        ax.axis('tight')
        ax.axis('off')

    if map_type == 'bivariate':
        table = ax.table \
            (cellText=cells, cellColours=legend_hex, cellLoc='center', loc="center",
             rowLabels=ticks[:len(xvar_bins)][::-1], colLabels=ticks[len(xvar_bins):],
             colWidths=[0.15 for x in range(len(xvar_bins))], fontsize=8)
        table.scale(1, plot_heigth / (nrows * 1.5))
    elif map_type == 'choropleth':
        table = ax.table \
            (cellText=cells, cellColours=legend_hex, cellLoc='center', loc="center",
             colLabels=ticks, colWidths=[0.15 for x in range(len(xvar_bins))], fontsize=8)
        table.scale(1, plot_heigth / (nrows * 1.5))

    if backend == 'pdf':
        plt.tight_layout()
        plt.savefig(path + "map_" + config.split('.')[0].split('_')[-1] + ".pdf", format="pdf", bbox_inches="tight")
        # plt.savefig("map.pdf", format="pdf", bbox_inches="tight")
    else:
        plt.show()
