# -*- coding: utf-8 -*-
import os
import numbers
import filecmp
import warnings
from collections import defaultdict

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.patches as patches
from matplotlib.lines import Line2D
from matplotlib_venn import venn3
import seaborn as sns
import numpy as np
import pandas as pd
import scipy as sc
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
from statsmodels.stats.proportion import proportions_ztest
from scipy.stats import linregress

from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Myriad Pro']


# If fonts don't work, move fonts to /usr/share/matplotlib/mpl-data/fonts/ttf/ 
# or .local/lib/python2.7/site-packages/matplotlib/mpl-data/fonts/ttf/ and run:
#from matplotlib import font_manager
#font_manager._rebuild()


class Plotter(object):
    
    linewidth = 1
    
    def __init__(self, output_path='results', page_size=7):
        self.output_path = output_path
        if not os.path.exists(output_path):
            os.mkdir(output_path)
            
        self.page_size = page_size
            
        plt.rc('xtick', labelsize=5)
        plt.rc('ytick', labelsize=5)
        plt.rc('axes', labelsize=6)
#        plt.rc('subplot', left=7)
        
        plt.rc('figure.subplot', left=0.15, bottom=0.15, top=0.95, right=0.95)
        
    
    def _legend(self, 
        ax, plots=[], labels=[], position='right', legend_title=None, reverse=True, handletextpad=None,
        shift_right=0, shift_top=0, ncol=1, columnspacing=None, markerscale=None, handlelength=None,
        handler_map=None
    ):
        if reverse:
            plots = plots[::-1]
            labels = labels[::-1]
        if ncol == 2:
            plots = plots[::2] + plots[1::2]
            labels = labels[::2] + labels[1::2]
        
        frameon = False
        borderpad = None
        labelspacing = 0.3
        bbox = None
        if position == 'right':
            bbox = (1+shift_right, 1+shift_top)
            loc='upper left'
        if position == 'above':
            bbox = (shift_right, 1+shift_top)
            loc='lower left'
        if position == 'floatright':
            loc='upper right'
            frameon = True
            borderpad = 0.5
            handletextpad = 0.5
            handlelength=0.8
            

        legend_props = {
            'bbox_to_anchor':bbox, 'loc': loc, 
            'prop': {'size': 5}, 'frameon': frameon,
            'title': legend_title, 
            'ncol': ncol,
            'borderpad': borderpad, 'handlelength': handlelength, 'handletextpad': handletextpad,
            'markerscale': markerscale,
            'columnspacing': columnspacing, 'labelspacing': labelspacing
        }

        if plots and labels:
            legend = ax.legend(plots, labels, handler_map=handler_map, **legend_props)
        else:
            legend = ax.legend(handler_map=handler_map, **legend_props)
            
        if legend_title:
            legend.set_title(legend_title, prop={'size': 5})
            legend._legend_box.align = "left"
        
        if position == 'floatright':
            legend.get_frame().set_linewidth(0.5)
        
        return legend


    def plot(
        self, plot_type, data, save='', show=True, size=0.5, colors=defaultdict(lambda: 'k'), save_scale=2,
        y_label='', x_label='', xlim=None, ylim=None, no_x=False, no_y=False,
        hline=None, vline=None, xticks=None, xticklabels=None, yticks=None, yticklabels=None,
        hlines=None,
        xpad=2, ypad=2, xtickpad=1, ytickpad=1, ylog=False, xlog=False,
        xgrid=False, ygrid=False, xintersect=None, yintersect=None, border=None,
        scalex=None, colorbar=None, margin={'left': 0.08, 'right': 0.02, 'top': 0.02, 'bottom': 0.08},
        *args, **kwargs
    ):

        plot_method = getattr(self, plot_type, None)
        assert callable(plot_method), 'Plot type ' + plot_type + ' does not exist'
        
        if show:
            plt.ion()
        else:
            plt.ioff()
            
        if isinstance(size, numbers.Number):
            size = (size, size * 0.75)
        size = np.array(size)
        self.figure_size = np.copy(size)
            

        size[0] += margin['left'] + margin['right']
        size[1] += margin['top'] + margin['bottom']
            
        fig = plt.figure(figsize=(size*self.page_size), dpi=150)
        ax = fig.gca()
        
 
        plt.subplots_adjust(
            left=margin['left']/size[0],
            right=1-margin['right']/size[0],
            bottom=margin['bottom']/size[1],
            top=1-margin['top']/size[1],
        )
        
        if not border:
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
        if no_x:
            ax.spines['bottom'].set_visible(False)
        if no_y:
            ax.spines['left'].set_visible(False)
        for d in ('top', 'left', 'right', 'bottom'):
            ax.spines[d].set_linewidth(border or Plotter.linewidth)
        
        ax.tick_params(width=Plotter.linewidth, length=3)
        if no_x and no_y:    
            ax.tick_params(width=Plotter.linewidth, length=0)
        ax.tick_params(axis='x', which='both', pad=xtickpad)
        ax.tick_params(axis='y', which='both', pad=ytickpad)
        
        if xgrid:
            plt.grid(which='major', axis='x', linewidth=Plotter.linewidth, color='#eeeeee', clip_on=False, zorder=0)
        if ygrid:
            plt.grid(which='major', axis='y', linewidth=Plotter.linewidth, color='#eeeeee', clip_on=False, zorder=0)

        if xlog:
            ax.set_xscale('log')
        if ylog:
            ax.set_yscale('log')
        if ylim:
            if type(ylim) == tuple:
                ax.set_ylim(ylim)
            else:
                ax.set_ylim(top=ylim) 
        if xlim:
            ax.set_xlim(xlim)
            
        # **********************************************************************
        
        plot_method(ax, data, colors, *args, **kwargs)
        
        # **********************************************************************
        
        if plot_type == 'hist' and xlim:
            ax.set_xlim(xlim)

        ax.set_ylabel(y_label, labelpad=ypad)
        ax.set_xlabel(x_label, labelpad=xpad)
        
        if xticks is not None:
            if str(xticks) == 'remove':
                ax.set_xticks([])
                ax.set_xticks([], minor=True)
            else:
                ax.set_xticks(xticks)
        if xticklabels is not None:
            ax.set_xticklabels(xticklabels)
        if yticks is not None:
            if str(yticks) == 'remove':
                ax.set_yticks([])
                ax.set_xticks([], minor=True)
            else:
                ax.set_yticks(list(yticks))
        if yticklabels is not None:
            ax.set_yticklabels(yticklabels)
            
        if hline is not None:
            ax.axhline(hline, c='k', linewidth=Plotter.linewidth, zorder=-1)
        if hlines:
            for y in hlines:
                ax.axhline(y, c='#eeeeee', linewidth=Plotter.linewidth/2.0, zorder=-1, clip_on=False)
                
        if vline is not None:
            ax.axvline(vline, color='#777777', linestyle='dotted', linewidth=Plotter.linewidth, zorder=-1)
        
        if xintersect is not None:
            ax.spines['bottom'].set_position(('data', xintersect))
        if yintersect is not None:
            ax.spines['left'].set_position(('data', yintersect))

        if colorbar:
            ax_bar = fig.add_axes([0.89, 0.41, 0.02, 0.4])
            (r1, g1, b1), (r2, g2, b2) = colorbar
            cmap = mpl.colors.ListedColormap(zip(*[np.linspace(r1, r2, 101), np.linspace(g1, g2, 101), np.linspace(b1, b2, 101)]))
            cb = mpl.colorbar.ColorbarBase(ax_bar, cmap=cmap, orientation='vertical', ticks=[0, 1])
            cb.set_ticklabels(['0%', '100%'])
            ax_bar.set_title('Convergent\nsimulations', size=5)

        # scale x-axis internally.
        if scalex:
            xticks = ax.get_xticks()
            ax.set_xticklabels([x*scalex for x in xticks])

        if save:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                png_path = os.path.join(self.output_path, save + '.png')
                png_path_old = os.path.join(self.output_path, save + '__old.png')
                pdf_path = os.path.join(self.output_path, save + '.pdf')
                # Only update pdf if the graph changed, to prevent confusing git.
                if os.path.exists(png_path):
                    os.rename(png_path, png_path_old)
                plt.savefig(png_path, dpi=fig.dpi*save_scale)

                graph_changed = True
                if os.path.exists(png_path_old):
                    graph_changed = not(filecmp.cmp(png_path, png_path_old))
                    os.remove(png_path_old)
                try:
                    if graph_changed:
                        plt.savefig(pdf_path) #, bbox_inches='tight')
                except Exception:
                    print('Error while trying to save pdf')

            print(f'Saved to `{os.path.join(self.output_path, save)}`:')

        if show:
            plt.show()
        else:
            plt.close()

    def empty(self, ax, data, colors):
        ax.plot([], [])

    def table(self, ax, columns, colors, width=1000, height=1000, header=None, row_names=None):

        if header is None:
            header = []
        if row_names is None:
            row_names = []

        ax.plot([], [])
        ax.set_xlim([0, width])
        ax.set_ylim([0, height])
        
        w = float(width)
        h = float(height)
        
        x_coords = list(range(0, width, width // len(columns)))
        for x in x_coords:
            ax.axvline(x, color='k', lw=0.5, clip_on=False)
        for y in [0, height]:
            ax.axhline(y, xmax=x_coords[-1]/w, color='k', lw=0.5, clip_on=False)
        
        row_height = 10
        column_width = x_coords[1]
        n_per_col = 3
        
        if header:
            height -= 20
            ax.axhline(height, xmax=x_coords[-1]/w, color='k', lw=0.5, clip_on=False)
            ax.add_patch(patches.Rectangle((0, height), column_width*len(columns), 20, linewidth=0, facecolor='#cccccc'))
                              
            for i, h in enumerate(header):
                ax.text(x_coords[i]+column_width*0.5, height+10, h, size=5,  
                        verticalalignment='center', horizontalalignment='center')
        
        margin = 5
        
        y_heights = []
        
        for col_i, column in enumerate(columns):
            
            y_height = []
            
            n_y = 0
            
            for row_i, row in enumerate(column):
                n_x = 0
                if row_i != 0:
                    n_y += 1 + (margin*2.0)/row_height
                    ax.axhline(
                        height - (n_y*row_height), 
                        xmin=x_coords[col_i]/w,
                        xmax=x_coords[col_i+1]/w,
                        color='k', lw=0.5, clip_on=False
                    )

                
                for n_i, n in enumerate(row):
                    
                    if n_x >= n_per_col:
                        n_x = 0
                        n_y += 1
                    
                    c = colors[n]
                    if n.startswith('BWM'):
                        n = n[4:]
                    
                    ax.text(
                        x_coords[col_i]+n_x*(column_width/n_per_col) + margin, 
                        height - (n_y*row_height + margin), 
                        n, size=5, color=c, verticalalignment='top')
                    
                    n_x += 1
                
                y_height.append(n_y+1)
            y_heights.append(y_height)
        
        for i in range(len(y_heights[-1])):
            end = y_heights[-1][i]
            start = y_heights[-1][i-1] if i else 0
            ax.text(
                x_coords[-1] + 10, 
                height - (np.mean([start, end])*row_height+margin*2), 
                row_names[i], size=5, verticalalignment='center'
            )
            
            
#        if row_names and col_i == len(columns)-1:
#            ax.text(
#                x_coords[-1], 
#                height - (n_y*row_height), 
#                row_names[row_i], size=5, 
#            )
        
        
        
    def pie(self, ax, data, colors):
        
        ax.pie(
            data, startangle=90, counterclock=False, colors=colors,
            wedgeprops={'edgecolor': (0.35, 0.25, 0), 'linewidth': 0.3}
        )
        
        
    def xy_graph(self, ax, data, colors, cumulative=False, err=[], dotted=True, rev_legend=False, no_legend=False,
                 stats=None, larval_stages=True, linkpoints=True, pvalues=None, smooth=False, alpha=1, legendpos=None,
                 legend_shift_right=0, legend_shift_top=0, markersize=6, legendcol=None, clipon=False, dots=True,
                 linewidth=None):

        """
        Generates a simple line graph, where the x-axis is the developmental
        timeline of C. elegans, with larval stages labeled instead of hours post-
        hatching.
        """
        xs, ys = data

        larval_stage_ends = [0, 16, 25, 34, 45]
        larval_stage_mids = [8, 20.5, 29.5, 39.5, 50]
        larval_stage_labels = ['L1', 'L2', 'L3', 'L4', ' Adult']


        # Set x-axis to larval stages.
        if larval_stages:
            ax.set_xlim([0, 55])
            ax.set_xticks([0])
    
            ax.tick_params(axis='x', labelbottom=False)
            ax.set_xticks(larval_stage_mids, minor=True)
            ax.set_xticklabels(larval_stage_labels, minor=True)
            ax.tick_params(axis='x', which='minor', bottom=False, pad=1)
            linestyle = 'dotted' if dotted else 'solid'
            for t in larval_stage_ends:
                plt.axvline(t, color='#999999', linestyle=linestyle, linewidth=Plotter.linewidth, zorder=-1)
        
        stats_values = []
        if stats:
            for i, (l, y) in enumerate(ys):
                reg = linregress(xs, y)
                a = reg.slope
                b = reg.intercept
                assert stats in ('regression', 'spearmanr', 'log', 'regression_only')
                if stats == 'regression':
                    p = reg.pvalue
                if stats == 'regression_only':
                    p = -1
                elif stats in ('spearmanr', 'log'):
                    if len(set(y)) == 1:
                        p = -1
                    else:
                        p = sc.stats.spearmanr(xs, y).pvalue
                    
                if stats == 'log':
                    xs[0] += 0.01
                    z = sc.optimize.curve_fit(lambda t, a, b: a+b*np.log(t), xs, y)
                    fit = lambda x: z[0][0] + z[0][1] * np.log(x)
                else:
                    fit = lambda x: b + a * x
                    
                x_fit = np.arange(0.01, 55, 0.01)
                y_fit = [fit(x) for x in x_fit]
                h = (ax.get_ylim()[1] - ax.get_ylim()[0]) / 50
                text_x = 55 + 2
                text_y = fit(text_x)
 
                if p < 0.05:
                    text_y -= h*0.5
                    
                stats_values.append((x_fit, y_fit, text_x, text_y, p))
                

        y_total = np.zeros(len(ys[0][1]))
        for i, (l, y) in enumerate(ys):
            
            color_i = i % 8
            color = ['k', 'grey', 'm', 'r', 'b', 'orange', 'g', 'y'][color_i]
            edge_color = color
            if colors and type(colors) == list:
                color = colors[i]
                edge_color = color
                
            if colors and type(colors) == dict:
                if l in colors:
                    color = colors[l]
                    edge_color = colors[l+'_edge'] if l+'_edge' in colors else color
                else:
                    edge_color = 'k'
            
            if err and sum(err[i]):
                err_xs, err_ys, err_hs = zip(*[(err_x, err_y, err_h) for err_x, err_y, err_h in zip(xs, y, err[i]) if err_h])
                e = ax.errorbar(err_xs, err_ys, yerr=err_hs, linestyle='None', elinewidth=Plotter.linewidth/1.5, color=color, capsize=0, capthick=Plotter.linewidth/1.5, clip_on=False)
                for b in e[2]:
                    b.set_clip_on(False)
            
            if cumulative:

                line_colors = {k[:-5]: v for k, v in colors.items() if k.endswith('_edge')}
                fill_colors = colors
                         
                
                colors = fill_colors
                
                x_raw = list(xs)
                y_raw = list(y)
                y_before = y_total
                y_total = y_total + y_raw
                
                if x_raw[-1] == x_raw[-2]:
                    x_average = list(xs[:-2]) + [np.mean(xs[-2:])]
                    y_average = list(y_total[:-2]) + [np.mean(y_total[-2:])]
                    y_before_average = list(y_before[:-2]) + [np.mean(y_before[-2:])]
                else:
                    x_average = x_raw
                    y_average = y_total
                    y_before_average = y_before

                ax.plot(
                    list(x_average) + [60], list(y_average) + [y_average[-1]], 
                    color=line_colors[l], label=l, linewidth=Plotter.linewidth, clip_on=clipon, 
                    markersize=0, alpha=(0.0 if line_colors[l] == fill_colors[l] else 1)
                )
                ax.plot(
                    x_raw, y_total, color=line_colors[l], label=l, linewidth=0, clip_on=clipon, 
                    marker='.', markeredgewidth=Plotter.linewidth/2, 
                    markeredgecolor=line_colors[l], markersize=markersize, 
                    alpha=(0.0 if line_colors[l] == fill_colors[l] else 1)
                )
                ax.fill_between(
                    list(x_average) + [60], 
                    list(y_average) + [y_average[-1]], 
                    list(y_before_average) + [y_before_average[-1]], 
                    facecolor=fill_colors[l], 
                    alpha=0.8
                )
                
                continue
            
            if stats:
                ax.plot(xs, y, markerfacecolor=color, markeredgecolor=edge_color, label=l, marker='.', linewidth=0, markersize=markersize, markeredgewidth=Plotter.linewidth/2.0, clip_on=clipon)

                x_fit, y_fit, text_x, text_y, p = stats_values[i]

                ps = [v[4] for v in stats_values]
                ps_corrected = fdrcorrection0(ps)[1]
                p = ps_corrected[i]

                print('Corrected p-value:', p)
                
                if p == -1:
                    t = ''
                elif p > 0.05:
                    t = 'ns'
                elif p > 0.01:
                    t = '*'
                elif p > 0.001:
                    t = '**'
                else:
                    t = '***'

                ax.plot(x_fit, y_fit, color=edge_color, marker=None, linewidth=Plotter.linewidth/2.0, linestyle='dashed', clip_on=clipon)  

                ax.text(text_x, text_y, t, ha='left', va='center', color='k', linespacing=20, size=5)

            if linkpoints:
                if dots:
                    ax.plot(xs, y, color=color, marker='.', linewidth=0, markersize=markersize, alpha=alpha, label='_nolegend_')
                
                # Take the mean of duplicate values.
                if type(y) == pd.Series:
                    y_line = y.groupby(y.index.name).mean()
                    xs_line = y_line.index
                else:
                    unique_y = defaultdict(list)
                    for j, x in enumerate(xs):
                        unique_y[x].append(y[j])
                    xs_line = sorted(set(xs))
                    y_line = [np.mean(unique_y[x]) for x in xs_line]

                linewidth = linewidth or Plotter.linewidth
                linestyle = 'solid'
                if smooth:
                    for j in range(len(y_line)):
                        if j == 0 or j == len(y_line)-1:
                            continue
                        y_line[j] = np.mean(y_line[j-1:j+2])

                ax.plot(xs_line, y_line, color=color, label=l, linewidth=linewidth, linestyle=linestyle, alpha=alpha)
            else:
                if not stats:
                    print(xs, y)
                    ax.plot(xs, y, marker='.', markerfacecolor=color, markeredgewidth=Plotter.linewidth/2, markeredgecolor=color, linewidth=0, markersize=markersize)
        
        if pvalues:
            for p, (text_x, text_y) in zip(*pvalues):
                if p > 0.05:
                    t = 'ns'
                elif p > 0.01:
                    t = '*'
                elif p > 0.001:
                    t = '**'
                else:
                    t = '***'
                ylim = ax.get_ylim()
                yrange = ylim[1]-ylim[0]
                if '*' in t:
                    text_y -= yrange*0.02
                text_y = min(ylim[1]-yrange*0.04, text_y)
                text_y = max(ylim[0]+yrange*0.04, text_y)
                ax.text(text_x, text_y, t, ha='center', va='center', color='k')
                
        
        if len(ys) > 1 and not no_legend:
            if cumulative:
                ncol = 3
                columnspacing = 1.5
                handletextpad = 0.4
                handlelength = 0.7
                labels = [l for l, y in ys]
                legend_elements = [Patch(facecolor=colors[l], edgecolor=colors[l+'_edge'], label=l) for l in labels]

            else:
                ncol = 1
                columnspacing = None
                handletextpad = 0
                handlelength = None
                position='right'
                legend_elements, labels = ax.get_legend_handles_labels()
                
            if err:
                ncol = 3
                position='above'
                columnspacing = 0.5
                handlelength = 1.3
            
            if legendcol:
                ncol = legendcol
                if ncol == 2:
                    columnspacing = -4
#            columnspacing = -5.5
#            ncol = 2
            position='above'
            
            if legendpos:
                position = legendpos
            
            self._legend(
                ax, legend_elements, labels, reverse=rev_legend, handletextpad=handletextpad,
                ncol=ncol, columnspacing=columnspacing, position=position,
                shift_right=legend_shift_right, shift_top=legend_shift_top,
                handlelength=handlelength
            )
            
            

    
    def skeleton1d(self, ax, data, color, scale=5000, datasets=[], show_post=True):
        
        ax.axis('equal')
        
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_xticks([])
        ax.tick_params(length=0, pad=10, labelsize=12)
        
        ax.axis([-1, 1, -1, 1])
        
        sep = 0.7
        
        max_dists = np.array(data[0])/scale
        pre_dists = [np.array(d)/scale for d in data[1]]
#        if show_post:
#            post_dists = [np.array(d)/scale for d in data[2]]

        ax.set_xlim((0, max(max_dists)+0.1))
        
        y_ticks = []
        
        for i, dist in enumerate(max_dists):
            
            y = -i * sep + (len(max_dists)/2.0)
            y_ticks.append(y)
        
            ax.plot((0, dist), (y, y), color='k', linewidth=Plotter.linewidth*2, solid_capstyle='round')
            ax.plot((-0.02, 0.02, -0.02, 0.02), (y-0.15, y-0.05, y+0.05, y+0.15), color='k', linewidth=Plotter.linewidth/2, clip_on=False)
            
            for pre_dist in pre_dists[i]:
                ax.arrow(pre_dist, y-0.08, 0, -0.1, head_length=0.12, width=0.01, head_width=0.1, overhang=0.05, color='r')
            

        ax.set_yticks(y_ticks)
        ax.set_yticklabels(datasets)
        
        scalebar_x = ((dist*scale-3000)/scale, dist)
        scalebar_y = (y_ticks[-1]-sep, y_ticks[-1]-sep)
        ax.plot(scalebar_x, scalebar_y, color='k', linewidth=Plotter.linewidth)
        ax.annotate(u'3 μm', xy=(np.mean(scalebar_x), np.mean(scalebar_y)-0.1), ha='center', va='top', fontsize=12)
        
        ax.arrow(dist-1, y_ticks[0]+0.11, 0, -0.1, head_length=0.12, width=0.01, head_width=0.1, overhang=0.05, color='r')
        ax.annotate(u'Synapse', xy=(dist-1+0.2, y_ticks[0]), ha='left', va='center', fontsize=12)
    
    
    def kde(self, ax, data, colors, fill=True, clear_x=True, linewidth=None):
        
        for i, (label, x) in enumerate(data):
            sns.kdeplot(x=x, label=label, color=colors[i], fill=fill, clip=(0,1), linewidth=linewidth or self.linewidth)
            
        if clear_x:
            ax.plot((0, 1), (0, 0), color='k')
            ax.tick_params(length=0)
        
        self._legend(ax, position='above')
        

    def hist(self, ax, data, colors, bins=20, fitbins=False, hist_range=None, 
             cumulative=False, relative=False, testweight=None, stats=None, 
             highlights=None):
        
        
        if hist_range and cumulative:
            hist_range = (hist_range[0], hist_range[1]*1.1)
        
        if isinstance(data, tuple):
            plots = []
            labels = []
            for i, (l, d) in enumerate(data):
                
                if relative:
                    weights = np.zeros_like(d) + 1. / len(d)
                    if testweight and len(testweight) == len(weights):
                        weights *= testweight
                else:
                    weights = np.zeros_like(d) + 1
                
                if cumulative:
                    facecolor = (0, 0, 0, 0)
                    edgecolor = colors[i]
                    plots.append(Line2D([0], [0], color=edgecolor, linewidth=Plotter.linewidth))
                else:
                    facecolor = colors[i]
                    edgecolor = np.array(colors[i][:3])/2
                    plots.append(Patch(facecolor=facecolor, edgecolor=edgecolor, linewidth=Plotter.linewidth, label=l))
                ax.hist(d, bins=(int((max(d)-min(d))) if fitbins else bins), weights=weights, histtype='stepfilled', density=cumulative, cumulative=cumulative, range=hist_range, facecolor=facecolor, linewidth=Plotter.linewidth, edgecolor=edgecolor, clip_on=cumulative)
                labels.append(l)
            
            self._legend(ax, plots[::-1], labels[::-1])
        else:
            edgecolor = np.array(colors[0][:3])/2
            ax.hist(data, bins=bins, histtype='stepfilled', density=cumulative, cumulative=cumulative, range=hist_range, facecolor=colors[0], linewidth=Plotter.linewidth, edgecolor=edgecolor, clip_on=False)
            
        if hist_range:
            r = hist_range[1] - hist_range[0]
            ax.set_xlim(hist_range[0]-r*0.02, hist_range[1]+r*0.05)
            
            if cumulative:
                ax.set_xlim(hist_range[0], hist_range[1]/1.1*1.05)

        if stats == 'wilcoxon':
            ylim = ax.get_ylim()
            h = (ylim[1] - ylim[0]) / 50
            y = ylim[1] + h
            
            p = sc.stats.wilcoxon(data).pvalue
            p = 'p < 10$^{' + str(int(np.log10(p))) + '}$'
                
            self._stats(ax, y, [(-0.25, 0.25)], [p])
            
        if highlights:
            
            for (label, (x, y)) in highlights:
                
                ax.arrow(x, y+500, 0, -500,  width=2, length_includes_head=True, head_width=6, head_length=200, linewidth=0, fc='k')
                ax.annotate(label, (x, y+500), ha='center', va='bottom')

    
    def simple_violin(self, ax, data, color, labels=[], x_positions=[], vert=True):
        
        x_positions = x_positions or range(1, len(labels)+1)
        
        violins = ax.violinplot(data, showextrema=False, showmeans=True, positions=x_positions, vert=vert)
        
        for i, violin in enumerate(violins['bodies']):
            violin.set_facecolor(color[i])
            violin.set_edgecolor('black')
            violin.set_alpha(1)
            
        violins['cmeans'].set_edgecolor('black')
        violins['cmeans'].set_linewidth(self.linewidth)
        if vert:
            ax.set_xticks(x_positions)
            ax.set_xticklabels([l.capitalize() for l in labels])
            ax.tick_params(axis='x', length=0)
        else:
            ax.set_yticks(x_positions)
            ax.set_yticklabels([l.capitalize() for l in labels])
            ax.tick_params(axis='y', length=0)
            
    
    
    
    def violin(self, ax, data, color, split=None, order=[], cut=2, inner=None):

        df = pd.DataFrame(columns=('x', 'y', 'hue'))
        df = df.append(data)
        hue = None
        if split:
            hue = split
            split = True
        sns.violinplot(ax=ax, x='x', y='y', hue=hue, data=df, order=order, cut=cut, inner=inner, split=split, palette=color)
        
        ax.set_xticklabels([l.capitalize() for l in order])

        means = df.groupby('x')['y'].mean()
        for i, item in enumerate(order):
            ax.plot((i-0.14, i+0.14), (means[item], means[item]), c='k')

        
#        plots = []
#        labels = []
#        for i, l in enumerate(color.keys()):
#            plots.append(Patch(facecolor=color[l], label='Synaptic '+l))
#            labels.append('Synaptic '+l)
#        
#        self._legend(ax, plots, labels)
    
    
    def swarm(self, ax, data, color, dotsize=3, stats=None):

        df = pd.DataFrame(columns=('x', 'y1', 'y2', 'y3'))
        df = df.append(data)
        
        x, y = 'x', 'y'

            
        sns.swarmplot(ax=ax, x=x, y=y, data=df, palette=color, size=dotsize, hue='x')
        
        sns.boxplot(ax=ax, x=x, y=y, data=df, showfliers=False, width=0.1, color='#444444', linewidth=Plotter.linewidth)
        
        # Set boxplot line color.
        for i,artist in enumerate(ax.artists):
            col = artist.get_facecolor()
            artist.set_edgecolor(col)
            artist.set_facecolor('None')
            for line in ax.lines:
                line.set_color(col)
                line.set_mfc(col)
                line.set_mec(col)
        
        # Remove legend title and border.
        self._legend(ax, markerscale=0.7)
        
        if stats == 'wilcoxon':
            ylim = ax.get_ylim()
            h = (ylim[1] - ylim[0]) / 50
            y = ylim[1] + h*3
            
            p = sc.stats.wilcoxon(df['y']).pvalue
            p = 'p < 10$^{' + str(int(np.log10(p))) + '}$'
                
            self._stats(ax, y, [(0, 0)], [p])
            
            
    
        
        
        
    def bar_graph(self, ax, data, colors, 
        linewidth=None, width=0.8, table=None, table_widths=None, table_rows=None, 
        legend=True, legend_title=None, rotate_ticks=False, dots=False, customstats=None,
        error=None, barerr=None, vlines=[], larval_stages=False
    ):
        
        if linewidth == None:
            linewidth = Plotter.linewidth/2,
        
        categories = []
        plots = []
        for i, (category, (xs, ys)) in enumerate(data):
            categories.append(category)
            if dots:
                if not error:
                    plots.append(plt.scatter(
                        xs, ys, color=colors[i], marker='o', zorder=1
                    ))
                if error:
                    for x, y, e, c in zip(xs, ys, error, colors[i]):
                        eb = plt.errorbar(x, y, e, lw=0, marker='o', markersize=3, color=c, elinewidth=2, clip_on=False)
                        for j, b in enumerate(eb[2]):
                            b.set_clip_on(False)
            else:
                plots.append(plt.bar(
                    xs, ys, width=width, color=colors, linewidth=linewidth, edgecolor='k', clip_on=False
                ))
                if barerr:
                    eb = plt.errorbar(xs, ys, barerr[i], lw=0, color='k', elinewidth=1)
                
        
        ax.tick_params(axis='x', pad=5, length=0)
        
        for vline in vlines:
            ax.axvline(vline, color='#EEEEEE', linestyle='dotted', linewidth=Plotter.linewidth, zorder=-1)

        if table:
                
            table = ax.table(
                cellText=table, cellLoc='center', edges='open',
                rowLabels=table_rows, colWidths=table_widths, 
                loc='bottom', fontsize=5
            )
            table.set_fontsize(5)
            table.scale(1, 0.5)
            for cell in table.properties()['children']:
                if cell._text.get_text() == table_rows[0]:
                    cell.set_text_props(fontweight='bold')
                    cell.set_fontsize(12)
                
                

        if legend:
            self._legend(ax, plots, categories, legend_title=legend_title)
    
        if customstats:
            
            xs, ys, text_x, text_y, t = customstats
            
            c = (0.5, 0.5, 0.5)
            c = 'k'
            if xs is not None:
                ax.plot(xs, ys, c=c, linewidth=Plotter.linewidth/2.0, clip_on=False, linestyle='dashed', zorder=-1)
            
            ax.text(text_x, text_y, t, ha='left', va='center', color='k')


        if larval_stages:
            larval_stage_ends = [0, 16, 25, 34, 45]
            larval_stage_mids = [8, 20.5, 29.5, 39.5, 50]
            larval_stage_labels = ['L1', 'L2', 'L3', 'L4', ' Adult']

            ax.set_xlim([0, 55])
            ax.set_xticks([0])
    
            ax.tick_params(axis='x', labelbottom=False)
            ax.set_xticks(larval_stage_mids, minor=True)
            ax.set_xticklabels(larval_stage_labels, minor=True)
            ax.tick_params(axis='x', which='minor', bottom=False, pad=1)
            linestyle = 'dotted'
            for t in larval_stage_ends:
                plt.axvline(t, color='#999999', linestyle=linestyle, linewidth=Plotter.linewidth, zorder=-1)
    
    
    
    
    
    
    def stacked_bar_graph(self, ax, pie, colors, adjacent_bars=False,
        stats=None, fancyborder=False, 
        horizontal=False, linewidth=None, relative=True, xlabels=[], 
        nospacing=False, directlegend=False, legend=True, legend_title=None, lines=[],
        legend_shift_top=0, legendpos='above', legendcol=2, legendreverse=False,
        width=0.7
    ):
                  
        categories = [n for n, ds in pie]
        pie = np.array([ds for n, ds in pie])
        data = np.copy(pie)
        
        n_datasets = len(pie[0])
    
        if relative:
            pie = pie.astype(float)/pie.sum(axis=0)
        
 
        if nospacing:
            width = 1
            
        cumulative = np.zeros(n_datasets)
        plots = []
        for i, ds in enumerate(pie):
            args = {
                'color': (colors[categories[i]] if colors else None), 
                'edgecolor': (colors[categories[i]+'_edge'] if colors and categories[i]+'_edge' in colors else 'k'),
                'linewidth': linewidth if linewidth != None else Plotter.linewidth/2.0, 'clip_on': False,
            }
            if horizontal:
                plots.append(plt.barh(
                    range(n_datasets), ds, left=cumulative, height=width, **args
                ))
                if fancyborder:
                    cs = [c[:3] for c in colors[categories[i]]]
                    plots.append(plt.barh(
                        [y+width/2 for y in range(n_datasets)], ds, left=cumulative, height=0.1, color=cs
                    ))
                    plots.append(plt.barh(
                        [y-width/2 for y in range(n_datasets)], ds, left=cumulative, height=0.1, color=cs
                    ))
                    if i == len(pie)-1:

                        plots.append(plt.barh(
                            range(n_datasets), [0]*n_datasets, left=cumulative+ds, height=width, linewidth=6, edgecolor=cs
                        ))
            else:
                plots.append(plt.bar(
                    range(n_datasets), ds, bottom=cumulative, width=width, **args
                ))
            cumulative += ds
            
        if lines:
            for (x1, y1), (x2, y2) in lines:
                ax.plot([x1, x2], [y1, y2], c='k', linewidth=Plotter.linewidth/2.0, alpha=0.3)
        
        if directlegend:
            ax.spines['left'].set_linewidth(0)
            ax.spines['bottom'].set_linewidth(0)
            ax.set_yticks([])
            ax.set_xlim(left=-width*0.5)
            y_bot = 0
            y_top = 0
            for t, vs in zip(categories, pie):
                y_top += vs[-1]
                ax.text(n_datasets - 0.3, np.mean([y_bot, y_top]), t, ha='left', va='center', color='k', size=5)
                y_bot = y_top
        elif legend:
            self._legend(ax, plots, categories, 
                position=legendpos, legend_title=legend_title, shift_top=legend_shift_top, 
                shift_right=0.014, 
                handlelength=0.8, handletextpad=0.5, ncol=legendcol, columnspacing=-4.7, reverse=legendreverse
            )
        
        if horizontal:
            ax.set_xlim(left=0)
            plt.xlim((0, 1))
        else:
            ax.set_ylim(bottom=0)
            if relative:
                plt.ylim((0, 1))
    
            if xlabels:
                ax.set_xticks(range(len(xlabels)))
                ax.set_xticklabels(xlabels)
            ax.tick_params(axis='x', length=0)
        
        if stats:
            
            y = 1.06
    
            comparisons = stats
            ps = [proportions_ztest([data[-1][x1-1], data[-1][x2-1]], [sum([ds[x1-1] for ds in data]), sum([ds[x2-1] for ds in data])])[1] for x1, x2 in comparisons]
            ps_corrected = fdrcorrection0(ps)[1]
            print(ps_corrected)
            
            self._stats(ax, y, [(x1-1, x2-1) for x1, x2 in comparisons], ps_corrected)
            


    def horizontal_bar_graph(self, ax, data, colors, pvalues=[], groups=[''], groupspacer=0.5, x_range=None):
        
        labels, correlations, error_lower, error_upper = zip(*data)
        y_pos = range(len(data))[::-1]
        
        group_patches = []
        if len(groups) > 1:
            y_pos = [y + int(y/len(groups))*groupspacer for y in y_pos]
            
            for i in range(len(groups)):
                group_patches.append(Patch(facecolor=colors[i], label=groups[i]))
            
        
        labels = labels[::len(groups)]
        labels_pos = np.array(y_pos[::len(groups)]) - 0.5*(len(groups)-1)
        
        error = [np.abs(np.array(error_lower) - correlations), np.abs(np.array(error_upper) - correlations)]
        
        
        ax.barh(y_pos, correlations, xerr=error, align='center', color=colors, ecolor='black', error_kw={'elinewidth': 1})
        ax.plot([0, 0], [min(y_pos)-1, max(y_pos)+1], color='k', linewidth=1)
        ax.set_yticks(labels_pos)
        ax.set_yticklabels(labels, size=5)
        
        ax.tick_params(length=0, axis='y')
        if x_range:
            ax.plot(x_range, (-1, -1), c='k', linewidth=Plotter.linewidth, clip_on=False)
        
        max_x = max(max(np.abs(error_lower)), max(np.abs(error_upper)))

        if len(groups) > 1:
            self._legend(ax, group_patches[::-1], groups[::-1], legend_title='Increase of:', handlelength=0.7, handletextpad=0.5, position='above')
            
            
#        ax, plots=[], labels=[], position='right', legend_title=None, reverse=True, handletextpad=None,
#        shift_right=0, shift_top=0, ncol=1, columnspacing=None, markerscale=None, handlelength=None,
#        handler_map=None
            
        if pvalues:
#            sig05 = np.where(fdrcorrection0(pvalues, alpha=0.05/2)[0])[0]
#            if len(sig05) > 0:
#                max_x += 0.05
#            for idx in sig05:
#                if correlations[idx] >= 0:
#                    x = error_upper[idx] + 0.05
#                else:
#                    x = error_lower[idx] - 0.05
#                ax.annotate(u'∗', xy=(x, y_pos[idx]), ha='center', va='center', fontname={'DejaVu Sans'}, fontsize=5)
            
            sig001 = np.where(fdrcorrection0(pvalues, alpha=0.001/2)[0])[0]
            if len(sig001) > 0:
                max_x += 0.15
            for idx in sig001:
                x = error_upper[idx] + 0.05
                ax.annotate(u'∗∗∗', xy=(x, y_pos[idx]), ha='left', va='center', fontname={'DejaVu Sans'}, fontsize=5)
                
        ax.set_ylim([-1, max(y_pos)+1])


    def _stats(self, ax, y, comparisons, pvalues, ticks=True):

        ylim = ax.get_ylim()
        h = (ylim[1] - ylim[0]) / 50.0
        
        if not ticks:
            h = 0

        for p, (x1, x2) in zip(pvalues, comparisons):
            
            if type(p) == str:
                t = p
            elif p > 0.05:
                t = 'ns'
            elif p > 0.01:
                t = '*'
            elif p > 0.001:
                t = '**'
            else:
                t = '***'

            c = 'k'
            if x1 != x2:
                ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], c=c, linewidth=Plotter.linewidth/2.0, clip_on=False)
            text_y = y if '*' in t else y + h
            
            ax.text(np.mean([x1, x2]), text_y, t, ha='center', va='bottom', color='k', fontsize=5)
            
            y += h*4
        
        ax.set_ylim(ylim)
    
          
            
    def box_plot(self, ax, box, colors, ylim=None, show_outliers=True, stats=None, darkmedian=True, vert=True):
                  

        boxprops = dict(linewidth=Plotter.linewidth/2.0, color='k')
        medianprops = dict(linewidth=Plotter.linewidth/2.0, color='#444444' if darkmedian else '#eeeeee')
        whiskerprops = dict(linewidth=Plotter.linewidth/2.0)
        
        bplot = ax.boxplot(
            box, patch_artist=True, showfliers=show_outliers,
            boxprops=boxprops, medianprops=medianprops, 
            whiskerprops=whiskerprops, capprops=whiskerprops,
            vert=vert
        )
        
        
        
        if ylim:
            if type(ylim) == tuple:
                ax.set_ylim(ylim)
            else:
                ax.set_ylim(top=ylim)
            
        for i, patch in enumerate(bplot['boxes']):
            patch.set_facecolor(colors[i])
            
        
        
        if stats:
                
            ylim = ax.get_ylim()
            h = (ylim[1] - ylim[0]) / 50.0
            y = max(max(cap.get_ydata()) for cap in bplot['caps']) + h*3
            
            if stats == 'anova':
                
                p = sc.stats.f_oneway(*box).pvalue
                self._stats(ax, y, [(1, len(box))], [p], ticks=False)
            
            else:
        
                comparisons = stats
                ps = [sc.stats.mannwhitneyu(box[x1-1], box[x2-1]).pvalue for x1, x2 in comparisons]
                ps_corrected = fdrcorrection0(ps)[1]

                print('Corrected p-values:', ps_corrected)
                print('n:', [len(b) for b in box])
                
                self._stats(ax, y, comparisons, ps_corrected)
            
        if vert:
            ax.tick_params(axis='x', length=0)
        else:
            ax.tick_params(axis='y', length=0)
        



    def scatter(self, ax, data, colors=None, line=False, crop=False, legend='', 
                legendpos='right', legendcol=1, marker='o', markersize=6, markeredgewidth=1,
                stats=None, alpha=0.8, 
                rev_legend=False, legend_shift_right=0, legend_shift_top=0, legend_title=None, 
                legendcolspace=1):
        
        
        plots = []
        if isinstance(data[0][0], (int, float, np.float32, np.float64, np.int64)):
            data = ([data[0]], [data[1]])
            colors = [colors]
            if legend:
                legend = [legend]
                
        for x, y, c in zip(data[0], data[1], colors):
            plot = ax.scatter(
                x, y, marker='o', facecolors=(c if marker == '.' else 'none'), edgecolors=c, 
                s=markersize, alpha=alpha, linewidth=(0 if marker == '.' else markeredgewidth), 
                clip_on=crop, 
            )
            plots.append(plot)
        

        xs = [x for xs in data[0] for x in xs]
        ys = [y for ys in data[1] for y in ys]
        

        if stats == 'spearman_combined':
            data = ([xs], [ys])

        if stats in ('spearman_combined', 'spearman_individual', 'regression_only'):

            for x, y, c in zip(data[0], data[1], colors):
                
                if stats == 'spearman_combined':
                    c = 'k'
                    

                reg = linregress(x, y)
                a, b = reg.slope, reg.intercept
                xlim = ax.get_xlim()
                min_x, max_x = min(xlim[0], min(xs)), max(xlim[1], max(xs))
                ax.plot([min_x, max_x], [min_x*a+b, max_x*a+b], marker=None, color=c, linewidth=Plotter.linewidth/2.0, linestyle='dashed')
                
                if stats == 'regression_only':
                    continue
                
                cor = sc.stats.spearmanr(x, y)

                r = cor.correlation
                p = cor.pvalue*len(data[0])

                text_x = min(ax.get_xlim()[1], max_x)
                text_x *= 1.02
                t = ''
                if stats == 'spearman_combined':
                    t += u'r = {0:.2f}\n'.format(r)
                if p > 0.05:
                    t += 'ns'
                elif p > 0.01:
                    t += 'p < 0.05'
                elif p > 0.001:
                    t += 'p < 0.01'
                else:
                    t += 'p < 10$^{' + str(int(np.log10(p)) if p != 0 else -9) + '}$'
                    
                print('p:', p, ', r:', r, ', n:', len(x))
                
                ax.text(text_x, text_x*a+b, t, ha='right', va='bottom', color='k', size=5)
        
        if legend:
            legend_obj = self._legend(
                ax, plots, legend, legend_title=legend_title, reverse=rev_legend, 
                shift_right=legend_shift_right, shift_top=legend_shift_top, position=legendpos, handletextpad=-0.4,
                ncol=legendcol, columnspacing=legendcolspace
            )
            
            for handles in legend_obj.legendHandles: 
                handles.set_alpha(1)

    
    def adjacency_matrix(self, ax, data, color, borders=None, post_spacers=[], pre_spacers=[], text=None):
        
        
        x = []
        y = []
        size = []
        c = []
        b = [] if borders else None
        
        data = np.array(data)
        
        cols, rows = data.shape
        
        col_size = (self.figure_size[0] * self.page_size * 64 / cols) ** 2
        
        text_xy = defaultdict(lambda: [[], []])
        
        for row in range(rows):
            for col in range(cols):
                x.append(col)
                y.append(row)
                size.append(data[col][row]*col_size+(0.5 if borders else 0))
                c.append(color[col][row])
                if borders is not None:
                    b.append(borders[col][row])
                if text is not None:
                    text_xy[text[col][row]][0].append(col)
                    text_xy[text[col][row]][1].append(row)

        ax.scatter(
            x, y, marker='s', clip_on=False, s=size, c=c, 
            linewidths=0.2 if borders else 0, edgecolors=b
        )
        
        for marker in text_xy:
            ax.scatter(*text_xy[marker], marker=marker, color='w', s=1, linewidths=0.3)
        
        
        xmin = -0.5
        xmax = cols-0.5
        hlines = [xmin, xmax]
        for s in pre_spacers:
            if s < hlines[-1]:
                hlines.insert(-1, s-0.5)
                hlines.insert(-1, s+0.5)
        for i in range(rows+1):
            for xfrom, xto in zip(hlines[::2], hlines[1::2]):
                ax.plot([xfrom, xto], [i-0.5]*2, c=(0.75, 0.75, 0.75), linewidth=0.3, clip_on=False)
        ymin = -0.5
        ymax = rows-0.5
        vlines = [ymin, ymax]
        for s in post_spacers:
            if s < vlines[-1]:
                vlines.insert(-1, s-0.5)
                vlines.insert(-1, s+0.5)
        for i in range(cols+1):
            for yfrom, yto in zip(vlines[::2], vlines[1::2]):
                ax.plot([i-0.5]*2, [yfrom, yto], c=(0.75, 0.75, 0.75), linewidth=0.3, clip_on=False)
        


        ax.set_xlim((xmin, xmax))
        ax.set_ylim((ymin, ymax))
        ax.xaxis.tick_top()
        ax.xaxis.set_tick_params(rotation=90)
        ax.xaxis.set_label_position('top')
        ax.xaxis.label.set_size(7)
        ax.yaxis.label.set_size(7)
        
 
    def venn(self, ax, data, color, labels=('Group1', 'Group2', 'Group3')):
        
        d1, d2, d3 = data
        
        v = venn3(subsets=(
            len(d1-d2-d3), 
            len(d2-d1-d3), 
            len(d1.intersection(d2)-d3), 
            len(d3-d1-d2), 
            len(d1.intersection(d3)-d2),
            len(d2.intersection(d3)-d1),
            len(d1.intersection(d2).intersection(d3)),
        ), set_labels=labels)
        
        v.get_patch_by_id('100').set_color('#aaaaaa')
        v.get_patch_by_id('010').set_color('#aaaaaa')
        v.get_patch_by_id('001').set_color('#aaaaaa')
        v.get_patch_by_id('110').set_color('#777777')
        v.get_patch_by_id('101').set_color('#777777')
        v.get_patch_by_id('011').set_color('#777777')
        v.get_patch_by_id('111').set_color('#222222')
        for pid in ('100', '010', '001', '110', '101', '011', '111'):
            v.get_patch_by_id(pid).set_linewidth(0)
            
        print(v)
        
