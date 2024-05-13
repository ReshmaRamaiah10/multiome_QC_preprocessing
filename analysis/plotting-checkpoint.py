import os
import matplotlib
import seaborn as sns
import numpy as np
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import pandas as pd

def set_colors(items, colormap_name):
    '''
    Assign colors to items based on colormap.

    Parameters:
    items (iterable) - list of items
    colormap_name (str) - a seaborn supported colormap

    Returns:
    color_dict (dict) - dict containing item:color pairs
    '''
    color_dict = {}
    sns_colormap = sns.color_palette(colormap_name, len(items))
    for i, item in enumerate(items):
        color_dict[item] = sns_colormap[i]
    return color_dict

def plt_umap_categorical(anndata, 
                         obs_col, 
                         ax,
                         umap_key = 'X_umap', 
                         palette = None,
                         title = None,
                         remove_legend = False,
                         legend_title = None,
                        size = 2):
    '''
    Plot a UMAP given where coordinates are stored and a categorical obs column to color cells by.
    
    Parameters:
    anndata (AnnData) - anndata to plot
    obs_col (str) - name of a categorical column in anndata.obs
    ax (Axes) - mpl ax to plot the figure on
    umap_key (str) - the key in anndata.obsm where 2-d coordinates are stored
    palette (dict) - a dictionary containing category:color pairs. If None, a palette will
                    be generated automatically.
    title (str) - a title for the plot
    remove_legend (bool) - if True, removes the legend from the plot
    legend_title (str) - title for the plot's legend

    '''
    # set title
    if title != None:
        ax.set(title=title)

    # get coords
    data = anndata.obs.copy()
    data['umap1'] = np.array(pd.DataFrame(anndata.obsm[umap_key]).iloc[:, 0])
    data['umap2'] = np.array(pd.DataFrame(anndata.obsm[umap_key]).iloc[:, 1])
    
    # gen palette if none
    if palette == None:
        palette = set_colors(data[obs_col].unique(), 'husl')
    
    # scatterplot
    sns.scatterplot(data,
                    x = 'umap1', 
                    y = 'umap2',
                    hue = obs_col,
                    palette = palette,
                    s = size,
                    linewidth = 0,
                    ax = ax)

    # remove axis lines
    sns.despine(ax=ax, left=True, bottom=True) 
    ax.tick_params(left=False, bottom=False)
    ax.set(xticklabels=[], yticklabels=[], xlabel = '', ylabel = '')

    # legend details
    if remove_legend:
        leg = ax.legend()
        leg.remove()
    else:
        leg = ax.legend(bbox_to_anchor=(1.02, 0.85), loc='upper left')
        leg.get_frame().set_linewidth(0)
        if legend_title != None:
            leg.set_title(legend_title)
            
def plt_umap_numerical(anndata, 
                         colorval, 
                         fig,
                         ax,
                         umap_key = 'X_umap', 
                         palette = 'viridis',
                         colorbar = True,
                         title = None,
                         remove_legend = False,
                         legend_title = None,
                         size = 2):
    '''
    Plot a UMAP given where coordinates are stored and a numerical obs column to color cells by.
    
    Parameters:
    anndata (AnnData) - anndata to plot
    colorval (str or iterable) - name of a numerical column in anndata.obs or list of values to color by
    fig (Figure) - mpl fig to plot the figure on
    ax (Axes) - mpl ax to plot the figure on
    umap_key (str) - the key in anndata.obsm where 2-d coordinates are stored
    palette (str) - a string indicating a matplotlib palette.
    colorbar (bool) - whether to include a colorbar key
    title (str) - a title for the plot
    legend_title (str) - title for the plot's legend
    '''
    # set title
    if title != None:
        ax.set(title=title)

    # get_coords
    if type(colorval) == str:
        if colorval in anndata.obs.columns:
            color = anndata.obs[colorval]
        else: 
            color = anndata[:,colorval].X.toarray().flatten()
    else: color = colorval
    umap1 = pd.DataFrame(anndata.obsm[umap_key]).iloc[:,0]
    umap2 = pd.DataFrame(anndata.obsm[umap_key]).iloc[:,1]
    
    #scatterplot
    scatterplt = ax.scatter(umap1,umap2,c=color,s=size)

    #add colorbar
    if colorbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cbar = fig.colorbar(scatterplt, cax,
                            ticks=[min(color),max(color)])
        for t in cbar.ax.get_yticklabels():
            t.set_fontsize(12)

    # remove axis lines
    ax.set_axis_off()
    fig.tight_layout()

def plt_scatter(adata,
                obs_col1,
                obs_col2,
                ax,
                hue = None,
                hue_order = None,
                palette = None,
                log = False):
    
    # get x and y values
    x = adata.obs[obs_col1]
    y = adata.obs[obs_col2]
    if log:
        x = np.log10(x+1)
        y = np.log10(y+1)
        
    # plot scatter
    if hue is None:
        sns.scatterplot(x=x,
                        y=y,
                        s=5,
                        linewidth = 0,
                        ax = ax)
    else:
        sns.scatterplot(x=x,
                        y=y,
                        s=5,
                        linewidth = 0,
                        ax = ax,
                       hue = hue,
                       hue_order = hue_order,
                       palette = palette)
    
    # make pretty
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    if log:
        ax.set(xlabel = 'log10(' + obs_col1 + ' + 1)',
               ylabel = 'log10(' + obs_col1 + ' + 1)')