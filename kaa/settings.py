"""
Simple settings file for the in-and-outs of Kaa.
"""
class KaaSettings:
    'Should we try to parallelize the generator calculations?'
    use_parallel = False

class PlotSettings:
    'Fonts for the indices on matplotlib plots'
    plot_font = 15

    'Toggle to save the figures to disk'
    save_fig = False

    'Path to save figures'
    fig_path = "/Users/edwardkim/Work/kaa-optimize/figures"
